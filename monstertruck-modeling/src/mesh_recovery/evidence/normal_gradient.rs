//! Two-direction normal-gradient evidence from actual source support positions
//! and normals. A one-dimensional chord strip cannot establish both local
//! principal curvatures, regardless of how accurately it fits a cylinder.

use super::{super::MeshRecoveryError, MeshEvidence};
use crate::{InnerSpace, Point3, Vector3};
use rustc_hash::FxHashSet as HashSet;

/// Explicit numerical bounds for discovering smooth source-support neighbors.
/// Sharp rejected interfaces remain intact in the original source graph.
#[derive(Clone, Copy, Debug)]
pub struct SourceNormalGradientOptions {
    /// Largest source facet-normal jump used as a local gradient observation.
    pub max_neighbor_normal_turn_radians: f64,
    /// Smallest allowed two-direction tangent sample covariance rank ratio.
    pub min_tangent_direction_rank_ratio: f64,
}

impl Default for SourceNormalGradientOptions {
    fn default() -> Self {
        Self {
            max_neighbor_normal_turn_radians: 0.5,
            min_tangent_direction_rank_ratio: 1.0e-6,
        }
    }
}

impl SourceNormalGradientOptions {
    fn valid(self) -> bool {
        self.max_neighbor_normal_turn_radians.is_finite()
            && (0.0..std::f64::consts::PI).contains(&self.max_neighbor_normal_turn_radians)
            && self.min_tangent_direction_rank_ratio.is_finite()
            && (0.0..=1.0).contains(&self.min_tangent_direction_rank_ratio)
    }
}

/// Symmetric source tangent normal derivative (S = -dn/dp), not a certified
/// source CAD curvature tensor. Principal values are signed and in 1/length.
#[derive(Clone, Debug)]
pub struct SourceNormalGradient {
    /// Canonical original source support ID.
    pub support_index: usize,
    /// Distinct source-neighbor supports with original two-sided interfaces.
    pub neighboring_supports: Vec<usize>,
    /// Mean original polygon vertices, calculated in a relative frame.
    pub source_centroid: Point3,
    /// Original support normal, not an analytic CAD source surface normal.
    pub normal: Vector3,
    /// Smallest/largest tangent sample direction covariance eigenvalue ratio.
    pub transverse_position_rank_ratio: f64,
    /// Signed principal eigenvalues of the approximate normal derivative.
    pub principal_normal_curvatures: [f64; 2],
    /// Maximum original normal-step residual of the first-order symmetric fit.
    pub max_normal_change_residual: f64,
    /// RMS original normal-step residual of the fit.
    pub rms_normal_change_residual: f64,
}

/// Lazily measure local normal gradients with source support centroids
/// computed once, independently of original native CAD primitive families.
pub struct SourceNormalGradientIndex<'a> {
    evidence: &'a MeshEvidence,
    centroids: Vec<Option<Point3>>,
    options: SourceNormalGradientOptions,
}

fn support_centroid(evidence: &MeshEvidence, id: usize) -> Option<Point3> {
    let support = evidence.planar_supports.get(id)?;
    let mut vertices = HashSet::default();
    for &polygon in &support.source_face_indices {
        vertices.extend(evidence.source_polygon_vertices(polygon)?.iter().copied());
    }
    let mut indices = vertices.into_iter().collect::<Vec<_>>();
    indices.sort_unstable();
    let reference = *evidence.positions.get(*indices.first()?)?;
    let mut sum = Vector3::new(0.0, 0.0, 0.0);
    for &index in &indices {
        sum += *evidence.positions.get(index)? - reference;
    }
    let point = reference + sum / crate::usize_to_f64(indices.len());
    [point.x, point.y, point.z]
        .into_iter()
        .all(f64::is_finite)
        .then_some(point)
}

fn tangent_basis(n: Vector3) -> (Vector3, Vector3) {
    let ref_axis = if n.z.abs() < 0.75 {
        Vector3::unit_z()
    } else {
        Vector3::unit_x()
    };
    let u = n.cross(ref_axis).normalize();
    (u, n.cross(u))
}

/// Each observation contains two tangent displacement components followed
/// by the opposing original support-normal change in the same coordinates.
/// Solve for the three independent symmetric shape-operator coefficients.
fn fit_symmetric_normal_gradient(
    data: &[[f64; 4]],
    min_tangent_rank: f64,
) -> Option<([f64; 3], f64, f64, f64)> {
    if data.len() < 2 {
        return None;
    }
    let (mut xx, mut xy, mut yy) = (0.0, 0.0, 0.0);
    let mut gram = [[0.0_f64; 3]; 3];
    let mut rhs = [0.0_f64; 3];
    for &[dx, dy, nx, ny] in data {
        let magnitude = dx.hypot(dy);
        if !magnitude.is_finite() || magnitude <= 0.0 {
            return None;
        }
        let x = dx / magnitude;
        let y = dy / magnitude;
        xx = x.mul_add(x, xx);
        xy = x.mul_add(y, xy);
        yy = y.mul_add(y, yy);
        for (row, target) in [([x, y, 0.0], nx / magnitude), ([0.0, x, y], ny / magnitude)] {
            for i in 0..3 {
                rhs[i] = row[i].mul_add(target, rhs[i]);
                for j in 0..3 {
                    gram[i][j] = row[i].mul_add(row[j], gram[i][j]);
                }
            }
        }
    }
    let maximum = 0.5 * (xx + yy + (xx - yy).hypot(2.0 * xy));
    let det = xx.mul_add(yy, -(xy * xy));
    if !maximum.is_finite() || maximum <= 0.0 || !det.is_finite() {
        return None;
    }
    let rank_ratio = det / (maximum * maximum);
    if rank_ratio <= min_tangent_rank {
        return None;
    }
    let scale = (0..3).map(|i| gram[i][i]).fold(0.0_f64, f64::max);
    for pivot in 0..3 {
        let best =
            (pivot..3).max_by(|&i, &j| gram[i][pivot].abs().total_cmp(&gram[j][pivot].abs()))?;
        gram.swap(pivot, best);
        rhs.swap(pivot, best);
        let diag = gram[pivot][pivot];
        if !diag.is_finite() || diag.abs() <= 1.0e-10 * scale {
            return None;
        }
        for lower in pivot + 1..3 {
            let factor = gram[lower][pivot] / diag;
            for j in pivot..3 {
                gram[lower][j] = factor.mul_add(-gram[pivot][j], gram[lower][j]);
            }
            rhs[lower] = factor.mul_add(-rhs[pivot], rhs[lower]);
        }
    }
    let mut fit = [0.0; 3];
    for i in (0..3).rev() {
        let rest = ((i + 1)..3).map(|j| gram[i][j] * fit[j]).sum::<f64>();
        fit[i] = (rhs[i] - rest) / gram[i][i];
    }
    if !fit.into_iter().all(f64::is_finite) {
        return None;
    }
    let (peak, squares) = data
        .iter()
        .fold((0.0_f64, 0.0_f64), |(peak, sum), &[x, y, nx, ny]| {
            let residual =
                (fit[0].mul_add(x, fit[1] * y) - nx).hypot(fit[1].mul_add(x, fit[2] * y) - ny);
            (peak.max(residual), residual.mul_add(residual, sum))
        });
    Some((
        fit,
        rank_ratio,
        peak,
        (squares / crate::usize_to_f64(data.len())).sqrt(),
    ))
}

impl<'a> SourceNormalGradientIndex<'a> {
    /// Index all original source vertices once for reusable local measurements.
    ///
    /// # Errors
    /// Invalid numerical options report a mesh evidence error.
    pub fn new(
        evidence: &'a MeshEvidence,
        options: SourceNormalGradientOptions,
    ) -> Result<Self, MeshRecoveryError> {
        if !options.valid() {
            return Err(MeshRecoveryError::InvalidEvidenceOptions);
        }
        let centroids = (0..evidence.planar_supports.len())
            .map(|id| support_centroid(evidence, id))
            .collect();
        Ok(Self {
            evidence,
            centroids,
            options,
        })
    }

    /// Measure local symmetric normal derivative only with two resolved source
    /// tangent directions. Underdetermined or invalid support returns None.
    #[must_use]
    pub fn measure(&self, support: usize) -> Option<SourceNormalGradient> {
        let anchor = self.evidence.planar_supports.get(support)?;
        let center = (*self.centroids.get(support)?)?;
        let (basis_u, basis_v) = tangent_basis(anchor.normal);
        let mut neighbors = HashSet::default();
        for &id in self.evidence.support_interface_indices.get(support)? {
            let edge = self.evidence.support_interfaces.get(id)?;
            let [a, b] = edge.supports;
            let other = if a == support {
                b
            } else if b == support {
                a
            } else {
                continue;
            };
            neighbors.insert(other);
        }
        let mut ids = neighbors.into_iter().collect::<Vec<_>>();
        ids.sort_unstable();
        let mut observations = Vec::new();
        let mut accepted = Vec::new();
        for id in ids {
            let Some(point) = self.centroids.get(id).copied().flatten() else {
                continue;
            };
            let other = self.evidence.planar_supports.get(id)?.normal;
            let normal_turn = anchor
                .normal
                .cross(other)
                .magnitude()
                .atan2(anchor.normal.dot(other));
            if !normal_turn.is_finite()
                || normal_turn > self.options.max_neighbor_normal_turn_radians
            {
                continue;
            }
            let displacement = point - center;
            let normal_change = anchor.normal - other;
            let obs = [
                displacement.dot(basis_u),
                displacement.dot(basis_v),
                normal_change.dot(basis_u),
                normal_change.dot(basis_v),
            ];
            if obs[0].hypot(obs[1]) <= 0.0 || !obs.into_iter().all(f64::is_finite) {
                continue;
            }
            observations.push(obs);
            accepted.push(id);
        }
        let ([first_diagonal, cross_component, second_diagonal], rank, maximum, rms) =
            fit_symmetric_normal_gradient(
                &observations,
                self.options.min_tangent_direction_rank_ratio,
            )?;
        let midpoint = first_diagonal.midpoint(second_diagonal);
        let half_gap = 0.5 * (first_diagonal - second_diagonal).hypot(2.0 * cross_component);
        let values = [midpoint - half_gap, midpoint + half_gap];
        if !values.into_iter().all(f64::is_finite) {
            return None;
        }
        Some(SourceNormalGradient {
            support_index: support,
            neighboring_supports: accepted,
            source_centroid: center,
            normal: anchor.normal,
            transverse_position_rank_ratio: rank,
            principal_normal_curvatures: values,
            max_normal_change_residual: maximum,
            rms_normal_change_residual: rms,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::fit_symmetric_normal_gradient;
    #[test]
    fn two_independent_source_normal_derivatives_resolve_symmetric_tensor() {
        let samples = [
            [1.0, 0.0, 0.5, 0.1],
            [0.0, 1.0, 0.1, 0.25],
            [-1.0, 0.0, -0.5, -0.1],
            [0.0, -1.0, -0.1, -0.25],
        ];
        let (values, rank, peak, rms) = fit_symmetric_normal_gradient(&samples, 1.0e-6).unwrap();
        assert!((values[0] - 0.5).abs() < 1.0e-12);
        assert!((values[1] - 0.1).abs() < 1.0e-12);
        assert!((values[2] - 0.25).abs() < 1.0e-12);
        assert!((rank - 1.0).abs() < 1.0e-12);
        assert!(peak < 1.0e-12 && rms < 1.0e-12);
    }

    #[test]
    fn one_resolved_tangent_direction_is_not_a_second_curvature() {
        let samples = [
            [1.0, 0.0, 0.5, 0.0],
            [-1.0, 0.0, -0.5, 0.0],
            [2.0, 0.0, 1.0, 0.0],
        ];
        assert!(fit_symmetric_normal_gradient(&samples, 1.0e-6).is_none());
    }
}

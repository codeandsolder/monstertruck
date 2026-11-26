//! Source-uncertainty-aware joint analytic reconciliation. The fitting
//! allowance is not a native incidence tolerance: the final shared vertex
//! must satisfy every incident carrier at kernel resolution.
use super::{
    MeshEvidence, ReconciliationOptions, SurfaceCarrier, SurfaceRegionPartition,
    reconcile::{preserves_polygon_orientation, reconcile_point},
};
use crate::{InnerSpace, KERNEL_EPSILON, Point3, Vector3};
use std::collections::BTreeSet;

#[derive(Clone, Copy)]
pub(in crate::mesh_recovery) struct JointFitOptions {
    pub maximum_adjustment: f64,
    pub declared_uncertainty: f64,
    pub angular_resolution: f64,
}

pub(in crate::mesh_recovery) struct JointFitCandidate {
    pub partition: SurfaceRegionPartition,
    pub positions: Vec<Point3>,
    pub allowances: Vec<f64>,
    pub maximum_displacement: f64,
    pub refitted_regions: usize,
    pub max_additional_carrier_residual: f64,
}

struct VertexSample {
    id: usize,
    regions: Vec<usize>,
    allowance: f64,
    boundary: bool,
}

struct Neighborhood<'a> {
    evidence: &'a MeshEvidence,
    partition: &'a SurfaceRegionPartition,
    samples: Vec<VertexSample>,
    supports: BTreeSet<usize>,
    options: JointFitOptions,
}

impl<'a> Neighborhood<'a> {
    fn new(
        evidence: &'a MeshEvidence,
        partition: &'a SurfaceRegionPartition,
        regions: &[usize],
        options: JointFitOptions,
    ) -> Option<Self> {
        let mut incidence = vec![Vec::new(); evidence.positions.len()];
        let mut active = vec![false; evidence.positions.len()];
        let mut allowances = vec![0.0_f64; evidence.positions.len()];
        let mut supports = BTreeSet::new();
        for (support, source) in evidence.planar_supports.iter().enumerate() {
            let region = *partition.region_by_support.get(support)?;
            for &polygon in &source.source_face_indices {
                for &id in evidence.source_polygon_vertices(polygon)? {
                    incidence[id].push(region);
                    active[id] |= regions.contains(&region);
                }
            }
        }
        // Only original edges incident to the affected neighborhood
        // supply uncertainty. Do not spread a coarse tessellation's
        // allowance to an unrelated small feature elsewhere in the mesh.
        for (support, source) in evidence.planar_supports.iter().enumerate() {
            let region = partition.region_by_support[support];
            let carrier = partition.carriers[region]?;
            for &polygon in &source.source_face_indices {
                let ids = evidence.source_polygon_vertices(polygon)?;
                if !ids.iter().any(|&id| active[id]) {
                    continue;
                }
                supports.insert(support);
                for (&a, &b) in ids.iter().zip(ids.iter().cycle().skip(1)) {
                    let (p, q) = (evidence.positions[a], evidence.positions[b]);
                    let Some(residual) = carrier.measure_segment(p, q) else {
                        continue;
                    };
                    // Use the measured sagitta/fit residual, never a loose
                    // analytic upper bound, as evidence for an allowance.
                    for id in [a, b] {
                        allowances[id] = allowances[id].max(residual.witnessed_max);
                    }
                }
            }
        }
        let mut boundary_vertices = vec![false; evidence.positions.len()];
        for path in &partition.contacts.paths {
            for edge in &path.edges {
                for id in edge.vertices {
                    boundary_vertices[id] = true;
                }
            }
        }
        let mut samples = Vec::new();
        for (id, mut regions) in incidence.into_iter().enumerate() {
            if !active[id] {
                continue;
            }
            regions.sort_unstable();
            regions.dedup();
            let p = evidence.positions[id];
            let roundoff = 32.0 * f64::EPSILON * p.x.abs().max(p.y.abs()).max(p.z.abs());
            let allowance = allowances[id]
                .max(options.declared_uncertainty)
                .max(roundoff)
                .max(KERNEL_EPSILON)
                .min(options.maximum_adjustment);
            let boundary = boundary_vertices[id] || regions.len() > 1;
            samples.push(VertexSample {
                id,
                regions,
                allowance,
                boundary,
            });
        }
        (!samples.is_empty() && samples.len() <= 4096).then_some(Self {
            evidence,
            partition,
            samples,
            supports,
            options,
        })
    }

    fn exact_positions(
        &self,
        carriers: &[Option<SurfaceCarrier>],
        prior: &[Point3],
    ) -> Option<(Vec<Point3>, Vec<f64>, f64)> {
        let mut positions = prior.to_vec();
        let mut allowances = vec![0.0; prior.len()];
        let tolerance = KERNEL_EPSILON * 0.01;
        let mut max_displacement = 0.0_f64;
        for sample in &self.samples {
            let source = self.evidence.positions[sample.id];
            let local = sample
                .regions
                .iter()
                .map(|&id| carriers[id])
                .collect::<Option<Vec<_>>>()?;
            // Interior source samples measure the carrier fit. They are
            // not native topological vertices and need not lie exactly
            // on the reconstructed surface. Only actual shared/open
            // boundary junctions impose exact incidence constraints.
            if !sample.boundary {
                if !local.iter().all(|carrier| {
                    carrier.signed_constraint(source).map_or_else(
                        || {
                            carrier
                                .measure_segment(source, source)
                                .is_some_and(|r| r.witnessed_max <= sample.allowance)
                        },
                        |(_, residual)| residual.abs() <= sample.allowance,
                    )
                }) {
                    return None;
                }
                allowances[sample.id] = sample.allowance;
                continue;
            }
            let point = if local.iter().all(|carrier| {
                carrier
                    .measure_segment(source, source)
                    .is_some_and(|r| r.upper_bound <= tolerance)
            }) {
                // At a cone apex the implicit normal is undefined, but
                // the existing point may already satisfy every carrier.
                source
            } else {
                reconcile_point(
                    source,
                    &local,
                    ReconciliationOptions {
                        max_displacement: sample.allowance,
                        parallel_angle_tolerance: self.options.angular_resolution,
                        max_normal_angle: self.options.angular_resolution * 2.0,
                        intersection_tolerance: tolerance,
                    },
                )?
            };
            let displacement = (point - source).magnitude();
            if !displacement.is_finite() || displacement > sample.allowance {
                return None;
            }
            positions[sample.id] = point;
            allowances[sample.id] = sample.allowance;
            max_displacement = max_displacement.max(displacement);
        }
        preserves_polygon_orientation(self.evidence, &positions, &self.supports).then_some((
            positions,
            allowances,
            max_displacement,
        ))
    }
}

// Rigid carrier placement plus independent meridian parameters. Rotations
// are scaled by the actual region span so all parameters have length units.
#[derive(Clone, Copy)]
struct CarrierModel {
    base: SurfaceCarrier,
    pivot: Point3,
    scale: f64,
    allowance: f64,
    parameters: [f64; 8],
}

impl CarrierModel {
    fn new(neighborhood: &Neighborhood<'_>, region: usize) -> Option<Self> {
        let mut samples = neighborhood
            .samples
            .iter()
            .filter(|s| s.regions.contains(&region));
        let first = samples.next()?;
        let pivot = neighborhood.evidence.positions[first.id];
        let (mut scale, mut allowance) = (0.0_f64, first.allowance);
        for sample in samples {
            scale = scale.max((neighborhood.evidence.positions[sample.id] - pivot).magnitude());
            allowance = allowance.max(sample.allowance);
        }
        let base = match neighborhood.partition.carriers[region]? {
            SurfaceCarrier::LinearRevolution {
                origin,
                axis,
                radius_at_origin,
                slope,
            } => {
                let z = (pivot - origin).dot(axis);
                SurfaceCarrier::LinearRevolution {
                    origin: origin + axis * z,
                    axis,
                    radius_at_origin: slope.mul_add(z, radius_at_origin),
                    slope,
                }
            }
            other => other,
        };
        Some(Self {
            base,
            pivot,
            scale: scale.max(allowance).max(KERNEL_EPSILON),
            allowance,
            parameters: [0.0; 8],
        })
    }

    fn carrier(self) -> SurfaceCarrier {
        let p = self.parameters;
        let translation = Vector3::new(p[0], p[1], p[2]);
        let rotation = Vector3::new(p[3], p[4], p[5]) / self.scale;
        let rotate = |v: Vector3| {
            let angle = rotation.magnitude();
            if angle == 0.0 {
                return v;
            }
            let axis = rotation / angle;
            let (sin, cos) = angle.sin_cos();
            v * cos + axis.cross(v) * sin + axis * (axis.dot(v) * (1.0 - cos))
        };
        let place = |v: Point3| self.pivot + rotate(v - self.pivot) + translation;
        match self.base {
            SurfaceCarrier::Plane { origin, normal } => SurfaceCarrier::Plane {
                origin: place(origin),
                normal: rotate(normal),
            },
            SurfaceCarrier::LinearRevolution {
                origin,
                axis,
                radius_at_origin,
                slope,
            } => SurfaceCarrier::LinearRevolution {
                origin: place(origin),
                axis: rotate(axis),
                radius_at_origin: radius_at_origin + p[6],
                slope: slope + p[7] / self.scale,
            },
            SurfaceCarrier::CircularRevolution {
                axial_center,
                axis,
                radial_center,
                profile_radius,
            } => SurfaceCarrier::CircularRevolution {
                axial_center: place(axial_center),
                axis: rotate(axis),
                radial_center: if radial_center == 0.0 {
                    0.0
                } else {
                    radial_center + p[7]
                },
                profile_radius: profile_radius + p[6],
            },
        }
    }

    // Numerical derivatives of the normalized signed distance. The
    // carrier family and sphere/torus distinction are held fixed.
    fn derivatives(self, point: Point3) -> Option<[f64; 8]> {
        let step = (self.scale * f64::EPSILON.sqrt()).max(KERNEL_EPSILON * 0.01);
        let mut row = [0.0; 8];
        for (i, value) in row.iter_mut().enumerate() {
            let (mut plus, mut minus) = (self, self);
            plus.parameters[i] += step;
            minus.parameters[i] -= step;
            let a = plus.carrier().signed_constraint(point)?.1;
            let b = minus.carrier().signed_constraint(point)?.1;
            *value = (a - b) / (2.0 * step);
        }
        Some(row)
    }

    // Bound the change in normalized carrier residual over the entire
    // convex hull of the original region samples, not just the vertices
    // used by the optimizer. Each signed profile distance is 1-Lipschitz
    // under rigid motion (radial norm and meridian norm are contractions).
    fn additional_residual_bound(self) -> f64 {
        let p = self.parameters;
        let translation = Vector3::new(p[0], p[1], p[2]).magnitude();
        let rotation = Vector3::new(p[3], p[4], p[5]).magnitude();
        // ||R p-p|| <= angle*||p||, with all source points in the
        // pivot-centered ball of radius scale.
        let motion = translation + rotation;
        match self.base {
            SurfaceCarrier::Plane { .. } => motion,
            SurfaceCarrier::CircularRevolution { radial_center, .. } => {
                motion
                    + p[6].abs()
                    + if radial_center == 0.0 {
                        0.0
                    } else {
                        p[7].abs()
                    }
            }
            SurfaceCarrier::LinearRevolution {
                origin,
                axis,
                radius_at_origin,
                slope,
            } => {
                let delta_slope = p[7] / self.scale;
                let old_norm = 1.0_f64.hypot(slope);
                let new_norm = 1.0_f64.hypot(slope + delta_slope);
                // Axis reference is centered at the pivot's axial
                // coordinate; include the actual coordinate roundoff.
                let delta = self.pivot - origin;
                let z = delta.dot(axis);
                let rho = (delta - axis * z).magnitude();
                let pivot_residual = (rho - slope.mul_add(z, radius_at_origin)).abs() / old_norm;
                let extent = self.scale + motion;
                let profile = (p[6].abs() + delta_slope.abs() * (extent + z.abs())) / new_norm;
                let normalization = (old_norm / new_norm - 1.0).abs() * (pivot_residual + extent);
                motion + profile + normalization
            }
        }
    }

    fn bounded(self) -> bool {
        self.additional_residual_bound() <= self.allowance
            && match self.carrier() {
                SurfaceCarrier::CircularRevolution {
                    radial_center,
                    profile_radius,
                    ..
                } => {
                    (radial_center == 0.0 || radial_center > KERNEL_EPSILON)
                        && profile_radius > KERNEL_EPSILON
                }
                _ => true,
            }
    }
}

// Cholesky for the small positive definite systems. All systems are
// explicitly damped; nonfinite or singular candidates are rejected.
fn solve(mut matrix: Vec<Vec<f64>>, mut rhs: Vec<f64>) -> Option<Vec<f64>> {
    let n = rhs.len();
    for i in 0..n {
        for j in 0..=i {
            let mut value = matrix[i][j];
            for k in 0..j {
                value = matrix[i][k].mul_add(-matrix[j][k], value);
            }
            if i == j {
                if value <= 0.0 || !value.is_finite() {
                    return None;
                }
                matrix[i][j] = value.sqrt();
            } else {
                matrix[i][j] = value / matrix[j][j];
            }
        }
        for j in 0..i {
            rhs[i] = matrix[i][j].mul_add(-rhs[j], rhs[i]);
        }
        rhs[i] /= matrix[i][i];
    }
    for i in (0..n).rev() {
        for j in i + 1..n {
            rhs[i] = matrix[j][i].mul_add(-rhs[j], rhs[i]);
        }
        rhs[i] /= matrix[i][i];
    }
    rhs.iter().all(|v| v.is_finite()).then_some(rhs)
}

// The per-vertex factor is reused for every active parameter column.
// Using fixed arrays avoids allocating and factoring a 3x3 matrix dozens
// of times per source vertex.
struct Cholesky3 {
    diagonal: [f64; 3],
    lower: [f64; 3],
}

impl Cholesky3 {
    fn new(a: [[f64; 3]; 3]) -> Option<Self> {
        let d0 = a[0][0].sqrt();
        let l10 = a[1][0] / d0;
        let l20 = a[2][0] / d0;
        let d1 = l10.mul_add(-l10, a[1][1]).sqrt();
        let l21 = l20.mul_add(-l10, a[2][1]) / d1;
        let d2 = l21.mul_add(-l21, l20.mul_add(-l20, a[2][2])).sqrt();
        let diagonal = [d0, d1, d2];
        diagonal
            .iter()
            .all(|&v| v.is_finite() && v > 0.0)
            .then_some(Self {
                diagonal,
                lower: [l10, l20, l21],
            })
    }

    fn solve(&self, rhs: [f64; 3]) -> [f64; 3] {
        let [d0, d1, d2] = self.diagonal;
        let [l10, l20, l21] = self.lower;
        let y0 = rhs[0] / d0;
        let y1 = l10.mul_add(-y0, rhs[1]) / d1;
        let y2 = l21.mul_add(-y1, l20.mul_add(-y0, rhs[2])) / d2;
        let x2 = y2 / d2;
        let x1 = l21.mul_add(-x2, y1) / d1;
        let x0 = l20.mul_add(-x2, l10.mul_add(-x1, y0)) / d0;
        [x0, x1, x2]
    }
}

fn accumulate(
    normal: &mut [Vec<f64>],
    rhs: &mut [f64],
    columns: &[usize],
    jacobian: &[f64],
    error: f64,
    weight: f64,
) {
    for (&i, &value) in columns.iter().zip(jacobian) {
        rhs[i] -= weight * value * error;
        for (&k, &other) in columns.iter().zip(jacobian) {
            normal[i][k] += weight * value * other;
        }
    }
}

struct ConstraintRow {
    normal: Vector3,
    residual: f64,
    jacobian: Vec<f64>,
}

fn sample_rows(
    neighborhood: &Neighborhood<'_>,
    sample: &VertexSample,
    models: &[(usize, CarrierModel)],
) -> Option<(Vec<usize>, Vec<ConstraintRow>)> {
    let point = neighborhood.evidence.positions[sample.id];
    let active = models
        .iter()
        .enumerate()
        .filter(|(_, m)| sample.regions.contains(&m.0))
        .map(|(i, _)| i)
        .collect::<Vec<_>>();
    let columns = active
        .iter()
        .flat_map(|&i| i * 8..(i + 1) * 8)
        .collect::<Vec<_>>();
    let mut rows = Vec::new();
    for &region in &sample.regions {
        let local = active.iter().position(|&i| models[i].0 == region);
        let carrier = local.map_or(neighborhood.partition.carriers[region]?, |i| {
            models[active[i]].1.carrier()
        });
        let Some((normal, residual)) = carrier.signed_constraint(point) else {
            continue;
        };
        let mut jacobian = vec![0.0; columns.len()];
        if let Some(i) = local {
            jacobian[i * 8..(i + 1) * 8].copy_from_slice(&models[active[i]].1.derivatives(point)?);
        }
        rows.push(ConstraintRow {
            normal,
            residual,
            jacobian,
        });
    }
    Some((columns, rows))
}

fn project_sample(
    normal: &mut [Vec<f64>],
    rhs: &mut [f64],
    columns: &[usize],
    rows: &[ConstraintRow],
    weight: f64,
    anchor_weight: f64,
) -> Option<()> {
    let mut local = [[0.0; 3]; 3];
    let mut residual = [0.0; 3];
    let mut cross = vec![[0.0; 3]; columns.len()];
    for (i, row) in local.iter_mut().enumerate() {
        row[i] = anchor_weight;
    }
    for row in rows {
        let n = row.normal;
        for (i, entries) in local.iter_mut().enumerate() {
            residual[i] = n[i].mul_add(row.residual, residual[i]);
            for (k, entry) in entries.iter_mut().enumerate() {
                *entry = n[i].mul_add(n[k], *entry);
            }
        }
        for (entry, &j) in cross.iter_mut().zip(&row.jacobian) {
            for (k, value) in entry.iter_mut().enumerate() {
                *value = j.mul_add(n[k], *value);
            }
        }
    }
    let factor = Cholesky3::new(local)?;
    let projected = factor.solve(residual);
    let inverse = cross
        .into_iter()
        .map(|row| factor.solve(row))
        .collect::<Vec<_>>();
    // Accumulate the projected Jacobian as a Gram matrix, avoiding
    // cancellation of nearly equal normal equations at constrained vertices.
    for row in rows {
        let project = |v: [f64; 3]| row.normal.dot(Vector3::new(v[0], v[1], v[2]));
        let error = row.residual - project(projected);
        let jacobian = row
            .jacobian
            .iter()
            .zip(&inverse)
            .map(|(j, &shift)| j - project(shift))
            .collect::<Vec<_>>();
        accumulate(normal, rhs, columns, &jacobian, error, weight);
    }
    for (k, &error) in projected.iter().enumerate() {
        let jacobian = inverse.iter().map(|row| row[k]).collect::<Vec<_>>();
        accumulate(
            normal,
            rhs,
            columns,
            &jacobian,
            error,
            weight * anchor_weight,
        );
    }
    Some(())
}

// Variable projection eliminates independent 3D vertex displacements;
// the dense system contains only carrier parameters. Each source vertex
// contributes only the columns of its incident adjustable carriers.
fn joint_step(
    neighborhood: &Neighborhood<'_>,
    models: &[(usize, CarrierModel)],
    anchor_weight: f64,
) -> Option<Vec<f64>> {
    let size = models.len() * 8;
    let mut normal = vec![vec![0.0; size]; size];
    let mut rhs = vec![0.0; size];
    for sample in &neighborhood.samples {
        let (columns, rows) = sample_rows(neighborhood, sample, models)?;
        let weight = (neighborhood.options.maximum_adjustment / sample.allowance).powi(2);
        if sample.boundary {
            project_sample(
                &mut normal,
                &mut rhs,
                &columns,
                &rows,
                weight,
                anchor_weight,
            )?;
        } else {
            // Source fit is a soft objective, not a native equality.
            for row in &rows {
                accumulate(
                    &mut normal,
                    &mut rhs,
                    &columns,
                    &row.jacobian,
                    row.residual,
                    weight * anchor_weight,
                );
            }
        }
    }
    for (i, (_, model)) in models.iter().enumerate() {
        for (j, p) in model.parameters.iter().enumerate() {
            let k = i * 8 + j;
            let damping = (anchor_weight * 0.001 + 1.0e-10) * (1.0 + normal[k][k].abs());
            normal[k][k] += damping;
            rhs[k] -= damping * p;
        }
    }
    solve(normal, rhs)
}

pub(in crate::mesh_recovery) fn fit_neighborhood(
    evidence: &MeshEvidence,
    partition: &SurfaceRegionPartition,
    prior: &[Point3],
    regions: &[usize],
    options: JointFitOptions,
    refit: bool,
) -> Option<JointFitCandidate> {
    let neighborhood = Neighborhood::new(evidence, partition, regions, options)?;
    let mut partition = partition.clone();
    let mut refitted_regions = 0;
    let mut max_additional_carrier_residual = 0.0_f64;
    if refit {
        // Keep a bounded local system, independent of the total model size.
        let mut models = regions
            .iter()
            .map(|&id| Some((id, CarrierModel::new(&neighborhood, id)?)))
            .collect::<Option<Vec<_>>>()?;
        for anchor_weight in [
            1.0e-2, 1.0e-3, 1.0e-4, 1.0e-5, 1.0e-6, 1.0e-7, 1.0e-8, 1.0e-9,
        ] {
            let step = joint_step(&neighborhood, &models, anchor_weight)?;
            if step.iter().all(|v| v.abs() <= KERNEL_EPSILON * 0.001) {
                continue;
            }
            let mut accepted = false;
            for scale in [1.0, 0.5, 0.25, 0.125] {
                let mut proposed = models.clone();
                for (i, (_, model)) in proposed.iter_mut().enumerate() {
                    for (j, p) in model.parameters.iter_mut().enumerate() {
                        *p = step[i * 8 + j].mul_add(scale, *p);
                    }
                }
                if proposed.iter().all(|(_, model)| model.bounded()) {
                    models = proposed;
                    accepted = true;
                    break;
                }
            }
            if !accepted {
                break;
            }
        }
        for (id, model) in models {
            if model
                .parameters
                .iter()
                .any(|p| p.abs() > KERNEL_EPSILON * 0.01)
            {
                partition.carriers[id] = Some(model.carrier());
                max_additional_carrier_residual =
                    max_additional_carrier_residual.max(model.additional_residual_bound());
                refitted_regions += 1;
            }
        }
        if refitted_regions == 0 {
            return None;
        }
    }
    let (positions, allowances, maximum_displacement) =
        neighborhood.exact_positions(&partition.carriers, prior)?;
    Some(JointFitCandidate {
        partition,
        positions,
        allowances,
        maximum_displacement,
        refitted_regions,
        max_additional_carrier_residual,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn jointly_refits_an_overconstrained_pyramid_apex() -> Result<(), Box<dyn std::error::Error>> {
        use crate::mesh_recovery::{EvidenceOptions, RegionBoundaryIndex, analyze_mesh_evidence};
        use monstertruck_mesh::{Faces, PolygonMesh, StandardAttributes};
        let mesh = PolygonMesh::new(
            StandardAttributes {
                positions: vec![
                    Point3::new(-1.0, -1.0, 0.0),
                    Point3::new(1.0, -1.0, 0.0),
                    Point3::new(1.0, 1.0, 0.0),
                    Point3::new(-1.0, 1.0, 0.0),
                    Point3::new(0.0, 0.0, 1.0),
                ],
                ..Default::default()
            },
            Faces::from_iter(vec![
                vec![0, 1, 4],
                vec![1, 2, 4],
                vec![2, 3, 4],
                vec![3, 0, 4],
                vec![3, 2, 1, 0],
            ]),
        )?;
        let evidence = analyze_mesh_evidence(&mesh, EvidenceOptions::default())?;
        let index = RegionBoundaryIndex::new(&evidence)?;
        let mut partition = index.select_surface_partition(&[], &evidence, 100)?;
        let mut regions = evidence.source_planar_support[..4].to_vec();
        regions.sort_unstable();
        regions.dedup();
        assert_eq!(regions.len(), 4);
        let bad = regions[0];
        let Some(SurfaceCarrier::Plane { origin, normal }) = partition.carriers[bad] else {
            return Err("expected source plane".into());
        };
        partition.carriers[bad] = Some(SurfaceCarrier::Plane {
            origin: origin + normal * 0.003,
            normal,
        });
        let options = JointFitOptions {
            maximum_adjustment: 0.01,
            declared_uncertainty: 0.01,
            angular_resolution: 1.0e-8,
        };
        assert!(
            fit_neighborhood(
                &evidence,
                &partition,
                &evidence.positions,
                &regions,
                options,
                false
            )
            .is_none()
        );
        let result = fit_neighborhood(
            &evidence,
            &partition,
            &evidence.positions,
            &regions,
            options,
            true,
        )
        .ok_or("joint fit should restore compatible planes")?;
        assert!(result.refitted_regions > 0);
        assert!(result.maximum_displacement <= 0.01);
        let trims = index.assemble_shared_trims(
            &result.partition,
            &result.positions,
            super::super::SharedTrimOptions {
                absolute_tolerance: 1.0e-7,
                relative_tolerance: 1.0e-4,
                max_interval_edges: 128,
            },
        )?;
        let native = index.construct_native_trims(
            &result.partition,
            &trims,
            &result.positions,
            super::super::NativeTrimOptions {
                max_vertex_displacement: 1.0e-7,
                intersection_tolerance: 1.0e-8,
                angular_resolution: 1.0e-8,
            },
        )?;
        assert!(native.trims.iter().all(Result::is_ok));
        let mut shell = crate::Shell::new();
        for region in &native.regions {
            assert!(!region.requires_chart_cut);
            shell.extend(region.faces.as_ref().map_err(|e| *e)?.iter().cloned());
        }
        assert_eq!(
            shell.shell_condition(),
            monstertruck_topology::shell::ShellCondition::Closed
        );
        assert!(shell.is_geometric_consistent());
        Ok(())
    }

    #[test]
    fn continuous_carrier_change_bound_covers_local_profile_variation() {
        let pivot = Point3::new(2.0, 0.0, 0.0);
        let scale = 3.0;
        let models = [
            SurfaceCarrier::Plane {
                origin: pivot,
                normal: Vector3::unit_x(),
            },
            SurfaceCarrier::LinearRevolution {
                origin: Point3::new(0.0, 0.0, 0.0),
                axis: Vector3::unit_z(),
                radius_at_origin: 2.0,
                slope: 0.3,
            },
            SurfaceCarrier::CircularRevolution {
                axial_center: Point3::new(0.0, 0.0, 0.0),
                axis: Vector3::unit_z(),
                radial_center: 2.0,
                profile_radius: 1.0,
            },
            SurfaceCarrier::CircularRevolution {
                axial_center: Point3::new(0.0, 0.0, 0.0),
                axis: Vector3::unit_z(),
                radial_center: 0.0,
                profile_radius: 2.0,
            },
        ];
        for base in models {
            let model = CarrierModel {
                base,
                pivot,
                scale,
                allowance: 0.01,
                parameters: [
                    0.001, -0.0004, 0.0005, 0.001, 0.0005, -0.0002, 0.0008, 0.0009,
                ],
            };
            assert!(model.bounded());
            for x in -2..=2 {
                for y in -2..=2 {
                    for z in -2..=2 {
                        let point =
                            pivot + Vector3::new(f64::from(x), f64::from(y), f64::from(z)) * 0.7;
                        let (Some((_, a)), Some((_, b))) = (
                            base.signed_constraint(point),
                            model.carrier().signed_constraint(point),
                        ) else {
                            continue;
                        };
                        assert!((a - b).abs() <= model.additional_residual_bound());
                    }
                }
            }
        }
    }

    #[test]
    fn damped_system_solves_coupled_constraints() {
        let result = solve(vec![vec![4.0, 2.0], vec![2.0, 3.0]], vec![6.0, 5.0]);
        assert!(
            result.is_some_and(|x| (x[0] - 1.0).abs() < 1.0e-12 && (x[1] - 1.0).abs() < 1.0e-12)
        );
        assert!(solve(vec![vec![0.0]], vec![0.0]).is_none());
    }
}

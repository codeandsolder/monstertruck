//! Recover placed axis observations from source-connected circular samples.
//! A circle through three points is only an interpolation. Repeated,
//! independently resolved circles may supply an axis even when local
//! surface-normal tracks do not expose the global rotational direction.

use super::super::sample_fit::{CircleCarrierFit, fit_circle_carrier};
use super::{MeshEvidence, source_edge_incidence};
use crate::{InnerProduct, Point3, Vector3};
use rustc_hash::{FxHashMap as HashMap, FxHashSet as HashSet};
use std::collections::VecDeque;

/// Accuracy budgets for original circular source samples and axis grouping.
#[derive(Clone, Copy, Debug)]
pub struct CircleAxisOptions {
    /// Absolute source sample fit error, in mesh units.
    pub vertex_tolerance: f64,
    /// Relative sample fit error, capped by observed source chord span.
    pub relative_tolerance: f64,
    /// Maximum separation of axis reference points, in mesh units.
    pub axis_position_tolerance: f64,
    /// Maximum undirected axis angular disagreement.
    pub axis_angle_tolerance_radians: f64,
}

impl CircleAxisOptions {
    fn valid(self) -> bool {
        [
            self.vertex_tolerance,
            self.relative_tolerance,
            self.axis_position_tolerance,
        ]
        .into_iter()
        .all(|v| v.is_finite() && v > 0.0)
            && (0.0..std::f64::consts::FRAC_PI_2).contains(&self.axis_angle_tolerance_radians)
    }
}

/// Independently resolved circular sample component.
#[derive(Clone, Debug)]
pub struct SourceCircleObservation {
    /// Original source-connected canonical vertex IDs.
    pub vertices: Vec<usize>,
    /// Placed circle center.
    pub center: Point3,
    /// Canonical unit normal.
    pub normal: Vector3,
    /// Fitted circle radius.
    pub radius: f64,
    /// Maximum absolute distance to the fitted plane.
    pub max_plane_residual: f64,
    /// Maximum projected radial residual.
    pub max_radial_residual: f64,
    /// Observed chord span times normalized noncollinearity.
    pub resolved_curvature_span: f64,
}

/// Repeated circles consistent with one placed axis.
#[derive(Clone, Debug)]
pub struct CircleAxisEvidence {
    /// Reference point on the seed axis, closest to a common mesh reference.
    pub origin: Point3,
    /// Fixed source circle normal, never a transitive average.
    pub axis: Vector3,
    /// Distinct axial/radial station observations on this axis.
    /// Separate arcs on a common ring contribute one station, not many.
    pub stations: Vec<CircleAxisStation>,
    /// Number of all contributing original circle components.
    pub circle_count: usize,
    /// Maximum observed angular discrepancy against the seed.
    pub max_axis_angle_radians: f64,
    /// Maximum observed reference-point separation against the seed.
    pub max_axis_position_error: f64,
}

/// Circle station with original connected sample provenance.
#[derive(Clone, Debug)]
pub struct CircleAxisStation {
    /// Axial coordinate in the seed reference frame.
    pub axial: f64,
    /// Observed circle radius.
    pub radius: f64,
    /// Original independently fitted circle component IDs.
    pub circle_indices: Vec<usize>,
}

/// Source-circle discovery keeps all observations for competing hypotheses.
#[derive(Debug)]
pub struct CircleAxisDiscovery {
    /// Distinct connected source triples considered.
    pub source_triples: usize,
    /// Original source-connected circles with at least four samples.
    pub circles: Vec<SourceCircleObservation>,
    /// Placed axes supported by multiple distinct circle stations.
    pub axes: Vec<CircleAxisEvidence>,
}

struct CircleSearch<'a> {
    evidence: &'a MeshEvidence,
    options: CircleAxisOptions,
    neighbors: Vec<Vec<usize>>,
    circles: Vec<SourceCircleObservation>,
    circles_by_vertex: Vec<Vec<usize>>,
    seeds: HashSet<[usize; 3]>,
}

impl<'a> CircleSearch<'a> {
    fn new(evidence: &'a MeshEvidence, options: CircleAxisOptions) -> Option<Self> {
        let mut neighbors = vec![Vec::new(); evidence.positions.len()];
        // A circular source path can cross a locally coplanar tangent
        // join. Feature-edge suppression is useful for normal tracks, but
        // it must not erase the source geometry used to observe circles.
        let incidence = source_edge_incidence(
            (0..evidence.source_planar_support.len())
                .map(|id| evidence.source_polygon_vertices(id))
                .collect::<Option<Vec<_>>>()?,
        );
        for ([a, b], uses) in incidence {
            if !matches!(uses.as_slice(), [_])
                && !matches!(uses.as_slice(),[left,right] if left.1!=right.1)
            {
                continue;
            }
            neighbors.get_mut(a)?.push(b);
            neighbors.get_mut(b)?.push(a);
        }
        for row in &mut neighbors {
            row.sort_unstable();
            row.dedup();
        }
        Some(Self {
            evidence,
            options,
            neighbors,
            circles: Vec::new(),
            circles_by_vertex: vec![Vec::new(); evidence.positions.len()],
            seeds: HashSet::default(),
        })
    }

    fn discover(mut self) -> CircleAxisDiscovery {
        for center in 0..self.neighbors.len() {
            for i in 0..self.neighbors[center].len() {
                let first = self.neighbors[center][i];
                for j in i + 1..self.neighbors[center].len() {
                    let last = self.neighbors[center][j];
                    let mut seed = [first, center, last];
                    seed.sort_unstable();
                    if self.seeds.insert(seed) && !self.already_resolved(seed) {
                        self.measure(seed);
                    }
                }
            }
        }
        let axes = group_axes(self.evidence, &self.circles, self.options);
        CircleAxisDiscovery {
            source_triples: self.seeds.len(),
            circles: self.circles,
            axes,
        }
    }

    fn already_resolved(&self, seed: [usize; 3]) -> bool {
        // Reuse a previously fitted exact circle only when all three
        // original samples already belong to it. Approximate circles do
        // not suppress alternative seed fits near their tolerance boundary.
        self.circles_by_vertex[seed[0]].iter().any(|&id| {
            let circle = &self.circles[id];
            let numerical =
                128.0 * f64::EPSILON * circle.radius.max(circle.resolved_curvature_span);
            circle.max_plane_residual.max(circle.max_radial_residual) <= numerical
                && seed[1..]
                    .iter()
                    .all(|vertex| circle.vertices.binary_search(vertex).is_ok())
        })
    }

    fn measure(&mut self, seed: [usize; 3]) {
        let Some(seed_fit) = fit_circle_carrier(&self.evidence.positions, &seed) else {
            return;
        };
        if !valid_fit(seed_fit, self.options) {
            return;
        }
        let mut members = HashSet::from_iter(seed);
        let mut tested = members.clone();
        let mut queue = VecDeque::from(seed);
        // Test each original neighbor once. A rejected point cannot connect
        // otherwise unrelated parts of the sampled circle.
        while let Some(id) = queue.pop_front() {
            for &other in &self.neighbors[id] {
                if tested.insert(other)
                    && sample_matches(
                        self.evidence.positions[other],
                        seed_fit,
                        self.options.vertex_tolerance,
                    )
                {
                    members.insert(other);
                    queue.push_back(other);
                }
            }
        }
        if members.len() < 4 {
            return;
        }
        let mut vertices = members.into_iter().collect::<Vec<_>>();
        vertices.sort_unstable();
        if self.circles_by_vertex[vertices[0]]
            .iter()
            .any(|&id| self.circles[id].vertices == vertices)
        {
            return;
        }
        let Some(fit) = fit_circle_carrier(&self.evidence.positions, &vertices) else {
            return;
        };
        if !valid_fit(fit, self.options) {
            return;
        }
        let id = self.circles.len();
        for &vertex in &vertices {
            self.circles_by_vertex[vertex].push(id);
        }
        self.circles.push(SourceCircleObservation {
            vertices,
            center: fit.center,
            normal: fit.normal,
            radius: fit.radius,
            max_plane_residual: fit.max_plane_residual,
            max_radial_residual: fit.max_radial_residual,
            resolved_curvature_span: fit.source_span * fit.plane_span_strength,
        });
    }
}

fn valid_fit(fit: CircleCarrierFit, options: CircleAxisOptions) -> bool {
    let tolerance = options
        .vertex_tolerance
        .min(options.relative_tolerance * fit.radius.min(fit.source_span));
    tolerance > 0.0
        && fit.source_span * fit.plane_span_strength > options.vertex_tolerance
        && fit.max_plane_residual.max(fit.max_radial_residual) <= tolerance
}

fn sample_matches(point: Point3, fit: CircleCarrierFit, tolerance: f64) -> bool {
    let delta = point - fit.center;
    let axial = delta.dot(&fit.normal);
    axial.abs() <= tolerance
        && ((delta - fit.normal * axial).norm() - fit.radius).abs() <= tolerance
}

fn reference_origin(center: Point3, axis: Vector3, reference: Point3) -> Point3 {
    center + axis * (reference - center).dot(&axis)
}

fn position_cell(origin: Point3, reference: Point3, width: f64) -> [f64; 3] {
    let delta = origin - reference;
    [delta.x, delta.y, delta.z].map(|value| (value / width).floor() + 0.0)
}

fn nearby_groups(bins: &HashMap<[u64; 3], Vec<usize>>, cell: [f64; 3]) -> Vec<usize> {
    let mut result = Vec::new();
    for x in [-1.0, 0.0, 1.0] {
        for y in [-1.0, 0.0, 1.0] {
            for z in [-1.0, 0.0, 1.0] {
                let key = [
                    (cell[0] + x + 0.0).to_bits(),
                    (cell[1] + y + 0.0).to_bits(),
                    (cell[2] + z + 0.0).to_bits(),
                ];
                if let Some(ids) = bins.get(&key) {
                    result.extend(ids);
                }
            }
        }
    }
    result
}

fn group_axes(
    evidence: &MeshEvidence,
    circles: &[SourceCircleObservation],
    options: CircleAxisOptions,
) -> Vec<CircleAxisEvidence> {
    let Some(&reference) = evidence.positions.first() else {
        return Vec::new();
    };
    let mut order = (0..circles.len()).collect::<Vec<_>>();
    order.sort_unstable_by(|&a, &b| {
        circles[b]
            .vertices
            .len()
            .cmp(&circles[a].vertices.len())
            .then_with(|| {
                circles[b]
                    .resolved_curvature_span
                    .total_cmp(&circles[a].resolved_curvature_span)
            })
            .then(a.cmp(&b))
    });
    let mut groups: Vec<CircleAxisEvidence> = Vec::new();
    let mut bins = HashMap::<[u64; 3], Vec<usize>>::default();
    for id in order {
        let circle = &circles[id];
        let origin = reference_origin(circle.center, circle.normal, reference);
        let cell = position_cell(origin, reference, options.axis_position_tolerance);
        let nearest = nearby_groups(&bins, cell)
            .into_iter()
            .filter_map(|index| {
                let group = &groups[index];
                let position = (origin - group.origin).norm();
                let angle = group
                    .axis
                    .cross(&circle.normal)
                    .norm()
                    .atan2(group.axis.dot(&circle.normal).abs());
                (position <= options.axis_position_tolerance
                    && angle <= options.axis_angle_tolerance_radians)
                    .then_some((index, position, angle))
            })
            .min_by(|a, b| {
                a.1.total_cmp(&b.1)
                    .then_with(|| a.2.total_cmp(&b.2))
                    .then(a.0.cmp(&b.0))
            });
        let index = if let Some((index, position, angle)) = nearest {
            groups[index].max_axis_position_error =
                groups[index].max_axis_position_error.max(position);
            groups[index].max_axis_angle_radians = groups[index].max_axis_angle_radians.max(angle);
            index
        } else {
            let index = groups.len();
            groups.push(CircleAxisEvidence {
                origin,
                axis: circle.normal,
                stations: Vec::new(),
                circle_count: 0,
                max_axis_angle_radians: 0.0,
                max_axis_position_error: 0.0,
            });
            bins.entry(cell.map(f64::to_bits)).or_default().push(index);
            index
        };
        let group = &mut groups[index];
        let axial = (circle.center - group.origin).dot(&group.axis);
        if let Some(station) = group.stations.iter_mut().find(|station| {
            (station.axial - axial).abs() <= options.axis_position_tolerance
                && (station.radius - circle.radius).abs() <= options.vertex_tolerance
        }) {
            station.circle_indices.push(id);
        } else {
            group.stations.push(CircleAxisStation {
                axial,
                radius: circle.radius,
                circle_indices: vec![id],
            });
        }
        group.circle_count += 1;
    }
    groups.retain(|group| group.stations.len() >= 2);
    groups.sort_unstable_by(|a, b| {
        b.stations
            .len()
            .cmp(&a.stations.len())
            .then(b.circle_count.cmp(&a.circle_count))
            .then_with(|| a.origin.x.total_cmp(&b.origin.x))
            .then_with(|| a.origin.y.total_cmp(&b.origin.y))
            .then_with(|| a.origin.z.total_cmp(&b.origin.z))
    });
    groups
}

/// Search regular source edges for independently sampled circles and placed axes.
///
/// An axis group is only a hypothesis. A planar or faceted feature can also
/// produce four co-circular vertices. Subsequent source-connected region
/// fitting must provide independent evidence of the actual surface.
#[must_use = "source circle observations do not identify a CAD surface"]
pub fn discover_circle_axes(
    evidence: &MeshEvidence,
    options: CircleAxisOptions,
) -> Option<CircleAxisDiscovery> {
    if !options.valid() {
        return None;
    }
    Some(CircleSearch::new(evidence, options)?.discover())
}

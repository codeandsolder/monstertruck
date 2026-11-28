//! Local-to-global normal-space rotation tracks over observed mesh supports.
//!
//! Three normals define a plane in normal space. A fourth tests whether
//! a fixed rotation axis still explains that plane. This is an evidence
//! stage, not a CAD surface or material certification.

use super::super::MeshRecoveryError;
use super::{MeshEvidence, canonical_unit};
use crate::{InnerProduct, Vector3};
use rustc_hash::{FxHashMap as HashMap, FxHashSet as HashSet};

/// Numerical hypothesis bounds. All original supports and interfaces remain
/// available in the evidence graph regardless of this track selection.
#[derive(Clone, Copy, Debug)]
pub struct CurvatureTrackOptions {
    /// Minimum rotation angle across an observed source support interface.
    pub min_normal_turn_radians: f64,
    /// Minimum sine between successive normal-difference vectors.
    pub min_normal_axis_strength: f64,
    /// Angular tolerance between two consecutive inferred axes.
    pub max_adjacent_axis_error_radians: f64,
    /// Angular tolerance against the first, fixed track axis.
    pub max_seed_axis_error_radians: f64,
    /// Maximum range of support-normal dot products with the fixed axis.
    pub max_normal_axis_projection_spread: f64,
    /// Minimum distinct supports needed to report a candidate track.
    pub min_distinct_supports: usize,
}

impl Default for CurvatureTrackOptions {
    fn default() -> Self {
        Self {
            min_normal_turn_radians: 1.0e-6,
            min_normal_axis_strength: 1.0e-5,
            max_adjacent_axis_error_radians: 0.15,
            max_seed_axis_error_radians: 0.15,
            max_normal_axis_projection_spread: 0.06,
            min_distinct_supports: 4,
        }
    }
}

/// Maximal unambiguous connected path with consistent normal-space axis.
#[derive(Clone, Debug)]
pub struct CurvatureTrack {
    /// Ordered source support IDs, repeating the first when closed.
    pub support_indices: Vec<usize>,
    /// One original support-interface ID between consecutive support IDs.
    pub interface_indices: Vec<usize>,
    /// Unique source support IDs. Non-simple source contact paths can revisit
    /// a region without returning to their starting directed contact.
    pub distinct_support_count: usize,
    /// The first triplet's canonical unit normal-space axis.
    pub seed_axis: Vector3,
    /// Maximum angular divergence of local axes from the fixed seed.
    pub max_seed_axis_error_radians: f64,
    /// Range of observed source-normal axial projections.
    pub normal_axis_projection_spread: f64,
    /// Mean axial component of the source support normals.
    pub mean_normal_axis_projection: f64,
    /// Sum of successive unsigned source-normal turns.
    pub accumulated_normal_turn_radians: f64,
    /// Maximum mismatch between consecutive normal turn magnitudes.
    pub max_adjacent_turn_mismatch_radians: f64,
    /// True when the support chain closes into a directed cycle.
    pub closed: bool,
}

#[derive(Clone, Copy)]
struct RotationStep {
    supports: [usize; 3],
    interfaces: [usize; 2],
    axis: Vector3,
}

impl RotationStep {
    const fn incoming(self) -> [usize; 3] {
        [self.interfaces[0], self.supports[0], self.supports[1]]
    }

    const fn outgoing(self) -> [usize; 3] {
        [self.interfaces[1], self.supports[1], self.supports[2]]
    }
}

fn axis_angle(a: Vector3, b: Vector3) -> f64 { a.cross(&b).norm().atan2(a.dot(&b).abs()) }

fn normal_angle(a: Vector3, b: Vector3) -> f64 { a.cross(&b).norm().atan2(a.dot(&b)) }

fn valid_options(options: CurvatureTrackOptions) -> bool {
    options.min_normal_turn_radians.is_finite()
        && (0.0..std::f64::consts::PI).contains(&options.min_normal_turn_radians)
        && options.min_normal_axis_strength.is_finite()
        && (0.0..=1.0).contains(&options.min_normal_axis_strength)
        && options.max_adjacent_axis_error_radians.is_finite()
        && (0.0..std::f64::consts::FRAC_PI_2).contains(&options.max_adjacent_axis_error_radians)
        && options.max_seed_axis_error_radians.is_finite()
        && (0.0..std::f64::consts::FRAC_PI_2).contains(&options.max_seed_axis_error_radians)
        && options.max_normal_axis_projection_spread.is_finite()
        && (0.0..=2.0).contains(&options.max_normal_axis_projection_spread)
        && options.min_distinct_supports >= 3
}

/// The local axis is inferred from differences between three source normals.
/// No manufactured primitive type or source STEP face ID is assumed here.
fn local_steps(evidence: &MeshEvidence, options: CurvatureTrackOptions) -> Vec<RotationStep> {
    let mut steps = Vec::new();
    for (center, contacts) in evidence.support_interface_indices.iter().enumerate() {
        for (i, &first) in contacts.iter().enumerate() {
            for &second in &contacts[i + 1..] {
                let Some(rotation) = evidence.support_normal_rotation(center, first, second) else {
                    continue;
                };
                let Some(axis) = rotation.normal_space_axis else {
                    continue;
                };
                if rotation.normal_space_axis_strength < options.min_normal_axis_strength
                    || rotation
                        .turn_radians
                        .iter()
                        .any(|&angle| angle < options.min_normal_turn_radians)
                {
                    continue;
                }
                let [left, right] = rotation.neighbor_supports;
                let [a, b, c] = [
                    evidence.planar_supports[left].normal,
                    evidence.planar_supports[center].normal,
                    evidence.planar_supports[right].normal,
                ];
                let sign_first = axis.dot(&a.cross(&b));
                let sign_second = axis.dot(&b.cross(&c));
                if sign_first * sign_second <= 0.0 {
                    continue;
                }
                let (supports, interfaces) = if sign_first > 0.0 {
                    ([left, center, right], [first, second])
                } else {
                    ([right, center, left], [second, first])
                };
                steps.push(RotationStep {
                    supports,
                    interfaces,
                    axis: canonical_unit(axis),
                });
            }
        }
    }
    steps
}

/// A link is selected only if each end has exactly one compatible partner.
/// Alternative branches remain as separate observations, never guessed through.
fn unambiguous_links(steps: &[RotationStep], options: CurvatureTrackOptions) -> Vec<Option<usize>> {
    let mut incoming = HashMap::<[usize; 3], Vec<usize>>::default();
    let mut outgoing = HashMap::<[usize; 3], Vec<usize>>::default();
    for (id, step) in steps.iter().enumerate() {
        incoming.entry(step.incoming()).or_default().push(id);
        outgoing.entry(step.outgoing()).or_default().push(id);
    }
    let mut result = vec![None; steps.len()];
    for (id, step) in steps.iter().enumerate() {
        let mut candidates = incoming
            .get(&step.outgoing())
            .into_iter()
            .flatten()
            .copied()
            .filter(|&next| {
                next != id
                    && axis_angle(step.axis, steps[next].axis)
                        <= options.max_adjacent_axis_error_radians
            });
        let Some(next) = candidates.next() else {
            continue;
        };
        if candidates.next().is_some() {
            continue;
        }
        let mut predecessors = outgoing
            .get(&steps[next].incoming())
            .into_iter()
            .flatten()
            .copied()
            .filter(|&previous| {
                previous != next
                    && axis_angle(steps[next].axis, steps[previous].axis)
                        <= options.max_adjacent_axis_error_radians
            });
        if predecessors.next() == Some(id) && predecessors.next().is_none() {
            result[id] = Some(next);
        }
    }
    result
}

fn track_from_path(
    evidence: &MeshEvidence,
    steps: &[RotationStep],
    indices: &[usize],
    closed: bool,
) -> Option<CurvatureTrack> {
    let first = steps.get(*indices.first()?)?;
    let mut supports = first.supports.to_vec();
    let mut interfaces = first.interfaces.to_vec();
    let axis = first.axis;
    let mut max_axis_error = 0.0_f64;
    for &index in &indices[1..] {
        let step = steps.get(index)?;
        max_axis_error = max_axis_error.max(axis_angle(axis, step.axis));
        supports.push(step.supports[2]);
        interfaces.push(step.interfaces[1]);
    }
    if closed {
        // The final step restates both initial vertices of a directed arc.
        supports.pop();
        interfaces.pop();
    }
    let projections = supports
        .iter()
        .map(|&support| evidence.planar_supports[support].normal.dot(&axis))
        .collect::<Vec<_>>();
    let min_projection = projections.iter().copied().fold(f64::INFINITY, f64::min);
    let max_projection = projections
        .iter()
        .copied()
        .fold(f64::NEG_INFINITY, f64::max);
    let turns = supports
        .windows(2)
        .map(|pair| {
            normal_angle(
                evidence.planar_supports[pair[0]].normal,
                evidence.planar_supports[pair[1]].normal,
            )
        })
        .collect::<Vec<_>>();
    let max_adjacent_turn_mismatch_radians = turns
        .windows(2)
        .map(|pair| (pair[0] - pair[1]).abs())
        .fold(0.0_f64, f64::max);
    let (projection_sum, projection_count) = projections
        .iter()
        .copied()
        .fold((0.0_f64, 0.0_f64), |(sum, count), value| {
            (sum + value, count + 1.0)
        });
    let distinct_support_count = supports.iter().copied().collect::<HashSet<_>>().len();
    Some(CurvatureTrack {
        support_indices: supports,
        interface_indices: interfaces,
        distinct_support_count,
        seed_axis: axis,
        max_seed_axis_error_radians: max_axis_error,
        normal_axis_projection_spread: max_projection - min_projection,
        mean_normal_axis_projection: projection_sum / projection_count,
        accumulated_normal_turn_radians: turns.iter().sum(),
        max_adjacent_turn_mismatch_radians,
        closed,
    })
}

/// Construct unambiguous normal-space-axis-coherent source support tracks.
///
/// Triplets of normals each define an axis by construction. A fourth support
/// is needed to test that axis's consistency, so the default output minimum
/// is four distinct supports. All original edges and ambiguous observations
/// remain available for alternative surface hypotheses.
///
/// # Errors
/// Returns `InvalidEvidenceOptions` for invalid hypothesis bounds.
pub fn trace_curvature_tracks(
    evidence: &MeshEvidence,
    options: CurvatureTrackOptions,
) -> Result<Vec<CurvatureTrack>, MeshRecoveryError> {
    if !valid_options(options) {
        return Err(MeshRecoveryError::InvalidEvidenceOptions);
    }
    let steps = local_steps(evidence, options);
    let next = unambiguous_links(&steps, options);
    let mut previous = vec![None; steps.len()];
    for (id, &successor) in next.iter().enumerate() {
        if let Some(successor) = successor {
            previous[successor] = Some(id);
        }
    }
    let mut visited = vec![false; steps.len()];
    let mut tracks = Vec::new();
    for origin in 0..steps.len() {
        if visited[origin] {
            continue;
        }
        let mut first = origin;
        for _ in 0..steps.len() {
            let Some(prev) = previous[first] else {
                break;
            };
            if prev == origin || visited[prev] {
                break;
            }
            first = prev;
        }
        let seed = steps[first].axis;
        let projections = steps[first]
            .supports
            .map(|support| evidence.planar_supports[support].normal.dot(&seed));
        let mut min_projection = projections.iter().copied().fold(f64::INFINITY, f64::min);
        let mut max_projection = projections
            .iter()
            .copied()
            .fold(f64::NEG_INFINITY, f64::max);
        let mut indices = Vec::new();
        let mut current = first;
        let mut closed = false;
        loop {
            if visited[current] {
                closed = current == first;
                break;
            }
            visited[current] = true;
            indices.push(current);
            let Some(successor) = next[current] else {
                break;
            };
            if visited[successor] {
                closed = successor == first;
                break;
            }
            if axis_angle(seed, steps[successor].axis) > options.max_seed_axis_error_radians {
                break;
            }
            let projection = evidence.planar_supports[steps[successor].supports[2]]
                .normal
                .dot(&seed);
            let new_min = min_projection.min(projection);
            let new_max = max_projection.max(projection);
            if new_max - new_min > options.max_normal_axis_projection_spread {
                break;
            }
            min_projection = new_min;
            max_projection = new_max;
            current = successor;
        }
        if let Some(track) = track_from_path(evidence, &steps, &indices, closed)
            && track.distinct_support_count >= options.min_distinct_supports
        {
            tracks.push(track);
        }
    }
    tracks.sort_unstable_by(|a, b| {
        b.distinct_support_count
            .cmp(&a.distinct_support_count)
            .then_with(|| a.support_indices.cmp(&b.support_indices))
    });
    Ok(tracks)
}

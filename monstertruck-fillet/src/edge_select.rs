use std::collections::{HashMap, VecDeque};
use std::env;

use crate::HashSet;

use monstertruck_core::tolerance::KERNEL_EPSILON;
use monstertruck_geometry::prelude::{NurbsSurface, Point3, Vector4};
use monstertruck_topology::shell::ShellCondition;
use monstertruck_traits::{
    BoundedCurve, ParametricCurve, ParametricSurface, SearchNearestParameter,
};
use smallvec::SmallVec;

use super::convert::{FilletableCurve, FilletableSurface, convert_shell_in, convert_shell_out};
use super::error::FilletError;
use super::ops;
use super::params::{FilletOptions, FilletProfile, RadiusSpec};
use super::topology::wire_edge_starts_and_spans;
use super::types::{
    Curve, Edge, EdgeId, Shell, VertexId, Wire, approximate_edge_length, squared_distance,
};

type Result<T> = std::result::Result<T, FilletError>;

/// Face list for an edge. Manifold edges always have exactly 2 faces,
/// so `SmallVec<[usize; 2]>` avoids heap allocation.
type FaceList = SmallVec<[usize; 2]>;

fn is_manifold_edge(edge_face_map: &HashMap<EdgeId, FaceList>, edge_id: EdgeId) -> bool {
    edge_face_map
        .get(&edge_id)
        .is_some_and(|faces| faces.len() == 2)
}

/// Builds a map from [`EdgeId`] to the face indices that contain it.
fn build_edge_face_map(shell: &Shell) -> HashMap<EdgeId, FaceList> {
    let mut map: HashMap<EdgeId, FaceList> = HashMap::new();
    shell.iter().enumerate().for_each(|(face_idx, face)| {
        face.edge_iter().for_each(|edge| {
            map.entry(edge.id())
                .and_modify(|v| {
                    if !v.contains(&face_idx) {
                        v.push(face_idx);
                    }
                })
                .or_insert_with(|| smallvec::smallvec![face_idx]);
        });
    });
    map
}

/// Project a surviving edge onto the original curve. The projected samples
/// must follow one monotone interval; a loop, reversal, or deviation must not
/// accidentally classify an unrelated edge as one of its descendants.
fn curve_interval_on_original(original: &Curve, candidate: &Curve) -> Option<(f64, f64)> {
    let (candidate_start, candidate_end) = candidate.range_tuple();
    let (original_start, original_end) = original.range_tuple();
    let interval_epsilon = (original_end - original_start).abs() * KERNEL_EPSILON;
    let mut previous = None;
    let mut projected = Vec::with_capacity(5);
    for fraction in [0.0, 0.25, 0.5, 0.75, 1.0] {
        let t = (candidate_end - candidate_start).mul_add(fraction, candidate_start);
        let point = candidate.evaluate(t);
        let u = original.search_nearest_parameter(point, None, 100)?;
        if !u.is_finite()
            || squared_distance(original.evaluate(u), point) >= KERNEL_EPSILON * KERNEL_EPSILON
        {
            return None;
        }
        if let Some(previous) = previous {
            projected.push(u - previous);
        }
        previous = Some(u);
    }
    let direction = projected.iter().sum::<f64>();
    if direction.abs() <= interval_epsilon
        || projected
            .iter()
            .any(|&step| step * direction < -interval_epsilon * direction.abs())
    {
        return None;
    }
    let first =
        original.search_nearest_parameter(candidate.evaluate(candidate_start), None, 100)?;
    let last = previous?;
    Some((first.min(last), first.max(last)))
}

/// Return every non-overlapping manifold segment of the original curve in
/// its original parameter order. Overlapping matches are ambiguous topology,
/// not a reason to pick whichever surviving edge happens to be longest.
fn rematch_selected_edge_ids(
    shell: &Shell,
    original_curve: &Curve,
    used_ids: &HashSet<EdgeId>,
    edge_face_map: &HashMap<EdgeId, FaceList>,
) -> Option<Vec<EdgeId>> {
    let mut seen_ids = HashSet::new();
    let mut candidates = shell
        .edge_iter()
        .filter(|edge| {
            seen_ids.insert(edge.id())
                && !used_ids.contains(&edge.id())
                && is_manifold_edge(edge_face_map, edge.id())
        })
        .filter_map(|edge| {
            curve_interval_on_original(original_curve, &edge.curve())
                .map(|interval| (edge.id(), interval))
        })
        .collect::<Vec<_>>();
    candidates.sort_by(|a, b| a.1.0.total_cmp(&b.1.0));
    let (start, end) = original_curve.range_tuple();
    let interval_epsilon = (end - start).abs() * KERNEL_EPSILON;
    if candidates.is_empty()
        || candidates
            .windows(2)
            .any(|pair| pair[1].1.0 < pair[0].1.1 - interval_epsilon)
    {
        return None;
    }
    Some(candidates.into_iter().map(|(id, _)| id).collect())
}

/// Exact identity takes precedence over geometry. A replaced edge may have
/// multiple non-overlapping successors; returning only the longest would omit
/// valid portions of the user's selection.
fn resolve_active_edge_ids(
    shell: &Shell,
    original_id: EdgeId,
    original_curve: &Curve,
    used_ids: &HashSet<EdgeId>,
    edge_face_map: &HashMap<EdgeId, FaceList>,
) -> Result<Vec<EdgeId>> {
    if is_manifold_edge(edge_face_map, original_id) && !used_ids.contains(&original_id) {
        return Ok(vec![original_id]);
    }
    rematch_selected_edge_ids(shell, original_curve, used_ids, edge_face_map)
        .ok_or(FilletError::UnresolvableSelectedEdge)
}

/// Visit every surviving fragment, including descendants created by a
/// preceding fillet in this same chain. Each pending fragment stores its own
/// curve so later rematching never falls back to an unrelated, longer portion
/// of the parent edge. The caller owns the surrounding shell transaction.
fn visit_surviving_fragments(
    shell: &mut Shell,
    segments: VecDeque<(EdgeId, EdgeId, Curve)>,
    mut apply: impl FnMut(&mut Shell, EdgeId, EdgeId) -> Result<()>,
) -> Result<()> {
    let mut pending = segments;
    let mut used = HashSet::<EdgeId>::new();
    // EdgeId is pointer-derived. Retain every pending and processed edge's
    // Arc-backed handle until the visit ends, preventing a freed ID from being
    // recycled for a new segment (the ABA problem).
    let mut pinned_edges: Vec<Edge> = pending
        .iter()
        .filter_map(|(_, id, _)| shell.edge_iter().find(|edge| edge.id() == *id))
        .collect();
    while let Some((original_id, active_id, hint_curve)) = pending.pop_front() {
        let edge_face_map = build_edge_face_map(shell);
        let successors =
            resolve_active_edge_ids(shell, active_id, &hint_curve, &used, &edge_face_map)?;
        let matched_edges: Vec<Edge> = successors
            .into_iter()
            .map(|id| {
                shell
                    .edge_iter()
                    .find(|edge| edge.id() == id)
                    .ok_or(FilletError::EdgeNotFound)
            })
            .collect::<Result<_>>()?;
        let first_id = matched_edges[0].id();
        // Descendants may be split again by the next topology edit. Queue
        // them with their actual geometry before touching the shell.
        for edge in matched_edges.iter().skip(1).rev() {
            pending.push_front((original_id, edge.id(), edge.oriented_curve()));
        }
        pinned_edges.extend(matched_edges);
        debug_assert!(pinned_edges.iter().any(|edge| edge.id() == first_id));
        used.insert(first_id);
        apply(shell, original_id, first_id)?;
    }
    Ok(())
}

fn sampled_curve_surface_max_distance(
    curve: &Curve,
    surface: &NurbsSurface<Vector4>,
    sample_count: usize,
) -> f64 {
    let (t0, t1) = curve.range_tuple();
    (0..=sample_count)
        .map(|i| t0 + (t1 - t0) * crate::count_as_f64(i) / crate::count_as_f64(sample_count))
        .map(|t| curve.evaluate(t))
        .filter_map(|point| {
            surface
                .search_nearest_parameter(point, None, 50)
                .map(|(u, v)| squared_distance(surface.evaluate(u, v), point).sqrt())
        })
        .fold(0.0, f64::max)
}

/// Scans all face boundaries and extracts contiguous runs of selected edges.
///
/// Returns `(face_idx, ordered edge IDs)` for each maximal run.
fn collect_boundary_runs(shell: &Shell, selected: &HashSet<EdgeId>) -> Vec<(usize, Vec<EdgeId>)> {
    let mut candidate_runs = Vec::new();

    for (face_idx, face) in shell.iter().enumerate() {
        for boundary in face.boundary_iters() {
            let boundary_edges: Vec<Edge> = boundary.collect();
            let n = boundary_edges.len();
            if n == 0 {
                continue;
            }

            let is_selected: Vec<bool> = boundary_edges
                .iter()
                .map(|e| selected.contains(&e.id()))
                .collect();

            let boundary_closed =
                n >= 2 && boundary_edges[0].front() == boundary_edges[n - 1].back();

            if boundary_closed {
                // Find a non-selected edge to start from so runs don't
                // get split at the wrap-around seam.
                if let Some(start) = is_selected.iter().position(|&s| !s) {
                    let mut run: Vec<EdgeId> = Vec::new();
                    for offset in 1..=n {
                        let idx = (start + offset) % n;
                        if is_selected[idx] {
                            run.push(boundary_edges[idx].id());
                        } else if !run.is_empty() {
                            candidate_runs.push((face_idx, std::mem::take(&mut run)));
                        }
                    }
                    if !run.is_empty() {
                        candidate_runs.push((face_idx, run));
                    }
                } else {
                    // ALL edges selected -> one closed run covering entire boundary.
                    let ids: Vec<EdgeId> = boundary_edges
                        .iter()
                        .map(monstertruck_topology::Edge::id)
                        .collect();
                    candidate_runs.push((face_idx, ids));
                }
            } else {
                // Open boundary: linear scan, split on non-selected gaps.
                let mut run: Vec<EdgeId> = Vec::new();
                for (i, &sel) in is_selected.iter().enumerate() {
                    if sel {
                        run.push(boundary_edges[i].id());
                    } else if !run.is_empty() {
                        candidate_runs.push((face_idx, std::mem::take(&mut run)));
                    }
                }
                if !run.is_empty() {
                    candidate_runs.push((face_idx, run));
                }
            }
        }
    }

    candidate_runs
}

/// An edge chain: a contiguous run of selected edges on a single face boundary.
///
/// `shared_face_idx` is the face whose boundary contains these edges.
/// The other face index is resolved at processing time (it may shift as
/// earlier chains mutate the shell).
#[derive(Debug)]
struct Chain {
    edge_ids: Vec<EdgeId>,
    shared_face_idx: usize,
}

/// Groups requested edges into chains by finding contiguous runs of selected
/// edges along each face boundary.
///
/// For each face boundary, maximal contiguous runs of selected edges are
/// extracted. Each selected edge appears in exactly two candidate runs (one per
/// adjacent face). The edge is assigned to the longer run (tiebreak: lower face
/// index). Final chains are sorted longest-first for processing stability.
fn group_edges_into_chains(
    shell: &Shell,
    edge_ids: &[EdgeId],
    edge_face_map: &HashMap<EdgeId, FaceList>,
) -> Result<Vec<Chain>> {
    // 1. Validate all edges are manifold.
    for &eid in edge_ids {
        let faces = edge_face_map.get(&eid).ok_or(FilletError::EdgeNotFound)?;
        if faces.len() != 2 {
            return Err(FilletError::NonManifoldEdge(faces.len()));
        }
    }

    let selected: HashSet<EdgeId> = edge_ids.iter().copied().collect();

    // 2. For each face boundary, find contiguous runs of selected edges.
    let candidate_runs = collect_boundary_runs(shell, &selected);

    // 3. Collect final chains: process runs longest-first, claiming unclaimed
    //    edges and splitting on gaps where edges were already claimed.
    let mut run_order: Vec<usize> = (0..candidate_runs.len()).collect();
    run_order.sort_by(|&a, &b| {
        candidate_runs[b]
            .1
            .len()
            .cmp(&candidate_runs[a].1.len())
            .then(candidate_runs[a].0.cmp(&candidate_runs[b].0))
    });

    let mut globally_claimed: HashSet<EdgeId> = HashSet::new();
    let mut chains: Vec<Chain> = Vec::new();

    for run_idx in run_order {
        let (face_idx, ref run_edges) = candidate_runs[run_idx];
        // Split this run on already-claimed edges, producing sub-chains.
        let mut current_run: Vec<EdgeId> = Vec::new();
        for &eid in run_edges {
            if globally_claimed.contains(&eid) {
                if !current_run.is_empty() {
                    chains.push(Chain {
                        edge_ids: std::mem::take(&mut current_run),
                        shared_face_idx: face_idx,
                    });
                }
            } else {
                current_run.push(eid);
            }
        }
        if !current_run.is_empty() {
            chains.push(Chain {
                edge_ids: current_run,
                shared_face_idx: face_idx,
            });
        }
        // Mark all original run edges as globally claimed.
        globally_claimed.extend(run_edges);
    }

    // 4. Sort final chains: longest first, then shared_face_idx for determinism.
    chains.sort_by(|a, b| {
        b.edge_ids
            .len()
            .cmp(&a.edge_ids.len())
            .then(a.shared_face_idx.cmp(&b.shared_face_idx))
    });

    Ok(chains)
}

/// Original parameterization of a selected edge within its initial wire.
/// Stored before modifying any chain, so later fragments retain the same
/// variable-radius function instead of independently restarting at t=0.
struct OriginalEdgeSpan {
    // Keep the selected topology's Arc-backed identity alive throughout the
    // whole request. Otherwise replacing an earlier chain can free a selected
    // edge, then pointer reuse may make its obsolete EdgeId match another edge.
    source_edge: Edge,
    start: f64,
    span: f64,
}

fn original_edge_spans(
    shell: &Shell,
    chains: &[Chain],
) -> Result<HashMap<EdgeId, OriginalEdgeSpan>> {
    let mut spans = HashMap::new();
    for chain in chains {
        let boundary_edges = chain
            .edge_ids
            .iter()
            .map(|&id| {
                shell[chain.shared_face_idx]
                    .edge_iter()
                    .find(|edge| edge.id() == id)
                    .ok_or(FilletError::EdgeNotFound)
            })
            .collect::<Result<Vec<_>>>()?;
        let wire: Wire = boundary_edges.into();
        let (starts, lengths) =
            wire_edge_starts_and_spans(&wire).ok_or(FilletError::DiscontinuousWire)?;
        for (edge, (start, span)) in wire.edge_iter().zip(starts.into_iter().zip(lengths)) {
            spans.insert(
                edge.id(),
                OriginalEdgeSpan {
                    source_edge: edge.clone(),
                    start,
                    span,
                },
            );
        }
    }
    Ok(spans)
}

/// Evaluates a fragment's radius by projecting each position back onto the
/// originally selected edge, then mapping to the original wire's normalized
/// parameter. Works for reversed and nonlinearly reparameterized descendants.
fn radius_on_original_wire<'a>(
    fragment_curve: Curve,
    original: &'a OriginalEdgeSpan,
    radius: &'a dyn Fn(f64) -> f64,
) -> impl Fn(f64) -> f64 + 'a {
    let (fragment_start, fragment_end) = fragment_curve.range_tuple();
    let original_curve = original.source_edge.oriented_curve();
    let (original_start, original_end) = original_curve.range_tuple();
    move |t| {
        let parameter = (fragment_end - fragment_start).mul_add(t, fragment_start);
        let point = fragment_curve.evaluate(parameter);
        let global_t = original_curve
            .search_nearest_parameter(point, None, 100)
            .map_or(f64::NAN, |original_parameter| {
                let edge_t =
                    (original_parameter - original_start) / (original_end - original_start);
                original
                    .span
                    .mul_add(edge_t, original.start)
                    .clamp(0.0, 1.0)
            });
        radius(global_t)
    }
}

#[derive(Clone, Copy)]
struct VariableRadiusMapping<'a> {
    original: &'a OriginalEdgeSpan,
    radius: &'a dyn Fn(f64) -> f64,
}

/// Finds the side face adjacent to a given vertex in a chain, if any.
///
/// A side face is one that shares the vertex and an edge with one of the two
/// main faces but is not itself one of the main faces.
fn find_side_face(
    shell: &Shell,
    shared_face_idx: usize,
    opposite_face_idx: usize,
    vertex_id: VertexId,
    edge_face_map: &HashMap<EdgeId, FaceList>,
) -> Option<usize> {
    // Find edges incident to vertex_id in face_a or face_b.
    [shared_face_idx, opposite_face_idx]
        .iter()
        .flat_map(|&fi| {
            shell[fi].edge_iter().filter_map(move |edge| {
                let (front, back) = edge.ends();
                if front.id() == vertex_id || back.id() == vertex_id {
                    edge_face_map.get(&edge.id()).and_then(|faces| {
                        faces
                            .iter()
                            .find(|&&f| f != shared_face_idx && f != opposite_face_idx)
                            .copied()
                    })
                } else {
                    None
                }
            })
        })
        .next()
}

fn apply_single_edge_fillet(
    shell: &mut Shell,
    edge_id: EdgeId,
    preferred_shared_face_idx: Option<usize>,
    options: &FilletOptions,
    radius_mapping: Option<VariableRadiusMapping<'_>>,
) -> Result<()> {
    let edge_face_map = build_edge_face_map(shell);
    let faces = edge_face_map
        .get(&edge_id)
        .ok_or(FilletError::EdgeNotFound)?;
    if faces.len() != 2 {
        return Err(FilletError::NonManifoldEdge(faces.len()));
    }

    // Asymmetric chamfer setbacks are tied to the two face roles: do not
    // arbitrarily swap them if a replaced fragment no longer touches the
    // originally selected shared face.
    if matches!(&options.profile, FilletProfile::AsymmetricChamfer { .. })
        && preferred_shared_face_idx.is_some_and(|idx| !faces.contains(&idx))
    {
        return Err(FilletError::SharedFaceNotFound);
    }
    let shared_face_idx = preferred_shared_face_idx
        .filter(|idx| faces.contains(idx))
        .unwrap_or(faces[0]);
    let opposite_face_idx = *faces
        .iter()
        .find(|&&f| f != shared_face_idx)
        .ok_or(FilletError::EdgeNotFound)?;

    let face_a = shell[shared_face_idx].clone();
    let face_b = shell[opposite_face_idx].clone();
    let edge = face_a
        .edge_iter()
        .find(|edge| edge.id() == edge_id)
        .ok_or(FilletError::EdgeNotFound)?;
    let (front_vertex_id, back_vertex_id) = {
        let (front, back) = edge.ends();
        (front.id(), back.id())
    };

    let side0_idx = find_side_face(
        shell,
        shared_face_idx,
        opposite_face_idx,
        front_vertex_id,
        &edge_face_map,
    );
    let side1_idx = find_side_face(
        shell,
        shared_face_idx,
        opposite_face_idx,
        back_vertex_id,
        &edge_face_map,
    );

    let side0 = side0_idx.map(|idx| shell[idx].clone());
    let side1 = side1_idx.map(|idx| shell[idx].clone());

    let mapped_radius = radius_mapping
        .map(|mapping| {
            let fragment = edge.oriented_curve();
            curve_interval_on_original(&mapping.original.source_edge.oriented_curve(), &fragment)
                .ok_or(FilletError::UnresolvableSelectedEdge)?;
            Ok(radius_on_original_wire(
                fragment,
                mapping.original,
                mapping.radius,
            ))
        })
        .transpose()?;
    let fillet_result = ops::fillet_with_side_radius(
        &face_a,
        &face_b,
        edge_id,
        side0.as_ref(),
        side1.as_ref(),
        options,
        mapped_radius
            .as_ref()
            .map(|radius| radius as &dyn Fn(f64) -> f64),
    );
    let (new_face_a, new_face_b, fillet, new_side0, new_side1) = match fillet_result {
        Ok(result) => result,
        Err(error) => {
            if env::var_os("MT_FILLET_DEBUG").is_some() {
                let curve = edge.oriented_curve();
                let surface_a = face_a.oriented_surface();
                let surface_b = face_b.oriented_surface();
                let drift_a = sampled_curve_surface_max_distance(&curve, &surface_a, 16);
                let drift_b = sampled_curve_surface_max_distance(&curve, &surface_b, 16);
                eprintln!(
                    "debug fillet edge {edge_id:?}: drift_a={drift_a:.3e}, drift_b={drift_b:.3e}, error={error}"
                );
            }
            return Err(error);
        }
    };

    shell[shared_face_idx] = new_face_a;
    shell[opposite_face_idx] = new_face_b;
    if let Some(new_side) = new_side0
        && let Some(idx) = side0_idx
    {
        shell[idx] = new_side;
    }
    if let Some(new_side) = new_side1
        && let Some(idx) = side1_idx
    {
        shell[idx] = new_side;
    }
    shell.push(fillet);
    Ok(())
}

/// Fillets the specified edges of a shell.
///
/// Resolves face adjacency automatically and dispatches to
/// [`fillet`](super::fillet)/[`fillet_with_side`](super::fillet_with_side)
/// for single edges or [`fillet_along_wire`](super::fillet_along_wire) for multi-edge chains.
///
/// # Errors
///
/// Returns [`FilletError`] when an edge is missing, duplicated, non-manifold,
/// or no longer uniquely resolvable, when radius/profile validation fails, or
/// when the required geometry/topology construction cannot be completed.
pub fn fillet_edges(
    shell: &mut Shell,
    edge_ids: &[EdgeId],
    params: Option<&FilletOptions>,
) -> Result<()> {
    // One request may contain several independent chains. Per-wire and
    // per-chain transactions are not enough: a hard failure in a later chain
    // must also roll back every chain successfully edited earlier.
    ops::edit_shell_transactionally(shell, |candidate| {
        fillet_edges_in_place(candidate, edge_ids, params)
    })
}

#[expect(
    clippy::too_many_lines,
    reason = "the transactional chain resolver deliberately keeps selection rematching and per-edge fallback in one mutation scope"
)]
fn fillet_edges_in_place(
    shell: &mut Shell,
    edge_ids: &[EdgeId],
    params: Option<&FilletOptions>,
) -> Result<()> {
    let default_options = FilletOptions::default();
    let options = params.unwrap_or(&default_options);
    if matches!(&options.profile, FilletProfile::AsymmetricChamfer { .. }) && edge_ids.len() != 1 {
        return Err(FilletError::AsymmetricChamferWireUnsupported);
    }

    // A repeated edge cannot be independently filleted twice, and a per-edge
    // radius list cannot assign two different radii to the same topology ID.
    let mut distinct_ids = HashSet::new();
    for &eid in edge_ids {
        if !distinct_ids.insert(eid) {
            return Err(FilletError::DuplicateSelectedEdge);
        }
    }

    // Validate all requested edges exist and are manifold up-front.
    {
        let edge_face_map = build_edge_face_map(shell);
        for &eid in edge_ids {
            let faces = edge_face_map.get(&eid).ok_or(FilletError::EdgeNotFound)?;
            if faces.len() != 2 {
                return Err(FilletError::NonManifoldEdge(faces.len()));
            }
        }
    }

    // Validate per-edge radius count.
    if options.profile.uses_radius()
        && let RadiusSpec::PerEdge(ref radii) = options.radius
        && radii.len() != edge_ids.len()
    {
        return Err(FilletError::PerEdgeRadiusMismatch {
            given: radii.len(),
            expected: edge_ids.len(),
        });
    }

    let variable_radius_upper_bound = options.validated_radius_upper_bound()?;

    // Reject edges that are too short for the requested rolling radius.
    // Asymmetric chamfer setbacks are perpendicular to the seam and do not
    // impose this along-edge radius bound.
    if options.profile.uses_radius() {
        for (i, &eid) in edge_ids.iter().enumerate() {
            let effective_radius = match &options.radius {
                RadiusSpec::Constant(r) => *r,
                RadiusSpec::Variable(_) => variable_radius_upper_bound,
                RadiusSpec::PerEdge(radii) => radii[i],
            };
            let edge = shell
                .edge_iter()
                .find(|e| e.id() == eid)
                .ok_or(FilletError::EdgeNotFound)?;
            let length = approximate_edge_length(&edge, 16);
            if length < 2.0 * effective_radius {
                return Err(FilletError::DegenerateEdge);
            }
        }
    }

    // Group edges into chains (by face pair topology). Chain stores only edge IDs.
    let chains = {
        let edge_face_map = build_edge_face_map(shell);
        group_edges_into_chains(shell, edge_ids, &edge_face_map)?
    };

    // Build EdgeId->index map for PerEdge radius lookup.
    let edge_id_to_idx: HashMap<EdgeId, usize> = edge_ids
        .iter()
        .enumerate()
        .map(|(i, &eid)| (eid, i))
        .collect();
    let original_spans = original_edge_spans(shell, &chains)?;

    // Process each chain with a fresh edge-face map, so that mutations from
    // earlier chains don't cause stale face indices.
    for chain in &chains {
        let edge_face_map = build_edge_face_map(shell);
        let mut used_ids = HashSet::<EdgeId>::new();
        let mut resolved_pairs: Vec<(EdgeId, EdgeId)> = Vec::new();
        for original_eid in chain.edge_ids.iter().copied() {
            let original_curve = original_spans
                .get(&original_eid)
                .ok_or(FilletError::EdgeNotFound)?
                .source_edge
                .oriented_curve();
            for resolved_id in resolve_active_edge_ids(
                shell,
                original_eid,
                &original_curve,
                &used_ids,
                &edge_face_map,
            )? {
                used_ids.insert(resolved_id);
                resolved_pairs.push((original_eid, resolved_id));
            }
        }
        let resolved_edge_ids: Vec<EdgeId> = resolved_pairs
            .iter()
            .map(|(_, resolved_eid)| *resolved_eid)
            .collect();
        let resolved_to_original: HashMap<EdgeId, EdgeId> = resolved_pairs
            .iter()
            .map(|(original_eid, resolved_eid)| (*resolved_eid, *original_eid))
            .collect();

        // Resolve face_a (shared face from grouping) and face_b from current map.
        let first_eid = resolved_edge_ids
            .first()
            .copied()
            .ok_or(FilletError::EdgeNotFound)?;
        let faces = edge_face_map
            .get(&first_eid)
            .ok_or(FilletError::EdgeNotFound)?;
        if faces.len() != 2 {
            return Err(FilletError::NonManifoldEdge(faces.len()));
        }
        // face_a is the shared face from boundary-run grouping; face_b is the other.
        let shared_face_idx = if faces.contains(&chain.shared_face_idx) {
            chain.shared_face_idx
        } else {
            faces[0]
        };

        if chain.edge_ids.len() == 1 && resolved_pairs.len() == 1 {
            let chain_opts;
            let opts = if options.profile.uses_radius()
                && let RadiusSpec::PerEdge(radii) = &options.radius
            {
                let original_eid = chain.edge_ids[0];
                chain_opts = FilletOptions {
                    radius: RadiusSpec::Constant(radii[edge_id_to_idx[&original_eid]]),
                    divisions: options.divisions,
                    profile: options.profile.clone(),
                };
                &chain_opts
            } else {
                options
            };

            match apply_single_edge_fillet(shell, first_eid, Some(shared_face_idx), opts, None) {
                Ok(()) | Err(FilletError::GeometryFailed { .. }) => {}
                Err(error) => return Err(error),
            }
        } else {
            // An original edge can have several surviving subsegments after
            // another fillet. Do not pass a larger, partially disconnected
            // selection to the original wire constructor; visit every segment
            // transactionally through the per-edge path below.
            let has_split_edge = resolved_pairs.len() != chain.edge_ids.len();
            let wire_fillet_failed = if has_split_edge {
                true
            } else {
                // Multi-edge chain: build a Wire from the shared face's boundary
                // to get edges in the correct orientation for that face.
                let chain_id_set: HashSet<EdgeId> = resolved_edge_ids.iter().copied().collect();
                let wire: Option<Wire> = shell[shared_face_idx]
                    .boundary_iters()
                    .into_iter()
                    .find_map(|boundary| {
                        let edges = boundary
                            .filter(|edge| chain_id_set.contains(&edge.id()))
                            .collect::<Vec<_>>();
                        (edges.len() == resolved_edge_ids.len()).then_some(edges.into())
                    });
                if let Some(wire) = wire {
                    let chain_opts;
                    let opts = if options.profile.uses_radius()
                        && let RadiusSpec::PerEdge(radii) = &options.radius
                    {
                        let chain_radii: Vec<f64> = wire
                            .edge_iter()
                            .map(|edge| {
                                resolved_to_original
                                    .get(&edge.id())
                                    .copied()
                                    .ok_or(FilletError::EdgeNotFound)
                                    .map(|eid| radii[edge_id_to_idx[&eid]])
                            })
                            .collect::<Result<_>>()?;
                        chain_opts = FilletOptions {
                            radius: RadiusSpec::PerEdge(chain_radii),
                            divisions: options.divisions,
                            profile: options.profile.clone(),
                        };
                        &chain_opts
                    } else {
                        options
                    };

                    match ops::fillet_along_wire(shell, &wire, opts) {
                        Ok(()) => false,
                        Err(
                            error @ (FilletError::FilletSurfaceComputationFailed
                            | FilletError::GeometryFailed { .. }
                            | FilletError::SharedFaceNotFound
                            | FilletError::AdjacentFacesNotFound
                            | FilletError::DiscontinuousWire),
                        ) => {
                            if env::var_os("MT_FILLET_DEBUG").is_some() {
                                eprintln!(
                                    "debug fillet chain {:?}: wire fillet failed ({error}), falling back per edge.",
                                    chain.edge_ids
                                );
                            }
                            true
                        }
                        Err(error) => return Err(error),
                    }
                } else {
                    // Carrier-face reconstruction may have split the original
                    // run over several boundaries. The per-edge visitor can
                    // still apply every surviving piece without a common wire.
                    if env::var_os("MT_FILLET_DEBUG").is_some() {
                        eprintln!(
                            "debug fillet chain {:?}: no shared carrier wire; retrying per-edge.",
                            chain.edge_ids
                        );
                    }
                    true
                }
            };

            if wire_fillet_failed {
                // The wire attempt is atomic. Run each fallback edge against
                // a second detached shell, and publish the chain only after
                // every edge has succeeded. Hard failures roll back too.
                let pending: VecDeque<(EdgeId, EdgeId, Curve)> = resolved_pairs
                    .iter()
                    .map(|&(original_id, active_id)| {
                        shell
                            .edge_iter()
                            .find(|edge| edge.id() == active_id)
                            .map(|edge| (original_id, active_id, edge.oriented_curve()))
                            .ok_or(FilletError::EdgeNotFound)
                    })
                    .collect::<Result<_>>()?;
                let fallback_result = ops::edit_shell_transactionally(shell, |candidate| {
                    visit_surviving_fragments(
                        candidate,
                        pending,
                        |current, original_id, active_id| {
                            let single_edge_opts;
                            let opts = if options.profile.uses_radius()
                                && let RadiusSpec::PerEdge(radii) = &options.radius
                            {
                                single_edge_opts = FilletOptions {
                                    radius: RadiusSpec::Constant(
                                        radii[edge_id_to_idx[&original_id]],
                                    ),
                                    divisions: options.divisions,
                                    profile: options.profile.clone(),
                                };
                                &single_edge_opts
                            } else {
                                options
                            };
                            let radius_mapping =
                                if let RadiusSpec::Variable(radius) = &options.radius {
                                    Some(VariableRadiusMapping {
                                        original: original_spans
                                            .get(&original_id)
                                            .ok_or(FilletError::EdgeNotFound)?,
                                        radius: radius.as_ref(),
                                    })
                                } else {
                                    None
                                };
                            apply_single_edge_fillet(
                                current,
                                active_id,
                                Some(chain.shared_face_idx),
                                opts,
                                radius_mapping,
                            )
                        },
                    )
                });
                match fallback_result {
                    Ok(()) | Err(FilletError::GeometryFailed { .. }) => {}
                    Err(error) => return Err(error),
                }
            }
        }
    }
    Ok(())
}

/// Fillets the specified edges of a shell with arbitrary curve/surface types.
///
/// Converts the shell to internal NURBS types, runs [`fillet_edges`], and
/// converts back. This is the main entry point for external callers whose
/// shells use types like `monstertruck_modeling::Curve` / `monstertruck_modeling::Surface`.
///
/// # Errors
///
/// Returns [`FilletError`] when conversion to/from the internal NURBS representation fails or
/// when [`fillet_edges`] rejects or cannot construct the requested fillet.
pub fn fillet_edges_generic<C, S>(
    shell: &mut monstertruck_topology::Shell<Point3, C, S>,
    edges: &[monstertruck_topology::Edge<Point3, C>],
    params: Option<&FilletOptions>,
) -> Result<()>
where
    C: FilletableCurve,
    S: FilletableSurface,
{
    let default_options = FilletOptions::default();
    let options = params.unwrap_or(&default_options);
    let (mut internal_shell, internal_edge_ids) = convert_shell_in(shell, edges)?;
    let original_shell = internal_shell.clone();
    // The NURBS shell is already a detached conversion. No additional
    // request-level clone is needed: the caller's shell is only updated after
    // both filleting and conversion back have succeeded.
    fillet_edges_in_place(&mut internal_shell, &internal_edge_ids, Some(options))?;
    if internal_shell.shell_condition() != ShellCondition::Closed
        && env::var_os("MT_FILLET_STRICT_CLOSED").is_some()
    {
        if env::var_os("MT_FILLET_DEBUG").is_some() {
            eprintln!(
                "debug fillet generic: rollback to original shell (non-closed result, strict mode)."
            );
        }
        internal_shell = original_shell;
    }
    *shell = convert_shell_out(&internal_shell)?;
    Ok(())
}

#[cfg(test)]
mod rematch_tests {
    use super::*;
    use monstertruck_geometry::prelude::{BsplineCurve, BsplineSurface, KnotVector, NurbsCurve};
    use monstertruck_topology::{Face, Vertex};

    fn line(start: Point3, end: Point3) -> Curve {
        Curve::NurbsCurve(NurbsCurve::from(BsplineCurve::new_unchecked(
            KnotVector::bezier_knot(1),
            vec![start, end],
        )))
    }

    /// Each requested segment has two incident faces and its own topology.
    /// Segments may share coordinates without sharing edge identities.
    fn manifold_shell(curves: &[(Point3, Point3, Curve)]) -> (Shell, Vec<EdgeId>) {
        let mut faces = Vec::new();
        let mut ids = Vec::new();
        for (front, back, curve) in curves {
            let offset = monstertruck_geometry::prelude::Vector3::new(0.0, 1.0, 0.0);
            let vertices = Vertex::from_points([*front, *back, *back + offset, *front + offset]);
            let selected = Edge::new_unchecked(&vertices[0], &vertices[1], curve.clone());
            ids.push(selected.id());
            let other_edges = [(1, 2), (2, 3), (3, 0)];
            let mut boundary = vec![selected];
            for (a, b) in other_edges {
                boundary.push(Edge::new_unchecked(
                    &vertices[a],
                    &vertices[b],
                    line(vertices[a].point(), vertices[b].point()),
                ));
            }
            let knots = KnotVector::bezier_knot(1);
            let surface = NurbsSurface::from(BsplineSurface::new_unchecked(
                (knots.clone(), knots),
                vec![
                    vec![vertices[0].point(), vertices[3].point()],
                    vec![vertices[1].point(), vertices[2].point()],
                ],
            ));
            let face = Face::new_unchecked(vec![boundary.into()], surface);
            faces.push(face.clone());
            faces.push(face.inverse());
        }
        (Shell::from(faces), ids)
    }

    #[test]
    fn indistinguishable_manifold_edges_are_not_rematched_arbitrarily() {
        let a = Point3::new(0.0, 0.0, 0.0);
        let b = Point3::new(1.0, 0.0, 0.0);
        let original = line(a, b);
        let (shell, ids) = manifold_shell(&[(a, b, original.clone()), (a, b, original.clone())]);
        let faces = build_edge_face_map(&shell);
        assert_eq!(faces.get(&ids[0]).map(|faces| faces.len()), Some(2));
        assert_eq!(faces.get(&ids[1]).map(|faces| faces.len()), Some(2));
        assert_eq!(
            rematch_selected_edge_ids(&shell, &original, &HashSet::new(), &faces),
            None
        );

        let used = HashSet::from([ids[0]]);
        assert_eq!(
            rematch_selected_edge_ids(&shell, &original, &used, &faces),
            Some(vec![ids[1]])
        );
    }

    #[test]
    fn ambiguous_rematch_is_an_error_instead_of_a_partial_fillet() {
        let a = Point3::new(0.0, 0.0, 0.0);
        let b = Point3::new(1.0, 0.0, 0.0);
        let original = line(a, b);
        let (shell, ids) = manifold_shell(&[(a, b, original.clone()), (a, b, original.clone())]);
        let face_map = build_edge_face_map(&shell);
        let missing_id =
            Edge::new_unchecked(&Vertex::new(a), &Vertex::new(b), original.clone()).id();
        assert!(!face_map.contains_key(&missing_id));
        assert!(matches!(
            resolve_active_edge_ids(&shell, missing_id, &original, &HashSet::new(), &face_map),
            Err(FilletError::UnresolvableSelectedEdge),
        ));
        // The exact identity remains authoritative even if coincident
        // alternatives would make geometric matching ambiguous.
        assert_eq!(
            resolve_active_edge_ids(&shell, ids[0], &original, &HashSet::new(), &face_map).unwrap(),
            vec![ids[0]],
        );
    }

    #[test]
    fn nonmanifold_candidate_is_an_error_instead_of_a_silent_success() {
        let a = Point3::new(0.0, 0.0, 0.0);
        let b = Point3::new(1.0, 0.0, 0.0);
        let (mut shell, ids) = manifold_shell(&[(a, b, line(a, b))]);
        shell.pop();
        assert!(matches!(
            apply_single_edge_fillet(&mut shell, ids[0], None, &FilletOptions::default(), None),
            Err(FilletError::NonManifoldEdge(1)),
        ));
    }

    #[test]
    fn rematch_accepts_non_grid_aligned_subdivision() {
        let a = Point3::new(0.0, 0.0, 0.0);
        let b = Point3::new(1.0, 0.0, 0.0);
        let cut = Point3::new(0.37, 0.0, 0.0);
        let original = line(a, b);
        let (shell, ids) = manifold_shell(&[(a, cut, line(a, cut))]);
        let faces = build_edge_face_map(&shell);
        assert_eq!(
            rematch_selected_edge_ids(&shell, &original, &HashSet::new(), &faces),
            Some(vec![ids[0]])
        );
    }

    #[test]
    fn rematch_rejects_overlapping_survivors_even_when_lengths_differ() {
        let a = Point3::new(0.0, 0.0, 0.0);
        let b = Point3::new(1.0, 0.0, 0.0);
        let short = Point3::new(0.37, 0.0, 0.0);
        let long = Point3::new(0.72, 0.0, 0.0);
        let original = line(a, b);
        let (shell, ids) = manifold_shell(&[(a, short, line(a, short)), (a, long, line(a, long))]);
        let faces = build_edge_face_map(&shell);
        assert_eq!(
            rematch_selected_edge_ids(&shell, &original, &HashSet::new(), &faces),
            None,
        );
        assert_eq!(ids.len(), 2);
    }

    #[test]
    fn rematch_returns_every_nonoverlapping_descendant_in_parameter_order() {
        let a = Point3::new(0.0, 0.0, 0.0);
        let b = Point3::new(1.0, 0.0, 0.0);
        let split = Point3::new(0.37, 0.0, 0.0);
        let original = line(a, b);
        // Store the later part first to check sorting, and reverse its curve
        // orientation. Neither must affect which pieces are visited.
        let (shell, ids) =
            manifold_shell(&[(b, split, line(b, split)), (a, split, line(a, split))]);
        let faces = build_edge_face_map(&shell);
        assert_eq!(
            rematch_selected_edge_ids(&shell, &original, &HashSet::new(), &faces),
            Some(vec![ids[1], ids[0]]),
        );
    }

    #[test]
    fn original_edge_snapshot_retains_selected_ids_after_topology_removal() {
        let a = Point3::new(0.0, 0.0, 0.0);
        let b = Point3::new(1.0, 0.0, 0.0);
        let (mut shell, selected) = manifold_shell(&[(a, b, line(a, b))]);
        let face_map = build_edge_face_map(&shell);
        let chains = group_edges_into_chains(&shell, &selected, &face_map).unwrap();
        let originals = original_edge_spans(&shell, &chains).unwrap();
        let source = &originals[&selected[0]].source_edge;
        assert_eq!(source.id(), selected[0]);
        shell.clear();
        assert_eq!(source.id(), selected[0]);
        assert!(curve_interval_on_original(&line(a, b), &source.oriented_curve()).is_some());

        // A selected ID represents an Arc allocation. As long as the original
        // snapshot owns a source handle, replacing the shell cannot recycle
        // that allocation for unrelated edges during a later chain.
        for _ in 0..128 {
            let candidate = Edge::new_unchecked(&Vertex::new(a), &Vertex::new(b), line(a, b));
            assert_ne!(candidate.id(), selected[0]);
        }
    }

    #[test]
    fn rematch_visits_both_survivors_after_an_intermediate_piece_was_removed() {
        let a = Point3::new(0.0, 0.0, 0.0);
        let b = Point3::new(1.0, 0.0, 0.0);
        let left = Point3::new(0.24, 0.0, 0.0);
        let right = Point3::new(0.81, 0.0, 0.0);
        let original = line(a, b);
        let (shell, ids) = manifold_shell(&[(right, b, line(right, b)), (a, left, line(a, left))]);
        let faces = build_edge_face_map(&shell);
        assert_eq!(
            rematch_selected_edge_ids(&shell, &original, &HashSet::new(), &faces),
            Some(vec![ids[1], ids[0]]),
        );
    }

    #[test]
    fn visiting_survivors_handles_fragments_split_again_by_an_earlier_edit() {
        let a = Point3::new(0.0, 0.0, 0.0);
        let b = Point3::new(1.0, 0.0, 0.0);
        let first_cut = Point3::new(0.37, 0.0, 0.0);
        let second_cut = Point3::new(0.68, 0.0, 0.0);
        let (mut shell, ids) = manifold_shell(&[
            (a, first_cut, line(a, first_cut)),
            (first_cut, b, line(first_cut, b)),
        ]);
        let original_curve = line(a, b);
        let original_id =
            Edge::new_unchecked(&Vertex::new(a), &Vertex::new(b), original_curve.clone()).id();
        let mut visited = Vec::new();
        let mut expected_new_ids = Vec::new();
        visit_surviving_fragments(
            &mut shell,
            VecDeque::from([(original_id, original_id, original_curve)]),
            |shell, parent_id, current_id| {
                assert_eq!(parent_id, original_id);
                visited.push(current_id);
                if current_id == ids[0] {
                    // A fillet on the first fragment modifies topology before
                    // its sibling is processed. Emulate the resulting split.
                    shell.pop();
                    shell.pop();
                    let (later, fragment_ids) = manifold_shell(&[
                        (first_cut, second_cut, line(first_cut, second_cut)),
                        (second_cut, b, line(second_cut, b)),
                    ]);
                    expected_new_ids = fragment_ids;
                    shell.extend(later);
                }
                Ok(())
            },
        )
        .unwrap();
        assert_eq!(visited, [ids[0], expected_new_ids[0], expected_new_ids[1]]);
        assert!(!visited.contains(&ids[1]));
    }

    #[test]
    fn variable_radius_uses_the_original_wire_position_for_reversed_fragments() {
        let a = Point3::new(0.0, 0.0, 0.0);
        let b = Point3::new(1.0, 0.0, 0.0);
        let original = OriginalEdgeSpan {
            source_edge: Edge::new_unchecked(&Vertex::new(a), &Vertex::new(b), line(a, b)),
            start: 0.2,
            span: 0.3,
        };
        let linear = line(Point3::new(0.75, 0.0, 0.0), Point3::new(0.25, 0.0, 0.0));
        let mapped = radius_on_original_wire(linear, &original, &|t| t);
        assert!((mapped(0.0) - 0.425).abs() < 1e-5);
        assert!((mapped(0.5) - 0.35).abs() < 1e-5);
        assert!((mapped(1.0) - 0.275).abs() < 1e-5);
    }

    #[test]
    fn variable_radius_tracks_nonlinear_reparameterized_fragment_geometry() {
        let a = Point3::new(0.0, 0.0, 0.0);
        let b = Point3::new(1.0, 0.0, 0.0);
        let original = OriginalEdgeSpan {
            source_edge: Edge::new_unchecked(&Vertex::new(a), &Vertex::new(b), line(a, b)),
            start: 0.2,
            span: 0.3,
        };
        let nonlinear = Curve::NurbsCurve(NurbsCurve::from(BsplineCurve::new_unchecked(
            KnotVector::bezier_knot(2),
            vec![
                Point3::new(0.25, 0.0, 0.0),
                Point3::new(0.25, 0.0, 0.0),
                Point3::new(0.75, 0.0, 0.0),
            ],
        )));
        let mapped = radius_on_original_wire(nonlinear, &original, &|t| t);
        assert!((mapped(0.0) - 0.275).abs() < 1e-5);
        // An endpoint-only affine map would incorrectly return 0.35 here.
        assert!((mapped(0.5) - 0.3125).abs() < 1e-5);
        assert!((mapped(1.0) - 0.425).abs() < 1e-5);
    }

    #[test]
    fn rematch_ignores_boundary_edges_even_with_identical_geometry() {
        let a = Point3::new(0.0, 0.0, 0.0);
        let b = Point3::new(1.0, 0.0, 0.0);
        let original = line(a, b);
        let (mut shell, _) = manifold_shell(&[(a, b, original.clone())]);
        shell.pop(); // The remaining edge is incident to only one face.
        let faces = build_edge_face_map(&shell);
        assert_eq!(
            rematch_selected_edge_ids(&shell, &original, &HashSet::new(), &faces),
            None,
        );
    }

    #[test]
    fn rematch_rejects_curve_with_same_endpoints_and_midpoint() {
        let a = Point3::new(0.0, 0.0, 0.0);
        let b = Point3::new(1.0, 0.0, 0.0);
        let original = line(a, b);
        let different = Curve::NurbsCurve(NurbsCurve::from(BsplineCurve::new_unchecked(
            KnotVector::bezier_knot(3),
            vec![
                a,
                Point3::new(1.0 / 3.0, 1.0, 0.0),
                Point3::new(2.0 / 3.0, -1.0, 0.0),
                b,
            ],
        )));
        let (shell, _) = manifold_shell(&[(a, b, different)]);
        let faces = build_edge_face_map(&shell);
        assert_eq!(
            rematch_selected_edge_ids(&shell, &original, &HashSet::new(), &faces),
            None
        );
    }
}

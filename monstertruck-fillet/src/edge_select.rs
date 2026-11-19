use std::collections::HashMap;
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

/// Checks a point against the actual curve rather than a fixed parameter grid.
/// Fixed-grid matching rejects valid subdivided edges whose cut parameters
/// do not coincide with the grid points.
fn curve_contains_point(curve: &Curve, point: Point3) -> bool {
    curve
        .search_nearest_parameter(point, None, 100)
        .is_some_and(|t| {
            squared_distance(curve.evaluate(t), point) < KERNEL_EPSILON * KERNEL_EPSILON
        })
}

fn sampled_variable_radius_upper_bound(radius: &dyn Fn(f64) -> f64, sample_count: usize) -> f64 {
    (0..=sample_count)
        .map(|i| radius(crate::count_as_f64(i) / crate::count_as_f64(sample_count)))
        .fold(f64::NEG_INFINITY, f64::max)
}

fn rematch_selected_edge_id(
    shell: &Shell,
    original_curve: &Curve,
    used_ids: &HashSet<EdgeId>,
    edge_face_map: &HashMap<EdgeId, FaceList>,
) -> Option<EdgeId> {
    // The same topological edge appears in both incident faces. Deduplicate it
    // before deciding whether several distinct candidates are indistinguishable.
    let mut seen_ids = HashSet::new();
    let mut candidates = shell
        .edge_iter()
        .filter(|edge| {
            seen_ids.insert(edge.id())
                && !used_ids.contains(&edge.id())
                && is_manifold_edge(edge_face_map, edge.id())
        })
        .filter(|edge| {
            let curve = edge.curve();
            let (start, end) = curve.range_tuple();
            // Interior samples distinguish different paths with the same ends.
            [0.0, 0.25, 0.5, 0.75, 1.0].into_iter().all(|fraction| {
                let t = (end - start).mul_add(fraction, start);
                curve_contains_point(original_curve, curve.evaluate(t))
            })
        })
        .map(|edge| (edge.id(), approximate_edge_length(&edge, 8)))
        .collect::<Vec<_>>();

    // When topology has changed there is no identity-preserving mapping. A
    // unique longest segment is a useful best-effort successor; equal-length
    // successors are ambiguous and must not be selected arbitrarily.
    candidates.sort_by(|a, b| b.1.total_cmp(&a.1));
    match candidates.as_slice() {
        [(id, _)] => Some(*id),
        [(id, longest), (_, runner_up), ..] if longest - runner_up > KERNEL_EPSILON => Some(*id),
        _ => None,
    }
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

    // 3. Each selected edge appears in exactly 2 candidate runs (one per face).
    //    Assign each edge to its longest run (tiebreak: lower face_idx).
    let mut edge_best_run: HashMap<EdgeId, usize> = HashMap::new();
    for (run_idx, (face_idx, run_edges)) in candidate_runs.iter().enumerate() {
        for &eid in run_edges {
            let replace = match edge_best_run.get(&eid) {
                None => true,
                Some(&prev_idx) => {
                    let prev_len = candidate_runs[prev_idx].1.len();
                    let this_len = run_edges.len();
                    this_len > prev_len
                        || (this_len == prev_len && *face_idx < candidate_runs[prev_idx].0)
                }
            };
            if replace {
                edge_best_run.insert(eid, run_idx);
            }
        }
    }

    // 4. Collect final chains: process runs longest-first, claiming unclaimed
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

    // 5. Sort final chains: longest first, then shared_face_idx for determinism.
    chains.sort_by(|a, b| {
        b.edge_ids
            .len()
            .cmp(&a.edge_ids.len())
            .then(a.shared_face_idx.cmp(&b.shared_face_idx))
    });

    Ok(chains)
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
) -> Result<()> {
    let edge_face_map = build_edge_face_map(shell);
    let faces = edge_face_map
        .get(&edge_id)
        .ok_or(FilletError::EdgeNotFound)?;
    if faces.len() != 2 {
        if env::var_os("MT_FILLET_DEBUG").is_some() {
            eprintln!(
                "debug fillet edge {edge_id:?}: skipped non-manifold candidate ({} faces).",
                faces.len()
            );
        }
        return Ok(());
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

    let fillet_result = ops::fillet_with_side(
        &face_a,
        &face_b,
        edge_id,
        side0.as_ref(),
        side1.as_ref(),
        options,
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
/// Returns [`FilletError`] when an edge is missing/non-manifold, radius/profile validation fails,
/// or the required geometry/topology construction cannot be completed.
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

    // Reject edges that are too short for the requested rolling radius.
    // Asymmetric chamfer setbacks are perpendicular to the seam and do not
    // impose this along-edge radius bound.
    if options.profile.uses_radius() {
        let variable_radius_upper_bound = if let RadiusSpec::Variable(f) = &options.radius {
            sampled_variable_radius_upper_bound(f.as_ref(), 32)
        } else {
            0.0
        };
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
    let original_edge_curves: HashMap<EdgeId, Curve> = edge_ids
        .iter()
        .copied()
        .map(|eid| {
            shell
                .edge_iter()
                .find(|edge| edge.id() == eid)
                .map(|edge| (eid, edge.oriented_curve()))
                .ok_or(FilletError::EdgeNotFound)
        })
        .collect::<Result<_>>()?;

    // Process each chain with a fresh edge-face map, so that mutations from
    // earlier chains don't cause stale face indices.
    for chain in &chains {
        let edge_face_map = build_edge_face_map(shell);
        let mut used_ids = HashSet::<EdgeId>::new();
        let mut resolved_pairs: Vec<(EdgeId, EdgeId)> = Vec::new();
        for original_eid in chain.edge_ids.iter().copied() {
            let resolved_eid = if is_manifold_edge(&edge_face_map, original_eid)
                && !used_ids.contains(&original_eid)
            {
                Some(original_eid)
            } else {
                original_edge_curves
                    .get(&original_eid)
                    .and_then(|original_curve| {
                        rematch_selected_edge_id(shell, original_curve, &used_ids, &edge_face_map)
                    })
            };
            if let Some(resolved_eid) = resolved_eid {
                used_ids.insert(resolved_eid);
                resolved_pairs.push((original_eid, resolved_eid));
            } else if env::var_os("MT_FILLET_DEBUG").is_some() {
                eprintln!(
                    "debug fillet chain {:?}: unresolved edge {original_eid:?}, skipping.",
                    chain.edge_ids
                );
            }
        }
        if resolved_pairs.is_empty() {
            continue;
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
        let Some(first_eid) = resolved_edge_ids.first().copied() else {
            continue;
        };
        let Some(faces) = edge_face_map.get(&first_eid) else {
            if env::var_os("MT_FILLET_DEBUG").is_some() {
                eprintln!(
                    "debug fillet chain {:?}: missing seed edge {first_eid:?}, skipping.",
                    chain.edge_ids
                );
            }
            continue;
        };
        if faces.len() != 2 {
            if env::var_os("MT_FILLET_DEBUG").is_some() {
                eprintln!(
                    "debug fillet chain {:?}: skipped non-manifold seed edge {first_eid:?} ({} faces).",
                    chain.edge_ids,
                    faces.len()
                );
            }
            continue;
        }
        // face_a is the shared face from boundary-run grouping; face_b is the other.
        let shared_face_idx = if faces.contains(&chain.shared_face_idx) {
            chain.shared_face_idx
        } else {
            faces[0]
        };

        if chain.edge_ids.len() == 1 {
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

            match apply_single_edge_fillet(shell, first_eid, Some(shared_face_idx), opts) {
                Ok(()) | Err(FilletError::GeometryFailed { .. }) => {}
                Err(error) => return Err(error),
            }
        } else {
            // Multi-edge chain: build a Wire from the shared face's boundary
            // to get edges in the correct orientation for that face.
            let chain_id_set: HashSet<EdgeId> = resolved_edge_ids.iter().copied().collect();
            let Some(wire): Option<Wire> = shell[shared_face_idx]
                .boundary_iters()
                .into_iter()
                .find_map(|boundary| {
                    let edges = boundary
                        .filter(|edge| chain_id_set.contains(&edge.id()))
                        .collect::<Vec<_>>();
                    (edges.len() == resolved_edge_ids.len()).then_some(edges.into())
                })
            else {
                if env::var_os("MT_FILLET_DEBUG").is_some() {
                    eprintln!(
                        "debug fillet chain {:?}: failed to construct wire, skipping.",
                        chain.edge_ids
                    );
                }
                continue;
            };

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

            let wire_fillet_failed = match ops::fillet_along_wire(shell, &wire, opts) {
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
            };

            if wire_fillet_failed {
                // The wire attempt is atomic. Run each fallback edge against
                // a second detached shell, and publish the chain only after
                // every edge has succeeded. Hard failures roll back too.
                let fallback_result = ops::edit_shell_transactionally(shell, |candidate| {
                    let original_to_resolved: HashMap<EdgeId, EdgeId> =
                        resolved_pairs.iter().copied().collect();
                    let mut fallback_used_ids = HashSet::<EdgeId>::new();
                    for original_eid in chain.edge_ids.iter().copied() {
                        let edge_face_map = build_edge_face_map(candidate);
                        let initial_eid = original_to_resolved
                            .get(&original_eid)
                            .copied()
                            .unwrap_or(original_eid);
                        let resolved_eid = if is_manifold_edge(&edge_face_map, initial_eid)
                            && !fallback_used_ids.contains(&initial_eid)
                        {
                            initial_eid
                        } else {
                            let original_curve = original_edge_curves.get(&original_eid).ok_or(
                                FilletError::GeometryFailed {
                                    context: "resolve per-edge fillet fallback curve",
                                },
                            )?;
                            rematch_selected_edge_id(
                                candidate,
                                original_curve,
                                &fallback_used_ids,
                                &edge_face_map,
                            )
                            .ok_or(FilletError::GeometryFailed {
                                context: "resolve per-edge fillet fallback edge",
                            })?
                        };
                        fallback_used_ids.insert(resolved_eid);

                        let single_edge_opts;
                        let opts = if options.profile.uses_radius()
                            && let RadiusSpec::PerEdge(radii) = &options.radius
                        {
                            single_edge_opts = FilletOptions {
                                radius: RadiusSpec::Constant(radii[edge_id_to_idx[&original_eid]]),
                                divisions: options.divisions,
                                profile: options.profile.clone(),
                            };
                            &single_edge_opts
                        } else {
                            options
                        };
                        apply_single_edge_fillet(
                            candidate,
                            resolved_eid,
                            Some(chain.shared_face_idx),
                            opts,
                        )?;
                    }
                    Ok(())
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
            rematch_selected_edge_id(&shell, &original, &HashSet::new(), &faces),
            None
        );

        let used = HashSet::from([ids[0]]);
        assert_eq!(
            rematch_selected_edge_id(&shell, &original, &used, &faces),
            Some(ids[1])
        );
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
            rematch_selected_edge_id(&shell, &original, &HashSet::new(), &faces),
            Some(ids[0])
        );
    }

    #[test]
    fn rematch_prefers_the_unique_longest_surviving_segment() {
        let a = Point3::new(0.0, 0.0, 0.0);
        let b = Point3::new(1.0, 0.0, 0.0);
        let short = Point3::new(0.37, 0.0, 0.0);
        let long = Point3::new(0.72, 0.0, 0.0);
        let original = line(a, b);
        let (shell, ids) = manifold_shell(&[(a, short, line(a, short)), (a, long, line(a, long))]);
        let faces = build_edge_face_map(&shell);
        assert_eq!(
            rematch_selected_edge_id(&shell, &original, &HashSet::new(), &faces),
            Some(ids[1]),
        );
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
            rematch_selected_edge_id(&shell, &original, &HashSet::new(), &faces),
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
            rematch_selected_edge_id(&shell, &original, &HashSet::new(), &faces),
            None
        );
    }
}

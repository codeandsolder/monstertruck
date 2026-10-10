//! Classic (0.3.2) boolean pipeline -- the self-contained default backend.
//!
//! This subtree grafts the proven upstream boolean assembly from the published
//! 0.3.2 crate. It is compiled only when the marching SSI is the active backend
//! (default features); a build with an external SSI backend never touches it
//! and runs the current `integrate` pipeline unchanged. The public entry points in
//! `transversal::integrate` cfg-dispatch `and`/`or`/`difference`/
//! `symmetric_difference` here.
//!
//! Graft-boundary adaptations vs. the 0.3.2 source (all documented inline):
//! - The Alternative intersection-curve arm is a raw
//!   `IntersectionCurve<PolylineCurve<Point3>, S, S>` (0.3.2 wrapped it with
//!   parameter polylines the pipeline never consumed).
//! - `altshell_to_shell` converts that arm back to `C` via
//!   `SurfaceCurve::with_boundaries(..).into()` (the conversion the current
//!   `ShapeOpsCurve` bound provides), not `IntersectionCurve::new(..).into()`.
//! - Failures map onto the current `ShapeOpsError` (`EmptyOutputShell` for a
//!   `None` pipeline result, `InvalidOutputShell` for a rejected solid), since
//!   the current error enum has no single `Internal` variant.

mod divide_face;
mod faces_classification;
mod intersection_curve;
mod loops_store;

use super::integrate::{ShapeOpsCurve, ShapeOpsError, ShapeOpsSurface};
use crate::alternative::Alternative;
use monstertruck_geometry::prelude::*;
use monstertruck_meshing::prelude::*;
use monstertruck_topology::*;
use std::{collections::VecDeque, iter};

type ClassicResult<T> = std::result::Result<T, ShapeOpsError>;

type AltCurve<C, S> = Alternative<C, IntersectionCurve<PolylineCurve<Point3>, S, S>>;
type AltCurveShell<C, S> = Shell<Point3, AltCurve<C, S>, S>;

fn weld_coincident_topology<C: ShapeOpsCurve<S>, S: ShapeOpsSurface>(
    shell: &AltCurveShell<C, S>,
) -> AltCurveShell<C, S> {
    type AltEdge<C, S> = Edge<Point3, AltCurve<C, S>>;

    let mut vertices = Vec::<Vertex<Point3>>::new();
    let mut edges = Vec::<AltEdge<C, S>>::new();
    let canonical_vertex = |point: Point3, vertices: &mut Vec<Vertex<Point3>>| {
        if let Some(vertex) = vertices.iter().find(|vertex| vertex.point().near(&point)) {
            vertex.clone()
        } else {
            let vertex = Vertex::new(point);
            vertices.push(vertex.clone());
            vertex
        }
    };
    let curves_coincident = |lhs: &AltEdge<C, S>, rhs: &AltEdge<C, S>| {
        let lhs_curve = lhs.oriented_curve();
        let rhs_curve = rhs.oriented_curve();
        let midpoint = |curve: &AltCurve<C, S>| {
            let (t0, t1) = curve.range_tuple();
            curve.evaluate((t0 + t1) * 0.5)
        };
        if midpoint(&lhs_curve).near(&midpoint(&rhs_curve)) {
            return true;
        }
        let sampled_on_other = |curve: &AltCurve<C, S>, other: &AltCurve<C, S>| {
            let (t0, t1) = curve.range_tuple();
            [0.25, 0.5, 0.75].into_iter().all(|fraction| {
                let point = curve.evaluate(t0 + (t1 - t0) * fraction);
                other
                    .search_nearest_parameter(point, None, 100)
                    .map(|parameter| other.evaluate(parameter).near(&point))
                    .unwrap_or(false)
            })
        };
        sampled_on_other(&lhs_curve, &rhs_curve) && sampled_on_other(&rhs_curve, &lhs_curve)
    };

    let faces = shell
        .iter()
        .map(|face| {
            let wires = face
                .absolute_boundaries()
                .iter()
                .map(|wire| {
                    let rebuilt = wire
                        .iter()
                        .map(|edge| {
                            let front =
                                canonical_vertex(edge.absolute_front().point(), &mut vertices);
                            let back =
                                canonical_vertex(edge.absolute_back().point(), &mut vertices);
                            let mut rebuilt = Edge::new(&front, &back, edge.curve());
                            if !edge.orientation() {
                                rebuilt.invert();
                            }
                            if let Some(canonical) = edges.iter().find(|candidate| {
                                let same_endpoints = candidate.front() == rebuilt.front()
                                    && candidate.back() == rebuilt.back();
                                let opposite_endpoints = candidate.front() == rebuilt.back()
                                    && candidate.back() == rebuilt.front();
                                (same_endpoints || opposite_endpoints)
                                    && curves_coincident(candidate, &rebuilt)
                            }) {
                                if canonical.front() == rebuilt.front() {
                                    canonical.clone()
                                } else {
                                    canonical.inverse()
                                }
                            } else {
                                edges.push(rebuilt.clone());
                                rebuilt
                            }
                        })
                        .collect::<Vec<_>>();
                    Wire::from(rebuilt)
                })
                .collect::<Vec<_>>();
            let mut rebuilt = Face::new_unchecked(wires, face.surface());
            if !face.orientation() {
                rebuilt.invert();
            }
            rebuilt
        })
        .collect::<Vec<_>>();
    Shell::from(faces)
}

/// Convert a polyline leader to an exactly equivalent degree-1 B-spline.
///
/// `PolylineCurve` uses integer segment parameters, so scale the uniform knot
/// vector to the same `[0, segment_count]` domain. This preserves both the
/// guide geometry and parameterization exactly.
fn polyline_leader_bspline(polyline: &PolylineCurve<Point3>) -> Option<BsplineCurve<Point3>> {
    let segment_count = polyline.len().checked_sub(1)?;
    if segment_count == 0 {
        return None;
    }
    let mut knots = KnotVector::uniform_knot(1, segment_count);
    knots.transform(segment_count as f64, 0.0);
    Some(BsplineCurve::new(knots, polyline.iter().copied().collect()))
}

/// Convert an Alternative-curve shell back to the target curve type.
///
/// Surface intersection curves already carry a cleaned polyline leader. Keep
/// that leader exactly as a degree-1 B-spline and wrap it in `SurfaceCurve`,
/// whose evaluation refines points against both supporting surfaces. Re-fitting
/// the refined curve here is both redundant and fragile at tight tolerances.
fn altshell_to_shell<C: ShapeOpsCurve<S>, S: ShapeOpsSurface>(
    altshell: &AltCurveShell<C, S>,
) -> Option<Shell<Point3, C, S>> {
    altshell.try_mapped(
        |p| Some(*p),
        |c| match c {
            Alternative::FirstType(c) => Some(c.clone()),
            Alternative::SecondType(ic) => {
                let surface0 = ic.surface0().clone();
                let surface1 = ic.surface1().clone();
                let bspline = polyline_leader_bspline(ic.leader())?;
                let boundary0: Option<ParameterCurve<BoundaryCurve2D, S>> = None;
                let boundary1: Option<ParameterCurve<BoundaryCurve2D, S>> = None;
                Some(
                    SurfaceCurve::with_boundaries(
                        surface0, surface1, bspline, boundary0, boundary1,
                    )
                    .into(),
                )
            }
        },
        |s| Some(s.clone()),
    )
}

/// Split one pair of shells into `[and_shell, or_shell]` (0.3.2 verbatim): build
/// intersection loops, divide faces, classify each divided face as `and`/`or`,
/// then ray-cast the remaining `unknown` faces against the other shell.
fn classify_inside_with_polyshell(
    poly_shell: &Shell<Point3, PolylineCurve<Point3>, Option<PolygonMesh>>,
    point: Point3,
) -> Option<bool> {
    let offsets = [
        Vector3::new(0.0, 0.0, 0.0),
        Vector3::new(0.613, -0.271, 0.149),
        Vector3::new(-0.347, 0.509, -0.221),
        Vector3::new(0.193, 0.401, 0.577),
        Vector3::new(-0.433, -0.127, 0.389),
    ];
    let (inside_votes, outside_votes) = offsets
        .into_iter()
        .map(|offset| hash::take_one_unit(point + offset))
        .try_fold((0usize, 0usize), |(inside, outside), direction| {
            let count = poly_shell.iter().try_fold(0isize, |count, face| {
                let polygon = face.surface()?;
                Some(count + polygon.signed_crossing_faces(point, direction))
            })?;
            if count == 0 {
                Some((inside, outside + 1))
            } else {
                Some((inside + 1, outside))
            }
        })?;
    Some(inside_votes > outside_votes)
}

fn sample_points_on_face<C, S>(face: &Face<Point3, C, S>) -> Option<Vec<Point3>> {
    let wire = face.absolute_boundaries().first()?;
    let vertices = wire
        .vertex_iter()
        .map(|vertex| vertex.point())
        .collect::<Vec<_>>();
    let (sum, count) = vertices
        .iter()
        .fold((Vector3::zero(), 0usize), |(sum, count), point| {
            (sum + point.to_vec(), count + 1)
        });
    if count == 0 {
        return None;
    }
    let centroid = Point3::from_vec(sum / count as f64);
    Some(
        iter::once(centroid)
            .chain(vertices.into_iter().map(|vertex| centroid.midpoint(vertex)))
            .collect(),
    )
}

fn classify_unknown_face<C, S>(
    poly_shell: &Shell<Point3, PolylineCurve<Point3>, Option<PolygonMesh>>,
    face: &Face<Point3, C, S>,
) -> Option<bool> {
    let points = sample_points_on_face(face)?;
    let (inside, outside) =
        points
            .into_iter()
            .try_fold((0usize, 0usize), |(inside, outside), point| {
                if classify_inside_with_polyshell(poly_shell, point)? {
                    Some((inside + 1, outside))
                } else {
                    Some((inside, outside + 1))
                }
            })?;
    Some(inside >= outside)
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum FaceRelation {
    Inside,
    Outside,
    CoincidentSameFacing,
}

fn classify_face_relation<C, S>(
    poly_shell: &Shell<Point3, PolylineCurve<Point3>, Option<PolygonMesh>>,
    face: &Face<Point3, C, S>,
    tol: f64,
) -> Option<FaceRelation>
where
    S: Clone
        + Invertible
        + ParametricSurface3D
        + SearchParameter<SurfaceParameter, Point = Point3>
        + SearchNearestParameter<SurfaceParameter, Point = Point3>,
{
    let points = sample_points_on_face(face)?;
    let surface = face.oriented_surface();
    let offset = tol.max(1.0e-4) * 8.0;
    let (inside, outside, coincident, total) = points.into_iter().try_fold(
        (0usize, 0usize, 0usize, 0usize),
        |(inside, outside, coincident, total), point| {
            let parameter = surface
                .search_parameter(point, None, 100)
                .or_else(|| surface.search_nearest_parameter(point, None, 100));
            let Some((u, v)) = parameter else {
                let is_inside = classify_inside_with_polyshell(poly_shell, point)?;
                return Some((
                    inside + usize::from(is_inside),
                    outside + usize::from(!is_inside),
                    coincident,
                    total + 1,
                ));
            };
            let normal = surface.normal(u, v);
            let outward_inside =
                classify_inside_with_polyshell(poly_shell, point + normal * offset)?;
            let inward_inside =
                classify_inside_with_polyshell(poly_shell, point - normal * offset)?;
            if outward_inside {
                Some((inside + 1, outside, coincident, total + 1))
            } else if inward_inside {
                Some((inside, outside, coincident + 1, total + 1))
            } else {
                Some((inside, outside + 1, coincident, total + 1))
            }
        },
    )?;
    let relation = if total != 0 && coincident > inside && coincident > outside {
        FaceRelation::CoincidentSameFacing
    } else if inside > outside {
        FaceRelation::Inside
    } else {
        FaceRelation::Outside
    };
    Some(relation)
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum PairMode {
    AndOr,
    Difference,
}

fn process_one_pair_of_shells<C: ShapeOpsCurve<S>, S: ShapeOpsSurface>(
    shell0: &Shell<Point3, C, S>,
    shell1: &Shell<Point3, C, S>,
    tol: f64,
    mode: PairMode,
) -> Option<[Shell<Point3, C, S>; 2]> {
    if tol <= 0.0 {
        return None;
    }
    // The triangulation is only an intersection seed; using the full geometric
    // tolerance here can massively over-refine curved shells and make SSI
    // classification less stable without improving the final curve tolerance.
    // Keep the actual divide/classification tolerance strict below.
    let seed_tol = tol.max(1.0e-3);
    let poly_shell0 = shell0.triangulation(seed_tol);
    let poly_shell1 = shell1.triangulation(seed_tol);
    let altshell0: AltCurveShell<C, S> =
        shell0.mapped(|x| *x, |c| Alternative::FirstType(c.clone()), Clone::clone);
    let altshell1: AltCurveShell<C, S> =
        shell1.mapped(|x| *x, |c| Alternative::FirstType(c.clone()), Clone::clone);
    let loops_store::LoopsStoreQuadruple {
        geom_loops_store0: loops_store0,
        geom_loops_store1: loops_store1,
        coplanar_overlap,
        ..
    } = loops_store::create_loops_stores(
        &altshell0,
        &poly_shell0,
        &altshell1,
        &poly_shell1,
        tol,
        mode == PairMode::AndOr,
    )?;
    let mut cls0 = divide_face::divide_faces(&altshell0, &loops_store0, tol)?;
    let mut cls1 = divide_face::divide_faces(&altshell1, &loops_store1, tol)?;

    if mode == PairMode::Difference {
        cls0.integrate_by_component();
        cls1.integrate_by_component();
        // Face division already classifies pieces adjacent to an intersection:
        // Or on A is outside B and survives A - B, while And on B is
        // inside A and becomes the inward-facing cut boundary. Prefer those
        // exact topological labels: a ray cast from a face lying exactly on
        // the other solid's boundary can choose the wrong side.
        let [_inside0, mut difference_faces0, unknown0] = cls0.and_or_unknown();
        unknown0.into_iter().try_for_each(|face| {
            if !classify_unknown_face(&poly_shell1, &face)? {
                difference_faces0.push(face);
            }
            Some(())
        })?;

        let [mut difference_faces1, _outside1, unknown1] = cls1.and_or_unknown();
        unknown1.into_iter().try_for_each(|face| {
            if classify_unknown_face(&poly_shell0, &face)? {
                difference_faces1.push(face);
            }
            Some(())
        })?;

        for face in difference_faces1.face_iter_mut() {
            // B contributes the inward-facing boundary of A - B.
            face.invert();
        }

        difference_faces0.append(&mut difference_faces1);
        let mut difference_shell = altshell_to_shell(&difference_faces0)?;
        let _ = orient_regular_shell(&mut difference_shell);
        let label_result_is_valid = difference_shell
            .connected_components()
            .iter()
            .all(|component| component.check_solid_boundary().is_ok());
        if label_result_is_valid {
            return Some([difference_shell, Shell::default()]);
        }

        // Axis-touching cuts can leave divider labels locally inconsistent even
        // though geometric classification still reconstructs a valid shell.
        // Fall back as a whole candidate rather than mixing the two schemes
        // face-by-face, so shared-edge closure remains the deciding invariant.
        let [and0, or0, unknown0] = cls0.and_or_unknown();
        let [and1, or1, unknown1] = cls1.and_or_unknown();

        let mut difference_faces0 = AltCurveShell::default();
        and0.into_iter()
            .chain(or0)
            .chain(unknown0)
            .try_for_each(|face| {
                if !classify_unknown_face(&poly_shell1, &face)? {
                    difference_faces0.push(face);
                }
                Some(())
            })?;

        let mut difference_faces1 = AltCurveShell::default();
        and1.into_iter()
            .chain(or1)
            .chain(unknown1)
            .try_for_each(|mut face| {
                if classify_unknown_face(&poly_shell0, &face)? {
                    face.invert();
                    difference_faces1.push(face);
                }
                Some(())
            })?;

        difference_faces0.append(&mut difference_faces1);
        let difference_shell = altshell_to_shell(&difference_faces0)?;
        return Some([difference_shell, Shell::default()]);
    }

    if !coplanar_overlap {
        cls0.integrate_by_component();
        cls1.integrate_by_component();
        let [mut and0, mut or0, unknown0] = cls0.and_or_unknown();
        unknown0.into_iter().try_for_each(|face| {
            if classify_unknown_face(&poly_shell1, &face)? {
                and0.push(face);
            } else {
                or0.push(face);
            }
            Some(())
        })?;
        let [mut and1, mut or1, unknown1] = cls1.and_or_unknown();
        unknown1.into_iter().try_for_each(|face| {
            if classify_unknown_face(&poly_shell0, &face)? {
                and1.push(face);
            } else {
                or1.push(face);
            }
            Some(())
        })?;
        and0.append(&mut and1);
        or0.append(&mut or1);
        let and_shell = altshell_to_shell(&and0)?;
        let or_shell = altshell_to_shell(&or0)?;
        return Some([and_shell, or_shell]);
    }

    // Exact divider labels are stronger evidence than a mesh ray cast,
    // especially for thin features where a point-in-solid vote can be
    // numerically unstable. Preserve labeled pieces and classify only regions
    // untouched by the divider. Coincident unknown ownership remains
    // asymmetric so a shared interface belongs to exactly one OR operand.
    let [mut and0, mut or0, unknown0] = cls0.and_or_unknown();
    unknown0.into_iter().try_for_each(|face| {
        let relation = classify_face_relation(&poly_shell1, &face, tol)?;
        match relation {
            FaceRelation::Inside => and0.push(face),
            FaceRelation::Outside | FaceRelation::CoincidentSameFacing => or0.push(face),
        }
        Some(())
    })?;

    let [mut and1, mut or1, unknown1] = cls1.and_or_unknown();
    unknown1.into_iter().try_for_each(|face| {
        let relation = classify_face_relation(&poly_shell0, &face, tol)?;
        match relation {
            FaceRelation::Inside | FaceRelation::CoincidentSameFacing => and1.push(face),
            FaceRelation::Outside => or1.push(face),
        }
        Some(())
    })?;

    and0.append(&mut and1);
    or0.append(&mut or1);
    and0 = weld_coincident_topology(&and0);
    or0 = weld_coincident_topology(&or0);
    let label_and_shell = altshell_to_shell(&and0)?;
    let label_or_shell = altshell_to_shell(&or0)?;
    let shell_is_valid = |shell: &Shell<Point3, C, S>| {
        shell.is_empty()
            || shell
                .connected_components()
                .iter()
                .all(|component| component.check_solid_boundary().is_ok())
    };
    if !label_or_shell.is_empty()
        && shell_is_valid(&label_and_shell)
        && shell_is_valid(&label_or_shell)
    {
        return Some([label_and_shell, label_or_shell]);
    }

    // Some coplanar configurations produce locally useful divider labels that
    // still cannot form complete output shells (for example overlapping boxes
    // with coplanar side planes). Fall back to geometric classification of all
    // divided pieces only after the topology-first candidate fails closure.
    let [and0_labeled, or0_labeled, unknown0] = cls0.and_or_unknown();
    let mut and0 = AltCurveShell::default();
    let mut or0 = AltCurveShell::default();
    and0_labeled
        .into_iter()
        .chain(or0_labeled)
        .chain(unknown0)
        .try_for_each(|face| {
            let relation = classify_face_relation(&poly_shell1, &face, tol)?;
            match relation {
                FaceRelation::Inside => and0.push(face),
                FaceRelation::Outside | FaceRelation::CoincidentSameFacing => or0.push(face),
            }
            Some(())
        })?;

    let [and1_labeled, or1_labeled, unknown1] = cls1.and_or_unknown();
    let mut and1 = AltCurveShell::default();
    let mut or1 = AltCurveShell::default();
    and1_labeled
        .into_iter()
        .chain(or1_labeled)
        .chain(unknown1)
        .try_for_each(|face| {
            let relation = classify_face_relation(&poly_shell0, &face, tol)?;
            match relation {
                FaceRelation::Inside | FaceRelation::CoincidentSameFacing => and1.push(face),
                FaceRelation::Outside => or1.push(face),
            }
            Some(())
        })?;

    and0.append(&mut and1);
    or0.append(&mut or1);
    and0 = weld_coincident_topology(&and0);
    or0 = weld_coincident_topology(&or0);
    let and_shell = altshell_to_shell(&and0)?;
    let or_shell = altshell_to_shell(&or0)?;
    Some([and_shell, or_shell])
}

/// Re-orient a regular manifold shell from shared-edge constraints.
fn orient_regular_shell<P, C, S>(shell: &mut Shell<P, C, S>) -> bool {
    if shell.shell_condition() != ShellCondition::Regular {
        return true;
    }

    let mut adjacency = vec![Vec::<(usize, bool)>::new(); shell.len()];
    for left in 0..shell.len() {
        for right in left + 1..shell.len() {
            let mut relation = None;
            for left_edge in shell[left].edge_iter() {
                for right_edge in shell[right].edge_iter() {
                    if left_edge.id() != right_edge.id() {
                        continue;
                    }
                    let must_differ = left_edge.orientation() == right_edge.orientation();
                    match relation {
                        Some(existing) if existing != must_differ => return false,
                        Some(_) => {}
                        None => relation = Some(must_differ),
                    }
                }
            }
            if let Some(must_differ) = relation {
                adjacency[left].push((right, must_differ));
                adjacency[right].push((left, must_differ));
            }
        }
    }

    let mut flips = vec![None; shell.len()];
    for root in 0..shell.len() {
        if flips[root].is_some() {
            continue;
        }
        flips[root] = Some(false);
        let mut queue = VecDeque::from([root]);
        while let Some(face) = queue.pop_front() {
            let Some(flip) = flips[face] else {
                return false;
            };
            for &(neighbor, must_differ) in &adjacency[face] {
                let required = flip ^ must_differ;
                match flips[neighbor] {
                    Some(existing) if existing != required => return false,
                    Some(_) => {}
                    None => {
                        flips[neighbor] = Some(required);
                        queue.push_back(neighbor);
                    }
                }
            }
        }
    }

    for (face, flip) in shell.iter_mut().zip(flips) {
        if flip == Some(true) {
            face.invert();
        }
    }
    shell.shell_condition() != ShellCondition::Regular
}

fn finalize<C: ShapeOpsCurve<S>, S: ShapeOpsSurface>(
    operation: &'static str,
    mut shell: Shell<Point3, C, S>,
) -> ClassicResult<Solid<Point3, C, S>> {
    let _ = orient_regular_shell(&mut shell);
    let boundaries = shell.connected_components();
    Solid::try_new(boundaries)
        .map_err(|source| ShapeOpsError::InvalidOutputShell { operation, source })
}

/// AND operation between two solids (classic backend).
fn and_shells<C: ShapeOpsCurve<S>, S: ShapeOpsSurface>(
    solid0: &Solid<Point3, C, S>,
    solid1: &Solid<Point3, C, S>,
    tol: f64,
) -> ClassicResult<Solid<Point3, C, S>> {
    let operation = "and";
    let pair = |a: &Shell<Point3, C, S>, b: &Shell<Point3, C, S>| {
        process_one_pair_of_shells(a, b, tol, PairMode::AndOr)
            .ok_or(ShapeOpsError::EmptyOutputShell { operation })
    };
    let mut iter0 = solid0.boundaries().iter();
    let mut iter1 = solid1.boundaries().iter();
    let shell0 = iter0.next().unwrap();
    let shell1 = iter1.next().unwrap();
    let [mut and_shell, _] = pair(shell0, shell1)?;
    for shell in iter0 {
        let [res, _] = pair(&and_shell, shell)?;
        and_shell = res;
    }
    for shell in iter1 {
        let [res, _] = pair(&and_shell, shell)?;
        and_shell = res;
    }
    finalize(operation, and_shell)
}

/// AND operation between two solids (classic backend).
pub(crate) fn and<C: ShapeOpsCurve<S>, S: ShapeOpsSurface>(
    solid0: &Solid<Point3, C, S>,
    solid1: &Solid<Point3, C, S>,
    tol: f64,
) -> ClassicResult<Solid<Point3, C, S>> {
    and_shells(solid0, solid1, tol)
}

/// OR operation between two solids (classic backend).
pub(crate) fn or<C: ShapeOpsCurve<S>, S: ShapeOpsSurface>(
    solid0: &Solid<Point3, C, S>,
    solid1: &Solid<Point3, C, S>,
    tol: f64,
) -> ClassicResult<Solid<Point3, C, S>> {
    let operation = "or";
    let pair = |a: &Shell<Point3, C, S>, b: &Shell<Point3, C, S>| {
        process_one_pair_of_shells(a, b, tol, PairMode::AndOr)
            .ok_or(ShapeOpsError::EmptyOutputShell { operation })
    };
    let mut iter0 = solid0.boundaries().iter();
    let mut iter1 = solid1.boundaries().iter();
    let shell0 = iter0.next().unwrap();
    let shell1 = iter1.next().unwrap();
    let [_, mut or_shell] = pair(shell0, shell1)?;
    for shell in iter0 {
        let [_, res] = pair(&or_shell, shell)?;
        or_shell = res;
    }
    for shell in iter1 {
        let [_, res] = pair(&or_shell, shell)?;
        or_shell = res;
    }
    finalize(operation, or_shell)
}

/// Difference: the region inside `solid0` but outside `solid1` (classic backend).
pub(crate) fn difference<C: ShapeOpsCurve<S>, S: ShapeOpsSurface>(
    solid0: &Solid<Point3, C, S>,
    solid1: &Solid<Point3, C, S>,
    tol: f64,
) -> ClassicResult<Solid<Point3, C, S>> {
    let operation = "difference";
    let difference_pair = |a: &Shell<Point3, C, S>, b: &Shell<Point3, C, S>| {
        process_one_pair_of_shells(a, b, tol, PairMode::Difference)
            .ok_or(ShapeOpsError::EmptyOutputShell { operation })
    };
    let and_pair = |a: &Shell<Point3, C, S>, b: &Shell<Point3, C, S>| {
        process_one_pair_of_shells(a, b, tol, PairMode::AndOr)
            .ok_or(ShapeOpsError::EmptyOutputShell { operation })
    };

    let mut iter0 = solid0.boundaries().iter();
    let mut iter1 = solid1.boundaries().iter();
    let shell0 = iter0.next().unwrap();
    let shell1 = iter1.next().unwrap();
    let [mut difference_shell, _] = difference_pair(shell0, shell1)?;

    for shell in iter0 {
        let [res, _] = and_pair(&difference_shell, shell)?;
        difference_shell = res;
    }
    for shell in iter1 {
        let [res, _] = difference_pair(&difference_shell, shell)?;
        difference_shell = res;
    }
    finalize(operation, difference_shell)
}

/// Symmetric difference (XOR): the region inside exactly one solid (classic backend).
pub(crate) fn symmetric_difference<C: ShapeOpsCurve<S>, S: ShapeOpsSurface>(
    solid0: &Solid<Point3, C, S>,
    solid1: &Solid<Point3, C, S>,
    tol: f64,
) -> ClassicResult<Solid<Point3, C, S>> {
    let d0 = difference(solid0, solid1, tol)?;
    let d1 = difference(solid1, solid0, tol)?;
    or(&d0, &d1, tol)
}

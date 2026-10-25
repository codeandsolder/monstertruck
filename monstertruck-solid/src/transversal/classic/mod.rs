//! Self-contained classic Boolean pipeline, descended from the published 0.3.2 backend.
//!
//! The historical backend remains the algorithmic base, but this subtree now has
//! explicit staged analysis/classification, typed failure boundaries, canonical
//! topology repair, and stricter geometric evidence semantics. It is compiled only
//! when the marching SSI is the active backend
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
//! - Internal stages return typed `ShapeOpsError` variants; only an assembled
//!   but topologically invalid candidate may select a documented fallback path.

mod classification;
mod conversion;
mod curve_equivalence;
mod divide_face;
mod faces_classification;
mod intersection_curve;
mod loops_store;
mod tolerance;
mod weld;

use super::integrate::{
    InvalidOutputShellConditionData, ShapeOpsCurve, ShapeOpsError, ShapeOpsSurface,
};
use crate::alternative::Alternative;
#[cfg(test)]
use classification::sample_points_on_face;
use classification::{
    FaceRelation, classify_face_relation, classify_inside_with_polyshell, classify_unknown_face,
    planar_partition_same_facing,
};
use conversion::altshell_to_shell;
use monstertruck_geometry::prelude::*;
use monstertruck_meshing::prelude::*;
use monstertruck_topology::*;
use std::collections::{HashMap, VecDeque};
use tolerance::BooleanTolerance;
use weld::weld_coincident_regions;

type ClassicResult<T> = std::result::Result<T, ShapeOpsError>;

type AltCurve<C, S> = Alternative<C, IntersectionCurve<PolylineCurve<Point3>, S, S>>;
type AltCurveShell<C, S> = Shell<Point3, AltCurve<C, S>, S>;
type AltRegion<C, S> = faces_classification::BooleanRegion<Point3, AltCurve<C, S>, S>;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum PairMode {
    AndOr,
    Difference,
}

type PolyShell = Shell<Point3, PolylineCurve<Point3>, Option<PolygonMesh>>;

#[derive(Clone, Debug)]
struct FullInterfaceAdjacency {
    shell0: Vec<usize>,
    shell1: Vec<usize>,
}

impl FullInterfaceAdjacency {
    fn contains(&self, shell_index: usize, source_face: usize) -> bool {
        match shell_index {
            0 => self.shell0.contains(&source_face),
            1 => self.shell1.contains(&source_face),
            _ => false,
        }
    }
}

#[derive(Clone, Debug)]
struct PairAnalysis<C, S> {
    tolerance: f64,
    poly_shell0: PolyShell,
    poly_shell1: PolyShell,
    cls0: faces_classification::FacesClassification<Point3, AltCurve<C, S>, S>,
    cls1: faces_classification::FacesClassification<Point3, AltCurve<C, S>, S>,
    coplanar_overlap: bool,
    full_interface_adjacency: Option<FullInterfaceAdjacency>,
}

fn analyze_pair<C: ShapeOpsCurve<S>, S: ShapeOpsSurface>(
    shell0: &Shell<Point3, C, S>,
    shell1: &Shell<Point3, C, S>,
    tol: f64,
    mode: PairMode,
) -> ClassicResult<PairAnalysis<C, S>> {
    if tol < TOLERANCE {
        return Err(ShapeOpsError::InvalidTolerance);
    }
    let seed_tol = BooleanTolerance::new(tol).mesh_seed;
    let poly_shell0 = shell0.triangulation(seed_tol);
    let poly_shell1 = shell1.triangulation(seed_tol);
    let altshell0: AltCurveShell<C, S> =
        shell0.mapped(|x| *x, |c| Alternative::FirstType(c.clone()), Clone::clone);
    let altshell1: AltCurveShell<C, S> =
        shell1.mapped(|x| *x, |c| Alternative::FirstType(c.clone()), Clone::clone);
    let loops_store::InterfaceLoops {
        geom: [loops_store0, loops_store1],
        coplanar_overlap,
        full_interface_adjacency,
    } = loops_store::create_loops_stores(
        &altshell0,
        &poly_shell0,
        &altshell1,
        &poly_shell1,
        tol,
        mode == PairMode::AndOr,
    )
    .ok_or(ShapeOpsError::InterfaceAnalysisFailed)?;
    let full_interface_adjacency = full_interface_adjacency.map(|matched| FullInterfaceAdjacency {
        shell0: matched.lhs_faces,
        shell1: matched.rhs_faces,
    });
    let cls0 = divide_face::divide_faces(&altshell0, &loops_store0, tol)
        .ok_or(ShapeOpsError::DivideFacesFailed { shell_index: 0 })?;
    let cls1 = divide_face::divide_faces(&altshell1, &loops_store1, tol)
        .ok_or(ShapeOpsError::DivideFacesFailed { shell_index: 1 })?;
    Ok(PairAnalysis {
        tolerance: tol,
        poly_shell0,
        poly_shell1,
        cls0,
        cls1,
        coplanar_overlap,
        full_interface_adjacency,
    })
}

fn classify_unknown<C, S>(
    poly_shell: &PolyShell,
    region: &faces_classification::BooleanRegion<Point3, C, S>,
    tol: f64,
    shell_index: usize,
) -> ClassicResult<bool>
where
    C: PolylineableCurve + ParameterDivision1D<Point = Point3> + Clone,
    S: MeshableSurface
        + Clone
        + ParametricSurface3D
        + SearchParameter<SurfaceParameter, Point = Point3>
        + SearchNearestParameter<SurfaceParameter, Point = Point3>,
{
    if let Some(face) = region.try_face() {
        return classify_unknown_face(poly_shell, &face, tol)
            .ok_or(ShapeOpsError::UnknownClassificationFailed { shell_index });
    }
    let points =
        divide_face::sample_points_on_region(region.absolute_boundaries(), &region.surface(), tol)
            .ok_or(ShapeOpsError::UnknownClassificationFailed { shell_index })?;
    let (inside, outside) = points
        .into_iter()
        .try_fold((0usize, 0usize), |(inside, outside), point| {
            let is_inside = classify_inside_with_polyshell(poly_shell, point)?;
            Some((
                inside + usize::from(is_inside),
                outside + usize::from(!is_inside),
            ))
        })
        .ok_or(ShapeOpsError::UnknownClassificationFailed { shell_index })?;
    Ok(inside >= outside)
}

fn classify_relation<C, S>(
    poly_shell: &PolyShell,
    region: &faces_classification::BooleanRegion<Point3, C, S>,
    tol: f64,
    shell_index: usize,
) -> ClassicResult<FaceRelation>
where
    C: PolylineableCurve + ParameterDivision1D<Point = Point3> + Clone,
    S: MeshableSurface
        + Clone
        + Invertible
        + ParametricSurface3D
        + SearchParameter<SurfaceParameter, Point = Point3>
        + SearchNearestParameter<SurfaceParameter, Point = Point3>,
{
    if let Some(face) = region.try_face() {
        return classify_face_relation(poly_shell, &face, tol)
            .ok_or(ShapeOpsError::UnknownClassificationFailed { shell_index });
    }
    let points =
        divide_face::sample_points_on_region(region.absolute_boundaries(), &region.surface(), tol)
            .ok_or(ShapeOpsError::UnknownClassificationFailed { shell_index })?;
    let mut surface = region.surface();
    if !region.orientation() {
        surface.invert();
    }
    let offset = BooleanTolerance::new(tol).classifier_offset;
    let (inside, outside, coincident, total) = points
        .into_iter()
        .try_fold(
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
        )
        .ok_or(ShapeOpsError::UnknownClassificationFailed { shell_index })?;
    Ok(
        if total != 0 && coincident > inside && coincident > outside {
            FaceRelation::CoincidentSameFacing
        } else if inside > outside {
            FaceRelation::Inside
        } else {
            FaceRelation::Outside
        },
    )
}

fn convert_alt_shell<C: ShapeOpsCurve<S>, S: ShapeOpsSurface>(
    shell: &AltCurveShell<C, S>,
    operation: &'static str,
) -> ClassicResult<Shell<Point3, C, S>> {
    altshell_to_shell(shell).ok_or(ShapeOpsError::AltShellConversionFailed { operation })
}

fn weld_and_convert<C: ShapeOpsCurve<S>, S: ShapeOpsSurface>(
    regions: &[AltRegion<C, S>],
    operation: &'static str,
) -> ClassicResult<Shell<Point3, C, S>> {
    let shell = weld_coincident_regions(regions)
        .ok_or(ShapeOpsError::RegionMaterializationFailed { operation })?;
    convert_alt_shell(&shell, operation)
}

fn shell_is_valid<P, C, S>(shell: &Shell<P, C, S>) -> bool {
    shell.is_empty()
        || shell
            .connected_components()
            .iter()
            .all(|component| component.check_solid_boundary().is_ok())
}

fn process_difference_analysis<C: ShapeOpsCurve<S>, S: ShapeOpsSurface>(
    mut analysis: PairAnalysis<C, S>,
) -> ClassicResult<[Shell<Point3, C, S>; 2]> {
    let tol = analysis.tolerance;
    analysis.cls0.integrate_by_component();
    analysis.cls1.integrate_by_component();

    let [_inside0, mut difference0, unknown0] = analysis.cls0.and_or_unknown();
    for region in unknown0 {
        if !classify_unknown(&analysis.poly_shell1, &region, tol, 0)? {
            difference0.push(region);
        }
    }
    let [mut difference1, _outside1, unknown1] = analysis.cls1.and_or_unknown();
    for region in unknown1 {
        if classify_unknown(&analysis.poly_shell0, &region, tol, 1)? {
            difference1.push(region);
        }
    }
    for region in &mut difference1 {
        region.invert();
    }
    difference0.append(&mut difference1);
    let mut label_shell = weld_and_convert(&difference0, "difference")?;
    let _ = orient_regular_shell(&mut label_shell);
    if label_shell
        .connected_components()
        .iter()
        .all(|component| component.check_solid_boundary().is_ok())
    {
        return Ok([label_shell, Shell::default()]);
    }

    let [and0, or0, unknown0] = analysis.cls0.and_or_unknown();
    let [and1, or1, unknown1] = analysis.cls1.and_or_unknown();
    let mut fallback = Vec::new();
    for region in and0.into_iter().chain(or0).chain(unknown0) {
        if !classify_unknown(&analysis.poly_shell1, &region, tol, 0)? {
            fallback.push(region);
        }
    }
    for mut region in and1.into_iter().chain(or1).chain(unknown1) {
        if classify_unknown(&analysis.poly_shell0, &region, tol, 1)? {
            region.invert();
            fallback.push(region);
        }
    }
    Ok([weld_and_convert(&fallback, "difference")?, Shell::default()])
}

fn process_non_coplanar_analysis<C: ShapeOpsCurve<S>, S: ShapeOpsSurface>(
    mut analysis: PairAnalysis<C, S>,
) -> ClassicResult<[Shell<Point3, C, S>; 2]> {
    let tol = analysis.tolerance;
    analysis.cls0.integrate_by_component();
    analysis.cls1.integrate_by_component();
    let [mut and0, mut or0, unknown0] = analysis.cls0.and_or_unknown();
    for region in unknown0 {
        if classify_unknown(&analysis.poly_shell1, &region, tol, 0)? {
            and0.push(region);
        } else {
            or0.push(region);
        }
    }
    let [mut and1, mut or1, unknown1] = analysis.cls1.and_or_unknown();
    for region in unknown1 {
        if classify_unknown(&analysis.poly_shell0, &region, tol, 1)? {
            and1.push(region);
        } else {
            or1.push(region);
        }
    }
    and0.append(&mut and1);
    or0.append(&mut or1);
    Ok([
        weld_and_convert(&and0, "and")?,
        weld_and_convert(&or0, "or")?,
    ])
}

fn partition_evidence<C: ShapeOpsCurve<S>, S: ShapeOpsSurface>(
    regions: &[AltRegion<C, S>],
    mesh_tol: f64,
) -> (Vec<AltRegion<C, S>>, PolyShell, HashMap<usize, PolygonMesh>) {
    let mut valid_regions = Vec::new();
    let mut poly_faces = Vec::new();
    let mut by_key = HashMap::new();
    for region in regions {
        let Some(face) = region.try_face() else {
            continue;
        };
        let triangulated = Shell::from(vec![face]).triangulation(mesh_tol);
        let Some(poly_face) = triangulated.first() else {
            continue;
        };
        let Some(polygon) = poly_face.surface() else {
            continue;
        };
        by_key.insert(region.key(), polygon.clone());
        valid_regions.push(region.clone());
        poly_faces.push(poly_face.clone());
    }
    (valid_regions, Shell::from(poly_faces), by_key)
}

fn process_coplanar_analysis<C: ShapeOpsCurve<S>, S: ShapeOpsSurface>(
    mut analysis: PairAnalysis<C, S>,
) -> ClassicResult<[Shell<Point3, C, S>; 2]> {
    let tol = analysis.tolerance;
    analysis.cls0.integrate_by_component();
    analysis.cls1.integrate_by_component();

    let [mut and0, mut or0, unknown0] = analysis.cls0.and_or_unknown();
    for region in unknown0 {
        if analysis
            .full_interface_adjacency
            .as_ref()
            .is_some_and(|interface| interface.contains(0, region.source_face()))
        {
            continue;
        }
        match classify_relation(&analysis.poly_shell1, &region, tol, 0)? {
            FaceRelation::Inside => and0.push(region),
            FaceRelation::Outside | FaceRelation::CoincidentSameFacing => or0.push(region),
        }
    }
    let [mut and1, mut or1, unknown1] = analysis.cls1.and_or_unknown();
    for region in unknown1 {
        if analysis
            .full_interface_adjacency
            .as_ref()
            .is_some_and(|interface| interface.contains(1, region.source_face()))
        {
            continue;
        }
        match classify_relation(&analysis.poly_shell0, &region, tol, 1)? {
            FaceRelation::Inside | FaceRelation::CoincidentSameFacing => and1.push(region),
            FaceRelation::Outside => or1.push(region),
        }
    }
    and0.append(&mut and1);
    or0.append(&mut or1);
    let label_and_shell = weld_and_convert(&and0, "and")?;
    let label_or_shell = weld_and_convert(&or0, "or")?;
    if !label_or_shell.is_empty()
        && shell_is_valid(&label_and_shell)
        && shell_is_valid(&label_or_shell)
    {
        return Ok([label_and_shell, label_or_shell]);
    }

    let all0 = analysis.cls0.all_regions();
    let all1 = analysis.cls1.all_regions();
    let partition_tolerance = BooleanTolerance::new(tol);
    let partition_geom_tol = partition_tolerance.length;
    let partition_mesh_tol = partition_tolerance.mesh_seed;
    let (valid0, poly0, candidate0) = partition_evidence(&all0, partition_mesh_tol);
    let (valid1, poly1, candidate1) = partition_evidence(&all1, partition_mesh_tol);

    let [and0_labeled, or0_labeled, unknown0] = analysis.cls0.and_or_unknown();
    let mut and0 = Vec::new();
    let mut or0 = Vec::new();
    for region in and0_labeled.into_iter().chain(or0_labeled).chain(unknown0) {
        let coincident = candidate0.get(&region.key()).is_some_and(|candidate| {
            planar_partition_same_facing(&region, candidate, &valid1, &poly1, partition_geom_tol)
                == Some(true)
        });
        let relation = if coincident {
            FaceRelation::CoincidentSameFacing
        } else {
            classify_relation(&analysis.poly_shell1, &region, tol, 0)?
        };
        match relation {
            FaceRelation::Inside => and0.push(region),
            FaceRelation::Outside | FaceRelation::CoincidentSameFacing => or0.push(region),
        }
    }

    let [and1_labeled, or1_labeled, unknown1] = analysis.cls1.and_or_unknown();
    let mut and1 = Vec::new();
    let mut or1 = Vec::new();
    for region in and1_labeled.into_iter().chain(or1_labeled).chain(unknown1) {
        let coincident = candidate1.get(&region.key()).is_some_and(|candidate| {
            planar_partition_same_facing(&region, candidate, &valid0, &poly0, partition_geom_tol)
                == Some(true)
        });
        let relation = if coincident {
            FaceRelation::CoincidentSameFacing
        } else {
            classify_relation(&analysis.poly_shell0, &region, tol, 1)?
        };
        match relation {
            FaceRelation::Inside | FaceRelation::CoincidentSameFacing => and1.push(region),
            FaceRelation::Outside => or1.push(region),
        }
    }
    and0.append(&mut and1);
    or0.append(&mut or1);
    Ok([
        weld_and_convert(&and0, "and")?,
        weld_and_convert(&or0, "or")?,
    ])
}

fn process_one_pair_of_shells<C: ShapeOpsCurve<S>, S: ShapeOpsSurface>(
    shell0: &Shell<Point3, C, S>,
    shell1: &Shell<Point3, C, S>,
    tol: f64,
    mode: PairMode,
) -> ClassicResult<[Shell<Point3, C, S>; 2]> {
    let analysis = analyze_pair(shell0, shell1, tol, mode)?;
    match mode {
        PairMode::Difference => process_difference_analysis(analysis),
        PairMode::AndOr if analysis.coplanar_overlap => process_coplanar_analysis(analysis),
        PairMode::AndOr => process_non_coplanar_analysis(analysis),
    }
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
    for (shell_index, boundary) in boundaries.iter().enumerate() {
        if boundary.check_solid_boundary().is_ok() {
            continue;
        }
        let open_boundaries = boundary.extract_boundaries();
        let first_boundary = open_boundaries.first();
        let singular = boundary.singular_vertices();
        let mut edge_uses = Vec::<(EdgeId<C>, Edge<Point3, C>, usize)>::new();
        for edge in boundary.edge_iter() {
            if let Some((_, _, uses)) = edge_uses.iter_mut().find(|(id, _, _)| *id == edge.id()) {
                *uses += 1;
            } else {
                edge_uses.push((edge.id(), edge.clone(), 1));
            }
        }
        let overused = edge_uses
            .iter()
            .filter(|(_, _, uses)| *uses > 2)
            .collect::<Vec<_>>();
        let first_overused = overused.first();
        return Err(ShapeOpsError::InvalidOutputShellCondition(Box::new(
            InvalidOutputShellConditionData {
                operation,
                shell_index,
                face_count: boundary.len(),
                empty: boundary.is_empty(),
                connected: boundary.is_connected(),
                condition: boundary.shell_condition(),
                boundary_loops: open_boundaries.len(),
                first_boundary_len: first_boundary.map(|wire| wire.len()),
                first_boundary_front: first_boundary
                    .and_then(Wire::front_vertex)
                    .map(Vertex::point),
                first_boundary_back: first_boundary
                    .and_then(Wire::back_vertex)
                    .map(Vertex::point),
                overused_edges: overused.len(),
                first_overused_uses: first_overused.map(|(_, _, uses)| *uses),
                first_overused_front: first_overused.map(|(_, edge, _)| edge.front().point()),
                first_overused_back: first_overused.map(|(_, edge, _)| edge.back().point()),
                singular_vertices: singular.len(),
                first_singular: singular.first().map(Vertex::point),
            },
        )));
    }
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
    };
    let and_pair = |a: &Shell<Point3, C, S>, b: &Shell<Point3, C, S>| {
        process_one_pair_of_shells(a, b, tol, PairMode::AndOr)
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

#[cfg(test)]
mod structural_classification_tests {
    use super::*;
    use monstertruck_modeling::Curve;

    fn tetrahedron_mesh(reverse: bool) -> PolygonMesh {
        let positions = vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
            Point3::new(0.0, 0.0, 1.0),
        ];
        let mut faces = vec![[0usize, 2, 1], [0, 1, 3], [0, 3, 2], [1, 2, 3]];
        if reverse {
            faces.iter_mut().for_each(|face| face.swap(1, 2));
        }
        PolygonMesh::new(
            StandardAttributes {
                positions,
                ..Default::default()
            },
            Faces::from_iter(faces),
        )
    }

    #[test]
    fn signed_winding_does_not_turn_negative_shell_into_inside() {
        let reversed = tetrahedron_mesh(true);
        let shell: Shell<Point3, PolylineCurve<Point3>, Option<PolygonMesh>> =
            Shell::from(vec![Face::new_unchecked(Vec::new(), Some(reversed))]);
        assert_eq!(
            classify_inside_with_polyshell(&shell, Point3::new(0.1, 0.1, 0.1)),
            Some(false)
        );
    }

    #[test]
    fn face_classification_samples_stay_inside_concave_trim() {
        let points = [
            (0.0, 0.0),
            (3.0, 0.0),
            (3.0, 3.0),
            (2.0, 3.0),
            (2.0, 1.0),
            (1.0, 1.0),
            (1.0, 3.0),
            (0.0, 3.0),
        ];
        let vertices = points
            .iter()
            .map(|&(x, y)| Vertex::new(Point3::new(x, y, 0.0)))
            .collect::<Vec<_>>();
        let edges = (0..vertices.len())
            .map(|index| {
                let next = (index + 1) % vertices.len();
                Edge::new(
                    &vertices[index],
                    &vertices[next],
                    Curve::from(Line(vertices[index].point(), vertices[next].point())),
                )
            })
            .collect::<Vec<_>>();
        let face = Face::new_unchecked(
            vec![Wire::from(edges)],
            monstertruck_modeling::Surface::Plane(Plane::xy()),
        );
        let boundary = PolylineCurve(
            points
                .iter()
                .map(|&(x, y)| Point2::new(x, y))
                .chain(std::iter::once(Point2::new(points[0].0, points[0].1)))
                .collect(),
        );
        let samples = sample_points_on_face(&face, 1.0e-6).unwrap_or_default();
        assert!(!samples.is_empty());
        assert!(samples.len() <= 9);
        assert!(
            samples
                .iter()
                .all(|point| boundary.include(Point2::new(point.x, point.y)))
        );
    }

    #[test]
    fn face_classification_sampling_keeps_regions_thinner_than_mesh_seed() {
        let points = [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0e-4), (0.0, 1.0e-4)];
        let vertices = points
            .iter()
            .map(|&(x, y)| Vertex::new(Point3::new(x, y, 0.0)))
            .collect::<Vec<_>>();
        let edges = (0..vertices.len())
            .map(|index| {
                let next = (index + 1) % vertices.len();
                Edge::new(
                    &vertices[index],
                    &vertices[next],
                    Curve::from(Line(vertices[index].point(), vertices[next].point())),
                )
            })
            .collect::<Vec<_>>();
        let face = Face::new_unchecked(
            vec![Wire::from(edges)],
            monstertruck_modeling::Surface::Plane(Plane::xy()),
        );

        let samples = sample_points_on_face(&face, 1.0e-6).unwrap_or_default();
        assert!(!samples.is_empty());
        assert!(
            samples
                .iter()
                .all(|point| point.x > 0.0 && point.x < 1.0 && point.y > 0.0 && point.y < 1.0e-4)
        );
    }

    fn scaled_tetrahedron_mesh(origin: Point3, scale: f64, reverse: bool) -> PolygonMesh {
        let positions = vec![
            origin,
            origin + Vector3::new(scale, 0.0, 0.0),
            origin + Vector3::new(0.0, scale, 0.0),
            origin + Vector3::new(0.0, 0.0, scale),
        ];
        let mut faces = vec![[0usize, 2, 1], [0, 1, 3], [0, 3, 2], [1, 2, 3]];
        if reverse {
            faces.iter_mut().for_each(|face| face.swap(1, 2));
        }
        PolygonMesh::new(
            StandardAttributes {
                positions,
                ..Default::default()
            },
            Faces::from_iter(faces),
        )
    }

    #[test]
    fn signed_winding_respects_cavity_orientation() {
        let outer = scaled_tetrahedron_mesh(Point3::origin(), 4.0, false);
        let cavity = scaled_tetrahedron_mesh(Point3::new(0.5, 0.5, 0.5), 1.0, true);
        let shell: Shell<Point3, PolylineCurve<Point3>, Option<PolygonMesh>> = Shell::from(vec![
            Face::new_unchecked(Vec::new(), Some(outer)),
            Face::new_unchecked(Vec::new(), Some(cavity)),
        ]);

        assert_eq!(
            classify_inside_with_polyshell(&shell, Point3::new(0.25, 0.25, 0.25)),
            Some(true)
        );
        assert_eq!(
            classify_inside_with_polyshell(&shell, Point3::new(0.6, 0.6, 0.6)),
            Some(false)
        );
    }

    #[test]
    fn face_classification_samples_avoid_holes() {
        fn wire(points: &[(f64, f64)]) -> Wire<Point3, Curve> {
            let vertices = points
                .iter()
                .map(|&(x, y)| Vertex::new(Point3::new(x, y, 0.0)))
                .collect::<Vec<_>>();
            (0..vertices.len())
                .map(|index| {
                    let next = (index + 1) % vertices.len();
                    Edge::new(
                        &vertices[index],
                        &vertices[next],
                        Curve::from(Line(vertices[index].point(), vertices[next].point())),
                    )
                })
                .collect::<Vec<_>>()
                .into()
        }

        let outer = [(0.0, 0.0), (4.0, 0.0), (4.0, 4.0), (0.0, 4.0)];
        let hole = [(1.0, 1.0), (1.0, 3.0), (3.0, 3.0), (3.0, 1.0)];
        let face = Face::new_unchecked(
            vec![wire(&outer), wire(&hole)],
            monstertruck_modeling::Surface::Plane(Plane::xy()),
        );

        let samples = sample_points_on_face(&face, 1.0e-6).unwrap_or_default();
        assert!(!samples.is_empty());
        assert!(samples.iter().all(|point| {
            let x = point.x;
            let y = point.y;
            x <= 1.0 + TOLERANCE
                || x >= 3.0 - TOLERANCE
                || y <= 1.0 + TOLERANCE
                || y >= 3.0 - TOLERANCE
        }));
    }
}

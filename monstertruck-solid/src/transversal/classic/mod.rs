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
mod curve_refinement;
mod divide_face;
mod faces_classification;
mod intersection_curve;
mod loops_store;
pub(super) mod planar_contact;
mod planar_union;
#[cfg(test)]
mod radial_slot_tests;
mod tolerance;
mod weld;

use super::integrate::{
    InvalidOutputShellConditionData, ShapeOpsCurve, ShapeOpsError, ShapeOpsSurface,
    validate_tolerance,
};
use crate::alternative::Alternative;
#[cfg(test)]
use classification::sample_points_on_face;
use classification::{
    FaceRelation, WindingClassifier, classify_face_relation, classify_points_inside,
    classify_points_relation, classify_unknown_face, planar_partition_same_facing,
};
use conversion::altshell_to_shell;
use monstertruck_geometry::prelude::*;
use monstertruck_meshing::prelude::*;
use monstertruck_topology::{Edge, EdgeId, Face, Shell, Solid, Vertex, VertexId, Wire};
use std::collections::{HashMap, VecDeque};
use tolerance::BooleanTolerance;
use weld::weld_coincident_regions;

type ClassicResult<T> = std::result::Result<T, ShapeOpsError>;

type AltCurve<C, S> = Alternative<C, IntersectionCurve<PolylineCurve<Point3>, S, S>>;
type AltCurveShell<C, S> = Shell<Point3, AltCurve<C, S>, S>;
type AltRegion<C, S> = faces_classification::BooleanRegion<Point3, AltCurve<C, S>, S>;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum PairMode {
    AndOr(PartitionOutput),
    Difference,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum PartitionOutput {
    Intersection,
    Union,
}

impl PartitionOutput {
    const fn operation(self) -> &'static str {
        match self {
            Self::Intersection => "and",
            Self::Union => "or",
        }
    }

    const fn includes_inside(self) -> bool { matches!(self, Self::Intersection) }

    fn select<T>(self, [intersection, union]: [T; 2]) -> T {
        match self {
            Self::Intersection => intersection,
            Self::Union => union,
        }
    }
}

type PolyShell = Shell<Point3, PolylineCurve<Point3>, Option<PolygonMesh>>;

#[derive(Clone, Debug)]
struct PairAnalysis<C, S> {
    tolerance: f64,
    poly_shell0: PolyShell,
    poly_shell1: PolyShell,
    winding0: WindingClassifier,
    winding1: WindingClassifier,
    cls0: faces_classification::FacesClassification<Point3, AltCurve<C, S>, S>,
    cls1: faces_classification::FacesClassification<Point3, AltCurve<C, S>, S>,
    coplanar_overlap: bool,
}

fn analyze_pair<C: ShapeOpsCurve<S>, S: ShapeOpsSurface>(
    shell0: &Shell<Point3, C, S>,
    shell1: &Shell<Point3, C, S>,
    tol: f64,
    mode: PairMode,
) -> ClassicResult<PairAnalysis<C, S>> {
    validate_tolerance(tol)?;
    let seed_tol = BooleanTolerance::new(tol).mesh_seed;
    let poly_shell0 = shell0.triangulation(seed_tol);
    let poly_shell1 = shell1.triangulation(seed_tol);
    let altshell0: AltCurveShell<C, S> =
        shell0.mapped(|x| *x, |c| Alternative::FirstType(c.clone()), Clone::clone);
    let altshell1: AltCurveShell<C, S> =
        shell1.mapped(|x| *x, |c| Alternative::FirstType(c.clone()), Clone::clone);
    let loops_store::InterfaceLoops {
        geom: [loops_store0, loops_store1],
        sampled: [sampled_store0, sampled_store1],
        coplanar_overlap,
    } = loops_store::create_loops_stores(
        &altshell0,
        &poly_shell0,
        &altshell1,
        &poly_shell1,
        tol,
        matches!(mode, PairMode::AndOr(_)),
    )
    .ok_or(ShapeOpsError::InterfaceAnalysisFailed)?;
    let cls0 = divide_face::divide_faces(&altshell0, &loops_store0, sampled_store0.as_ref(), tol)
        .ok_or(ShapeOpsError::DivideFacesFailed { shell_index: 0 })?;
    let cls1 = divide_face::divide_faces(&altshell1, &loops_store1, sampled_store1.as_ref(), tol)
        .ok_or(ShapeOpsError::DivideFacesFailed { shell_index: 1 })?;
    Ok(PairAnalysis {
        tolerance: tol,
        winding0: WindingClassifier::new(&poly_shell0),
        winding1: WindingClassifier::new(&poly_shell1),
        poly_shell0,
        poly_shell1,
        cls0,
        cls1,
        coplanar_overlap,
    })
}

fn classify_unknown<C, S>(
    winding: &WindingClassifier,
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
        return classify_unknown_face(winding, &face, tol)
            .ok_or(ShapeOpsError::UnknownClassificationFailed { shell_index });
    }
    let points =
        divide_face::sample_points_on_region(region.absolute_boundaries(), &region.surface(), tol)
            .ok_or(ShapeOpsError::UnknownClassificationFailed { shell_index })?;
    classify_points_inside(winding, points)
        .ok_or(ShapeOpsError::UnknownClassificationFailed { shell_index })
}

#[derive(Clone, Copy)]
enum DifferenceContactSide {
    /// Source material survives only where the cutter is absent immediately
    /// on the interior side of its original face.
    SourceInterior,
    /// An inverted cutter face is exposed where original material occupies
    /// the exterior side of the cutter.
    CutterExterior,
}

fn classify_difference_unknown<C, S>(
    winding: &WindingClassifier,
    region: &faces_classification::BooleanRegion<Point3, C, S>,
    tol: f64,
    shell_index: usize,
    side: DifferenceContactSide,
) -> ClassicResult<bool>
where
    C: PolylineableCurve + ParameterDivision1D<Point = Point3> + Clone,
    S: MeshableSurface
        + Clone
        + ParametricSurface3D
        + SearchParameter<SurfaceParameter, Point = Point3>
        + SearchNearestParameter<SurfaceParameter, Point = Point3>,
{
    if let Ok(label) = classify_unknown(winding, region, tol, shell_index) {
        return Ok(label);
    }
    let points = region
        .try_face()
        .map_or_else(
            || {
                divide_face::sample_points_on_region(
                    region.absolute_boundaries(),
                    &region.surface(),
                    tol,
                )
            },
            |face| classification::sample_points_on_face(&face, tol),
        )
        .ok_or(ShapeOpsError::UnknownClassificationFailed { shell_index })?;
    let [exterior, interior] = classification::unanimous_boundary_contact_sides(
        winding,
        &points,
        &region.surface(),
        region.orientation(),
        tol,
    )
    .ok_or(ShapeOpsError::UnknownClassificationFailed { shell_index })?;
    Ok(match side {
        DifferenceContactSide::SourceInterior => interior,
        DifferenceContactSide::CutterExterior => exterior,
    })
}

/// Corroborate a preclassified region against actual opposing-solid volume.
///
/// A boundary-local And/Or label does not establish occupancy throughout the
/// face, especially after SSI imprinting. We correct it only if every sampled
/// interior point is classifiable, agrees, and has the same inside/outside
/// relation after equal-magnitude offsets to both sides of its source surface.
/// Mixed samples, boundary contacts, and zero-area regions are NOT evidence
/// for either label and must never be force-classified.
fn unanimous_region_containment<C, S>(
    winding: &WindingClassifier,
    region: &faces_classification::BooleanRegion<Point3, C, S>,
    tol: f64,
) -> Option<bool>
where
    C: PolylineableCurve + ParameterDivision1D<Point = Point3> + Clone,
    S: MeshableSurface
        + Clone
        + ParametricSurface3D
        + SearchParameter<SurfaceParameter, Point = Point3>
        + SearchNearestParameter<SurfaceParameter, Point = Point3>,
{
    let samples = if let Some(face) = region.try_face() {
        classification::sample_points_on_face(&face, tol)?
    } else {
        divide_face::sample_points_on_region(region.absolute_boundaries(), &region.surface(), tol)?
    };
    let surface = region.surface();
    let classify = |point: Point3| -> Option<bool> {
        let verdict = winding.inside(point)?;
        // A point on the other solid's boundary can receive a misleading
        // signed-winding verdict. Both small physical displacements normal to
        // this face must agree with the point's label; a coincident/tangent
        // contact therefore cannot override a topological boundary label.
        let (u, v) = surface
            .search_parameter(point, None, 100)
            .or_else(|| surface.search_nearest_parameter(point, None, 100))?;
        let normal = surface.normal(u, v);
        let length = normal.magnitude();
        if !length.is_finite() || length <= f64::EPSILON {
            return None;
        }
        let offset = normal / length * BooleanTolerance::new(tol).classifier_offset;
        (winding.inside(point + offset)? == verdict && winding.inside(point - offset)? == verdict)
            .then_some(verdict)
    };
    let mut samples = samples.into_iter();
    let first = classify(samples.next()?)?;
    samples.try_fold(first, |label, point| {
        (classify(point)? == label).then_some(label)
    })
}

fn classify_relation<C, S>(
    winding: &WindingClassifier,
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
        return classify_face_relation(winding, &face, tol)
            .ok_or(ShapeOpsError::UnknownClassificationFailed { shell_index });
    }
    let points =
        divide_face::sample_points_on_region(region.absolute_boundaries(), &region.surface(), tol)
            .ok_or(ShapeOpsError::UnknownClassificationFailed { shell_index })?;
    let mut surface = region.surface();
    if !region.orientation() {
        surface.invert();
    }
    classify_points_relation(winding, points, &surface, tol)
        .ok_or(ShapeOpsError::UnknownClassificationFailed { shell_index })
}

fn convert_alt_shell<C: ShapeOpsCurve<S>, S: ShapeOpsSurface>(
    shell: &AltCurveShell<C, S>,
    operation: &'static str,
    tolerances: BooleanTolerance,
) -> ClassicResult<Shell<Point3, C, S>> {
    altshell_to_shell(shell, tolerances)
        .ok_or(ShapeOpsError::AltShellConversionFailed { operation })
}

fn weld_and_convert<C: ShapeOpsCurve<S>, S: ShapeOpsSurface>(
    regions: &[AltRegion<C, S>],
    operation: &'static str,
    tolerances: BooleanTolerance,
) -> ClassicResult<Shell<Point3, C, S>> {
    let shell = weld_coincident_regions(regions, tolerances)
        .ok_or(ShapeOpsError::RegionMaterializationFailed { operation })?;
    convert_alt_shell(&shell, operation, tolerances)
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
) -> ClassicResult<Shell<Point3, C, S>> {
    let tol = analysis.tolerance;
    analysis.cls0.integrate_by_component();
    analysis.cls1.integrate_by_component();
    let full_positive_regions = analysis.cls0.all_regions();
    let full_cutter_regions = analysis.cls1.all_regions();
    // Bind the material-side selection once; both the labeled overlay and
    // the fallback reconstruction must make identical contact decisions.
    let source_inside_cutter = |region: &AltRegion<C, S>| {
        classify_difference_unknown(
            &analysis.winding1,
            region,
            tol,
            0,
            DifferenceContactSide::SourceInterior,
        )
    };
    let cutter_exterior_inside_source = |region: &AltRegion<C, S>| {
        classify_difference_unknown(
            &analysis.winding0,
            region,
            tol,
            1,
            DifferenceContactSide::CutterExterior,
        )
    };

    let [known_and0, known_or0, unknown0] = analysis.cls0.and_or_unknown();
    // An SSI boundary label describes the winding at one interface; it is
    // not evidence that the entire resulting face lies inside/outside the
    // other solid. In particular, an imprinted, retraced edge can mark a
    // large region And although all its interior witnesses are outside.
    // Only contradict labels with UNANIMOUS geometric witnesses. If the
    // region cannot be sampled or has mixed classifications, retain its
    // existing status; never fabricate an interior for zero-area regions.
    let mut difference0 = Vec::new();
    for region in known_and0 {
        if unanimous_region_containment(&analysis.winding1, &region, tol) == Some(false) {
            difference0.push(region);
        }
    }
    for region in known_or0 {
        if unanimous_region_containment(&analysis.winding1, &region, tol) != Some(true) {
            difference0.push(region);
        }
    }
    for region in unknown0 {
        if !source_inside_cutter(&region)? {
            difference0.push(region);
        }
    }
    let [mut difference1, _outside1, unknown1] = analysis.cls1.and_or_unknown();
    for region in unknown1 {
        if cutter_exterior_inside_source(&region)? {
            difference1.push(region);
        }
    }
    for region in &mut difference1 {
        region.invert();
    }
    let label_positive_len = difference0.len();
    difference0.append(&mut difference1);
    // A failed arrangement or weld is not evidence that the alternative
    // geometric classification is invalid. Try it before reporting failure,
    // just as when the labeled shell materializes but has invalid topology.
    if let Some(label_regions) = planar_union::overlay_difference_regions(
        difference0,
        label_positive_len,
        &full_positive_regions,
        &full_cutter_regions,
        tol,
    ) && let Ok(mut label_shell) =
        weld_and_convert(&label_regions, "difference", BooleanTolerance::new(tol))
    {
        let _ = orient_regular_shell(&mut label_shell);
        if label_shell
            .connected_components()
            .iter()
            .all(|component| component.check_solid_boundary().is_ok())
        {
            return Ok(label_shell);
        }
    }

    let [and0, or0, unknown0] = analysis.cls0.and_or_unknown();
    let [and1, or1, unknown1] = analysis.cls1.and_or_unknown();
    let mut fallback = Vec::new();
    for region in and0.into_iter().chain(or0).chain(unknown0) {
        if !source_inside_cutter(&region)? {
            fallback.push(region);
        }
    }
    let fallback_positive_len = fallback.len();
    for mut region in and1.into_iter().chain(or1).chain(unknown1) {
        if cutter_exterior_inside_source(&region)? {
            region.invert();
            fallback.push(region);
        }
    }
    let fallback_regions = planar_union::overlay_difference_regions(
        fallback,
        fallback_positive_len,
        &full_positive_regions,
        &full_cutter_regions,
        tol,
    )
    .ok_or(ShapeOpsError::RegionMaterializationFailed {
        operation: "difference_planar_overlay",
    })?;
    weld_and_convert(&fallback_regions, "difference", BooleanTolerance::new(tol))
}

fn process_non_coplanar_analysis<C: ShapeOpsCurve<S>, S: ShapeOpsSurface>(
    mut analysis: PairAnalysis<C, S>,
    output: PartitionOutput,
) -> ClassicResult<Shell<Point3, C, S>> {
    let tol = analysis.tolerance;
    analysis.cls0.integrate_by_component();
    analysis.cls1.integrate_by_component();
    let mut selected = Vec::new();
    for (operand, (classification, opposing_winding)) in [
        (analysis.cls0, &analysis.winding1),
        (analysis.cls1, &analysis.winding0),
    ]
    .into_iter()
    .enumerate()
    {
        let [known, unknown] = classification.into_selected(output);
        // Match the original known-first, then unknown, operand traversal.
        selected.extend(known);
        for region in unknown {
            if classify_unknown(opposing_winding, &region, tol, operand)?
                == output.includes_inside()
            {
                selected.push(region);
            }
        }
    }
    // The complementary partition is not a prerequisite for materializing
    // this operation. Coplanar label-only proofs are handled separately.
    weld_and_convert(&selected, output.operation(), BooleanTolerance::new(tol))
}

struct PartitionEvidence<C, S> {
    regions: Vec<AltRegion<C, S>>,
    polygons: Vec<Option<PolygonMesh>>,
    by_key: HashMap<usize, usize>,
}

fn partition_evidence<C: ShapeOpsCurve<S>, S: ShapeOpsSurface>(
    regions: &[AltRegion<C, S>],
    mesh_tol: f64,
) -> PartitionEvidence<C, S> {
    let mut valid_regions = Vec::new();
    let mut polygons = Vec::new();
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
        by_key.insert(region.key(), polygons.len());
        valid_regions.push(region.clone());
        polygons.push(Some(polygon));
    }
    PartitionEvidence {
        regions: valid_regions,
        polygons,
        by_key,
    }
}

/// Bounding containment is necessary (but not sufficient) for a Boolean
/// union. In particular it rejects a closed output containing only one
/// operand, which topological checks alone cannot identify.
fn covers_both_input_bounding_boxes<C: ShapeOpsCurve<S>, S: ShapeOpsSurface>(
    candidate: &Shell<Point3, C, S>,
    lhs: &PolyShell,
    rhs: &PolyShell,
    mesh_tol: f64,
) -> bool {
    fn bounds(shell: &PolyShell) -> Option<(Point3, Point3)> {
        let mut points = shell
            .iter()
            .filter_map(Face::surface)
            .flat_map(|polygon| polygon.positions().to_vec());
        let first = points.next()?;
        Some(points.fold((first, first), |(mut lo, mut hi), p| {
            for axis in 0..3 {
                lo[axis] = lo[axis].min(p[axis]);
                hi[axis] = hi[axis].max(p[axis]);
            }
            (lo, hi)
        }))
    }
    let sampled = candidate.triangulation(mesh_tol);
    let (Some((actual_lo, actual_hi)), Some((lhs_lo, lhs_hi)), Some((rhs_lo, rhs_hi))) =
        (bounds(&sampled), bounds(lhs), bounds(rhs))
    else {
        return false;
    };
    (0..3).all(|axis| {
        actual_lo[axis] <= lhs_lo[axis].min(rhs_lo[axis]) + mesh_tol
            && actual_hi[axis] >= lhs_hi[axis].max(rhs_hi[axis]) - mesh_tol
    })
}

/// A non-materializable arrangement region can contain several *positive*
/// planar patches meeting at vertices. Classifying the entire region averages
/// unrelated cells (e.g. one outside and one coincident), creating a tie.
/// Split only when each boundary is independently a positive, valid face.
/// A negative contour is a hole and must stay with its parent, not be filled.
fn positive_region_components<C: ShapeOpsCurve<S>, S: ShapeOpsSurface>(
    region: &AltRegion<C, S>,
    tol: f64,
) -> Vec<AltRegion<C, S>> {
    let boundaries = region.absolute_boundaries();
    if boundaries.len() < 2 || region.try_face().is_some() {
        return vec![region.clone()];
    }
    let mut parts = Vec::with_capacity(boundaries.len());
    for wire in boundaries {
        if Face::new(vec![wire.clone()], region.surface()).is_err()
            || divide_face::sample_points_on_region(
                std::slice::from_ref(wire),
                &region.surface(),
                tol,
            )
            .is_none()
        {
            return vec![region.clone()];
        }
        parts.push(AltRegion::new(
            vec![wire.clone()],
            region.surface(),
            region.orientation(),
            region.key(),
            region.source_face(),
        ));
    }
    parts
}

/// Accept the fast label-only partition only when it is closed *and* covers
/// the two input bounds; otherwise the full coplanar partition must run.
fn classify_complete_label_shells<C: ShapeOpsCurve<S>, S: ShapeOpsSurface>(
    analysis: &PairAnalysis<C, S>,
) -> ClassicResult<Option<[Shell<Point3, C, S>; 2]>> {
    let tol = analysis.tolerance;
    let [mut and0, mut or0, unknown0] = analysis.cls0.and_or_unknown();
    for region in unknown0 {
        for region in positive_region_components(&region, tol) {
            match classify_relation(&analysis.winding1, &region, tol, 0)? {
                FaceRelation::Inside => and0.push(region),
                FaceRelation::Outside | FaceRelation::CoincidentSameFacing => or0.push(region),
            }
        }
    }
    let [mut and1, mut or1, unknown1] = analysis.cls1.and_or_unknown();
    for region in unknown1 {
        for region in positive_region_components(&region, tol) {
            match classify_relation(&analysis.winding0, &region, tol, 1)? {
                FaceRelation::Inside | FaceRelation::CoincidentSameFacing => and1.push(region),
                FaceRelation::Outside => or1.push(region),
            }
        }
    }
    and0.append(&mut and1);
    or0.append(&mut or1);
    let label_and_shell = weld_and_convert(&and0, "and", BooleanTolerance::new(tol))?;
    let label_or_shell = weld_and_convert(&or0, "or", BooleanTolerance::new(tol))?;
    if !label_or_shell.is_empty()
        && shell_is_valid(&label_and_shell)
        && shell_is_valid(&label_or_shell)
        && covers_both_input_bounding_boxes(
            &label_or_shell,
            &analysis.poly_shell0,
            &analysis.poly_shell1,
            BooleanTolerance::new(tol).mesh_seed,
        )
    {
        return Ok(Some([label_and_shell, label_or_shell]));
    }
    Ok(None)
}

fn process_coplanar_analysis<C: ShapeOpsCurve<S>, S: ShapeOpsSurface>(
    mut analysis: PairAnalysis<C, S>,
    output: PartitionOutput,
) -> ClassicResult<Shell<Point3, C, S>> {
    let tol = analysis.tolerance;
    analysis.cls0.integrate_by_component();
    analysis.cls1.integrate_by_component();

    if let Some(shells) = classify_complete_label_shells(&analysis)? {
        // The complete-label shortcut requires both complementary shells as
        // evidence that no input material was lost. Retain that proof.
        return Ok(output.select(shells));
    }

    // A topologically valid set of classified faces need not cover the input
    // solids. With partially contacting coplanar faces, the label-only path
    // may form a valid closed subset (e.g. the attached hook but not the
    // annular host), silently discarding source material. Always partition
    // every source region and overlay the surviving coplanar cells.
    let all0 = analysis.cls0.all_regions();
    let all1 = analysis.cls1.all_regions();
    let partition_tolerance = BooleanTolerance::new(tol);
    let partition_geom_tol = partition_tolerance.length;
    let partition_mesh_tol = partition_tolerance.mesh_seed;
    let evidence0 = partition_evidence(&all0, partition_mesh_tol);
    let evidence1 = partition_evidence(&all1, partition_mesh_tol);
    let planar0 = loops_store::PlanarMesh::from_polygons(&evidence0.polygons, partition_geom_tol);
    let planar1 = loops_store::PlanarMesh::from_polygons(&evidence1.polygons, partition_geom_tol);

    let [and0_labeled, or0_labeled, unknown0] = analysis.cls0.and_or_unknown();
    let mut and0 = Vec::new();
    let mut or0 = Vec::new();
    for region in and0_labeled.into_iter().chain(or0_labeled).chain(unknown0) {
        for region in positive_region_components(&region, tol) {
            let coincident = evidence0
                .by_key
                .get(&region.key())
                .and_then(|&index| planar0[index].as_ref())
                .is_some_and(|candidate| {
                    planar_partition_same_facing(
                        &region,
                        candidate,
                        &evidence1.regions,
                        &planar1,
                        partition_geom_tol,
                    ) == Some(true)
                });
            let relation = if coincident {
                FaceRelation::CoincidentSameFacing
            } else {
                classify_relation(&analysis.winding1, &region, tol, 0)?
            };
            match relation {
                FaceRelation::Inside => and0.push(region),
                FaceRelation::Outside | FaceRelation::CoincidentSameFacing => or0.push(region),
            }
        }
    }

    let [and1_labeled, or1_labeled, unknown1] = analysis.cls1.and_or_unknown();
    let mut and1 = Vec::new();
    let mut or1 = Vec::new();
    for region in and1_labeled.into_iter().chain(or1_labeled).chain(unknown1) {
        for region in positive_region_components(&region, tol) {
            let coincident = evidence1
                .by_key
                .get(&region.key())
                .and_then(|&index| planar1[index].as_ref())
                .is_some_and(|candidate| {
                    planar_partition_same_facing(
                        &region,
                        candidate,
                        &evidence0.regions,
                        &planar0,
                        partition_geom_tol,
                    ) == Some(true)
                });
            let relation = if coincident {
                FaceRelation::CoincidentSameFacing
            } else {
                classify_relation(&analysis.winding0, &region, tol, 1)?
            };
            match relation {
                FaceRelation::Inside | FaceRelation::CoincidentSameFacing => and1.push(region),
                FaceRelation::Outside => or1.push(region),
            }
        }
    }
    let mut selected = planar_union::overlay_source_regions(output.select([and0, or0]), tol);
    selected.extend(planar_union::overlay_source_regions(
        output.select([and1, or1]),
        tol,
    ));
    weld_and_convert(&selected, output.operation(), BooleanTolerance::new(tol))
}

fn process_one_pair_of_shells<C: ShapeOpsCurve<S>, S: ShapeOpsSurface>(
    shell0: &Shell<Point3, C, S>,
    shell1: &Shell<Point3, C, S>,
    tol: f64,
    mode: PairMode,
) -> ClassicResult<Shell<Point3, C, S>> {
    let analysis = analyze_pair(shell0, shell1, tol, mode)?;
    match mode {
        PairMode::Difference => process_difference_analysis(analysis),
        PairMode::AndOr(output) if analysis.coplanar_overlap => {
            process_coplanar_analysis(analysis, output)
        }
        PairMode::AndOr(output) => process_non_coplanar_analysis(analysis, output),
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

/// Separate closed boundary sheets that share only a non-manifold edge.
///
/// A Boolean union of two solids touching along an edge is not a single
/// two-manifold boundary, but it can be represented by separate closed
/// boundary shells. Only edges with exactly two incident face uses may
/// connect faces within a manifold sheet. Every resulting sheet must still
/// pass the full canonical solid-boundary validator before being accepted.
fn partition_manifold_sheets<C: ShapeOpsCurve<S>, S: ShapeOpsSurface>(
    shell: &Shell<Point3, C, S>,
) -> Option<Vec<Shell<Point3, C, S>>> {
    let faces: Vec<_> = shell.iter().cloned().collect();
    let mut edge_users = HashMap::<EdgeId<C>, Vec<usize>>::new();
    for (face_index, face) in faces.iter().enumerate() {
        for edge in face.edge_iter() {
            edge_users.entry(edge.id()).or_default().push(face_index);
        }
    }
    if !edge_users.values().any(|uses| uses.len() > 2) {
        return None;
    }
    let mut adjacency = vec![Vec::<usize>::new(); faces.len()];
    for uses in edge_users.values() {
        if let [a, b] = uses.as_slice()
            && a != b
        {
            adjacency[*a].push(*b);
            adjacency[*b].push(*a);
        }
    }
    let mut seen = vec![false; faces.len()];
    let mut components = Vec::new();
    for start in 0..faces.len() {
        if seen[start] {
            continue;
        }
        seen[start] = true;
        let mut queue = VecDeque::from([start]);
        let mut members = Vec::new();
        while let Some(index) = queue.pop_front() {
            members.push(index);
            for &neighbor in &adjacency[index] {
                if !seen[neighbor] {
                    seen[neighbor] = true;
                    queue.push_back(neighbor);
                }
            }
        }
        // Preserve input face order independent of hash-table traversal.
        members.sort_unstable();
        let mut component: Shell<_, _, _> = members
            .into_iter()
            .map(|index| faces[index].clone())
            .collect();
        let _ = orient_regular_shell(&mut component);
        if component.check_solid_boundary().is_err() {
            return None;
        }
        components.push(component);
    }
    (components.len() > 1).then_some(components)
}

fn finalize<C: ShapeOpsCurve<S>, S: ShapeOpsSurface>(
    operation: &'static str,
    mut shell: Shell<Point3, C, S>,
    tol: f64,
) -> ClassicResult<Solid<Point3, C, S>> {
    let _ = orient_regular_shell(&mut shell);
    if !shell_is_valid(&shell)
        && let Some(mut repaired) =
            weld::weld_converted_duplicate_edges(&shell, BooleanTolerance::new(tol))
    {
        let _ = orient_regular_shell(&mut repaired);
        if shell_is_valid(&repaired)
            && repaired
                .edge_iter()
                .all(|edge| edge.is_geometric_consistent())
        {
            shell = repaired;
        }
    }
    let shell = planar_union::split_disconnected_planar_faces(&shell, tol);
    let mut boundaries = shell.connected_components();
    if boundaries
        .iter()
        .any(|boundary| boundary.check_solid_boundary().is_err())
        && let Some(manifold_sheets) = partition_manifold_sheets(&shell)
    {
        boundaries = manifold_sheets;
    }
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
                edge_uses.push((edge.id(), edge, 1));
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
    Solid::new(boundaries).map_err(|source| ShapeOpsError::InvalidOutputShell { operation, source })
}

/// Fold shells in source order, retaining the left-hand cavity semantics of difference.
fn combine_solids<C: ShapeOpsCurve<S>, S: ShapeOpsSurface>(
    solid0: &Solid<Point3, C, S>,
    solid1: &Solid<Point3, C, S>,
    tol: f64,
    mode: PairMode,
) -> ClassicResult<Solid<Point3, C, S>> {
    validate_tolerance(tol)?;
    let mut left = solid0.boundaries().iter();
    let mut right = solid1.boundaries().iter();
    let Some(first_left) = left.next() else {
        return Ok(match mode {
            PairMode::AndOr(PartitionOutput::Union) => solid1.clone(),
            PairMode::AndOr(PartitionOutput::Intersection) | PairMode::Difference => Solid::empty(),
        });
    };
    let Some(first_right) = right.next() else {
        return Ok(match mode {
            PairMode::AndOr(PartitionOutput::Intersection) => Solid::empty(),
            PairMode::AndOr(PartitionOutput::Union) | PairMode::Difference => solid0.clone(),
        });
    };
    let left_mode = match mode {
        PairMode::Difference => PairMode::AndOr(PartitionOutput::Intersection),
        PairMode::AndOr(_) => mode,
    };
    let initial = process_one_pair_of_shells(first_left, first_right, tol, mode)?;
    let combined = left
        .map(|shell| (shell, left_mode))
        .chain(right.map(|shell| (shell, mode)))
        .try_fold(initial, |acc, (shell, pair_mode)| {
            process_one_pair_of_shells(&acc, shell, tol, pair_mode)
        })?;
    let operation = match mode {
        PairMode::AndOr(output) => output.operation(),
        PairMode::Difference => "difference",
    };
    finalize(operation, combined, tol)
}

/// AND operation between two solids (classic backend).
pub fn and<C: ShapeOpsCurve<S>, S: ShapeOpsSurface>(
    solid0: &Solid<Point3, C, S>,
    solid1: &Solid<Point3, C, S>,
    tol: f64,
) -> ClassicResult<Solid<Point3, C, S>> {
    combine_solids(
        solid0,
        solid1,
        tol,
        PairMode::AndOr(PartitionOutput::Intersection),
    )
}

/// OR operation between two solids (classic backend).
pub fn or<C: ShapeOpsCurve<S>, S: ShapeOpsSurface>(
    solid0: &Solid<Point3, C, S>,
    solid1: &Solid<Point3, C, S>,
    tol: f64,
) -> ClassicResult<Solid<Point3, C, S>> {
    combine_solids(solid0, solid1, tol, PairMode::AndOr(PartitionOutput::Union))
}

/// Difference: the region inside `solid0` but outside `solid1` (classic backend).
pub fn difference<C: ShapeOpsCurve<S>, S: ShapeOpsSurface>(
    solid0: &Solid<Point3, C, S>,
    solid1: &Solid<Point3, C, S>,
    tol: f64,
) -> ClassicResult<Solid<Point3, C, S>> {
    combine_solids(solid0, solid1, tol, PairMode::Difference)
}

/// Symmetric difference (XOR): the region inside exactly one solid (classic backend).
pub fn symmetric_difference<C: ShapeOpsCurve<S>, S: ShapeOpsSurface>(
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
    use monstertruck_modeling::{Curve, Surface, builder};

    fn unit_cube() -> Solid<Point3, Curve, Surface> {
        let vertex = builder::vertex(Point3::origin());
        let edge = builder::extrude(&vertex, Vector3::unit_x());
        let face = builder::extrude(&edge, Vector3::unit_y());
        builder::extrude(&face, Vector3::unit_z())
    }

    #[test]
    fn preclassified_difference_regions_require_unanimous_geometric_corroboration() {
        let solid = unit_cube();
        let poly_shell = solid.boundaries()[0].triangulation(0.01);
        let winding = WindingClassifier::new(&poly_shell);
        let flat_region = |xlo: f64, xhi: f64, z: f64, key: usize| {
            let points = [
                Point3::new(xlo, 0.25, z),
                Point3::new(xhi, 0.25, z),
                Point3::new(xhi, 0.75, z),
                Point3::new(xlo, 0.75, z),
            ];
            let vertices = points.map(Vertex::new);
            let wire: Wire<Point3, AltCurve<Curve, Surface>> = (0..4)
                .map(|index| {
                    let next = (index + 1) % 4;
                    Edge::new_unchecked(
                        &vertices[index],
                        &vertices[next],
                        Alternative::FirstType(Curve::from(Line(points[index], points[next]))),
                    )
                })
                .collect();
            AltRegion::new(
                vec![wire],
                Surface::Plane(Plane::new(
                    Point3::new(0.0, 0.0, z),
                    Point3::new(1.0, 0.0, z),
                    Point3::new(0.0, 1.0, z),
                )),
                true,
                key,
                key,
            )
        };
        let outside = flat_region(0.25, 0.75, 1.5, 1);
        let inside = flat_region(0.25, 0.75, 0.5, 2);
        let crossing = flat_region(-0.3, 0.3, 0.5, 3);
        assert_eq!(
            unanimous_region_containment(&winding, &outside, 1.0e-5),
            Some(false),
            "confident outside evidence must recover a region incorrectly marked And"
        );
        assert_eq!(
            unanimous_region_containment(&winding, &inside, 1.0e-5),
            Some(true),
            "confident inside evidence must reject a region incorrectly marked Or"
        );
        assert_eq!(
            unanimous_region_containment(&winding, &crossing, 1.0e-5),
            None,
            "a region straddling the other shell cannot inherit either label"
        );

        let boundary_contact = {
            let points = [
                Point3::new(1.0, 0.25, 0.25),
                Point3::new(1.0, 0.75, 0.25),
                Point3::new(1.0, 0.75, 0.75),
                Point3::new(1.0, 0.25, 0.75),
            ];
            let vertices = points.map(Vertex::new);
            let wire: Wire<Point3, AltCurve<Curve, Surface>> = (0..4)
                .map(|index| {
                    let next = (index + 1) % 4;
                    Edge::new_unchecked(
                        &vertices[index],
                        &vertices[next],
                        Alternative::FirstType(Curve::from(Line(points[index], points[next]))),
                    )
                })
                .collect();
            AltRegion::new(
                vec![wire],
                Surface::Plane(Plane::new(
                    Point3::new(1.0, 0.0, 0.0),
                    Point3::new(1.0, 1.0, 0.0),
                    Point3::new(1.0, 0.0, 1.0),
                )),
                true,
                5,
                5,
            )
        };
        assert_eq!(
            unanimous_region_containment(&winding, &boundary_contact, 1.0e-5),
            None,
            "an exact coincident boundary cannot certify inside or outside from the surface itself"
        );

        // This is the small topology pattern seen in the source-derived stage6:
        // two forward collinear edges and one retraced inverse edge. All exact
        // curve samples are collinear: there is no real 2-D area to classify.
        let points = [
            Point3::new(0.2, 0.2, 0.5),
            Point3::new(0.4, 0.2, 0.5),
            Point3::new(0.8, 0.2, 0.5),
        ];
        let vertices = points.map(Vertex::new);
        let retraced: Wire<Point3, AltCurve<Curve, Surface>> = (0..3)
            .map(|index| {
                let next = (index + 1) % 3;
                Edge::new_unchecked(
                    &vertices[index],
                    &vertices[next],
                    Alternative::FirstType(Curve::from(Line(points[index], points[next]))),
                )
            })
            .collect();
        let zero_area = AltRegion::new(vec![retraced], Surface::Plane(Plane::xy()), true, 4, 4);
        assert_eq!(
            unanimous_region_containment(&winding, &zero_area, 1.0e-5),
            None,
            "zero-area retraces must not gain an invented interior witness"
        );
    }

    #[test]
    fn classify_touching_positive_face_patches_independently_but_keep_holes() {
        let coords = [
            (0.0, 0.0),
            (2.0, 0.0),
            (2.0, 2.0),
            (0.0, 2.0),
            (4.0, 2.0),
            (4.0, 4.0),
            (2.0, 4.0),
            (0.5, 0.5),
            (0.5, 1.5),
            (1.5, 1.5),
            (1.5, 0.5),
        ];
        let vertices = coords.map(|(x, y)| Vertex::new(Point3::new(x, y, 0.0)));
        let wire = |path: &[usize]| -> Wire<Point3, AltCurve<Curve, Surface>> {
            (0..path.len())
                .map(|index| {
                    let front = &vertices[path[index]];
                    let back = &vertices[path[(index + 1) % path.len()]];
                    Edge::new_unchecked(
                        front,
                        back,
                        Alternative::FirstType(Curve::from(Line(front.point(), back.point()))),
                    )
                })
                .collect()
        };
        let left = wire(&[0, 1, 2, 3]);
        let right = wire(&[2, 4, 5, 6]);
        let support = Surface::Plane(Plane::xy());
        let touching = AltRegion::new(vec![left.clone(), right], support.clone(), true, 7, 11);
        assert!(touching.try_face().is_none());
        let parts = positive_region_components(&touching, 1.0e-6);
        assert_eq!(
            parts.len(),
            2,
            "two additive regions must not share one classifier vote"
        );
        assert!(parts.iter().all(|part| part.try_face().is_some()));
        assert!(parts.iter().all(|part| part.source_face() == 11));

        let hole = wire(&[7, 8, 9, 10]);
        let annulus = AltRegion::new(vec![left, hole], support, true, 8, 12);
        let parts = positive_region_components(&annulus, 1.0e-6);
        assert_eq!(
            parts.len(),
            1,
            "negative hole must remain attached to its host"
        );
        assert_eq!(parts[0].absolute_boundaries().len(), 2);
    }

    #[test]
    fn empty_set_boolean_identities_are_total() {
        let empty = Solid::<Point3, Curve, Surface>::empty();
        let cube = unit_cube();
        let cube_faces = cube.face_iter().count();

        assert!(matches!(and(&empty, &cube, 1.0e-6), Ok(solid) if solid.is_empty()));
        assert!(matches!(and(&cube, &empty, 1.0e-6), Ok(solid) if solid.is_empty()));
        assert!(
            matches!(or(&empty, &cube, 1.0e-6), Ok(solid) if solid.face_iter().count() == cube_faces)
        );
        assert!(
            matches!(or(&cube, &empty, 1.0e-6), Ok(solid) if solid.face_iter().count() == cube_faces)
        );
        assert!(matches!(difference(&empty, &cube, 1.0e-6), Ok(solid) if solid.is_empty()));
        assert!(
            matches!(difference(&cube, &empty, 1.0e-6), Ok(solid) if solid.face_iter().count() == cube_faces)
        );
        assert!(
            matches!(symmetric_difference(&empty, &cube, 1.0e-6), Ok(solid) if solid.face_iter().count() == cube_faces)
        );
    }

    #[test]
    fn invalid_tolerance_is_rejected_before_empty_set_shortcuts() {
        let empty = Solid::<Point3, Curve, Surface>::empty();
        let cube = unit_cube();
        for tol in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            assert!(matches!(
                and(&empty, &cube, tol),
                Err(ShapeOpsError::InvalidTolerance)
            ));
            assert!(matches!(
                or(&empty, &cube, tol),
                Err(ShapeOpsError::InvalidTolerance)
            ));
            assert!(matches!(
                difference(&cube, &empty, tol),
                Err(ShapeOpsError::InvalidTolerance)
            ));
        }
    }

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
        .expect("fixture mesh is valid")
    }

    #[test]
    fn signed_winding_does_not_turn_negative_shell_into_inside() {
        let reversed = tetrahedron_mesh(true);
        let shell: Shell<Point3, PolylineCurve<Point3>, Option<PolygonMesh>> =
            Shell::from(vec![Face::new_unchecked(Vec::new(), Some(reversed))]);
        assert_eq!(
            WindingClassifier::new(&shell).inside(Point3::new(0.1, 0.1, 0.1)),
            Some(false)
        );
    }

    #[test]
    fn signed_winding_honors_topological_face_orientation() {
        let mut face = Face::new_unchecked(Vec::new(), Some(tetrahedron_mesh(false)));
        face.invert();
        let shell: Shell<Point3, PolylineCurve<Point3>, Option<PolygonMesh>> =
            Shell::from(vec![face]);
        assert_eq!(
            WindingClassifier::new(&shell).inside(Point3::new(0.1, 0.1, 0.1)),
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
                Edge::new_unchecked(
                    &vertices[index],
                    &vertices[next],
                    Curve::from(Line(vertices[index].point(), vertices[next].point())),
                )
            })
            .collect::<Vec<_>>();
        let face = Face::new_unchecked(vec![Wire::from(edges)], Surface::Plane(Plane::xy()));
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
                Edge::new_unchecked(
                    &vertices[index],
                    &vertices[next],
                    Curve::from(Line(vertices[index].point(), vertices[next].point())),
                )
            })
            .collect::<Vec<_>>();
        let face = Face::new_unchecked(vec![Wire::from(edges)], Surface::Plane(Plane::xy()));

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
        .expect("fixture mesh is valid")
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
            WindingClassifier::new(&shell).inside(Point3::new(0.25, 0.25, 0.25)),
            Some(true)
        );
        assert_eq!(
            WindingClassifier::new(&shell).inside(Point3::new(0.6, 0.6, 0.6)),
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
                    Edge::new_unchecked(
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
        let face =
            Face::new_unchecked(vec![wire(&outer), wire(&hole)], Surface::Plane(Plane::xy()));

        let samples = sample_points_on_face(&face, 1.0e-6).unwrap_or_default();
        assert!(!samples.is_empty());
        assert!(samples.iter().all(|point| {
            let x = point.x;
            let y = point.y;
            x <= 1.0 + KERNEL_EPSILON
                || x >= 3.0 - KERNEL_EPSILON
                || y <= 1.0 + KERNEL_EPSILON
                || y >= 3.0 - KERNEL_EPSILON
        }));
    }
}

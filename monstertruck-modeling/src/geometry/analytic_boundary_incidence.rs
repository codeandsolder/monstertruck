//! Validate all native analytic carrier/edge incidences in a B-rep solid.
//!
//! This is *local* boundary evidence, not an equality certificate for
//! oriented occupied sets. A complete solid comparator must additionally
//! certify the same support-patch partitions, complete physical intersection
//! loci, loop coverage, seam identity and shell nesting.

use super::analytic_curve_incidence::certify_analytic_curve_on_surface;
use super::surface::analytic_carrier;
use super::{BoundedCurve, Curve, MetricSpace, Point3, Solid, Surface};
use crate::Wire;
use monstertruck_geometry::analytic_carrier::{AnalyticCarrier, CarrierDomainRelation};
use monstertruck_geometry::analytic_curve_incidence::CurveIncidence;
use monstertruck_topology::{Edge, EdgeId, Face, Shell};
use rustc_hash::FxHashMap as HashMap;

/// A specific reason why a native B-rep's complete local incidence could
/// not be proven. Neither invalidity nor lack of a certificate proves that
/// two occupied sets differ.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum BoundaryIncidenceFailure {
    InvalidTolerance,
    InvalidShell,
    InvalidFaceWire,
    UnsupportedFaceCarrier,
    UnsupportedEdgeCurve,
    CurveNotOnFaceCarrier,
    NonFiniteEdgeEndpoint,
    CurveEndpointMismatch,
    InvalidEdgeUseCount,
    InconsistentEdgeOrientation,
    OpposingCoincidentFaces,
}

/// Location of a failed proof. Face/edge indices are local to the shell; an
/// edge index refers to its position in a face's oriented wire traversal.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct BoundaryIncidenceError {
    pub shell: usize,
    pub face: Option<usize>,
    pub edge: Option<usize>,
    pub reason: BoundaryIncidenceFailure,
}

/// A face's globally certified carrier and its B-rep orientation. Both are
/// necessary to derive the physical oriented surface normal.
#[derive(Clone, Copy, Debug)]
pub struct NativeFaceIncidence {
    pub carrier: AnalyticCarrier,
    pub face_orientation: bool,
}

/// A directed edge use in one face's ordered boundary traversal.
///
/// Uses the face and edge ordinals, not process-local pointer identity.
/// The positions remain stable for the same stored B-rep topology.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub struct NativeEdgeUse {
    pub face_index: usize,
    pub face_edge_index: usize,
}

/// One topological edge's two boundary uses and their support relation.
///
/// The indices refer to `NativeShellIncidence::faces`. The edge's geometry
/// has been checked against BOTH face carriers, including all nonzero spline
/// knot intervals and the topological endpoints. This record does not itself
/// prove the edge is the unique complete curve of carrier intersection.
#[derive(Clone, Copy, Debug)]
pub struct NativeEdgeIncidence {
    /// Both oriented uses in ascending face/edge traversal order.
    pub uses: [NativeEdgeUse; 2],
    /// True for one edge used in opposite directions by a single periodic face.
    pub seam: bool,
    /// True when adjacent face normals and geometric supports coincide.
    pub same_oriented_support: bool,
}

/// Complete native *incidence* evidence for one manifold boundary shell.
#[derive(Clone, Debug)]
pub struct NativeShellIncidence {
    pub faces: Vec<NativeFaceIncidence>,
    pub edges: Vec<NativeEdgeIncidence>,
}

/// Per-shell native incidence evidence. Zero shells is a valid empty B-rep;
/// this is not a certificate that another solid has the same occupied set.
#[derive(Clone, Debug)]
pub struct NativeAnalyticBoundaryIncidence {
    pub shells: Vec<NativeShellIncidence>,
}

#[derive(Clone, Copy)]
struct EdgeUse {
    face: usize,
    edge_index: usize,
    orientation: bool,
}

impl EdgeUse {
    const fn position(self) -> NativeEdgeUse {
        NativeEdgeUse {
            face_index: self.face,
            face_edge_index: self.edge_index,
        }
    }
}

const fn failure(
    shell: usize,
    face: Option<usize>,
    edge: Option<usize>,
    reason: BoundaryIncidenceFailure,
) -> BoundaryIncidenceError {
    BoundaryIncidenceError {
        shell,
        face,
        edge,
        reason,
    }
}

fn finite(point: Point3) -> bool { [point.x, point.y, point.z].into_iter().all(f64::is_finite) }

fn certify_edge(
    edge: &Edge<Point3, Curve>,
    surface: &Surface,
    tolerance: f64,
    location: (usize, usize, usize),
) -> Result<(), BoundaryIncidenceError> {
    let (shell, face, index) = location;
    let location_error = |reason| failure(shell, Some(face), Some(index), reason);
    let front = edge.front().point();
    let back = edge.back().point();
    if !finite(front) || !finite(back) {
        return Err(location_error(
            BoundaryIncidenceFailure::NonFiniteEdgeEndpoint,
        ));
    }
    let curve = edge.oriented_curve();
    // Certify the curve representation before evaluating its endpoints:
    // unsupported composed/intersection curves may run a Newton projection
    // during evaluation, which is not admissible geometric evidence.
    match certify_analytic_curve_on_surface(&curve, surface, tolerance) {
        CurveIncidence::Proven => {}
        CurveIncidence::NotProven => {
            return Err(location_error(
                BoundaryIncidenceFailure::CurveNotOnFaceCarrier,
            ));
        }
        CurveIncidence::Unsupported => {
            return Err(location_error(
                BoundaryIncidenceFailure::UnsupportedEdgeCurve,
            ));
        }
    }
    let (start, end) = (curve.front(), curve.back());
    if !finite(start)
        || !finite(end)
        || front.distance(start) > tolerance
        || back.distance(end) > tolerance
    {
        return Err(location_error(
            BoundaryIncidenceFailure::CurveEndpointMismatch,
        ));
    }
    Ok(())
}

fn certify_face(
    face: &Face<Point3, Curve, Surface>,
    shell_index: usize,
    face_index: usize,
    tolerance: f64,
    edges: &mut HashMap<EdgeId<Curve>, Vec<EdgeUse>>,
) -> Result<NativeFaceIncidence, BoundaryIncidenceError> {
    let surface = face.surface();
    let Some(carrier) = analytic_carrier(&surface) else {
        return Err(failure(
            shell_index,
            Some(face_index),
            None,
            BoundaryIncidenceFailure::UnsupportedFaceCarrier,
        ));
    };
    let wires = face.boundaries();
    if wires.is_empty() || !Wire::disjoint_wires(&wires) {
        return Err(failure(
            shell_index,
            Some(face_index),
            None,
            BoundaryIncidenceFailure::InvalidFaceWire,
        ));
    }
    let mut edge_index = 0;
    for wire in &wires {
        if wire.is_empty() || !wire.is_closed() || !wire.is_simple() {
            return Err(failure(
                shell_index,
                Some(face_index),
                None,
                BoundaryIncidenceFailure::InvalidFaceWire,
            ));
        }
        for edge in wire.edge_iter() {
            certify_edge(
                edge,
                &surface,
                tolerance,
                (shell_index, face_index, edge_index),
            )?;
            edges.entry(edge.id()).or_default().push(EdgeUse {
                face: face_index,
                edge_index,
                orientation: edge.orientation(),
            });
            edge_index += 1;
        }
    }
    Ok(NativeFaceIncidence {
        carrier,
        face_orientation: face.orientation(),
    })
}

fn shared_edge_relation(
    first: NativeFaceIncidence,
    second: NativeFaceIncidence,
    linear_tolerance: f64,
    angular_tolerance: f64,
) -> Result<bool, BoundaryIncidenceFailure> {
    let relation = first
        .carrier
        .same_domain(second.carrier, linear_tolerance, angular_tolerance);
    let same_normal = match relation {
        None => return Ok(false),
        Some(CarrierDomainRelation::SameOrientation) => {
            first.face_orientation == second.face_orientation
        }
        Some(CarrierDomainRelation::ReversedOrientation) => {
            first.face_orientation != second.face_orientation
        }
    };
    if same_normal {
        Ok(true)
    } else {
        Err(BoundaryIncidenceFailure::OpposingCoincidentFaces)
    }
}

fn collect_edge_incidence(
    shell_index: usize,
    faces: &[NativeFaceIncidence],
    edges: HashMap<EdgeId<Curve>, Vec<EdgeUse>>,
    linear_tolerance: f64,
    angular_tolerance: f64,
) -> Result<Vec<NativeEdgeIncidence>, BoundaryIncidenceError> {
    let mut results = Vec::with_capacity(edges.len());
    // FxHashMap is deliberately not an ordering source for certificates or
    // first-failure diagnostics. Order edge groups by their first B-rep use.
    let mut grouped = edges.into_values().collect::<Vec<_>>();
    grouped.sort_unstable_by_key(|uses| uses.iter().map(|edge_use| edge_use.position()).min());
    for pair in grouped {
        let [first, second] = pair.as_slice() else {
            return Err(failure(
                shell_index,
                None,
                None,
                BoundaryIncidenceFailure::InvalidEdgeUseCount,
            ));
        };
        if first.orientation == second.orientation {
            return Err(failure(
                shell_index,
                Some(first.face),
                None,
                BoundaryIncidenceFailure::InconsistentEdgeOrientation,
            ));
        }
        let same = shared_edge_relation(
            faces[first.face],
            faces[second.face],
            linear_tolerance,
            angular_tolerance,
        )
        .map_err(|reason| failure(shell_index, Some(first.face), None, reason))?;
        let mut uses = [first.position(), second.position()];
        uses.sort_unstable();
        results.push(NativeEdgeIncidence {
            uses,
            seam: first.face == second.face,
            same_oriented_support: same,
        });
    }
    Ok(results)
}

fn certify_shell(
    shell: &Shell<Point3, Curve, Surface>,
    shell_index: usize,
    linear_tolerance: f64,
    angular_tolerance: f64,
) -> Result<NativeShellIncidence, BoundaryIncidenceError> {
    shell.check_solid_boundary().map_err(|_| {
        failure(
            shell_index,
            None,
            None,
            BoundaryIncidenceFailure::InvalidShell,
        )
    })?;
    let mut edges = HashMap::<EdgeId<Curve>, Vec<EdgeUse>>::default();
    let mut faces = Vec::with_capacity(shell.len());
    for (face_index, face) in shell.face_iter().enumerate() {
        faces.push(certify_face(
            face,
            shell_index,
            face_index,
            linear_tolerance,
            &mut edges,
        )?);
    }
    let edges = collect_edge_incidence(
        shell_index,
        &faces,
        edges,
        linear_tolerance,
        angular_tolerance,
    )?;
    Ok(NativeShellIncidence { faces, edges })
}

/// Certifies complete native local B-rep topology/geometry incidence.
///
/// Every boundary shell must pass the kernel's manifold validation; every
/// oriented face must have a structurally certified analytic carrier; every
/// wire must be closed and simple; every edge curve must have matched finite
/// endpoints and a global nonsampling carrier proof on BOTH adjacent faces;
/// every topological edge must have exactly two opposing uses. Same-support
/// adjacencies and periodic self-seams are recorded, not discarded.
///
/// This does not prove that nonadjacent face interiors do not intersect,
/// that all intersections are represented, that every loop is correctly
/// nested, or that two solids enclose the same volume. The forthcoming
/// global boundary certificate must establish those claims independently.
///
/// # Errors
/// Returns the first typed topology, curve or carrier refusal. In particular,
/// a curve with only sampled coincidence never counts as an incident edge.
pub fn certify_native_analytic_boundary_incidence(
    solid: &Solid,
    linear_tolerance: f64,
    angular_tolerance: f64,
) -> Result<NativeAnalyticBoundaryIncidence, BoundaryIncidenceError> {
    if !linear_tolerance.is_finite()
        || linear_tolerance <= 0.0
        || !angular_tolerance.is_finite()
        || angular_tolerance <= 0.0
    {
        return Err(failure(
            0,
            None,
            None,
            BoundaryIncidenceFailure::InvalidTolerance,
        ));
    }
    let shells = solid
        .boundaries()
        .iter()
        .enumerate()
        .map(|(index, shell)| certify_shell(shell, index, linear_tolerance, angular_tolerance))
        .collect::<Result<Vec<_>, _>>()?;
    Ok(NativeAnalyticBoundaryIncidence { shells })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{BsplineCurve, KnotVector, Line, Vector3, builder};

    fn cube() -> Solid {
        let vertex = builder::vertex(Point3::new(-0.5, -0.5, -0.5));
        let edge = builder::extrude(&vertex, Vector3::unit_x());
        let face = builder::extrude(&edge, Vector3::unit_y());
        builder::extrude(&face, Vector3::unit_z())
    }

    const LINEAR: f64 = 1.0e-7;
    const ANGULAR: f64 = 1.0e-9;

    #[test]
    fn cube_proves_all_six_faces_and_twelve_incident_edges() {
        let proof = certify_native_analytic_boundary_incidence(&cube(), LINEAR, ANGULAR)
            .expect("regular cube should certify");
        assert_eq!(proof.shells.len(), 1);
        assert_eq!(proof.shells[0].faces.len(), 6);
        assert_eq!(proof.shells[0].edges.len(), 12);
        assert!(
            proof.shells[0]
                .edges
                .iter()
                .all(|e| !e.seam && !e.same_oriented_support)
        );
    }

    #[test]
    fn empty_solid_has_valid_vacuous_incidence() {
        let proof = certify_native_analytic_boundary_incidence(&Solid::empty(), LINEAR, ANGULAR)
            .expect("empty native solid");
        assert!(proof.shells.is_empty());
    }

    #[test]
    fn a_topologically_valid_bumped_curve_cannot_claim_incidence() {
        let solid = cube();
        let edge = solid.edge_iter().next().expect("cube has edges");
        let start = edge.front().point();
        let end = edge.back().point();
        let bump = Vector3::new(0.1, 0.1, 0.1);
        edge.set_curve(Curve::BsplineCurve(BsplineCurve::new_unchecked(
            KnotVector::bezier_knot(3),
            vec![start, start + bump, end + bump, end],
        )));
        let refused = certify_native_analytic_boundary_incidence(&solid, LINEAR, ANGULAR)
            .expect_err("between-endpoint excursion must refuse");
        assert_eq!(
            refused.reason,
            BoundaryIncidenceFailure::CurveNotOnFaceCarrier
        );
    }

    #[test]
    fn endpoint_mismatch_is_not_masked_by_surface_incidence() {
        let solid = cube();
        let edge = solid.edge_iter().next().expect("cube has edges");
        let front = edge.front().point();
        edge.set_curve(Curve::Line(Line(front, front)));
        let refused = certify_native_analytic_boundary_incidence(&solid, LINEAR, ANGULAR)
            .expect_err("changed topological endpoint must refuse");
        assert_eq!(
            refused.reason,
            BoundaryIncidenceFailure::CurveEndpointMismatch
        );
    }

    #[test]
    fn reparameterized_collinear_edge_remains_incident() {
        let solid = cube();
        let edge = solid.edge_iter().next().expect("cube has edges");
        let start = edge.front().point();
        let end = edge.back().point();
        let delta = end - start;
        edge.set_curve(Curve::BsplineCurve(BsplineCurve::new_unchecked(
            KnotVector::bezier_knot(3),
            vec![start, start + delta * 0.2, start + delta * 0.8, end],
        )));
        let proof = certify_native_analytic_boundary_incidence(&solid, LINEAR, ANGULAR)
            .expect("a polynomial line cannot lose incidence through reparameterization");
        assert_eq!(proof.shells[0].edges.len(), 12);
    }

    #[test]
    fn flipped_face_breaks_manifold_shell_before_geometry_proof() {
        let mut shells = cube().into_boundaries();
        shells[0][0].invert();
        let invalid = Solid::new_unchecked(shells);
        let refused = certify_native_analytic_boundary_incidence(&invalid, LINEAR, ANGULAR)
            .expect_err("one flipped face breaks the oriented shell");
        assert_eq!(refused.reason, BoundaryIncidenceFailure::InvalidShell);
    }

    #[test]
    fn output_edge_use_order_is_stable_across_new_topology_instances() {
        let first = certify_native_analytic_boundary_incidence(&cube(), LINEAR, ANGULAR)
            .expect("first cube");
        let second = certify_native_analytic_boundary_incidence(&cube(), LINEAR, ANGULAR)
            .expect("second independent cube");
        let first_edges = &first.shells[0].edges;
        let second_edges = &second.shells[0].edges;
        assert_eq!(
            first_edges.iter().map(|e| e.uses).collect::<Vec<_>>(),
            second_edges.iter().map(|e| e.uses).collect::<Vec<_>>(),
        );
        assert!(
            first_edges
                .windows(2)
                .all(|pair| pair[0].uses[0] < pair[1].uses[0])
        );
        for edge in first_edges {
            assert!(edge.uses[0] < edge.uses[1]);
        }
    }

    #[test]
    fn globally_reversed_shell_keeps_incidence_but_not_face_orientation() {
        let ordinary = cube();
        let original = certify_native_analytic_boundary_incidence(&ordinary, LINEAR, ANGULAR)
            .expect("regular cube");
        let mut reverse = cube();
        reverse.not();
        let reversed = certify_native_analytic_boundary_incidence(&reverse, LINEAR, ANGULAR)
            .expect("reversed complete shell still has local incidence");
        assert_eq!(
            original.shells[0].edges.len(),
            reversed.shells[0].edges.len()
        );
        for (first, second) in original.shells[0]
            .faces
            .iter()
            .zip(&reversed.shells[0].faces)
        {
            assert_ne!(first.face_orientation, second.face_orientation);
        }
    }

    #[test]
    fn malformed_shell_is_explicitly_rejected() {
        let invalid = Solid::new_unchecked(vec![Shell::new()]);
        let refused = certify_native_analytic_boundary_incidence(&invalid, LINEAR, ANGULAR)
            .expect_err("empty boundary shell invalid");
        assert_eq!(refused.reason, BoundaryIncidenceFailure::InvalidShell);
    }

    #[test]
    fn invalid_tolerances_refuse() {
        let cube = cube();
        for (linear, angular) in [
            (f64::NAN, ANGULAR),
            (0.0, ANGULAR),
            (LINEAR, 0.0),
            (LINEAR, f64::INFINITY),
        ] {
            let refused = certify_native_analytic_boundary_incidence(&cube, linear, angular)
                .expect_err("invalid tolerances");
            assert_eq!(refused.reason, BoundaryIncidenceFailure::InvalidTolerance);
        }
    }
}

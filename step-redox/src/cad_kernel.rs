use crate::cad_ir::{CadModel, NodeId};
use anyhow::Result;
use monstertruck_core::tolerance::DEFAULT_STEP_RECOVERY_TOLERANCES;

/// User-facing geometry accuracy contract at the unit-aware Step Redox boundary.
///
/// Step Redox CAD IR stores physical lengths in millimetres, so Standard policy
/// is expressed in millimetres here and translated to the unitless kernel's
/// current model units only inside a backend adapter.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum GeometryPolicy {
    /// Geometry-changing operations must establish equivalence or refuse.
    Equivalent,
    /// Bounded geometric change under an explicit millimetre budget.
    Standard { tolerance_mm: f64 },
}

impl GeometryPolicy {
    /// Project-wide Standard preset: 0.01 mm.
    pub const PROJECT_STANDARD_TOLERANCE_MM: f64 =
        DEFAULT_STEP_RECOVERY_TOLERANCES.cad_project_standard_mm;

    #[must_use]
    pub const fn standard_mm(tolerance_mm: f64) -> Self { Self::Standard { tolerance_mm } }

    #[must_use]
    pub const fn project_standard() -> Self {
        Self::standard_mm(Self::PROJECT_STANDARD_TOLERANCE_MM)
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct KernelSummary {
    pub geometrically_consistent: bool,
    pub shells: usize,
    pub faces: usize,
}

impl KernelSummary {
    /// One non-empty, geometrically valid closed boundary component.
    ///
    /// The evaluated value already owns a single native Solid; re-exporting
    /// it to STEP and reparsing entities adds no geometric information.
    #[must_use]
    pub const fn is_single_solid(self) -> bool {
        self.geometrically_consistent && self.shells == 1 && self.faces != 0
    }
}

/// Narrow execution boundary for constructive CAD.
///
/// The associated evaluated value is intentionally backend-defined and opaque to
/// recovery code, but always represents one native solid (which may contain
/// multiple boundary shells). Its in-memory summary is authoritative for shell
/// count and topology validity. Serialization is requested only when an external
/// STEP representation or independent boundary comparison is actually needed.
pub trait CadKernel {
    type Evaluated;

    /// Evaluate one canonical CAD root in the backend.
    ///
    /// # Errors
    /// Returns an error when the model cannot be represented or evaluated by this kernel.
    fn evaluate(&self, model: &CadModel, root: NodeId) -> Result<Self::Evaluated>;

    fn summarize(&self, evaluated: &Self::Evaluated) -> KernelSummary;

    /// Whether this kernel supports exporting evaluated CAD as STEP.
    ///
    /// Native topology editing and replay need no serialization. Independent
    /// cross-STEP source-boundary admission requires it, and is skipped
    /// entirely if the adapter cannot provide it.
    fn supports_step_export(&self) -> bool;

    /// Serialize an evaluated shape as a STEP exchange.
    ///
    /// # Errors
    /// Returns an error when the evaluated topology cannot be serialized as STEP.
    fn to_step(&self, evaluated: &Self::Evaluated) -> Result<String>;
}

#[cfg(feature = "cad-kernel-monstertruck")]
pub mod monstertruck {
    use super::{CadKernel, GeometryPolicy, KernelSummary};
    use crate::cad_ir::{
        Axis3, BooleanOp, CadModel, CadNode, Curve2d, EdgeSelector, EdgeTreatmentEdge,
        EdgeTreatmentProfile, LoftSection, NodeId, Primitive, Profile2d, RigidTransform,
        SweepPath3d, SweepSegment3d,
    };
    use anyhow::{Context, Result, bail, ensure};
    use monstertruck_core::tolerance::{
        DEFAULT_NUMERICAL_TOLERANCES, DEFAULT_STEP_RECOVERY_TOLERANCES,
    };
    use monstertruck_fillet::{FilletOptions, FilletProfile, RadiusSpec};
    use monstertruck_geometry::analytic_curve_incidence::CurveIncidence;
    use monstertruck_io::step::save::{self, CompleteStepDisplay};
    use monstertruck_modeling::{
        BoundedCurve, BsplineCurve, Curve, Edge, EuclideanSpace, InnerSpace, KnotVector, Matrix4,
        NurbsCurve, ParametricCurve, ParametricSurface, ParametricSurface3D, Point3, Rad,
        SearchNearestParameter, SearchParameter, Shell, Solid, Surface, SurfaceCurve, Vector3,
        Vector4, Wire, algo, builder, certify_analytic_curve_on_surface, primitive, profile,
    };
    use std::collections::{HashMap, HashSet};

    mod circular_selection;
    mod planar_loft;

    #[derive(Debug, Clone, Copy)]
    pub struct MonstertruckKernel {
        policy: GeometryPolicy,
    }

    impl MonstertruckKernel {
        #[must_use]
        pub const fn new(policy: GeometryPolicy) -> Self { Self { policy } }

        #[must_use]
        pub const fn equivalent() -> Self { Self::new(GeometryPolicy::Equivalent) }

        #[must_use]
        pub const fn standard_mm(tolerance_mm: f64) -> Self {
            Self::new(GeometryPolicy::standard_mm(tolerance_mm))
        }

        #[must_use]
        pub const fn project_standard() -> Self { Self::new(GeometryPolicy::project_standard()) }
    }

    /// STEP is a sampled representation of a non-analytic surface intersection.
    /// The numerical export precision is separate from the modeling policy:
    /// Standard provides a physical upper bound; Equivalent uses the kernel's
    /// numerical representation threshold, without altering any CAD geometry.
    fn step_intersection_accuracy(policy: GeometryPolicy) -> Result<f64> {
        let tolerance = match policy {
            GeometryPolicy::Equivalent => DEFAULT_NUMERICAL_TOLERANCES.kernel_f64,
            GeometryPolicy::Standard { tolerance_mm } => {
                ensure!(
                    tolerance_mm.is_finite()
                        && tolerance_mm >= DEFAULT_NUMERICAL_TOLERANCES.kernel_f64,
                    "STEP modeling tolerance must be finite and >= {} mm; requested {tolerance_mm}",
                    DEFAULT_NUMERICAL_TOLERANCES.kernel_f64
                );
                // The modeling budget limits acceptable changes to CAD geometry;
                // it does not need to coarsen the independent STEP serialization
                // of already-evaluated intersection curves. Smaller user-requested
                // tolerances remain authoritative.
                tolerance_mm.min(DEFAULT_STEP_RECOVERY_TOLERANCES.cad_step_export_accuracy_mm)
            }
        };
        ensure!(
            tolerance.is_finite() && tolerance >= DEFAULT_NUMERICAL_TOLERANCES.kernel_f64,
            "STEP intersection accuracy must be finite and >= {} mm; requested {tolerance}",
            DEFAULT_NUMERICAL_TOLERANCES.kernel_f64,
        );
        Ok(tolerance)
    }

    fn polyline_closest_distance(point: Point3, polyline: &[Point3]) -> Option<f64> {
        polyline
            .windows(2)
            .map(|segment| {
                let [a, b] = [segment[0], segment[1]];
                let ab = b - a;
                let squared = ab.magnitude2();
                let t = if squared > 0.0 {
                    ((point - a).dot(ab) / squared).clamp(0.0, 1.0)
                } else {
                    0.0
                };
                (point - (a + ab * t)).magnitude()
            })
            .reduce(f64::min)
    }

    /// Evaluating an `IntersectionCurve` can fall back to its marching leader
    /// if the two-surface solver cannot find a corrected common point.
    /// Certify every proposed STEP guide vertex against *both* carrier
    /// surfaces rather than serializing potentially unsupported fallback data.
    fn validated_step_surface_support(
        surface: &Surface,
        points: &[Point3],
        tolerance: f64,
    ) -> bool {
        let mut hint = None;
        points.iter().copied().all(|point| {
            let uv = surface
                .search_parameter(point, hint, 100)
                .or_else(|| surface.search_nearest_parameter(point, hint, 100))
                .or_else(|| surface.search_nearest_parameter(point, None, 100));
            let Some((u, v)) = uv else { return false };
            if !u.is_finite()
                || !v.is_finite()
                || (surface.evaluate(u, v) - point).magnitude() > tolerance
            {
                return false;
            }
            hint = Some((u, v));
            true
        })
    }

    /// A STEP PCURVE has its *own* chart and parameter range, unrelated to
    /// the 3-D curve parameter values. Validate the geometric locus in both
    /// directions; never compare equal normalized t-values.
    fn validated_step_parameter_curve(
        boundary: &monstertruck_modeling::ParameterCurve<
            monstertruck_modeling::Curve2D,
            Box<Surface>,
        >,
        intersection_points: &[Point3],
        tolerance: f64,
    ) -> bool {
        let (_, points) =
            algo::curve::parameter_division(boundary, boundary.range_tuple(), tolerance);
        !points.is_empty()
            && points.iter().all(|point| {
                point.x.is_finite()
                    && point.y.is_finite()
                    && point.z.is_finite()
                    && polyline_closest_distance(*point, intersection_points)
                        .is_some_and(|dist| dist <= tolerance)
            })
            && intersection_points.iter().all(|point| {
                polyline_closest_distance(*point, &points).is_some_and(|dist| dist <= tolerance)
            })
    }

    /// A degree-one two-control-point rational spline has an exact straight
    /// point-set locus (including nonuniform positive homogeneous weights).
    /// Keep that structural information instead of introducing a spurious
    /// transverse component by projecting its endpoints independently.
    fn straight_step_guide(guide: &Curve) -> bool {
        match guide {
            Curve::Line(_) => true,
            Curve::BsplineCurve(spline) => {
                spline.degree() == 1 && spline.control_points().len() == 2
            }
            Curve::NurbsCurve(spline) => {
                spline.degree() == 1
                    && spline.control_points().len() == 2
                    && spline
                        .control_points()
                        .iter()
                        .all(|point| point.w.is_finite() && point.w > 0.0)
            }
            Curve::ParameterCurve(_) | Curve::Intersection(_) => false,
        }
    }

    /// Preserve a source NURBS guide only with independent whole-curve
    /// Bernstein incidence certificates for both analytic carriers. For straight
    /// guides, the exact point-set remains certified by their simple structure.
    /// The evaluated SSI still controls both the endpoints and acceptance.
    fn certified_step_guide(
        guide: &Curve,
        surfaces: [&Surface; 2],
        parameters: &[f64],
        evaluated: &[Point3],
        tolerance: f64,
    ) -> Option<Vec<Point3>> {
        let proof_tolerance = DEFAULT_STEP_RECOVERY_TOLERANCES
            .solid_geometry_mm
            .min(tolerance);
        let globally_proven = matches!(guide, Curve::NurbsCurve(_))
            && surfaces.iter().all(|&carrier| {
                certify_analytic_curve_on_surface(guide, carrier, proof_tolerance)
                    == CurveIncidence::Proven
            });
        if !straight_step_guide(guide) && !globally_proven {
            return None;
        }
        let guide_points = parameters
            .iter()
            .map(|&parameter| guide.evaluate(parameter))
            .collect::<Vec<_>>();
        let accuracy = if globally_proven {
            proof_tolerance
        } else {
            tolerance
        };
        (guide_points.len() == evaluated.len()
            && guide_points
                .iter()
                .zip(evaluated)
                .all(|(guide, point)| (*guide - *point).magnitude() <= accuracy)
            && surfaces
                .iter()
                .all(|carrier| validated_step_surface_support(carrier, &guide_points, accuracy)))
        .then_some(guide_points)
    }

    /// Construct a STEP-exchange guide from the *evaluated* two-surface
    /// intersection, not the marching leader used to find it. This is export
    /// only: replacing Boolean/SSI sampling globally breaks interface
    /// construction in real CAD models.
    fn corrected_step_curve(curve: &Curve, tolerance: f64) -> Result<Curve> {
        let Curve::Intersection(intersection) = curve else {
            return Ok(curve.clone());
        };
        let (start, end) = intersection.range_tuple();
        ensure!(
            start.is_finite() && end.is_finite() && end > start,
            "invalid STEP intersection parameter range {start}..{end}"
        );
        let (parameters, points) =
            algo::curve::parameter_division(intersection, (start, end), tolerance);
        ensure!(
            points.len() >= 2,
            "STEP intersection has fewer than two vertices"
        );
        ensure!(
            points
                .iter()
                .all(|point| point.x.is_finite() && point.y.is_finite() && point.z.is_finite()),
            "non-finite evaluated STEP intersection guide"
        );
        ensure!(
            validated_step_surface_support(intersection.surface0(), &points, tolerance)
                && validated_step_surface_support(intersection.surface1(), &points, tolerance),
            "evaluated STEP intersection guide does not lie on both carriers within {tolerance} mm"
        );
        for (index, boundary) in [intersection.boundary0(), intersection.boundary1()]
            .into_iter()
            .enumerate()
        {
            if let Some(boundary) = boundary {
                ensure!(
                    validated_step_parameter_curve(boundary, &points, tolerance),
                    "STEP intersection PCURVE {index} does not match evaluated 3-D curve within {tolerance} mm"
                );
            }
        }
        if let Some(guide_points) = certified_step_guide(
            intersection.leader(),
            [intersection.surface0(), intersection.surface1()],
            &parameters,
            &points,
            tolerance,
        ) && [intersection.boundary0(), intersection.boundary1()]
            .into_iter()
            .flatten()
            .all(|boundary| validated_step_parameter_curve(boundary, &guide_points, tolerance))
        {
            return Ok(curve.clone());
        }
        // Preserve each native parameter, including nonuniform adaptive
        // subdivisions and shifted input domains. A uniform 0..N knot
        // vector changes the curve's parameterization after STEP export.
        ensure!(
            parameters.len() == points.len()
                && parameters.first().copied() == Some(start)
                && parameters.last().copied() == Some(end)
                && parameters
                    .windows(2)
                    .all(|pair| pair[0].is_finite() && pair[1] > pair[0]),
            "STEP intersection sampling must cover a strictly increasing native parameter range"
        );
        let knots = KnotVector::from(
            std::iter::once(start)
                .chain(parameters.iter().copied())
                .chain(std::iter::once(end))
                .collect::<Vec<_>>(),
        );
        let leader = Curve::BsplineCurve(BsplineCurve::new_unchecked(knots, points));
        let mut exported = SurfaceCurve::with_boundaries(
            intersection.surface0().clone(),
            intersection.surface1().clone(),
            Box::new(leader),
            intersection.boundary0().cloned(),
            intersection.boundary1().cloned(),
        );
        let (front, back) = intersection.endpoint_anchors();
        exported.set_endpoint_anchors(front, back);
        Ok(Curve::Intersection(exported))
    }

    /// No topology, edge IDs or evaluated CAD curves are modified during
    /// exchange. A malformed 3-D/PCURVE association is an export error, not
    /// permission to serialize a misleading curve.
    fn step_export_solid(source: &Solid, policy: GeometryPolicy) -> Result<Solid> {
        let tolerance = step_intersection_accuracy(policy)?;
        let first_error = std::cell::RefCell::new(None);
        let exported = source.try_mapped(
            |point| Some(*point),
            |curve| match corrected_step_curve(curve, tolerance) {
                Ok(converted) => Some(converted),
                Err(error) => {
                    if first_error.borrow().is_none() {
                        *first_error.borrow_mut() = Some(error);
                    }
                    None
                }
            },
            |surface| Some(surface.clone()),
        );
        if let Some(error) = first_error.into_inner() {
            return Err(error).context("STEP intersection geometry conversion failed");
        }
        exported.context("STEP export could not reconstruct original solid topology")
    }

    /// Opaque wrapper: Monstertruck topology does not escape the backend module.
    pub struct EvaluatedShape {
        solid: Solid,
    }

    // Evaluating a node under an override must not reuse an object evaluated
    // at a different numerical policy.
    #[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
    enum PolicyKey {
        Equivalent,
        Standard(u64),
    }

    impl From<GeometryPolicy> for PolicyKey {
        fn from(policy: GeometryPolicy) -> Self {
            match policy {
                GeometryPolicy::Equivalent => Self::Equivalent,
                GeometryPolicy::Standard { tolerance_mm } => Self::Standard(tolerance_mm.to_bits()),
            }
        }
    }

    type EvaluationKey = (NodeId, PolicyKey);

    /// A reusable session retains independent snapshots for future roots.
    /// A one-shot evaluator retains only shared dependencies that still have
    /// a possible future use, and moves the owned snapshot on the final use.
    #[derive(Default)]
    struct EvaluationCache {
        solids: HashMap<EvaluationKey, Solid>,
        remaining: Option<HashMap<EvaluationKey, usize>>,
    }

    impl EvaluationCache {
        fn single_root(model: &CadModel, root: NodeId, policy: GeometryPolicy) -> Result<Self> {
            let mut remaining = HashMap::new();
            let mut pending = vec![(root, policy)];
            while let Some((id, policy)) = pending.pop() {
                let uses = remaining.entry((id, PolicyKey::from(policy))).or_insert(0);
                *uses += 1;
                // The first use builds this node; later uses read its
                // independent snapshot. Count its evaluated dependencies once.
                if *uses > 1 {
                    continue;
                }
                match model.node(id)? {
                    CadNode::SolverTolerance {
                        child,
                        tolerance_mm,
                    } => {
                        pending.push((*child, GeometryPolicy::standard_mm(*tolerance_mm)));
                    }
                    node => {
                        pending.extend(node.children().into_iter().map(|child| (child, policy)));
                    }
                }
            }
            Ok(Self {
                solids: HashMap::new(),
                remaining: Some(remaining),
            })
        }

        /// Decrement the bounded use count before evaluating the key. Native
        /// fused fast paths may bypass some dependencies, so counts can
        /// conservatively overestimate but never underestimate actual calls.
        fn begin(&mut self, key: EvaluationKey) -> Result<(Option<Solid>, bool)> {
            let retain = if let Some(remaining) = self.remaining.as_mut() {
                let uses = remaining
                    .get_mut(&key)
                    .context("unplanned CAD evaluation dependency")?;
                *uses = uses
                    .checked_sub(1)
                    .context("CAD evaluation dependency used too often")?;
                *uses != 0
            } else {
                true
            };
            let value = if retain {
                self.solids
                    .get(&key)
                    .map(|solid| solid.mapped(Clone::clone, Clone::clone, Clone::clone))
            } else {
                // Last consumer owns the cached mutable geometry exclusively.
                self.solids.remove(&key)
            };
            Ok((value, retain))
        }

        fn remember(&mut self, key: EvaluationKey, solid: &Solid, retain: bool) {
            if retain {
                self.solids
                    .insert(key, solid.mapped(Clone::clone, Clone::clone, Clone::clone));
            }
        }
    }

    /// Evaluate several CAD roots without recomputing common geometry.
    ///
    /// Cached topology and returned shapes own independent mutable geometry;
    /// a plain `Solid::clone()` would share the underlying
    /// `Arc<Mutex<_>>` vertices, curves, and surfaces across evaluations.
    pub struct EvaluationSession<'a> {
        model: &'a CadModel,
        policy: GeometryPolicy,
        cache: EvaluationCache,
    }

    impl EvaluationSession<'_> {
        /// Evaluate a root, memoizing successful DAG nodes by ID and policy.
        ///
        /// # Errors
        /// Returns any kernel evaluation error; failed results are not cached.
        pub fn evaluate(&mut self, root: NodeId) -> Result<EvaluatedShape> {
            Ok(EvaluatedShape {
                solid: evaluate_node(self.model, root, self.policy, &mut self.cache)?,
            })
        }
    }

    impl MonstertruckKernel {
        /// Validate a model once and create an incremental evaluation session.
        ///
        /// # Errors
        /// Returns a validation error for malformed CAD IR.
        pub fn session<'a>(&self, model: &'a CadModel) -> Result<EvaluationSession<'a>> {
            model.validate()?;
            Ok(EvaluationSession {
                model,
                policy: self.policy,
                cache: EvaluationCache::default(),
            })
        }
    }

    fn profile_wires(profile: &Profile2d) -> Result<Vec<Wire>> {
        if profile.loops.is_empty() {
            bail!("profile must contain at least one loop");
        }
        profile
            .loops
            .iter()
            .map(|loop_| profile_loop_wire(&loop_.curves))
            .collect::<Result<Vec<_>>>()
    }

    fn extrude_profile_z(profile: &Profile2d, length_mm: f64) -> Result<Solid> {
        if !length_mm.is_finite() || length_mm == 0.0 {
            bail!("extrusion length must be finite and nonzero");
        }
        Ok(profile::solid_from_planar_profile(
            profile_wires(profile)?,
            Vector3::new(0.0, 0.0, length_mm),
        )?)
    }

    fn revolve_profile(profile: &Profile2d, axis: Axis3, angle_rad: f64) -> Result<Solid> {
        revolve_wires(profile_wires(profile)?, axis, angle_rad)
    }

    fn sweep_profile(profile: &Profile2d, path: &SweepPath3d) -> Result<Solid> {
        if path.segments.is_empty() {
            bail!("sweep path must contain at least one segment");
        }
        let face = profile::attach_plane_normalized(profile_wires(profile)?)?;
        let face = builder::transformed(&face, rigid_matrix(&path.initial_transform));
        let segments = path.segments.iter().map(|segment| match segment {
            SweepSegment3d::Translation { vector_mm } => {
                builder::CompositeSweepSegment::Translation(Vector3::new(
                    vector_mm[0],
                    vector_mm[1],
                    vector_mm[2],
                ))
            }
            SweepSegment3d::Rotation { axis, angle_rad } => {
                builder::CompositeSweepSegment::Rotation {
                    origin: Point3::new(axis.origin_mm[0], axis.origin_mm[1], axis.origin_mm[2]),
                    axis: Vector3::new(axis.direction[0], axis.direction[1], axis.direction[2]),
                    angle: Rad(*angle_rad),
                    division: 1,
                }
            }
        });
        if path.closed {
            Ok(builder::composite_closed_sweep(&face, segments)?)
        } else {
            Ok(builder::composite_sweep(&face, segments)?)
        }
    }

    /// The loft skin inherits each section wire's traversal direction.
    /// A reversed sketch produces inward walls even when the start/end caps
    /// are normalized separately. Orient the *wires themselves* along the
    /// actual section progression before building any common coedges.
    fn orient_loft_section_wires(wires: &mut [Wire]) -> Result<()> {
        let centroids = wires
            .iter()
            .map(|wire| {
                let vertices = wire
                    .vertex_iter()
                    .map(|vertex| vertex.point())
                    .collect::<Vec<_>>();
                if vertices.is_empty() {
                    bail!("loft section has no vertices");
                }
                let sum = vertices
                    .iter()
                    .fold(Vector3::new(0.0, 0.0, 0.0), |sum, point| {
                        sum + point.to_vec()
                    });
                Ok(Point3::from_vec(
                    sum / f64::from(u32::try_from(vertices.len())?),
                ))
            })
            .collect::<Result<Vec<_>>>()?;
        for index in 0..wires.len() {
            let before = centroids[index.saturating_sub(1)];
            let after = centroids[(index + 1).min(centroids.len() - 1)];
            let progression = after - before;
            if progression.magnitude() <= 0.0 {
                continue;
            }
            let face = builder::try_attach_plane::<Curve, Surface>(vec![wires[index].clone()])?;
            let normal = face.oriented_surface().normal(0.0, 0.0);
            let projected = normal.dot(progression);
            // A nearly tangential section does not offer reliable winding
            // evidence. Preserve it rather than inventing a face orientation.
            let threshold = DEFAULT_STEP_RECOVERY_TOLERANCES.source_cad_frame_orthogonality
                * normal.magnitude()
                * progression.magnitude();
            if projected < -threshold {
                wires[index].invert();
            }
        }
        Ok(())
    }

    fn loft_sections(sections: &[LoftSection]) -> Result<Solid> {
        if sections.len() < 2 {
            bail!("loft needs at least two sections");
        }
        if let Some(solid) = planar_loft::build_if_needed(sections)? {
            return Ok(solid);
        }
        let mut wires = sections
            .iter()
            .map(|section| {
                let mut wires = profile_wires(&section.profile)?;
                if wires.len() != 1 {
                    bail!("Monstertruck loft currently supports one loop per section");
                }
                let wire = wires
                    .pop()
                    .ok_or_else(|| anyhow::anyhow!("missing loft section wire"))?;
                Ok(builder::transformed(
                    &wire,
                    rigid_matrix(&section.transform),
                ))
            })
            .collect::<Result<Vec<_>>>()?;
        orient_loft_section_wires(&mut wires)?;
        let first = wires
            .first()
            .ok_or_else(|| anyhow::anyhow!("missing first loft wire"))?
            .clone();
        let last = wires
            .last()
            .ok_or_else(|| anyhow::anyhow!("missing last loft wire"))?
            .clone();
        let start_cap = profile::attach_plane_normalized(vec![first])?;
        let end_cap = profile::attach_plane_normalized(vec![last])?;
        // A ruled wall between straight source edges may be an exact plane.
        // Retain its analytic carrier and coedges when provably coplanar;
        // twisted and curved walls keep the original homotopy surface.
        let mut shell = Shell::new();
        for face in builder::try_skin_wires(&wires)? {
            let boundaries = face.boundaries();
            if let [wire] = boundaries.as_slice()
                && wire.len() == 4
                && wire
                    .edge_iter()
                    .all(|edge| matches!(edge.curve(), Curve::Line(_)))
                && let Ok(plane) = builder::try_attach_plane(boundaries)
            {
                shell.push(plane);
            } else {
                shell.push(face);
            }
        }
        shell.push(start_cap.inverse());
        shell.push(end_cap);
        Ok(Solid::new(vec![shell])?)
    }

    fn revolve_profile_transformed(
        profile: &Profile2d,
        axis: Axis3,
        angle_rad: f64,
        transform: &RigidTransform,
    ) -> Result<Solid> {
        let matrix = rigid_matrix(transform);
        let wires = profile_wires(profile)?
            .into_iter()
            .map(|wire| builder::transformed(&wire, matrix))
            .collect::<Vec<_>>();
        revolve_wires(wires, transform_axis(axis, transform), angle_rad)
    }

    /// A generated straight radial profile sweeps an affine planar cap.
    /// Retain its exact, already connected boundary curves, but give the
    /// newly built face a regular affine chart rather than the polar chart
    /// (which is singular at the rotation axis and unsuitable for overlay).
    /// Only native swept structural lines are considered; imported pcurves
    /// and general revolutions are never reparameterized here.
    fn affine_revolution_cap(
        face: &monstertruck_modeling::Face,
    ) -> Option<monstertruck_modeling::Face> {
        use monstertruck_modeling::{Line, Plane};
        let Surface::Revolution(surface) = face.oriented_surface() else {
            return None;
        };
        let revolution = surface.entity();
        let Curve::Line(Line(start, end)) = revolution.entity_curve() else {
            return None;
        };
        let generator = *end - *start;
        let axis = revolution.axis();
        let length = generator.magnitude();
        if !length.is_finite()
            || length <= DEFAULT_STEP_RECOVERY_TOLERANCES.cad_axis_locus_mm
            || generator.dot(axis).abs()
                > length * DEFAULT_STEP_RECOVERY_TOLERANCES.source_cad_frame_orthogonality
        {
            return None;
        }
        // The revolved profile must lie on one radial half-line. Its
        // supporting planar face is regular away from the optional pole.
        let radial = |p: Point3| {
            let relative = p - revolution.origin();
            relative - axis * relative.dot(axis)
        };
        let a = radial(*start);
        let b = radial(*end);
        if a.cross(b).magnitude() > length * DEFAULT_STEP_RECOVERY_TOLERANCES.cad_axis_locus_mm
            || a.dot(b) < -DEFAULT_STEP_RECOVERY_TOLERANCES.cad_axis_locus_mm.powi(2)
        {
            return None;
        }
        let (u_range, v_range) = surface.try_range_tuple();
        let (u0, u1) = u_range?;
        let (v0, v1) = v_range?;
        let u = f64::midpoint(u0, u1);
        let v = f64::midpoint(v0, v1);
        let point = surface.evaluate(u, v);
        let normal = surface.normal(u, v);
        let tangent = surface.derivative_u(u, v).normalize();
        let bitangent = normal.cross(tangent).normalize();
        if ![
            normal.x,
            normal.y,
            normal.z,
            tangent.x,
            tangent.y,
            tangent.z,
            bitangent.x,
            bitangent.y,
            bitangent.z,
        ]
        .into_iter()
        .all(f64::is_finite)
        {
            return None;
        }
        let plane = Plane::new(point, point + tangent, point + bitangent);
        let cap =
            monstertruck_modeling::Face::new(face.boundaries(), Surface::Plane(plane)).ok()?;
        cap.is_geometric_consistent().then_some(cap)
    }

    fn regularize_generated_revolution_caps(solid: Solid) -> Result<Solid> {
        let shells = solid
            .into_boundaries()
            .into_iter()
            .map(|shell| {
                shell
                    .into_iter()
                    .map(|face| affine_revolution_cap(&face).unwrap_or(face))
                    .collect::<Shell>()
            })
            .collect::<Vec<_>>();
        Ok(Solid::new(shells)?)
    }

    fn revolve_wires(mut wires: Vec<Wire>, axis: Axis3, angle_rad: f64) -> Result<Solid> {
        if !angle_rad.is_finite() || angle_rad == 0.0 {
            bail!("revolution angle must be finite and nonzero");
        }
        if angle_rad.abs()
            > std::f64::consts::TAU
                + DEFAULT_STEP_RECOVERY_TOLERANCES.source_cad_frame_orthogonality
        {
            bail!("revolution angle must not exceed one full turn");
        }
        if axis.origin_mm.iter().any(|value| !value.is_finite())
            || axis.direction.iter().any(|value| !value.is_finite())
        {
            bail!("revolution axis must be finite");
        }

        let direction = Vector3::new(axis.direction[0], axis.direction[1], axis.direction[2]);
        let magnitude = direction.magnitude();
        if !magnitude.is_finite() || magnitude <= DEFAULT_STEP_RECOVERY_TOLERANCES.cad_normal_unit {
            bail!("revolution axis direction must be nonzero");
        }
        let direction = direction / magnitude;
        let origin = Point3::new(axis.origin_mm[0], axis.origin_mm[1], axis.origin_mm[2]);
        let closed = (angle_rad.abs() - std::f64::consts::TAU).abs()
            <= DEFAULT_STEP_RECOVERY_TOLERANCES.source_cad_frame_orthogonality;
        let sweep = if closed {
            builder::SweepAngle::Closed
        } else {
            builder::SweepAngle::Partial(Rad(angle_rad))
        };

        let touches_axis = wires.iter().any(|wire| {
            wire.vertex_iter().any(|vertex| {
                point_axis_distance(vertex.point(), origin, direction)
                    <= DEFAULT_STEP_RECOVERY_TOLERANCES.cad_axis_locus_mm
            })
        });

        if touches_axis {
            if !closed {
                bail!("partial revolution with a profile touching the axis is not yet supported");
            }
            if wires.len() != 1 {
                bail!("axis-touching revolution with profile holes is not yet supported");
            }

            let mut wire = wires
                .pop()
                .ok_or_else(|| anyhow::anyhow!("missing profile wire"))?;
            let axis_edges = wire
                .iter()
                .enumerate()
                .filter_map(|(index, edge)| {
                    let curve = edge.oriented_curve();
                    let (start, end) = curve.range_tuple();
                    let midpoint = curve.evaluate(f64::midpoint(start, end));
                    (point_axis_distance(edge.front().point(), origin, direction)
                        <= DEFAULT_STEP_RECOVERY_TOLERANCES.cad_axis_locus_mm
                        && point_axis_distance(edge.back().point(), origin, direction)
                            <= DEFAULT_STEP_RECOVERY_TOLERANCES.cad_axis_locus_mm
                        && point_axis_distance(midpoint, origin, direction)
                            <= DEFAULT_STEP_RECOVERY_TOLERANCES.cad_axis_locus_mm)
                        .then_some(index)
                })
                .collect::<Vec<_>>();
            let [axis_edge_index] = axis_edges.as_slice() else {
                bail!("axis-touching closed profile needs exactly one line segment on the axis");
            };

            wire.rotate_left(*axis_edge_index);
            wire.pop_front()
                .ok_or_else(|| anyhow::anyhow!("missing axis profile edge"))?;
            let Some(front) = wire.front_vertex() else {
                bail!("axis-touching revolution has no generatrix");
            };
            let Some(back) = wire.back_vertex() else {
                bail!("axis-touching revolution has no generatrix");
            };
            if point_axis_distance(front.point(), origin, direction)
                > DEFAULT_STEP_RECOVERY_TOLERANCES.cad_axis_locus_mm
                || point_axis_distance(back.point(), origin, direction)
                    > DEFAULT_STEP_RECOVERY_TOLERANCES.cad_axis_locus_mm
            {
                bail!("revolution generatrix endpoints must lie on the axis");
            }

            // A generatrix running along the positive sweep axis creates
            // an inward-oriented boundary. The missing axis segment closes
            // the profile but does not contribute a swept face: choose the
            // outward material orientation before constructing any coedges.
            // This is axis-relative and independent of the world frame.
            if (back.point() - front.point()).dot(direction) > 0.0 {
                wire.invert();
            }
            let shell =
                builder::revolve_wire(&wire, origin, direction, builder::SweepAngle::Closed, 2)?;
            return regularize_generated_revolution_caps(Solid::new(vec![shell])?);
        }

        let face = profile::attach_plane_normalized(wires)?;
        regularize_generated_revolution_caps(builder::revolve(
            &face,
            origin,
            direction,
            sweep,
            if closed { 4 } else { 1 },
        )?)
    }

    const fn rigid_matrix(transform: &RigidTransform) -> Matrix4 {
        let m = transform.matrix;
        Matrix4::from_cols(
            Vector4::new(m[0][0], m[1][0], m[2][0], m[3][0]),
            Vector4::new(m[0][1], m[1][1], m[2][1], m[3][1]),
            Vector4::new(m[0][2], m[1][2], m[2][2], m[3][2]),
            Vector4::new(m[0][3], m[1][3], m[2][3], m[3][3]),
        )
    }

    fn transform_axis(axis: Axis3, transform: &RigidTransform) -> Axis3 {
        let m = transform.matrix;
        let [x, y, z] = axis.origin_mm;
        let [dx, dy, dz] = axis.direction;
        Axis3 {
            origin_mm: [
                m[0][2].mul_add(z, m[0][1].mul_add(y, m[0][0].mul_add(x, m[0][3]))),
                m[1][2].mul_add(z, m[1][1].mul_add(y, m[1][0].mul_add(x, m[1][3]))),
                m[2][2].mul_add(z, m[2][1].mul_add(y, m[2][0].mul_add(x, m[2][3]))),
            ],
            direction: [
                m[0][2].mul_add(dz, m[0][1].mul_add(dy, m[0][0] * dx)),
                m[1][2].mul_add(dz, m[1][1].mul_add(dy, m[1][0] * dx)),
                m[2][2].mul_add(dz, m[2][1].mul_add(dy, m[2][0] * dx)),
            ],
        }
    }

    fn point_axis_distance(point: Point3, origin: Point3, direction: Vector3) -> f64 {
        (point - origin).cross(direction).magnitude()
    }

    fn profile_loop_wire(curves: &[Curve2d]) -> Result<Wire> {
        if curves.is_empty() {
            bail!("profile loop must contain at least one curve");
        }
        if curves.len() == 1
            && let Curve2d::CircleArc {
                center_mm,
                radius_mm,
                start_angle_rad,
                end_angle_rad,
            } = &curves[0]
            && ((end_angle_rad - start_angle_rad).abs() - std::f64::consts::TAU).abs()
                <= DEFAULT_STEP_RECOVERY_TOLERANCES.cad_axis_locus_mm
        {
            let start = Point3::new(
                center_mm[0] + radius_mm * start_angle_rad.cos(),
                center_mm[1] + radius_mm * start_angle_rad.sin(),
                0.0,
            );
            return Ok(primitive::circle(
                start,
                Point3::new(center_mm[0], center_mm[1], 0.0),
                Vector3::new(0.0, 0.0, 1.0),
                2,
            ));
        }
        profile_wire(curves)
    }

    #[expect(
        clippy::too_many_lines,
        reason = "topology recovery keeps sequential fail-closed proof stages explicit for auditability"
    )]
    fn profile_wire(curves: &[Curve2d]) -> Result<Wire> {
        let starts = curves
            .iter()
            .map(curve_start_point)
            .collect::<Result<Vec<_>>>()?;
        let vertices = starts
            .iter()
            .map(|[x, y]| builder::vertex(Point3::new(*x, *y, 0.0)))
            .collect::<Vec<_>>();
        let mut edges = Vec::with_capacity(curves.len());

        for (index, curve) in curves.iter().enumerate() {
            let next = (index + 1) % curves.len();
            let expected_end = curve_end_point(curve)?;
            if point2_distance(expected_end, starts[next])
                > DEFAULT_STEP_RECOVERY_TOLERANCES.cad_profile_join_mm
            {
                bail!("profile curve endpoints are not topologically continuous");
            }
            match curve {
                Curve2d::Line { .. } => {
                    edges.push(builder::line(&vertices[index], &vertices[next]));
                }
                Curve2d::CircleArc {
                    center_mm,
                    radius_mm,
                    start_angle_rad,
                    end_angle_rad,
                } => {
                    let sweep = end_angle_rad - start_angle_rad;
                    if !sweep.is_finite()
                        || sweep.abs() <= DEFAULT_STEP_RECOVERY_TOLERANCES.cad_normal_unit
                    {
                        bail!("circle arc sweep must be finite and nonzero");
                    }
                    if sweep.abs()
                        >= std::f64::consts::TAU
                            - DEFAULT_STEP_RECOVERY_TOLERANCES.cad_axis_locus_mm
                    {
                        bail!("full circles must be represented by one profile curve");
                    }
                    let midpoint = f64::midpoint(*start_angle_rad, *end_angle_rad);
                    let transit = Point3::new(
                        center_mm[0] + radius_mm * midpoint.cos(),
                        center_mm[1] + radius_mm * midpoint.sin(),
                        0.0,
                    );
                    edges.push(builder::circle_arc(
                        &vertices[index],
                        &vertices[next],
                        transit,
                    )?);
                }
                Curve2d::Bezier { control_points_mm } => {
                    if control_points_mm.len() < 2 {
                        bail!("Bezier profile curve needs at least two control points");
                    }
                    let inter_points = control_points_mm[1..control_points_mm.len() - 1]
                        .iter()
                        .map(|[x, y]| Point3::new(*x, *y, 0.0))
                        .collect::<Vec<_>>();
                    edges.push(builder::bezier(
                        &vertices[index],
                        &vertices[next],
                        inter_points,
                    ));
                }
                Curve2d::BSpline {
                    degree,
                    control_points_mm,
                    knots,
                    weights,
                } => {
                    let Some(min_control_points) = degree.checked_add(1) else {
                        bail!("B-spline profile degree overflows dimensions");
                    };
                    let Some(expected_knots) =
                        control_points_mm.len().checked_add(min_control_points)
                    else {
                        bail!("B-spline profile knot count overflows dimensions");
                    };
                    let finite_points = control_points_mm
                        .iter()
                        .flatten()
                        .all(|value| value.is_finite());
                    let endpoint_clamped =
                        min_control_points.checked_mul(2).is_some_and(|min_knots| {
                            knots.len() >= min_knots
                                && knots[..min_control_points]
                                    .iter()
                                    .all(|knot| knot.to_bits() == knots[0].to_bits())
                                && knots[knots.len() - min_control_points..]
                                    .iter()
                                    .all(|knot| knot.to_bits() == knots[knots.len() - 1].to_bits())
                        });
                    if *degree == 0
                        || control_points_mm.len() < min_control_points
                        || knots.len() != expected_knots
                        || !finite_points
                        || !knots.iter().all(|knot| knot.is_finite())
                        || knots.windows(2).any(|pair| pair[0] > pair[1])
                        || knots.first().map(|value| value.to_bits())
                            == knots.last().map(|value| value.to_bits())
                        || !endpoint_clamped
                    {
                        bail!("invalid B-spline profile geometry");
                    }
                    if let Some(weights) = weights
                        && (weights.len() != control_points_mm.len()
                            || weights
                                .iter()
                                .any(|weight| !weight.is_finite() || *weight <= 0.0))
                    {
                        bail!("invalid B-spline profile weights");
                    }
                    let control_points = control_points_mm
                        .iter()
                        .map(|[x, y]| Point3::new(*x, *y, 0.0))
                        .collect::<Vec<_>>();
                    let curve = BsplineCurve::new(KnotVector::from(knots.clone()), control_points)?;
                    let geometry = if let Some(weights) = weights {
                        Curve::NurbsCurve(NurbsCurve::try_from_bspline_and_weights(
                            curve,
                            weights.clone(),
                        )?)
                    } else {
                        Curve::BsplineCurve(curve)
                    };
                    edges.push(Edge::new(&vertices[index], &vertices[next], geometry)?);
                }
            }
        }

        Ok(edges.into())
    }

    fn curve_start_point(curve: &Curve2d) -> Result<[f64; 2]> {
        match curve {
            Curve2d::Line { start_mm, .. } => Ok(*start_mm),
            Curve2d::CircleArc {
                center_mm,
                radius_mm,
                start_angle_rad,
                ..
            } => Ok([
                center_mm[0] + radius_mm * start_angle_rad.cos(),
                center_mm[1] + radius_mm * start_angle_rad.sin(),
            ]),
            Curve2d::Bezier { control_points_mm }
            | Curve2d::BSpline {
                control_points_mm, ..
            } => control_points_mm
                .first()
                .copied()
                .ok_or_else(|| anyhow::anyhow!("spline profile has no control points")),
        }
    }

    fn curve_end_point(curve: &Curve2d) -> Result<[f64; 2]> {
        match curve {
            Curve2d::Line { end_mm, .. } => Ok(*end_mm),
            Curve2d::CircleArc {
                center_mm,
                radius_mm,
                end_angle_rad,
                ..
            } => Ok([
                center_mm[0] + radius_mm * end_angle_rad.cos(),
                center_mm[1] + radius_mm * end_angle_rad.sin(),
            ]),
            Curve2d::Bezier { control_points_mm }
            | Curve2d::BSpline {
                control_points_mm, ..
            } => control_points_mm
                .last()
                .copied()
                .ok_or_else(|| anyhow::anyhow!("spline profile has no control points")),
        }
    }

    fn point2_distance(a: [f64; 2], b: [f64; 2]) -> f64 { (a[0] - b[0]).hypot(a[1] - b[1]) }

    const fn monstertruck_policy(policy: GeometryPolicy) -> monstertruck_solid::GeometryPolicy {
        match policy {
            GeometryPolicy::Equivalent => monstertruck_solid::GeometryPolicy::Equivalent,
            GeometryPolicy::Standard { tolerance_mm } => {
                monstertruck_solid::GeometryPolicy::standard(tolerance_mm)
            }
        }
    }

    /// Extract a prism's generating profile and exact placement, without
    /// evaluating a solid. A union of coextensive prisms distributively
    /// equals extrusion of the union of their planar profiles.
    fn placed_prism(model: &CadModel, mut id: NodeId) -> Option<(&Profile2d, f64, Matrix4)> {
        let mut transform = rigid_matrix(&RigidTransform::identity());
        let node = loop {
            match model.node(id).ok()? {
                CadNode::Transform {
                    child,
                    transform: local,
                } => {
                    transform = transform * rigid_matrix(local);
                    id = *child;
                }
                node => break node,
            }
        };
        let CadNode::Extrude {
            profile,
            vector_mm: [x, y, z],
        } = node
        else {
            return None;
        };
        (*x == 0.0 && *y == 0.0 && z.is_finite() && *z != 0.0).then_some((profile, *z, transform))
    }

    fn coextensive_prism_union(
        model: &CadModel,
        children: &[NodeId],
        policy: GeometryPolicy,
    ) -> Result<Option<Solid>> {
        if children.len() < 2 {
            return Ok(None);
        }
        let mut pending = children.to_vec();
        let mut seen = HashSet::new();
        let mut prisms = Vec::new();
        while let Some(id) = pending.pop() {
            if !seen.insert(id) {
                continue;
            }
            if let CadNode::Boolean {
                op: BooleanOp::Union,
                children,
            } = model.node(id)?
            {
                pending.extend(children.iter().copied());
            } else if let Some(prism) = placed_prism(model, id) {
                prisms.push(prism);
            } else {
                return Ok(None);
            }
        }
        let Some(&(_, length, transform)) = prisms.first() else {
            return Ok(None);
        };
        if prisms
            .iter()
            .any(|(_, l, t)| l.to_bits() != length.to_bits() || *t != transform)
        {
            return Ok(None);
        }
        let mut faces = Vec::with_capacity(prisms.len());
        for (profile, _, _) in prisms {
            let face: monstertruck_modeling::Face =
                profile::attach_plane_normalized(profile_wires(profile)?)?;
            let mut face = face;
            let Surface::Plane(plane) = face.oriented_surface() else {
                return Ok(None);
            };
            // A shared plane chart must have an unambiguous common sense.
            // The final extrusion uses the same local positive/negative Z
            // orientation contract as an individual prism.
            if plane.normal().z * length < 0.0 {
                face.invert();
            }
            faces.push(face);
        }
        let Some(united) =
            monstertruck_solid::try_planar_profile_union(&faces, monstertruck_policy(policy))?
        else {
            return Ok(None);
        };
        let mut boundaries = Vec::new();
        for face in &united {
            let solid: Solid = builder::extrude(face, Vector3::new(0.0, 0.0, length));
            boundaries.extend(solid.into_boundaries());
        }
        let joined = Solid::new(boundaries)?;
        Ok(Some(builder::transformed(&joined, transform)))
    }

    fn evaluate_boolean(
        model: &CadModel,
        op: BooleanOp,
        children: &[NodeId],
        policy: GeometryPolicy,
        cache: &mut EvaluationCache,
    ) -> Result<Solid> {
        if op == BooleanOp::Union
            && let Some(result) = coextensive_prism_union(model, children, policy)?
        {
            return Ok(result);
        }
        if matches!(op, BooleanOp::Union | BooleanOp::Difference) && children.len() > 2 {
            let solids = children
                .iter()
                .map(|&id| evaluate_node(model, id, policy, cache))
                .collect::<Result<Vec<_>>>()?;
            let tools = solids[1..].iter().collect::<Vec<_>>();
            let contact_op = match op {
                BooleanOp::Union => monstertruck_solid::PlanarContactOp::Union,
                BooleanOp::Difference => monstertruck_solid::PlanarContactOp::Difference,
                BooleanOp::Intersection => {
                    unreachable!("contact arrangement requires union or difference")
                }
            };
            let joined = monstertruck_solid::try_planar_contact_many(
                &solids[0],
                &tools,
                contact_op,
                monstertruck_policy(policy),
            )?;
            // Contact certification borrows, rather than consumes, the
            // operands. If it cannot prove contact, fold those same owned
            // solids, certifying each pair independently before general
            // intersection. One non-contact tool must not disable exact
            // contact for every other tool in the batch.
            return joined.map_or_else(|| fold_boolean(op, solids.into_iter().map(Ok), policy), Ok);
        }
        fold_boolean(
            op,
            children
                .iter()
                .map(|&id| evaluate_node(model, id, policy, cache)),
            policy,
        )
    }

    fn fold_boolean(
        op: BooleanOp,
        mut operands: impl Iterator<Item = Result<Solid>>,
        policy: GeometryPolicy,
    ) -> Result<Solid> {
        let mut result = operands
            .next()
            .context("boolean node needs at least one child")??;
        let kernel_policy = monstertruck_policy(policy);
        for rhs in operands {
            let rhs = rhs?;
            result = match op {
                BooleanOp::Union => {
                    if let Some(joined) =
                        monstertruck_solid::try_planar_contact_union(&result, &rhs, kernel_policy)?
                    {
                        joined
                    } else {
                        monstertruck_solid::or(&result, &rhs, kernel_policy)?
                    }
                }
                BooleanOp::Intersection => monstertruck_solid::and(&result, &rhs, kernel_policy)?,
                BooleanOp::Difference => {
                    if let Some(joined) = monstertruck_solid::try_planar_contact_many(
                        &result,
                        &[&rhs],
                        monstertruck_solid::PlanarContactOp::Difference,
                        kernel_policy,
                    )? {
                        joined
                    } else {
                        monstertruck_solid::difference(&result, &rhs, kernel_policy)?
                    }
                }
            };
        }
        Ok(result)
    }

    fn primitive_solid(primitive: &Primitive) -> Result<Solid> {
        match primitive {
            Primitive::Box { size_mm } => {
                let [x, y, z] = *size_mm;
                if [x, y, z]
                    .into_iter()
                    .any(|value| !value.is_finite() || value <= 0.0)
                {
                    bail!("box dimensions must be finite and positive");
                }
                let vertex = builder::vertex(Point3::origin());
                let edge = builder::extrude(&vertex, Vector3::new(x, 0.0, 0.0));
                let face = builder::extrude(&edge, Vector3::new(0.0, y, 0.0));
                Ok(builder::extrude(&face, Vector3::new(0.0, 0.0, z)))
            }
            Primitive::Cylinder {
                radius_mm,
                height_mm,
            } => {
                if !radius_mm.is_finite()
                    || *radius_mm <= 0.0
                    || !height_mm.is_finite()
                    || *height_mm <= 0.0
                {
                    bail!("cylinder dimensions must be finite and positive");
                }
                let vertex = builder::vertex(Point3::new(0.0, 0.0, *radius_mm));
                let circle = builder::revolve(
                    &vertex,
                    Point3::origin(),
                    Vector3::unit_y(),
                    builder::SweepAngle::Closed,
                    2,
                )?;
                let disk = profile::attach_plane_normalized(vec![circle])?;
                Ok(builder::extrude(&disk, Vector3::new(0.0, *height_mm, 0.0)))
            }
            Primitive::Sphere { .. } => {
                bail!("Monstertruck backend does not yet evaluate full-sphere primitives")
            }
        }
    }

    const EDGE_SELECTOR_TOLERANCE_MM: f64 = DEFAULT_STEP_RECOVERY_TOLERANCES.cad_edge_selector_mm;

    const fn segment_point(point_mm: [f64; 3]) -> Point3 {
        Point3::new(point_mm[0], point_mm[1], point_mm[2])
    }

    fn points_match(left: Point3, right: Point3) -> bool {
        (left - right).magnitude() <= EDGE_SELECTOR_TOLERANCE_MM
    }

    const CIRCLE_SELECTOR_SAMPLES: usize = 16;

    fn circle_point_matches(
        point: Point3,
        center: Point3,
        normal: Vector3,
        radius_mm: f64,
    ) -> bool {
        let offset = point - center;
        offset.dot(normal).abs() <= EDGE_SELECTOR_TOLERANCE_MM
            && (offset.magnitude() - radius_mm).abs() <= EDGE_SELECTOR_TOLERANCE_MM
    }

    fn edge_matches_circle(edge: &Edge, center: Point3, normal: Vector3, radius_mm: f64) -> bool {
        let curve = edge.oriented_curve();
        let (start, end) = curve.range_tuple();
        if !start.is_finite()
            || !end.is_finite()
            || (start - end).abs() <= f64::EPSILON * start.abs().max(end.abs()).max(1.0)
        {
            return false;
        }
        (0..=CIRCLE_SELECTOR_SAMPLES).all(|index| {
            let Ok(index) = u32::try_from(index) else {
                return false;
            };
            let fraction = f64::from(index) / 16.0;
            let parameter = f64::mul_add(end - start, fraction, start);
            circle_point_matches(curve.evaluate(parameter), center, normal, radius_mm)
        })
    }

    fn endpoint_cluster(points: &mut Vec<Point3>, point: Point3) -> usize {
        points
            .iter()
            .position(|candidate| points_match(*candidate, point))
            .unwrap_or_else(|| {
                points.push(point);
                points.len() - 1
            })
    }

    fn circle_edges_form_one_cycle(edges: &[Edge]) -> bool {
        if edges.is_empty() {
            return false;
        }
        let mut points = Vec::<Point3>::new();
        let endpoints = edges
            .iter()
            .map(|edge| {
                let front = endpoint_cluster(&mut points, edge.absolute_front().point());
                let back = endpoint_cluster(&mut points, edge.absolute_back().point());
                (front, back)
            })
            .collect::<Vec<_>>();
        let mut degrees = vec![0_usize; points.len()];
        let mut adjacency = vec![Vec::<usize>::new(); points.len()];
        for &(front, back) in &endpoints {
            if front == back {
                degrees[front] += 2;
                adjacency[front].push(front);
            } else {
                degrees[front] += 1;
                degrees[back] += 1;
                adjacency[front].push(back);
                adjacency[back].push(front);
            }
        }
        if degrees.iter().any(|degree| *degree != 2) {
            return false;
        }
        let mut visited = vec![false; points.len()];
        let mut stack = vec![0_usize];
        while let Some(index) = stack.pop() {
            if visited[index] {
                continue;
            }
            visited[index] = true;
            stack.extend(adjacency[index].iter().copied());
        }
        visited.into_iter().all(std::convert::identity)
    }

    fn select_edges(shell: &Shell, selector: &EdgeSelector) -> Result<Vec<Edge>> {
        match selector {
            EdgeSelector::Segment { start_mm, end_mm } => {
                let start = segment_point(*start_mm);
                let end = segment_point(*end_mm);
                let mut matched = None;
                let mut matched_id = None;
                for edge in shell.edge_iter() {
                    let front = edge.absolute_front().point();
                    let back = edge.absolute_back().point();
                    let endpoints_match = (points_match(front, start) && points_match(back, end))
                        || (points_match(front, end) && points_match(back, start));
                    if !endpoints_match || matched_id == Some(edge.id()) {
                        continue;
                    }
                    if matched_id.is_some() {
                        bail!("edge-treatment selector matched more than one topological edge");
                    }
                    matched_id = Some(edge.id());
                    matched = Some(edge);
                }
                Ok(vec![matched.ok_or_else(|| {
                    anyhow::anyhow!("edge-treatment selector matched no edge")
                })?])
            }
            EdgeSelector::CircleArc {
                center_mm,
                normal,
                radius_mm,
                start_mm,
                sweep_rad,
            } => circular_selection::select_arc(
                shell, *center_mm, *normal, *radius_mm, *start_mm, *sweep_rad,
            ),
            EdgeSelector::Circle {
                center_mm,
                normal,
                radius_mm,
            } => {
                let center = segment_point(*center_mm);
                let normal = Vector3::new(normal[0], normal[1], normal[2]);
                let normal_magnitude = normal.magnitude();
                if !normal_magnitude.is_finite()
                    || normal_magnitude <= DEFAULT_STEP_RECOVERY_TOLERANCES.cad_normal_unit
                {
                    bail!("edge-treatment circle selector has invalid normal");
                }
                let normal = normal / normal_magnitude;
                let mut matched = Vec::<Edge>::new();
                for edge in shell.edge_iter() {
                    if matched.iter().any(|candidate| candidate.id() == edge.id())
                        || !edge_matches_circle(&edge, center, normal, *radius_mm)
                    {
                        continue;
                    }
                    matched.push(edge);
                }
                if matched.is_empty() {
                    bail!("edge-treatment circle selector matched no edge");
                }
                if !circle_edges_form_one_cycle(&matched) {
                    bail!(
                        "edge-treatment circle selector did not resolve to one closed edge cycle"
                    );
                }
                Ok(matched)
            }
        }
    }

    fn resolve_treatment_edges(
        shell: &Shell,
        edges: &[EdgeTreatmentEdge],
    ) -> Result<(Vec<Edge>, Vec<f64>)> {
        let mut selected = Vec::<Edge>::new();
        let mut radii = Vec::new();
        for treatment in edges {
            for edge in select_edges(shell, &treatment.selector).with_context(|| {
                format!("resolve edge treatment selector {:?}", treatment.selector)
            })? {
                if selected.iter().any(|existing| existing.id() == edge.id()) {
                    bail!("edge-treatment selectors overlap on one topological edge");
                }
                selected.push(edge);
                radii.push(treatment.radius_mm);
            }
        }
        Ok((selected, radii))
    }

    fn evaluate_edge_treatment(
        model: &CadModel,
        child: NodeId,
        profile: EdgeTreatmentProfile,
        edges: &[EdgeTreatmentEdge],
        policy: GeometryPolicy,
        cache: &mut EvaluationCache,
    ) -> Result<Solid> {
        let child = evaluate_node(model, child, policy, cache)?;
        let mut boundaries = child.into_boundaries();
        let [shell] = boundaries.as_mut_slice() else {
            bail!("edge treatment currently requires exactly one boundary shell");
        };
        let (selected, radii) = match resolve_treatment_edges(shell, edges) {
            Ok(resolved) => resolved,
            Err(original) => {
                // An exact constructive seam may be partitioned by sectional
                // Boolean interfaces. Normalize the same geometric shell,
                // then retry the unchanged selector (not a wider tolerance).
                let normalized = monstertruck_modeling::normalize_planar_partitions(shell)
                    .with_context(|| {
                        format!("normalize planar partitions after selector failure: {original:#}")
                    })?;
                let resolved = resolve_treatment_edges(&normalized, edges).with_context(|| {
                    format!("selector remained unresolved after planar normalization: {original:#}")
                })?;
                *shell = normalized;
                resolved
            }
        };
        let fillet_profile = match profile {
            EdgeTreatmentProfile::Chamfer => {
                monstertruck_modeling::chamfer_planar_edges(shell, &selected, &radii)?;
                return Ok(Solid::new(boundaries)?);
            }
            EdgeTreatmentProfile::CarrierChamfer {
                first_carrier_normal: [x, y, z],
                second_setback_mm,
            } => {
                monstertruck_modeling::chamfer_carrier_edges(
                    shell,
                    &selected,
                    &radii,
                    Vector3::new(x, y, z),
                    second_setback_mm,
                )?;
                return Ok(Solid::new(boundaries)?);
            }
            EdgeTreatmentProfile::Round => FilletProfile::Round,
            EdgeTreatmentProfile::AsymmetricChamfer {
                first_setback_mm,
                second_setback_mm,
            } => FilletProfile::AsymmetricChamfer {
                first_setback_mm,
                second_setback_mm,
            },
        };
        let options = FilletOptions {
            radius: RadiusSpec::PerEdge(radii),
            profile: fillet_profile,
            ..FilletOptions::default()
        };
        monstertruck_modeling::fillet_edges(shell, &selected, Some(&options))
            .context("native edge treatment construction")?;
        Solid::new(boundaries).context("validate treated solid topology")
    }

    fn evaluate_node(
        model: &CadModel,
        root: NodeId,
        policy: GeometryPolicy,
        cache: &mut EvaluationCache,
    ) -> Result<Solid> {
        let key = (root, PolicyKey::from(policy));
        let (cached, retain) = cache.begin(key)?;
        if let Some(solid) = cached {
            return Ok(solid);
        }
        let evaluated = match model.node(root)? {
            CadNode::Primitive(primitive) => primitive_solid(primitive),
            CadNode::Extrude { profile, vector_mm } => {
                if vector_mm[0] != 0.0 || vector_mm[1] != 0.0 {
                    bail!("Monstertruck backend currently supports local Z extrusion only");
                }
                extrude_profile_z(profile, vector_mm[2])
            }
            CadNode::Revolve {
                profile,
                axis,
                angle_rad,
            } => revolve_profile(profile, *axis, *angle_rad),
            CadNode::Sweep { profile, path } => sweep_profile(profile, path),
            CadNode::Loft { sections } => loft_sections(sections),
            CadNode::Boolean { op, children } => {
                evaluate_boolean(model, *op, children, policy, cache)
            }
            CadNode::EdgeTreatment {
                child,
                profile,
                edges,
            } => evaluate_edge_treatment(model, *child, *profile, edges, policy, cache),
            CadNode::SolverTolerance {
                child,
                tolerance_mm,
            } => evaluate_node(
                model,
                *child,
                GeometryPolicy::standard_mm(*tolerance_mm),
                cache,
            ),
            CadNode::Transform { transform, child } => {
                if let CadNode::Revolve {
                    profile,
                    axis,
                    angle_rad,
                } = model.node(*child)?
                {
                    revolve_profile_transformed(profile, *axis, *angle_rad, transform)
                } else {
                    let child = evaluate_node(model, *child, policy, cache)?;
                    Ok(builder::transformed(&child, rigid_matrix(transform)))
                }
            }
            unsupported => {
                bail!("Monstertruck backend does not yet evaluate node {root:?}: {unsupported:?}")
            }
        }?;
        cache.remember(key, &evaluated, retain);
        Ok(evaluated)
    }

    impl CadKernel for MonstertruckKernel {
        type Evaluated = EvaluatedShape;

        fn evaluate(&self, model: &CadModel, root: NodeId) -> Result<Self::Evaluated> {
            model.validate()?;
            let mut cache = EvaluationCache::single_root(model, root, self.policy)?;
            Ok(EvaluatedShape {
                solid: evaluate_node(model, root, self.policy, &mut cache)?,
            })
        }

        fn summarize(&self, evaluated: &Self::Evaluated) -> KernelSummary {
            let boundaries = evaluated.solid.boundaries();
            KernelSummary {
                geometrically_consistent: evaluated.solid.is_geometric_consistent(),
                shells: boundaries.len(),
                faces: boundaries.iter().map(|shell| shell.len()).sum(),
            }
        }

        fn supports_step_export(&self) -> bool { true }

        fn to_step(&self, evaluated: &Self::Evaluated) -> Result<String> {
            let export_solid = step_export_solid(&evaluated.solid, self.policy)?;
            let compressed = export_solid.compress();
            Ok(CompleteStepDisplay::new(
                save::StepModel::from(&compressed),
                save::StepHeaderDescriptor {
                    organization_system: "step-redox Monstertruck backend".into(),
                    ..Default::default()
                },
            )
            .to_string())
        }
    }

    #[cfg(test)]
    mod tests {
        #![expect(
            clippy::manual_let_else,
            reason = "legacy test fixture code is allowed narrow lint exceptions; production paths remain strict"
        )]
        use super::*;
        use monstertruck_modeling::{Curve2D, Line, ParameterCurve, Plane, Point2};

        fn bowed_intersection(boundary_y: Option<f64>) -> Curve {
            let first = Surface::Plane(Plane::xy());
            let second = Surface::Plane(Plane::new(
                Point3::new(0.0, 0.0, 0.0),
                Point3::new(1.0, 0.0, 0.0),
                Point3::new(0.0, 0.0, 1.0),
            ));
            let guide = Curve::BsplineCurve(BsplineCurve::new_unchecked(
                KnotVector::uniform_knot(2, 1),
                vec![
                    Point3::new(0.0, 0.0, 0.0),
                    Point3::new(0.5, 0.0, 0.01),
                    Point3::new(1.0, 0.0, 0.0),
                ],
            ));
            Curve::Intersection(SurfaceCurve::with_boundaries(
                Box::new(first.clone()),
                Box::new(second),
                Box::new(guide),
                boundary_y.map(|offset| {
                    ParameterCurve::new(
                        Curve2D::Line(Line(Point2::new(0.0, offset), Point2::new(1.0, offset))),
                        Box::new(first),
                    )
                }),
                None,
            ))
        }

        #[test]
        fn step_export_uses_evaluated_intersection_not_marching_guide() -> Result<()> {
            let original = bowed_intersection(Some(0.0));
            let old_at_middle = original.evaluate(0.5);
            assert!(old_at_middle.z.abs() < 1e-7);
            let Curve::Intersection(raw) = &original else {
                bail!("test fixture is not an intersection");
            };
            assert!(raw.leader().evaluate(0.5).z > 0.004);

            let corrected = corrected_step_curve(&original, 1e-5)?;
            let Curve::Intersection(intersection) = corrected else {
                bail!("STEP correction lost intersection support");
            };
            assert!(intersection.boundary0().is_some());
            assert!(intersection.boundary1().is_none());
            let (lo, hi) = intersection.range_tuple();
            for i in 0..=32_u32 {
                let fraction = f64::from(i) / 32.0;
                let point = intersection
                    .leader()
                    .evaluate(f64::mul_add(hi - lo, fraction, lo));
                assert!(
                    point.y.abs() < 1e-5 && point.z.abs() < 1e-5,
                    "STEP leader must follow the true shared line: {point:?}"
                );
            }
            // Original evaluated CAD geometry is immutable under STEP export.
            let Curve::Intersection(source) = original else {
                bail!("original intersection was replaced");
            };
            assert!(source.leader().evaluate(0.5).z > 0.004);
            Ok(())
        }

        #[test]
        fn step_export_preserves_nonuniform_native_intersection_parameters() -> Result<()> {
            let Curve::Intersection(mut intersection) = bowed_intersection(None) else {
                bail!("fixture must be an intersection");
            };
            let mut knots = KnotVector::bezier_knot(2);
            knots.transform(2.5, 7.0);
            **intersection.leader_mut() = Curve::BsplineCurve(BsplineCurve::new_unchecked(
                knots,
                vec![
                    Point3::new(0.0, 0.0, 0.0),
                    Point3::new(0.5, 0.0, 0.01),
                    Point3::new(1.0, 0.0, 0.0),
                ],
            ));
            let original = Curve::Intersection(intersection);
            let range = (7.0, 9.5);
            assert_eq!(original.range_tuple(), range);
            let (parameters, points) = algo::curve::parameter_division(&original, range, 1e-5);
            let Curve::Intersection(exported) = corrected_step_curve(&original, 1e-5)? else {
                bail!("STEP correction lost intersection support");
            };
            assert_eq!(exported.range_tuple(), range);
            for (parameter, point) in parameters.into_iter().zip(points) {
                let rebuilt = exported.leader().evaluate(parameter);
                assert!(
                    (rebuilt - point).magnitude() < 1e-9,
                    "STEP remapped native sample {parameter}: {rebuilt:?} versus {point:?}"
                );
            }
            Ok(())
        }

        #[test]
        fn step_export_retains_certified_straight_rational_guide() -> Result<()> {
            let Curve::Intersection(mut intersection) = bowed_intersection(None) else {
                bail!("fixture must be an intersection");
            };
            let exact_line = NurbsCurve::new(BsplineCurve::new_unchecked(
                KnotVector::bezier_knot(1),
                vec![
                    Vector4::new(0.0, 0.0, 0.0, 1.0),
                    Vector4::new(1.0, 0.0, 0.0, 1.0),
                ],
            ));
            **intersection.leader_mut() = Curve::NurbsCurve(exact_line);
            let curve = Curve::Intersection(intersection);
            let corrected = corrected_step_curve(&curve, 1e-5)?;
            let Curve::Intersection(exported) = corrected else {
                bail!("STEP conversion lost the intersection");
            };
            assert!(matches!(exported.leader().as_ref(), Curve::NurbsCurve(_)));
            let (lo, hi) = exported.range_tuple();
            for fraction in [0.0, 0.25, 0.5, 0.75, 1.0] {
                let p = exported
                    .leader()
                    .evaluate(f64::mul_add(hi - lo, fraction, lo));
                assert!(p.y.abs() < 1e-10 && p.z.abs() < 1e-10);
            }
            Ok(())
        }

        #[test]
        fn step_export_preserves_globally_proven_quadratic_circle_not_its_chord() -> Result<()> {
            use monstertruck_geometry::analytic_curve_incidence::CurveIncidence;
            use monstertruck_modeling::{Line, Plane, Processor, RevolutionSurface};
            let cylinder = Surface::Revolution(Processor::new(RevolutionSurface::by_revolution(
                Curve::Line(Line(Point3::new(2.0, 0.0, 0.0), Point3::new(2.0, 0.0, 1.0))),
                Point3::origin(),
                Vector3::unit_z(),
            )));
            let plane = Surface::Plane(Plane::xy());
            let w = f64::sqrt(0.5);
            let original = Curve::NurbsCurve(NurbsCurve::new(BsplineCurve::new_unchecked(
                KnotVector::bezier_knot(2),
                vec![
                    Vector4::new(2.0, 0.0, 0.0, 1.0),
                    Vector4::new(2.0 * w, 2.0 * w, 0.0, w),
                    Vector4::new(0.0, 2.0, 0.0, 1.0),
                ],
            )));
            for surface in [&cylinder, &plane] {
                assert_eq!(
                    certify_analytic_curve_on_surface(&original, surface, 1.0e-7),
                    CurveIncidence::Proven,
                );
            }
            let intersection = Curve::Intersection(SurfaceCurve::with_boundaries(
                Box::new(cylinder),
                Box::new(plane),
                Box::new(original),
                None,
                None,
            ));
            let Curve::Intersection(exported) = corrected_step_curve(&intersection, 1.0e-5)? else {
                bail!("STEP conversion lost the intersection");
            };
            let Curve::NurbsCurve(exact) = exported.leader().as_ref() else {
                bail!("globally certified quadratic rational arc degraded into a polyline");
            };
            assert_eq!(exact.degree(), 2);
            for i in 0..=64 {
                let t = f64::from(i) / 64.0;
                let p = exact.evaluate(t);
                assert!((p.x.hypot(p.y) - 2.0).abs() < 1.0e-9);
                assert!(p.z.abs() < 1.0e-9);
            }
            Ok(())
        }

        #[test]
        fn step_export_rejects_off_carrier_straight_guide() -> Result<()> {
            let Curve::Intersection(mut intersection) = bowed_intersection(None) else {
                bail!("fixture must be an intersection");
            };
            **intersection.leader_mut() = Curve::Line(monstertruck_modeling::Line(
                Point3::new(0.0, 0.0, 0.01),
                Point3::new(1.0, 0.0, 0.01),
            ));
            let error = corrected_step_curve(&Curve::Intersection(intersection), 1e-5)
                .err()
                .ok_or_else(|| {
                    anyhow::anyhow!("unsupported straight guide must not be exported")
                })?;
            assert!(error.to_string().contains("both carriers"));
            Ok(())
        }

        #[test]
        fn step_export_rejects_inconsistent_associated_parameter_curve() {
            let wrong = bowed_intersection(Some(0.1));
            assert!(
                corrected_step_curve(&wrong, 1e-5)
                    .is_err_and(|error| error.to_string().contains("PCURVE 0"))
            );
        }

        #[test]
        fn step_export_rejects_unresolved_nonintersecting_carriers() -> Result<()> {
            let original = bowed_intersection(None);
            let Curve::Intersection(mut intersection) = original else {
                return Err(anyhow::anyhow!("fixture must be an intersection"));
            };
            **intersection.surface1_mut() = Surface::Plane(Plane::new(
                Point3::new(0.0, 0.0, 0.1),
                Point3::new(1.0, 0.0, 0.1),
                Point3::new(0.0, 1.0, 0.1),
            ));
            let wrong = Curve::Intersection(intersection);
            assert!(
                corrected_step_curve(&wrong, 1e-5)
                    .is_err_and(|error| error.to_string().contains("both carriers"))
            );
            Ok(())
        }

        #[test]
        fn step_export_accuracy_never_silently_widens_requested_tolerance() {
            let export_accuracy = DEFAULT_STEP_RECOVERY_TOLERANCES.cad_step_export_accuracy_mm;
            assert_eq!(export_accuracy, 1.0e-5);
            assert_eq!(
                step_intersection_accuracy(GeometryPolicy::project_standard()).ok(),
                Some(export_accuracy),
                "the default modeling allowance must not coarsen STEP edge fidelity"
            );
            assert_eq!(
                step_intersection_accuracy(GeometryPolicy::standard_mm(0.1)).ok(),
                Some(export_accuracy)
            );
            assert_eq!(
                step_intersection_accuracy(GeometryPolicy::standard_mm(5e-6)).ok(),
                Some(5e-6),
                "a stricter caller-requested accuracy must remain authoritative"
            );
            assert_eq!(
                step_intersection_accuracy(GeometryPolicy::standard_mm(1e-6)).ok(),
                Some(DEFAULT_NUMERICAL_TOLERANCES.kernel_f64)
            );
            for bad in [0.0, -1.0, 1e-8, f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
                assert!(
                    step_intersection_accuracy(GeometryPolicy::standard_mm(bad)).is_err(),
                    "invalid modeling tolerance must not become valid after an export cap: {bad}"
                );
            }
            assert_eq!(
                step_intersection_accuracy(GeometryPolicy::Equivalent).ok(),
                Some(DEFAULT_NUMERICAL_TOLERANCES.kernel_f64)
            );
        }

        use crate::cad_ir::{
            Axis3, CadNode, Curve2d, LoftSection, Profile2d, ProfileLoop, RigidTransform,
            SweepPath3d, SweepSegment3d,
        };

        fn box_model() -> Result<(CadModel, NodeId)> {
            let mut model = CadModel::new();
            let profile =
                Profile2d::polygon(vec![[0.0, 0.0], [10.0, 0.0], [10.0, 6.0], [0.0, 6.0]])?;
            let root = model.add_node(CadNode::Extrude {
                profile,
                vector_mm: [0.0, 0.0, 2.0],
            });
            model.add_root(root)?;
            Ok((model, root))
        }

        /// Kernel feasibility probe for the general bored-extrusion and
        /// paired ruled-lobe construction. Dimensions in this test are taken
        /// from one real STEP source; production recognition must derive them
        /// from topology rather than reuse this fixture.
        #[test]
        fn real_bored_extrusion_with_two_ruled_lobes() -> Result<()> {
            use crate::cad_ir::BooleanOp;
            use std::f64::consts::TAU;

            let y_back = -3.227_218_606_110_5;
            let y_stage = -0.677_218_606_110_5;
            let y_front = -0.277_218_606_110_5;
            let z_upper = 1.407_981_393_889_5;
            let z_lower = -1.132_018_606_110_5;
            let x_center = 71.12;
            let xy_frame = |y| RigidTransform {
                matrix: [
                    [1.0, 0.0, 0.0, 0.0],
                    [0.0, 0.0, 1.0, y],
                    [0.0, -1.0, 0.0, 0.0],
                    [0.0, 0.0, 0.0, 1.0],
                ],
            };
            let back_profile = [
                [70.131_985_117_132, z_lower],
                [72.108_014_882_868, z_lower],
                [72.29, -0.632_018_606_110_5],
                [72.39, -0.632_018_606_110_5],
                [72.39, 0.907_981_393_889_5],
                [72.29, 0.907_981_393_889_5],
                [72.108_014_882_868, z_upper],
                [70.131_985_117_132, z_upper],
                [69.95, 0.907_981_393_889_5],
                [69.85, 0.907_981_393_889_5],
                [69.85, -0.632_018_606_110_5],
                [69.95, -0.632_018_606_110_5],
            ];
            let rear_polygon =
                Profile2d::polygon(back_profile.iter().map(|[x, z]| [*x, -*z]).collect())?;
            let mut rear_with_bore = rear_polygon;
            rear_with_bore.loops.push(ProfileLoop {
                curves: vec![Curve2d::CircleArc {
                    center_mm: [x_center, -0.137_981_393_889_5],
                    radius_mm: 0.69,
                    start_angle_rad: 0.0,
                    end_angle_rad: TAU,
                }],
            });

            let mut model = CadModel::new();
            let rear_extrude = model.add_node(CadNode::Extrude {
                profile: rear_with_bore,
                vector_mm: [0.0, 0.0, y_stage - y_back],
            });
            let rear = model.add_node(CadNode::Transform {
                transform: xy_frame(y_back),
                child: rear_extrude,
            });
            let lobes = [
                (
                    [70.006_201_337_707_14, 72.233_798_662_292_86],
                    1.062_393_300_183,
                    [70.059_191_070_279_85, 72.180_808_929_720_14],
                    1.207_981_393_89,
                    z_upper,
                ),
                (
                    [70.006_201_337_707_14, 72.233_798_662_292_86],
                    -0.786_430_512_404,
                    [70.059_191_070_279_85, 72.180_808_929_720_14],
                    -0.932_018_606_11,
                    z_lower,
                ),
            ];
            let mut components = vec![rear];
            for (back, back_z, front, front_z, clip) in lobes {
                let profile = |ends: [f64; 2], base_z: f64| -> Result<Profile2d> {
                    let points = vec![
                        [ends[0], -base_z],
                        [ends[1], -base_z],
                        [72.108_014_882_868, -clip],
                        [70.131_985_117_132, -clip],
                    ];
                    Profile2d::polygon(points)
                };
                components.push(model.add_node(CadNode::Loft {
                    sections: vec![
                        LoftSection {
                            profile: profile(back, back_z)?,
                            transform: xy_frame(y_stage),
                        },
                        LoftSection {
                            profile: profile(front, front_z)?,
                            transform: xy_frame(y_front),
                        },
                    ],
                }));
            }
            let root = model.add_node(CadNode::Boolean {
                op: BooleanOp::Union,
                children: components.clone(),
            });
            model.add_root(root)?;
            model.validate()?;

            let kernel = MonstertruckKernel::project_standard();
            let evaluated = kernel.evaluate(&model, root)?;
            assert!(kernel.summarize(&evaluated).is_single_solid());
            let step = kernel.to_step(&evaluated)?;
            let exchange =
                monstertruck_io::step::load::ParsedStepSource::parse(step.as_bytes())?.exchange;
            let source_solid_ids = exchange
                .data
                .iter()
                .flat_map(|section| &section.entities)
                .filter_map(|entity| {
                    let record = crate::step_graph::simple_record(entity)?;
                    (record.name == "MANIFOLD_SOLID_BREP")
                        .then(|| crate::step_graph::entity_id(entity))
                })
                .collect::<Vec<_>>();
            let [rebuilt_id] = source_solid_ids.as_slice() else {
                bail!("native three-part reconstruction must export exactly one STEP solid");
            };
            let comparison = crate::compare_analytic_solid_boundaries_bytes(
                include_bytes!("../validation/fixtures/native_mixed_carrier_chamfer_06d018ae.step"),
                19,
                step.as_bytes(),
                *rebuilt_id,
            )?;
            assert!(
                comparison.equivalent,
                "source-derived bored extrusion and two ruled lofts must have identical oriented analytic boundaries",
            );
            Ok(())
        }

        #[test]
        fn evaluates_tapered_rectangular_loft() -> Result<()> {
            let mut model = CadModel::new();
            let base =
                Profile2d::polygon(vec![[-1.0, -1.0], [1.0, -1.0], [1.0, 1.0], [-1.0, 1.0]])?;
            let end = Profile2d::polygon(vec![[-0.5, -0.5], [0.5, -0.5], [0.5, 0.5], [-0.5, 0.5]])?;
            let root = model.add_node(CadNode::Loft {
                sections: vec![
                    LoftSection {
                        profile: base,
                        transform: RigidTransform::identity(),
                    },
                    LoftSection {
                        profile: end,
                        transform: RigidTransform::translation_mm([0.0, 0.0, 2.0]),
                    },
                ],
            });
            model.add_root(root)?;
            let kernel = MonstertruckKernel::project_standard();
            let evaluated = kernel.evaluate(&model, root)?;
            let summary = kernel.summarize(&evaluated);
            assert!(summary.geometrically_consistent);
            assert_eq!(summary.shells, 1);
            assert_eq!(summary.faces, 6);
            let step = kernel.to_step(&evaluated)?;
            step_p21::parser::parse(&step)?;
            assert!(
                !step.contains("B_SPLINE_SURFACE"),
                "coplanar straight-edge loft walls should remain analytic planes"
            );
            Ok(())
        }

        fn curved_profile_model(curve: Curve2d) -> Result<(CadModel, NodeId)> {
            let mut model = CadModel::new();
            let profile = Profile2d {
                loops: vec![ProfileLoop {
                    curves: vec![
                        curve,
                        Curve2d::Line {
                            start_mm: [2.0, 0.0],
                            end_mm: [0.0, 0.0],
                        },
                    ],
                }],
            };
            let root = model.add_node(CadNode::Extrude {
                profile,
                vector_mm: [0.0, 0.0, 2.0],
            });
            model.add_root(root)?;
            Ok((model, root))
        }

        fn assert_curved_profile_evaluates(
            curve: Curve2d,
            expected_step_fragment: &str,
        ) -> Result<()> {
            let (model, root) = curved_profile_model(curve)?;
            let kernel = MonstertruckKernel::project_standard();
            let evaluated = kernel.evaluate(&model, root)?;
            let summary = kernel.summarize(&evaluated);
            assert!(summary.geometrically_consistent);
            assert_eq!(summary.shells, 1);
            assert_eq!(summary.faces, 4);

            let step = kernel.to_step(&evaluated)?;
            assert!(step.contains(expected_step_fragment));
            step_p21::parser::parse(&step)?;
            Ok(())
        }

        #[test]
        fn evaluates_exact_line_arc_line_sweep() -> Result<()> {
            let mut model = CadModel::new();
            let profile = Profile2d::polygon(vec![
                [-0.25, -0.25],
                [0.25, -0.25],
                [0.25, 0.25],
                [-0.25, 0.25],
            ])?;
            let root = model.add_node(CadNode::Sweep {
                profile,
                path: SweepPath3d {
                    initial_transform: RigidTransform::translation_mm([2.0, 0.0, 4.0]),
                    closed: false,
                    segments: vec![
                        SweepSegment3d::Translation {
                            vector_mm: [0.0, 0.0, -4.0],
                        },
                        SweepSegment3d::Rotation {
                            axis: Axis3 {
                                origin_mm: [0.0, 0.0, 0.0],
                                direction: [0.0, 1.0, 0.0],
                            },
                            angle_rad: std::f64::consts::FRAC_PI_2,
                        },
                        SweepSegment3d::Translation {
                            vector_mm: [-4.0, 0.0, 0.0],
                        },
                    ],
                },
            });
            model.add_root(root)?;

            let kernel = MonstertruckKernel::project_standard();
            let evaluated = kernel.evaluate(&model, root)?;
            let summary = kernel.summarize(&evaluated);
            assert!(summary.geometrically_consistent);
            assert_eq!(summary.shells, 1);
            assert_eq!(summary.faces, 14);

            let step = kernel.to_step(&evaluated)?;
            step_p21::parser::parse(&step)?;
            Ok(())
        }

        fn box_edge_treatment(profile: EdgeTreatmentProfile) -> Result<()> {
            let mut model = CadModel::new();
            let box_node = model.add_node(CadNode::Primitive(Primitive::Box {
                size_mm: [4.0, 3.0, 2.0],
            }));
            let root = model.add_node(CadNode::EdgeTreatment {
                child: box_node,
                profile,
                edges: vec![EdgeTreatmentEdge {
                    selector: EdgeSelector::Segment {
                        start_mm: [0.0, 0.0, 0.0],
                        end_mm: [4.0, 0.0, 0.0],
                    },
                    radius_mm: 0.25,
                }],
            });
            model.add_root(root)?;
            let kernel = MonstertruckKernel::project_standard();
            let evaluated = kernel.evaluate(&model, root)?;
            let summary = kernel.summarize(&evaluated);
            assert!(summary.geometrically_consistent);
            assert_eq!(summary.shells, 1);
            assert!(summary.faces > 6);
            step_p21::parser::parse(&kernel.to_step(&evaluated)?)?;
            Ok(())
        }

        #[test]
        fn evaluates_round_edge_treatment() -> Result<()> {
            box_edge_treatment(EdgeTreatmentProfile::Round)
        }

        #[test]
        fn evaluates_chamfer_edge_treatment() -> Result<()> {
            box_edge_treatment(EdgeTreatmentProfile::Chamfer)
        }

        #[test]
        fn evaluates_transverse_chamfer_on_thin_prism() -> Result<()> {
            let mut model = CadModel::new();
            let child = model.add_node(CadNode::Primitive(Primitive::Box {
                size_mm: [4.0, 3.0, 0.1],
            }));
            let root = model.add_node(CadNode::EdgeTreatment {
                child,
                profile: EdgeTreatmentProfile::Chamfer,
                edges: vec![EdgeTreatmentEdge {
                    selector: EdgeSelector::Segment {
                        start_mm: [0.0, 0.0, 0.0],
                        end_mm: [0.0, 0.0, 0.1],
                    },
                    radius_mm: 0.137,
                }],
            });
            model.add_root(root)?;
            let kernel = MonstertruckKernel::project_standard();
            let evaluated = kernel.evaluate(&model, root)?;
            let summary = kernel.summarize(&evaluated);
            assert!(summary.is_single_solid());
            assert_eq!(summary.faces, 7);
            step_p21::parser::parse(&kernel.to_step(&evaluated)?)?;
            Ok(())
        }

        #[test]
        fn evaluates_chamfers_at_reentrant_shoulder() -> Result<()> {
            let mut model = CadModel::new();
            let base = model.add_node(CadNode::Primitive(Primitive::Box {
                size_mm: [6.0, 5.0, 1.0],
            }));
            let upper = model.add_node(CadNode::Extrude {
                profile: Profile2d::polygon(vec![[1.0, 1.0], [5.0, 1.0], [5.0, 4.0], [1.0, 4.0]])?,
                vector_mm: [0.0, 0.0, 0.1],
            });
            let upper = model.add_node(CadNode::Transform {
                transform: RigidTransform::translation_mm([0.0, 0.0, 1.0]),
                child: upper,
            });
            let mut root = model.add_node(CadNode::Boolean {
                op: BooleanOp::Union,
                children: vec![base, upper],
            });
            // A later chamfer selects the exact seam generated by an unequal
            // earlier bevel; the unoriented carrier normal fixes its role.
            for ([x, y], profile, setback) in [
                (
                    [1.0, 1.0],
                    EdgeTreatmentProfile::CarrierChamfer {
                        first_carrier_normal: [0.0, -1.0, 0.0],
                        second_setback_mm: 0.4,
                    },
                    0.25,
                ),
                ([1.25, 1.0], EdgeTreatmentProfile::Chamfer, 0.05),
                ([5.0, 4.0], EdgeTreatmentProfile::Chamfer, 0.25),
            ] {
                root = model.add_node(CadNode::EdgeTreatment {
                    child: root,
                    profile,
                    edges: vec![EdgeTreatmentEdge {
                        selector: EdgeSelector::Segment {
                            start_mm: [x, y, 1.0],
                            end_mm: [x, y, 1.1],
                        },
                        radius_mm: setback,
                    }],
                });
            }
            model.add_root(root)?;
            let kernel = MonstertruckKernel::project_standard();
            let evaluated = kernel.evaluate(&model, root)?;
            assert!(kernel.summarize(&evaluated).is_single_solid());
            assert!(
                evaluated
                    .solid
                    .face_iter()
                    .all(|f| matches!(f.surface(), Surface::Plane(_)))
            );
            assert!(
                evaluated
                    .solid
                    .edge_iter()
                    .all(|e| matches!(e.curve(), Curve::Line(_)))
            );
            let mut reference = CadModel::new();
            let base = reference.add_node(CadNode::Primitive(Primitive::Box {
                size_mm: [6.0, 5.0, 1.0],
            }));
            let upper = reference.add_node(CadNode::Extrude {
                profile: Profile2d::polygon(vec![
                    [1.30, 1.0],
                    [5.0, 1.0],
                    [5.0, 3.75],
                    [4.75, 4.0],
                    [1.0, 4.0],
                    [1.0, 1.4],
                    [
                        1.25 - 0.05 * 0.25 / 0.25_f64.hypot(0.4),
                        1.0 + 0.05 * 0.4 / 0.25_f64.hypot(0.4),
                    ],
                ])?,
                vector_mm: [0.0, 0.0, 0.1],
            });
            let upper = reference.add_node(CadNode::Transform {
                transform: RigidTransform::translation_mm([0.0, 0.0, 1.0]),
                child: upper,
            });
            let reference_root = reference.add_node(CadNode::Boolean {
                op: BooleanOp::Union,
                children: vec![base, upper],
            });
            reference.add_root(reference_root)?;
            let expected = kernel.evaluate(&reference, reference_root)?;
            assert_same_analytic_boundary(&kernel, &evaluated, &expected)
        }

        fn assert_same_analytic_boundary(
            kernel: &MonstertruckKernel,
            evaluated: &EvaluatedShape,
            expected: &EvaluatedShape,
        ) -> Result<()> {
            let actual_step = kernel.to_step(evaluated)?;
            let expected_step = kernel.to_step(expected)?;
            let actual_exchange =
                monstertruck_io::step::load::ParsedStepSource::parse(actual_step.as_bytes())?
                    .exchange;
            let expected_exchange =
                monstertruck_io::step::load::ParsedStepSource::parse(expected_step.as_bytes())?
                    .exchange;
            let actual = &actual_exchange.data[0].entities;
            let expected = &expected_exchange.data[0].entities;
            let solid_id = |entities: &[step_p21::ast::EntityInstance]| -> Result<u64> {
                let ids = entities
                    .iter()
                    .filter(|entity| {
                        crate::step_graph::simple_record(entity)
                            .is_some_and(|record| record.name == "MANIFOLD_SOLID_BREP")
                    })
                    .map(crate::step_graph::entity_id)
                    .collect::<Vec<_>>();
                let [id] = ids.as_slice() else {
                    bail!("expected one solid")
                };
                Ok(*id)
            };
            let comparison = crate::boundary_proof::compare_analytic_boundaries(
                solid_id(actual)?,
                actual,
                solid_id(expected)?,
                expected,
            );
            assert!(
                comparison.equivalent,
                "replayed solid differs from explicit reference: {comparison:?}"
            );
            Ok(())
        }

        #[test]
        fn evaluates_partitioned_longitudinal_round() -> Result<()> {
            let mut model = CadModel::new();
            let mut children = Vec::new();
            for (z, height) in [(0.0, 0.125), (0.125, 0.125), (0.25, 0.25)] {
                let child = model.add_node(CadNode::Primitive(Primitive::Box {
                    size_mm: [4.0, 3.0, height],
                }));
                children.push(model.add_node(CadNode::Transform {
                    transform: RigidTransform::translation_mm([0.0, 0.0, z]),
                    child,
                }));
            }
            let child = model.add_node(CadNode::Boolean {
                op: BooleanOp::Union,
                children,
            });
            let kernel = MonstertruckKernel::project_standard();
            let original = kernel.evaluate(&model, child)?;
            let selector = EdgeSelector::Segment {
                start_mm: [0.0, 0.0, 0.0],
                end_mm: [0.0, 0.0, 0.5],
            };
            assert!(
                select_edges(&original.solid.boundaries()[0], &selector).is_err(),
                "fixture must contain partitioned seam"
            );
            let root = model.add_node(CadNode::EdgeTreatment {
                child,
                profile: EdgeTreatmentProfile::Round,
                edges: vec![EdgeTreatmentEdge {
                    selector,
                    radius_mm: 0.75,
                }],
            });
            model.add_root(root)?;
            let evaluated = kernel.evaluate(&model, root)?;
            assert!(kernel.summarize(&evaluated).is_single_solid());
            assert_eq!(kernel.summarize(&evaluated).faces, 7);

            let points = [[0.75, 0.0], [4.0, 0.0], [4.0, 3.0], [0.0, 3.0], [0.0, 0.75]];
            let mut curves = points
                .windows(2)
                .map(|p| Curve2d::Line {
                    start_mm: p[0],
                    end_mm: p[1],
                })
                .collect::<Vec<_>>();
            curves.push(Curve2d::CircleArc {
                center_mm: [0.75, 0.75],
                radius_mm: 0.75,
                start_angle_rad: std::f64::consts::PI,
                end_angle_rad: 3.0 * std::f64::consts::FRAC_PI_2,
            });
            let mut reference = CadModel::new();
            let root = reference.add_node(CadNode::Extrude {
                profile: Profile2d {
                    loops: vec![crate::cad_ir::ProfileLoop { curves }],
                },
                vector_mm: [0.0, 0.0, 0.5],
            });
            reference.add_root(root)?;
            let expected = kernel.evaluate(&reference, root)?;
            assert_same_analytic_boundary(&kernel, &evaluated, &expected)
        }

        #[test]
        fn evaluates_asymmetric_chamfer_edge_treatment() -> Result<()> {
            box_edge_treatment(EdgeTreatmentProfile::AsymmetricChamfer {
                first_setback_mm: 0.2,
                second_setback_mm: 0.3,
            })
        }

        #[test]
        fn bounded_circle_edge_selector_distinguishes_complements_and_split_arcs() -> Result<()> {
            use std::f64::consts::{FRAC_PI_2, PI, TAU};
            let profile = Profile2d {
                loops: vec![ProfileLoop {
                    curves: (0..4)
                        .map(|i| Curve2d::CircleArc {
                            center_mm: [0.0, 0.0],
                            radius_mm: 2.0,
                            start_angle_rad: f64::from(i) * FRAC_PI_2,
                            end_angle_rad: f64::from(i + 1) * FRAC_PI_2,
                        })
                        .collect(),
                }],
            };
            let solid = extrude_profile_z(&profile, 3.0)?;
            let shell = &solid.boundaries()[0];
            let arc = |normal, sweep_rad| EdgeSelector::CircleArc {
                center_mm: [0.0, 0.0, 3.0],
                normal,
                radius_mm: 2.0,
                start_mm: [2.0, 0.0, 3.0],
                sweep_rad,
            };
            let short = select_edges(shell, &arc([0.0, 0.0, 1.0], FRAC_PI_2))?;
            let long = select_edges(shell, &arc([0.0, 0.0, -1.0], 3.0 * FRAC_PI_2))?;
            assert_eq!(short.len(), 1);
            assert_eq!(long.len(), 3);
            assert!(long.iter().all(|edge| edge.id() != short[0].id()));
            assert_eq!(select_edges(shell, &arc([0.0, 0.0, 1.0], PI))?.len(), 2);
            assert!(select_edges(shell, &arc([0.0, 0.0, 1.0], FRAC_PI_2 / 2.0)).is_err());
            assert!(select_edges(shell, &arc([0.0, 0.0, 1.0], TAU)).is_err());
            Ok(())
        }

        #[test]
        fn evaluates_coupled_bounded_arc_and_line_perimeter() -> Result<()> {
            use std::f64::consts::{FRAC_PI_2, PI};
            let arc = |center_mm, start_angle_rad, end_angle_rad| Curve2d::CircleArc {
                center_mm,
                radius_mm: 1.0,
                start_angle_rad,
                end_angle_rad,
            };
            let profile = Profile2d {
                loops: vec![ProfileLoop {
                    curves: vec![
                        Curve2d::Line {
                            start_mm: [0.0, -1.0],
                            end_mm: [4.0, -1.0],
                        },
                        arc([4.0, 0.0], -FRAC_PI_2, 0.0),
                        arc([4.0, 0.0], 0.0, FRAC_PI_2),
                        Curve2d::Line {
                            start_mm: [4.0, 1.0],
                            end_mm: [0.0, 1.0],
                        },
                        arc([0.0, 0.0], FRAC_PI_2, PI),
                        arc([0.0, 0.0], PI, 3.0 * FRAC_PI_2),
                    ],
                }],
            };
            let end_arc = |x, y| EdgeSelector::CircleArc {
                center_mm: [x, 0.0, 3.0],
                normal: [0.0, 0.0, 1.0],
                radius_mm: 1.0,
                start_mm: [x, y, 3.0],
                sweep_rad: PI,
            };
            let selectors = [
                EdgeSelector::Segment {
                    start_mm: [0.0, -1.0, 3.0],
                    end_mm: [4.0, -1.0, 3.0],
                },
                end_arc(4.0, -1.0),
                EdgeSelector::Segment {
                    start_mm: [4.0, 1.0, 3.0],
                    end_mm: [0.0, 1.0, 3.0],
                },
                end_arc(0.0, 1.0),
            ];
            let solid = extrude_profile_z(&profile, 3.0)?;
            let shell = &solid.boundaries()[0];
            let mut selected = Vec::new();
            for selector in &selectors {
                selected.extend(select_edges(shell, selector)?);
            }
            assert_eq!(selected.len(), 6);
            let mut model = CadModel::new();
            let child = model.add_node(CadNode::Extrude {
                profile,
                vector_mm: [0.0, 0.0, 3.0],
            });
            let root = model.add_node(CadNode::EdgeTreatment {
                child,
                profile: EdgeTreatmentProfile::Round,
                edges: selectors
                    .into_iter()
                    .map(|selector| EdgeTreatmentEdge {
                        selector,
                        radius_mm: 0.15,
                    })
                    .collect(),
            });
            model.add_root(root)?;
            let kernel = MonstertruckKernel::project_standard();
            let evaluated = kernel.evaluate(&model, root)?;
            let summary = kernel.summarize(&evaluated);
            assert!(summary.geometrically_consistent);
            assert_eq!(summary.shells, 1);
            assert_eq!(summary.faces, 14);
            step_p21::parser::parse(&kernel.to_step(&evaluated)?)?;
            Ok(())
        }

        fn tangent_rim_box(model: &mut CadModel) -> NodeId {
            let child = model.add_node(CadNode::Primitive(Primitive::Box {
                size_mm: [6.0, 5.0, 4.0],
            }));
            model.add_node(CadNode::EdgeTreatment {
                child,
                profile: EdgeTreatmentProfile::Round,
                edges: [[0.0, 0.0], [6.0, 0.0], [6.0, 5.0], [0.0, 5.0]]
                    .into_iter()
                    .map(|[x, y]| EdgeTreatmentEdge {
                        selector: EdgeSelector::Segment {
                            start_mm: [x, y, 0.0],
                            end_mm: [x, y, 4.0],
                        },
                        radius_mm: 1.0,
                    })
                    .collect(),
            })
        }

        fn tangent_rim_selectors(z: f64, offset: [f64; 2]) -> Vec<EdgeTreatmentEdge> {
            use std::f64::consts::FRAC_PI_2;
            let points = [
                [1.0, 0.0],
                [5.0, 0.0],
                [6.0, 1.0],
                [6.0, 4.0],
                [5.0, 5.0],
                [1.0, 5.0],
                [0.0, 4.0],
                [0.0, 1.0],
            ];
            let centers = [[5.0, 1.0], [5.0, 4.0], [1.0, 4.0], [1.0, 1.0]];
            let at_top = |[x, y]: [f64; 2]| [x + offset[0], y + offset[1], z];
            let mut edges = Vec::new();
            for i in 0..4 {
                edges.push(EdgeTreatmentEdge {
                    selector: EdgeSelector::Segment {
                        start_mm: at_top(points[2 * i]),
                        end_mm: at_top(points[2 * i + 1]),
                    },
                    radius_mm: 0.25,
                });
                edges.push(EdgeTreatmentEdge {
                    selector: EdgeSelector::CircleArc {
                        center_mm: at_top(centers[i]),
                        normal: [0.0, 0.0, 1.0],
                        radius_mm: 1.0,
                        start_mm: at_top(points[2 * i + 1]),
                        sweep_rad: FRAC_PI_2,
                    },
                    radius_mm: 0.25,
                });
            }
            edges
        }

        #[test]
        fn evaluates_native_mixed_tangent_rim() -> Result<()> {
            let mut model = CadModel::new();
            let core = tangent_rim_box(&mut model);
            let mut edges = tangent_rim_selectors(4.0, [0.0; 2]);
            let kernel = MonstertruckKernel::project_standard();
            let root = model.add_node(CadNode::EdgeTreatment {
                child: core,
                profile: EdgeTreatmentProfile::Round,
                edges: edges.clone(),
            });
            let evaluated = kernel.evaluate(&model, root)?;
            let summary = kernel.summarize(&evaluated);
            assert!(summary.is_single_solid());
            assert_eq!(summary.faces, 18);
            let torus_count = evaluated
                .solid
                .face_iter()
                .filter(|face| {
                    if let Surface::Torus(torus) = face.surface() {
                        assert!(
                            (torus.entity().large_radius() - 0.75).abs()
                                < crate::numeric::NUMERICAL_GEOMETRY_TOL_MM
                        );
                        assert!(
                            (torus.entity().small_radius() - 0.25).abs()
                                < crate::numeric::NUMERICAL_GEOMETRY_TOL_MM
                        );
                        true
                    } else {
                        false
                    }
                })
                .count();
            assert_eq!(torus_count, 4, "native analytic tori must survive");
            edges.reverse();
            let reversed = model.add_node(CadNode::EdgeTreatment {
                child: core,
                profile: EdgeTreatmentProfile::Round,
                edges,
            });
            let expected = kernel.evaluate(&model, reversed)?;
            // Exact native carrier incidence and deterministic complete
            // exported topology are the order-invariance contract here.
            // This also checks deterministic native serialization.
            monstertruck_modeling::certify_native_analytic_boundary_incidence(
                &evaluated.solid,
                crate::numeric::NUMERICAL_GEOMETRY_TOL_MM,
                DEFAULT_STEP_RECOVERY_TOLERANCES.edge_treatment_direction,
            )
            .map_err(|error| anyhow::anyhow!("native rim incidence: {error:?}"))?;
            let first = monstertruck_io::step::load::ParsedStepSource::parse(
                kernel.to_step(&evaluated)?.as_bytes(),
            )?
            .exchange;
            let second = monstertruck_io::step::load::ParsedStepSource::parse(
                kernel.to_step(&expected)?.as_bytes(),
            )?
            .exchange;
            assert_eq!(first.data[0].entities, second.data[0].entities);
            Ok(())
        }

        #[test]
        fn evaluates_translated_reentrant_tangent_rim() -> Result<()> {
            let mut model = CadModel::new();
            let core = tangent_rim_box(&mut model);
            let kernel = MonstertruckKernel::project_standard();
            let base = model.add_node(CadNode::Primitive(Primitive::Box {
                size_mm: [10.0, 9.0, 1.0],
            }));
            let boss = model.add_node(CadNode::Transform {
                transform: RigidTransform::translation_mm([2.0, 2.0, 1.0]),
                child: core,
            });
            let shoulder = model.add_node(CadNode::Boolean {
                op: BooleanOp::Union,
                children: vec![base, boss],
            });
            let round = model.add_node(CadNode::EdgeTreatment {
                child: shoulder,
                profile: EdgeTreatmentProfile::Round,
                edges: tangent_rim_selectors(1.0, [2.0; 2]),
            });
            let reentrant = kernel.evaluate(&model, round)?;
            assert!(
                kernel.summarize(&reentrant).is_single_solid(),
                "{:?}",
                kernel.summarize(&reentrant)
            );
            let tori = reentrant
                .solid
                .face_iter()
                .filter(|face| {
                    if let Surface::Torus(torus) = face.oriented_surface() {
                        assert!(
                            (torus.entity().large_radius() - 1.25).abs()
                                < crate::numeric::NUMERICAL_GEOMETRY_TOL_MM
                        );
                        assert!(
                            (torus.entity().small_radius() - 0.25).abs()
                                < crate::numeric::NUMERICAL_GEOMETRY_TOL_MM
                        );
                        assert!(!torus.orientation());
                        true
                    } else {
                        false
                    }
                })
                .count();
            assert_eq!(tori, 4, "reentrant rim must retain inverted analytic tori");
            monstertruck_modeling::certify_native_analytic_boundary_incidence(
                &reentrant.solid,
                crate::numeric::NUMERICAL_GEOMETRY_TOL_MM,
                DEFAULT_STEP_RECOVERY_TOLERANCES.edge_treatment_direction,
            )
            .map_err(|error| anyhow::anyhow!("native reentrant rim incidence: {error:?}"))?;
            Ok(())
        }

        #[test]
        fn evaluates_three_way_spherical_round_network() -> Result<()> {
            use crate::brep::SurfaceSupport;
            let mut model = CadModel::new();
            let child = model.add_node(CadNode::Primitive(Primitive::Box {
                size_mm: [6.0, 5.0, 4.0],
            }));
            let corners = [
                [0.0, 0.0, 4.0],
                [6.0, 0.0, 4.0],
                [6.0, 5.0, 4.0],
                [0.0, 5.0, 4.0],
            ];
            let mut edges = Vec::new();
            for (index, point) in corners.iter().copied().enumerate() {
                for end_mm in [corners[(index + 1) % 4], [point[0], point[1], 0.0]] {
                    edges.push(EdgeTreatmentEdge {
                        selector: EdgeSelector::Segment {
                            start_mm: point,
                            end_mm,
                        },
                        radius_mm: 0.5,
                    });
                }
            }
            let root = model.add_node(CadNode::EdgeTreatment {
                child,
                profile: EdgeTreatmentProfile::Round,
                edges,
            });
            model.add_root(root)?;
            let kernel = MonstertruckKernel::project_standard();
            let evaluated = kernel.evaluate(&model, root)?;
            let summary = kernel.summarize(&evaluated);
            assert!(summary.is_single_solid(), "{summary:?}");
            assert_eq!(summary.faces, 18);
            let step = kernel.to_step(&evaluated)?;
            let exchange = step_p21::parser::parse(&step)?;
            let section = &exchange.data[0];
            let index = crate::step_graph::build_index(&section.entities);
            let solid_id = section
                .entities
                .iter()
                .find_map(|entity| {
                    let record = crate::step_graph::simple_record(entity)?;
                    (record.name == "MANIFOLD_SOLID_BREP")
                        .then(|| crate::step_graph::entity_id(entity))
                })
                .context("missing rounded solid")?;
            let solid =
                crate::reduced_brep::ReducedSolid::from_step(solid_id, &section.entities, &index)
                    .context("cannot resolve analytic rounded solid")?;
            assert_eq!(
                solid
                    .faces()
                    .filter(|face| matches!(face.surface, SurfaceSupport::Sphere(_)))
                    .count(),
                4
            );
            assert_eq!(
                solid
                    .faces()
                    .filter(|face| matches!(face.surface, SurfaceSupport::Cylinder(_)))
                    .count(),
                8
            );
            assert!(solid.validate_closed_topology().is_some());
            Ok(())
        }

        #[test]
        fn evaluates_round_circle_edge_treatment() -> Result<()> {
            let mut model = CadModel::new();
            let cylinder = model.add_node(CadNode::Primitive(Primitive::Cylinder {
                radius_mm: 1.0,
                height_mm: 2.0,
            }));
            let root = model.add_node(CadNode::EdgeTreatment {
                child: cylinder,
                profile: EdgeTreatmentProfile::Round,
                edges: vec![EdgeTreatmentEdge {
                    selector: EdgeSelector::Circle {
                        center_mm: [0.0, 2.0, 0.0],
                        normal: [0.0, 1.0, 0.0],
                        radius_mm: 1.0,
                    },
                    radius_mm: 0.1,
                }],
            });
            model.add_root(root)?;
            let kernel = MonstertruckKernel::project_standard();
            let evaluated = kernel.evaluate(&model, root)?;
            let summary = kernel.summarize(&evaluated);
            assert!(summary.geometrically_consistent);
            assert_eq!(summary.shells, 1);
            assert!(summary.faces > 3);
            step_p21::parser::parse(&kernel.to_step(&evaluated)?)?;
            Ok(())
        }

        #[test]
        fn project_standard_maps_exactly_to_point_zero_one_mm() {
            assert_eq!(
                monstertruck_policy(GeometryPolicy::project_standard()).tolerance(),
                Some(0.01)
            );
        }

        fn shared_tolerance_graph(
            first_tolerance_mm: f64,
            second_tolerance_mm: f64,
        ) -> Result<(CadModel, NodeId, NodeId, NodeId)> {
            let mut model = CadModel::new();
            let shared = model.add_node(CadNode::Primitive(Primitive::Box {
                size_mm: [1.0, 2.0, 3.0],
            }));
            let first = model.add_node(CadNode::SolverTolerance {
                child: shared,
                tolerance_mm: first_tolerance_mm,
            });
            let second = model.add_node(CadNode::SolverTolerance {
                child: shared,
                tolerance_mm: second_tolerance_mm,
            });
            let root = model.add_node(CadNode::Boolean {
                op: BooleanOp::Union,
                children: vec![first, second],
            });
            model.add_root(root)?;
            model.validate()?;
            Ok((model, first, second, root))
        }

        #[test]
        fn one_shot_cache_only_retains_shared_geometry_until_last_use() -> Result<()> {
            let (model, first, second, root) = shared_tolerance_graph(0.01, 0.01)?;
            let policy = GeometryPolicy::project_standard();
            let mut cache = EvaluationCache::single_root(&model, root, policy)?;

            let first_shape = evaluate_node(&model, first, policy, &mut cache)?;
            assert_eq!(
                cache.solids.len(),
                1,
                "only the shared descendant is retained"
            );
            let vertex = first_shape
                .vertex_iter()
                .next()
                .context("box has no vertices")?;
            let original = vertex.point();
            vertex.set_point(Point3::new(100.0, 100.0, 100.0));
            let second_shape = evaluate_node(&model, second, policy, &mut cache)?;
            assert_eq!(
                second_shape.vertex_iter().next().map(|v| v.point()),
                Some(original)
            );
            assert!(
                cache.solids.is_empty(),
                "last consumer moves owned cached geometry"
            );
            Ok(())
        }

        #[test]
        fn one_shot_cache_distinguishes_numerical_policy_and_drops_linear_nodes() -> Result<()> {
            let (model, first, second, root) = shared_tolerance_graph(0.001, 0.01)?;
            let policy = GeometryPolicy::project_standard();
            let mut cache = EvaluationCache::single_root(&model, root, policy)?;
            for id in [first, second] {
                let evaluated = evaluate_node(&model, id, policy, &mut cache)?;
                assert_eq!(evaluated.boundaries()[0].len(), 6);
                assert!(
                    cache.solids.is_empty(),
                    "different policy keys are not shared"
                );
            }
            Ok(())
        }

        #[test]
        fn evaluation_session_memoizes_shared_nodes_per_numerical_policy() -> Result<()> {
            let mut model = CadModel::new();
            let box_node = model.add_node(CadNode::Primitive(Primitive::Box {
                size_mm: [2.0, 3.0, 4.0],
            }));
            let tight = model.add_node(CadNode::SolverTolerance {
                child: box_node,
                tolerance_mm: 1.0e-5,
            });
            let loose = model.add_node(CadNode::SolverTolerance {
                child: box_node,
                tolerance_mm: 0.01,
            });
            model.add_root(tight)?;
            model.add_root(loose)?;

            let kernel = MonstertruckKernel::project_standard();
            let mut session = kernel.session(&model)?;
            assert_eq!(kernel.summarize(&session.evaluate(tight)?).faces, 6);
            assert_eq!(session.cache.solids.len(), 2);
            assert_eq!(kernel.summarize(&session.evaluate(tight)?).faces, 6);
            assert_eq!(
                session.cache.solids.len(),
                2,
                "repeat must reuse both nodes"
            );
            assert_eq!(kernel.summarize(&session.evaluate(loose)?).faces, 6);
            assert_eq!(
                session.cache.solids.len(),
                4,
                "same child under a distinct SolverTolerance must be evaluated separately"
            );
            Ok(())
        }

        #[test]
        fn evaluation_session_does_not_share_mutable_geometry_with_cached_shapes() -> Result<()> {
            let mut model = CadModel::new();
            let root = model.add_node(CadNode::Primitive(Primitive::Box {
                size_mm: [1.0, 2.0, 3.0],
            }));
            model.add_root(root)?;
            let mut session = MonstertruckKernel::project_standard().session(&model)?;

            let initial = session.evaluate(root)?;
            let first_vertex = initial
                .solid
                .vertex_iter()
                .next()
                .ok_or_else(|| anyhow::anyhow!("box must have vertices"))?;
            let expected = first_vertex.point();
            first_vertex.set_point(Point3::new(99.0, 98.0, 97.0));

            let replay = session.evaluate(root)?;
            let replay_vertex = replay
                .solid
                .vertex_iter()
                .next()
                .ok_or_else(|| anyhow::anyhow!("cached box must have vertices"))?;
            assert_eq!(
                replay_vertex.point(),
                expected,
                "first evaluation may mutate its geometry without corrupting the cache"
            );
            replay_vertex.set_point(Point3::new(96.0, 95.0, 94.0));
            let replay_again = session.evaluate(root)?;
            assert_eq!(
                replay_again
                    .solid
                    .vertex_iter()
                    .next()
                    .map(|vertex| vertex.point()),
                Some(expected),
                "cache hits must also own independent mutable geometry"
            );
            Ok(())
        }

        #[test]
        fn equivalent_refuses_approximate_boolean_backend() -> Result<()> {
            let mut model = CadModel::new();
            let first = model.add_node(CadNode::Primitive(Primitive::Box {
                size_mm: [2.0, 2.0, 2.0],
            }));
            let second_body = model.add_node(CadNode::Primitive(Primitive::Box {
                size_mm: [2.0, 2.0, 2.0],
            }));
            let second = model.add_node(CadNode::Transform {
                transform: RigidTransform::translation_mm([1.0, 0.0, 0.0]),
                child: second_body,
            });
            let root = model.add_node(CadNode::Boolean {
                op: BooleanOp::Union,
                children: vec![first, second],
            });
            model.add_root(root)?;

            let error = match MonstertruckKernel::equivalent().evaluate(&model, root) {
                Ok(_) => bail!("Equivalent unexpectedly admitted the approximate Boolean backend"),
                Err(error) => error,
            };
            assert!(
                error.to_string().contains("no equivalence-proving backend"),
                "unexpected refusal: {error:#}"
            );
            Ok(())
        }

        #[test]
        fn evaluates_adjacent_cylinders_union() -> Result<()> {
            let mut model = CadModel::new();
            let first = model.add_node(CadNode::Primitive(Primitive::Cylinder {
                radius_mm: 1.0,
                height_mm: 1.0,
            }));
            let second_body = model.add_node(CadNode::Primitive(Primitive::Cylinder {
                radius_mm: 1.0,
                height_mm: 1.0,
            }));
            let second = model.add_node(CadNode::Transform {
                transform: RigidTransform::translation_mm([0.0, 1.0, 0.0]),
                child: second_body,
            });
            let root = model.add_node(CadNode::Boolean {
                op: BooleanOp::Union,
                children: vec![first, second],
            });
            model.add_root(root)?;
            let kernel = MonstertruckKernel::project_standard();
            let evaluated = kernel.evaluate(&model, root)?;
            let summary = kernel.summarize(&evaluated);
            assert!(summary.geometrically_consistent);
            assert_eq!(summary.shells, 1);
            Ok(())
        }

        #[test]
        fn evaluates_boolean_difference() -> Result<()> {
            use crate::cad_ir::BooleanOp;

            let mut model = CadModel::new();
            let base = model.add_node(CadNode::Extrude {
                profile: Profile2d::polygon(vec![
                    [0.0, 0.0],
                    [10.0, 0.0],
                    [10.0, 6.0],
                    [0.0, 6.0],
                ])?,
                vector_mm: [0.0, 0.0, 2.0],
            });
            let cutter_body = model.add_node(CadNode::Extrude {
                profile: Profile2d::polygon(vec![[3.0, 2.0], [7.0, 2.0], [7.0, 4.0], [3.0, 4.0]])?,
                vector_mm: [0.0, 0.0, 4.0],
            });
            let cutter = model.add_node(CadNode::Transform {
                transform: RigidTransform::translation_mm([0.0, 0.0, -1.0]),
                child: cutter_body,
            });
            let root = model.add_node(CadNode::Boolean {
                op: BooleanOp::Difference,
                children: vec![base, cutter],
            });
            model.add_root(root)?;

            let kernel = MonstertruckKernel::project_standard();
            let evaluated = kernel.evaluate(&model, root)?;
            let summary = kernel.summarize(&evaluated);
            assert!(summary.geometrically_consistent);
            assert!(summary.faces >= 8);
            step_p21::parser::parse(&kernel.to_step(&evaluated)?)?;
            Ok(())
        }

        #[test]
        fn mixed_difference_batch_preserves_certified_contact() -> Result<()> {
            let kernel = MonstertruckKernel::project_standard();
            let host = primitive_solid(&Primitive::Box {
                size_mm: [10.0, 6.0, 2.0],
            })?;
            let blind = builder::translated(
                &primitive_solid(&Primitive::Box {
                    size_mm: [1.0, 1.0, 0.01],
                })?,
                Vector3::new(1.0, 1.0, 1.99),
            );
            let through = builder::translated(
                &primitive_solid(&Primitive::Box {
                    size_mm: [1.0, 1.0, 4.0],
                })?,
                Vector3::new(6.0, 2.0, -1.0),
            );
            let policy = monstertruck_policy(kernel.policy);
            let op = monstertruck_solid::PlanarContactOp::Difference;
            assert!(
                monstertruck_solid::try_planar_contact_many(
                    &host,
                    &[&blind, &through],
                    op,
                    policy
                )?
                .is_none()
            );
            assert!(
                monstertruck_solid::try_planar_contact_many(&host, &[&blind], op, policy)?
                    .is_some()
            );
            let solid = fold_boolean(
                BooleanOp::Difference,
                [Ok(host), Ok(blind), Ok(through)].into_iter(),
                kernel.policy,
            )?;
            let actual = EvaluatedShape { solid };
            assert!(kernel.summarize(&actual).is_single_solid());

            // Independent sectional construction: the lower section has a
            // through-hole; the upper section has the same hole plus a pocket.
            let outer = Profile2d::polygon(vec![[0.0, 0.0], [10.0, 0.0], [10.0, 6.0], [0.0, 6.0]])?;
            let hole = Profile2d::polygon(vec![[6.0, 2.0], [6.0, 3.0], [7.0, 3.0], [7.0, 2.0]])?;
            let pocket = Profile2d::polygon(vec![[1.0, 1.0], [1.0, 2.0], [2.0, 2.0], [2.0, 1.0]])?;
            let mut lower = outer;
            lower.loops.extend(hole.loops);
            let mut upper = lower.clone();
            upper.loops.extend(pocket.loops);
            let bottom = extrude_profile_z(&lower, 1.99)?;
            let top = builder::translated(
                &extrude_profile_z(&upper, 0.01)?,
                Vector3::new(0.0, 0.0, 1.99),
            );
            let expected = EvaluatedShape {
                solid: monstertruck_solid::try_planar_contact_union(&bottom, &top, policy)?
                    .context("sectional reference must have certified planar contact")?,
            };
            assert_same_analytic_boundary(&kernel, &actual, &expected)
        }

        #[test]
        fn evaluates_full_profile_revolution() -> Result<()> {
            let mut model = CadModel::new();
            let profile =
                Profile2d::polygon(vec![[1.0, -1.0], [2.0, -1.0], [2.0, 1.0], [1.0, 1.0]])?;
            let root = model.add_node(CadNode::Revolve {
                profile,
                axis: Axis3 {
                    origin_mm: [0.0, 0.0, 0.0],
                    direction: [0.0, 1.0, 0.0],
                },
                angle_rad: std::f64::consts::TAU,
            });
            model.add_root(root)?;

            let kernel = MonstertruckKernel::project_standard();
            let evaluated = kernel.evaluate(&model, root)?;
            let summary = kernel.summarize(&evaluated);
            assert!(summary.geometrically_consistent);
            assert_eq!(summary.shells, 1);
            let step = kernel.to_step(&evaluated)?;
            step_p21::parser::parse(&step)?;
            let recovered = crate::detect_solid_revolutions_bytes(step.as_bytes())?;
            assert_eq!(recovered.len(), 1);
            let fragment = crate::cad_recovery::recover_solid_revolution_fragment(&recovered[0])?;
            let rebuilt = kernel.evaluate(&fragment.model, fragment.root)?;
            assert!(kernel.summarize(&rebuilt).geometrically_consistent);
            step_p21::parser::parse(&kernel.to_step(&rebuilt)?)?;
            Ok(())
        }

        #[test]
        fn evaluates_partial_profile_revolution() -> Result<()> {
            let mut model = CadModel::new();
            let profile = Profile2d::polygon(vec![[2.0, 0.0], [3.0, 0.0], [3.0, 1.0], [2.0, 1.0]])?;
            let root = model.add_node(CadNode::Revolve {
                profile,
                axis: Axis3 {
                    origin_mm: [0.0, 0.0, 0.0],
                    direction: [0.0, 1.0, 0.0],
                },
                angle_rad: std::f64::consts::FRAC_PI_2,
            });
            model.add_root(root)?;

            let kernel = MonstertruckKernel::project_standard();
            let evaluated = kernel.evaluate(&model, root)?;
            let summary = kernel.summarize(&evaluated);
            assert!(summary.geometrically_consistent);
            assert_eq!(summary.shells, 1);
            let step = kernel.to_step(&evaluated)?;
            step_p21::parser::parse(&step)?;
            Ok(())
        }

        #[test]
        fn evaluates_on_axis_profile_revolution() -> Result<()> {
            let mut model = CadModel::new();
            let profile =
                Profile2d::polygon(vec![[0.0, -1.0], [2.0, -1.0], [2.0, 1.0], [0.0, 1.0]])?;
            let root = model.add_node(CadNode::Revolve {
                profile,
                axis: Axis3 {
                    origin_mm: [0.0, 0.0, 0.0],
                    direction: [0.0, 1.0, 0.0],
                },
                angle_rad: std::f64::consts::TAU,
            });
            model.add_root(root)?;

            let kernel = MonstertruckKernel::project_standard();
            let evaluated = kernel.evaluate(&model, root)?;
            let summary = kernel.summarize(&evaluated);
            assert!(summary.geometrically_consistent);
            assert_eq!(summary.shells, 1);
            let step = kernel.to_step(&evaluated)?;
            step_p21::parser::parse(&step)?;
            let recovered = crate::detect_solid_revolutions_bytes(step.as_bytes())?;
            assert_eq!(recovered.len(), 1);
            let fragment = crate::cad_recovery::recover_solid_revolution_fragment(&recovered[0])?;
            let rebuilt = kernel.evaluate(&fragment.model, fragment.root)?;
            assert!(kernel.summarize(&rebuilt).geometrically_consistent);
            step_p21::parser::parse(&kernel.to_step(&rebuilt)?)?;
            Ok(())
        }

        #[test]
        fn evaluates_apex_cone_profile_revolution() -> Result<()> {
            let mut model = CadModel::new();
            let profile = Profile2d::polygon(vec![[0.0, -1.0], [2.0, 1.0], [0.0, 1.0]])?;
            let root = model.add_node(CadNode::Revolve {
                profile,
                axis: Axis3 {
                    origin_mm: [0.0, 0.0, 0.0],
                    direction: [0.0, 1.0, 0.0],
                },
                angle_rad: std::f64::consts::TAU,
            });
            model.add_root(root)?;

            let kernel = MonstertruckKernel::project_standard();
            let evaluated = kernel.evaluate(&model, root)?;
            assert!(kernel.summarize(&evaluated).geometrically_consistent);
            let step = kernel.to_step(&evaluated)?;
            step_p21::parser::parse(&step)?;

            // Re-recognition of whichever STEP surface grammar Monstertruck chooses
            // is a separate detector concern. This test proves exact construction and
            // exchange serialization without pinning the writer to one representation.
            Ok(())
        }

        #[test]
        fn evaluates_transformed_profile_revolution_without_post_transform_breakage() -> Result<()>
        {
            let mut model = CadModel::new();
            let profile = Profile2d::polygon(vec![[0.0, 0.0], [2.0, 0.0], [2.0, 3.0], [0.0, 3.0]])?;
            let body = model.add_node(CadNode::Revolve {
                profile,
                axis: Axis3 {
                    origin_mm: [0.0, 0.0, 0.0],
                    direction: [0.0, 1.0, 0.0],
                },
                angle_rad: std::f64::consts::TAU,
            });
            let root = model.add_node(CadNode::Transform {
                transform: RigidTransform {
                    matrix: [
                        [0.0, 0.0, 1.0, 10.0],
                        [1.0, 0.0, 0.0, 20.0],
                        [0.0, 1.0, 0.0, 30.0],
                        [0.0, 0.0, 0.0, 1.0],
                    ],
                },
                child: body,
            });
            model.add_root(root)?;

            let kernel = MonstertruckKernel::project_standard();
            let evaluated = kernel.evaluate(&model, root)?;
            let summary = kernel.summarize(&evaluated);
            assert!(summary.geometrically_consistent);
            assert_eq!(summary.shells, 1);
            step_p21::parser::parse(&kernel.to_step(&evaluated)?)?;
            Ok(())
        }

        #[test]
        fn evaluates_quadratic_bezier_profile_extrusion() -> Result<()> {
            assert_curved_profile_evaluates(
                Curve2d::Bezier {
                    control_points_mm: vec![[0.0, 0.0], [1.0, 1.0], [2.0, 0.0]],
                },
                "B_SPLINE_CURVE",
            )
        }

        #[test]
        fn evaluates_quadratic_bspline_profile_extrusion() -> Result<()> {
            assert_curved_profile_evaluates(
                Curve2d::BSpline {
                    degree: 2,
                    control_points_mm: vec![[0.0, 0.0], [1.0, 1.0], [2.0, 0.0]],
                    knots: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
                    weights: None,
                },
                "B_SPLINE_CURVE",
            )
        }

        #[test]
        fn evaluates_rational_bspline_profile_extrusion() -> Result<()> {
            assert_curved_profile_evaluates(
                Curve2d::BSpline {
                    degree: 2,
                    control_points_mm: vec![[0.0, 0.0], [1.0, 1.0], [2.0, 0.0]],
                    knots: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
                    weights: Some(vec![1.0, 0.75, 1.0]),
                },
                "RATIONAL_B_SPLINE_CURVE",
            )
        }

        #[test]
        fn roundtrips_holed_extrusion_through_recovery() -> Result<()> {
            let mut model = CadModel::new();
            let profile = Profile2d {
                loops: vec![
                    ProfileLoop {
                        curves: vec![
                            Curve2d::Line {
                                start_mm: [0.0, 0.0],
                                end_mm: [8.0, 0.0],
                            },
                            Curve2d::Line {
                                start_mm: [8.0, 0.0],
                                end_mm: [8.0, 6.0],
                            },
                            Curve2d::Line {
                                start_mm: [8.0, 6.0],
                                end_mm: [0.0, 6.0],
                            },
                            Curve2d::Line {
                                start_mm: [0.0, 6.0],
                                end_mm: [0.0, 0.0],
                            },
                        ],
                    },
                    ProfileLoop {
                        curves: vec![Curve2d::CircleArc {
                            center_mm: [4.0, 3.0],
                            radius_mm: 1.25,
                            start_angle_rad: 0.0,
                            end_angle_rad: -std::f64::consts::TAU,
                        }],
                    },
                ],
            };
            let root = model.add_node(CadNode::Extrude {
                profile,
                vector_mm: [0.0, 0.0, 2.5],
            });
            model.add_root(root)?;

            let kernel = MonstertruckKernel::project_standard();
            let evaluated = kernel.evaluate(&model, root)?;
            let summary = kernel.summarize(&evaluated);
            assert!(summary.geometrically_consistent);
            assert_eq!(summary.shells, 1);
            assert_eq!(summary.faces, 8);

            let step = kernel.to_step(&evaluated)?;
            step_p21::parser::parse(&step)?;
            let recovered = crate::detect_solid_extrusions_bytes(step.as_bytes())?;
            assert_eq!(recovered.len(), 1);
            assert_eq!(recovered[0].inner_profile_loops.len(), 1);
            assert_eq!(recovered[0].inner_profile_loops[0].len(), 2);
            assert!(
                recovered[0].inner_profile_loops[0]
                    .iter()
                    .all(|curve| matches!(
                        curve,
                        crate::solid_extrusions::RecoveredProfileCurve::BSpline {
                            weights: Some(_),
                            ..
                        }
                    ))
            );

            let fragment = crate::cad_recovery::recover_solid_extrusion_fragment(&recovered[0])?;
            let CadNode::Transform { child, .. } = fragment.model.node(fragment.root)? else {
                bail!("expected transform root");
            };
            let CadNode::Extrude { profile, .. } = fragment.model.node(*child)? else {
                bail!("expected extrusion child");
            };
            assert_eq!(profile.loops.len(), 2);

            let rebuilt = kernel.evaluate(&fragment.model, fragment.root)?;
            let rebuilt_summary = kernel.summarize(&rebuilt);
            assert!(rebuilt_summary.geometrically_consistent);
            assert_eq!(rebuilt_summary.faces, 8);
            step_p21::parser::parse(&kernel.to_step(&rebuilt)?)?;
            Ok(())
        }

        #[test]
        fn roundtrips_rational_bspline_extrusion_through_recovery() -> Result<()> {
            let (model, root) = curved_profile_model(Curve2d::BSpline {
                degree: 2,
                control_points_mm: vec![[0.0, 0.0], [1.0, 1.0], [2.0, 0.0]],
                knots: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
                weights: Some(vec![1.0, 0.75, 1.0]),
            })?;
            let kernel = MonstertruckKernel::project_standard();
            let evaluated = kernel.evaluate(&model, root)?;
            let step = kernel.to_step(&evaluated)?;
            let recovered = crate::detect_solid_extrusions_bytes(step.as_bytes())?;
            assert_eq!(recovered.len(), 1);
            assert!(recovered[0].profile_curves.iter().any(|curve| matches!(
                curve,
                crate::solid_extrusions::RecoveredProfileCurve::BSpline {
                    weights: Some(_),
                    ..
                }
            )));

            let fragment = crate::cad_recovery::recover_solid_extrusion_fragment(&recovered[0])?;
            let rebuilt = kernel.evaluate(&fragment.model, fragment.root)?;
            assert!(kernel.summarize(&rebuilt).geometrically_consistent);
            step_p21::parser::parse(&kernel.to_step(&rebuilt)?)?;
            Ok(())
        }

        #[test]
        fn malformed_bspline_profile_fails_closed() -> Result<()> {
            let mut model = CadModel::new();
            let profile = Profile2d {
                loops: vec![ProfileLoop {
                    curves: vec![
                        Curve2d::Line {
                            start_mm: [0.0, 0.0],
                            end_mm: [2.0, 0.0],
                        },
                        Curve2d::BSpline {
                            degree: usize::MAX,
                            control_points_mm: vec![[2.0, 0.0], [1.0, 1.0], [0.0, 0.0]],
                            knots: Vec::new(),
                            weights: None,
                        },
                    ],
                }],
            };
            let root = model.add_node(CadNode::Extrude {
                profile,
                vector_mm: [0.0, 0.0, 1.0],
            });
            model.add_root(root)?;
            let kernel = MonstertruckKernel::project_standard();
            assert!(kernel.evaluate(&model, root).is_err());
            Ok(())
        }

        #[test]
        fn evaluates_semicircular_profile_extrusion() -> Result<()> {
            let mut model = CadModel::new();
            let profile = Profile2d {
                loops: vec![ProfileLoop {
                    curves: vec![
                        Curve2d::Line {
                            start_mm: [-1.0, 0.0],
                            end_mm: [1.0, 0.0],
                        },
                        Curve2d::CircleArc {
                            center_mm: [0.0, 0.0],
                            radius_mm: 1.0,
                            start_angle_rad: 0.0,
                            end_angle_rad: std::f64::consts::PI,
                        },
                    ],
                }],
            };
            let root = model.add_node(CadNode::Extrude {
                profile,
                vector_mm: [0.0, 0.0, 2.0],
            });
            model.add_root(root)?;

            let kernel = MonstertruckKernel::project_standard();
            let evaluated = kernel.evaluate(&model, root)?;
            let summary = kernel.summarize(&evaluated);
            assert!(summary.geometrically_consistent);
            assert_eq!(summary.shells, 1);
            let step = kernel.to_step(&evaluated)?;
            step_p21::parser::parse(&step)?;
            Ok(())
        }

        #[test]
        fn evaluates_full_circle_profile_extrusion() -> Result<()> {
            let mut model = CadModel::new();
            let profile = Profile2d {
                loops: vec![ProfileLoop {
                    curves: vec![Curve2d::CircleArc {
                        center_mm: [0.0, 0.0],
                        radius_mm: 2.0,
                        start_angle_rad: 0.0,
                        end_angle_rad: std::f64::consts::TAU,
                    }],
                }],
            };
            let root = model.add_node(CadNode::Extrude {
                profile,
                vector_mm: [0.0, 0.0, 5.0],
            });
            model.add_root(root)?;

            let kernel = MonstertruckKernel::project_standard();
            let evaluated = kernel.evaluate(&model, root)?;
            let summary = kernel.summarize(&evaluated);
            assert!(summary.geometrically_consistent);
            assert_eq!(summary.shells, 1);
            let step = kernel.to_step(&evaluated)?;
            assert!(step.starts_with("ISO-10303-21;"));
            step_p21::parser::parse(&step)?;
            Ok(())
        }

        #[test]
        fn refuses_unimplemented_pattern_execution() -> Result<()> {
            use crate::cad_ir::PatternSpec;

            let (mut model, body) = box_model()?;
            model.roots.clear();
            let patterned = model.add_node(CadNode::Pattern {
                pattern: PatternSpec::Linear {
                    count: 4,
                    step_mm: [2.54, 0.0, 0.0],
                },
                child: body,
            });
            model.add_root(patterned)?;

            let error = match MonstertruckKernel::project_standard().evaluate(&model, patterned) {
                Ok(_) => bail!("unimplemented pattern unexpectedly evaluated"),
                Err(error) => error,
            };
            assert!(error.to_string().contains("does not yet evaluate"));
            Ok(())
        }

        /// The independent OCC constructive decomposition of a curved lip
        /// combines two exact circular boundary arcs with planar tapered
        /// side limits. Native intersection must retain both types of carrier.
        ///
        /// The residual treatment is not an opaque B-rep: both rounded
        /// terminals are executable analytic arc extrusions and intersections.
        #[test]
        fn evaluates_source_rounding_caps_as_editable_cad() -> Result<()> {
            let (model, upper, lower) = crate::cad_recovery::mixed_fixture_models::rounding_caps()?;
            let kernel = MonstertruckKernel::project_standard();
            model.validate()?;
            for (name, root) in [("upper", upper), ("lower", lower)] {
                let evaluated = kernel.evaluate(&model, root)?;
                assert!(kernel.summarize(&evaluated).is_single_solid());
                if let Ok(dir) = std::env::var("MIXED_ROUND_CAP_OUTPUT_DIR") {
                    std::fs::write(
                        std::path::Path::new(&dir).join(format!("native-rounding-{name}.step")),
                        kernel.to_step(&evaluated)?,
                    )?;
                }
            }
            Ok(())
        }

        #[test]
        fn evaluates_curved_lip_intersected_with_tapered_envelope() -> Result<()> {
            let (model, root) = crate::cad_recovery::mixed_fixture_models::curved_lip()?;
            let kernel = MonstertruckKernel::project_standard();
            let evaluated = kernel.evaluate(&model, root)?;
            assert!(kernel.summarize(&evaluated).is_single_solid());
            if let Ok(path) = std::env::var("CURVED_LIP_STEP_OUTPUT") {
                std::fs::write(path, kernel.to_step(&evaluated)?)?;
            }
            Ok(())
        }

        #[test]
        fn evaluates_box_without_exposing_monstertruck_topology() -> Result<()> {
            let (model, root) = box_model()?;
            let kernel = MonstertruckKernel::project_standard();
            let evaluated = kernel.evaluate(&model, root)?;
            let summary = kernel.summarize(&evaluated);
            assert!(summary.geometrically_consistent);
            assert_eq!(summary.shells, 1);
            assert_eq!(summary.faces, 6);

            let step = kernel.to_step(&evaluated)?;
            assert!(step.starts_with("ISO-10303-21;"));
            step_p21::parser::parse(&step)?;
            Ok(())
        }
    }
}

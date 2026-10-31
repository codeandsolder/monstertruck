//! Recovery of malformed STEP `EDGE_CURVE`s whose `edge_geometry` is `$`.
//!
//! The topology and face supports are still useful evidence when an exporter
//! omits the curve itself.  Recovery here is deliberately fail-closed:
//!
//! * a face contributes path/winding evidence only when exactly one edge is unresolved,
//! * every other edge in such a complementary chain must have geometry available,
//!   either from the source or from an earlier recovery pass,
//! * a face with multiple unresolved edges contributes only its native support and
//!   endpoint incidence -- never an invented path,
//! * complementary chains are sampled from their actual curves and projected
//!   continuously into the incident face's parameter domain,
//! * periodic coordinates must occupy less than one complete period, so the
//!   missing edge has an unambiguous local zero-wrap closure,
//! * at least two incident face supports must independently describe the edge,
//! * at least one support must carry complementary-path evidence,
//! * the transverse surface/surface intersection must continue from both edge
//!   endpoints through the same local branch and reach the opposite endpoint,
//! * the final native-support intersection curve is independently re-solved at
//!   inter-sample points rather than trusting its numerical guide.
//!
//! The recovered representation is an [`IntersectionCurve`] carrying the two
//! native support surfaces and a polyline leader.  It is therefore not a
//! tessellated replacement for the edge: the leader is only the numerical guide
//! used by Monstertruck's existing intersection-curve implementation.

use super::{
    BoundResolution, Curve3D, EdgeCurve, EuclideanSpace, FaceBoundHolder, HashMap, InnerSpace,
    IntersectionCurve, MetricSpace, Name, ParametricSurface3D, PlaceHolder, Point2, Point3,
    PolylineCurve, SearchNearestParameter, SearchParameter, ShellHolder, Surface, TOLERANCE, Table,
    step_geometry, step_p21,
};
use monstertruck_geometry::prelude::{
    BoundedCurve, Invertible, ParameterDivision1D, ParametricCurve, ParametricSurface,
};
use step_p21::tables::IntoOwned;

const PROJECTION_TRIALS: usize = 100;
const RECOVERY_SAMPLES: usize = 17;
const GEOMETRY_TOLERANCE: f64 = 100.0 * TOLERANCE;

#[inline]
fn trace_recovery(edge_id: u64, stage: &str) {
    if std::env::var_os("MT_STEP_MISSING_EDGE_TRACE").is_some() {
        eprintln!("trace step missing_edge edge=#{edge_id} {stage}");
    }
}

#[derive(Clone, Debug)]
struct RecoveryUse {
    edge_id: u64,
    orientation: bool,
    edge: EdgeCurve,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum SurfaceHintEvidence {
    /// The known complementary face boundary proves the missing edge's local UV path.
    ComplementPath,
    /// The face proves only that this edge is incident to the native support and
    /// that its endpoints lie on that support. It contributes no path/winding evidence.
    SupportOnly,
}

#[derive(Clone, Debug)]
struct SurfaceArcHint {
    surface: Surface,
    evidence: SurfaceHintEvidence,
    start_uv: Point2,
    end_uv: Point2,
    period_u: Option<f64>,
    period_v: Option<f64>,
}

type RecoveryIntersection = IntersectionCurve<Box<Curve3D>, Box<Surface>, Box<Surface>>;

#[derive(Clone, Copy, Debug)]
struct MarchState {
    point: Point3,
    uv0: Point2,
    uv1: Point2,
    guide: Point3,
}

#[inline]
fn unwrap_near(value: f64, reference: f64, period: Option<f64>) -> f64 {
    let Some(period) = period.filter(|period| period.is_finite() && *period > TOLERANCE) else {
        return value;
    };
    ((reference - value) / period)
        .round()
        .mul_add(period, value)
}

#[inline]
fn unwrap_uv_near(surface: &Surface, mut value: Point2, reference: Point2) -> Point2 {
    value.x = unwrap_near(value.x, reference.x, surface.period_u());
    value.y = unwrap_near(value.y, reference.y, surface.period_v());
    value
}

#[inline]
fn point_close(left: Point3, right: Point3) -> bool {
    left.distance2(right) <= GEOMETRY_TOLERANCE * GEOMETRY_TOLERANCE
}

fn project_on_surface(
    surface: &Surface,
    point: Point3,
    previous: Option<Point2>,
) -> Option<Point2> {
    let hint = previous.map(|uv| (uv.x, uv.y));
    let parameter = surface
        .search_parameter(point, hint, PROJECTION_TRIALS)
        .or_else(|| surface.search_nearest_parameter(point, hint, PROJECTION_TRIALS))?;
    let mut uv = Point2::new(parameter.0, parameter.1);
    if let Some(previous) = previous {
        uv = unwrap_uv_near(surface, uv, previous);
    }
    point_close(surface.evaluate(uv.x, uv.y), point).then_some(uv)
}

fn sample_oriented_curve(mut curve: Curve3D, orientation: bool) -> Option<Vec<Point3>> {
    if !orientation {
        curve.invert();
    }
    let range = curve.range_tuple();
    if !range.0.is_finite() || !range.1.is_finite() {
        return None;
    }

    // Follow the curve's actual 3D shape rather than sampling a magic number of
    // evenly-spaced parameters.  This matters on periodic surfaces: sparse
    // samples can alias a multi-turn edge back into a falsely local UV walk.
    // `try_parameter_division` is already the kernel's bounded-error curve
    // sampler and refuses curves it cannot subdivide reliably.
    let (_, points) = curve.try_parameter_division(range, GEOMETRY_TOLERANCE * 0.25)?;
    (points.len() >= 2).then_some(points)
}

#[inline]
fn span_is_local(
    values: &[Point2],
    period: Option<f64>,
    coordinate: impl Fn(Point2) -> f64,
) -> bool {
    let Some(period) = period.filter(|period| period.is_finite() && *period > TOLERANCE) else {
        return true;
    };
    let (minimum, maximum) = values.iter().copied().fold(
        (f64::INFINITY, f64::NEG_INFINITY),
        |(minimum, maximum), uv| {
            let value = coordinate(uv);
            (minimum.min(value), maximum.max(value))
        },
    );
    let numerical_margin = 1024.0 * f64::EPSILON * period.abs().max(1.0);
    maximum - minimum < period - numerical_margin
}

impl Table {
    fn ordered_recovery_uses(
        &self,
        bound: &FaceBoundHolder,
        face_sense: bool,
    ) -> Option<Vec<RecoveryUse>> {
        use PlaceHolder::Ref;

        let absolute_orientation = Self::absolute_bound_orientation(bound, face_sense);
        let edge_loop = match bound.resolve_bound(self) {
            BoundResolution::EdgeLoop(edge_loop) => edge_loop,
            BoundResolution::VertexLoop(_) | BoundResolution::Unresolved(_) => return None,
        };

        let mut uses = edge_loop
            .edge_list
            .into_iter()
            .map(|edge| {
                let Ref(Name::Entity(reference_id)) = edge else {
                    return None;
                };
                let (edge_id, edge_holder) = self.place_holder_edge_any_to_index_and_edge_curve(
                    &Ref(Name::Entity(reference_id)),
                )?;
                let orientation = self
                    .oriented_edge
                    .get(&reference_id)
                    .map_or(absolute_orientation, |oriented_edge| {
                        oriented_edge.orientation == absolute_orientation
                    });
                Some(RecoveryUse {
                    edge_id,
                    orientation,
                    edge: edge_holder.into_owned(self).ok()?,
                })
            })
            .collect::<Option<Vec<_>>>()?;

        if !absolute_orientation {
            uses.reverse();
        }
        Some(uses)
    }

    fn arc_hint_from_complement(
        surface: &Surface,
        recovery_uses: &[RecoveryUse],
        missing_index: usize,
        recovered: &HashMap<u64, Curve3D>,
    ) -> Option<SurfaceArcHint> {
        if recovery_uses.len() < 2 {
            return None;
        }

        let missing = &recovery_uses[missing_index];
        let edge_id = missing.edge_id;
        let intrinsic_start = Point3::from(&missing.edge.edge_start.vertex_geometry);
        let intrinsic_end = Point3::from(&missing.edge.edge_end.vertex_geometry);
        let (face_start, face_end) = if missing.orientation {
            (intrinsic_start, intrinsic_end)
        } else {
            (intrinsic_end, intrinsic_start)
        };

        // Walk the known complementary chain from the missing edge's end back
        // to its start.  Sampling the real curves, rather than just vertices,
        // preserves long periodic windings carried by valid neighboring edges.
        let mut projected = Vec::<Point2>::new();
        let mut previous = None;
        for offset in 1..recovery_uses.len() {
            let recovery_use = &recovery_uses[(missing_index + offset) % recovery_uses.len()];
            let curve = if recovery_use.edge.edge_geometry.is_none() {
                let Some(curve) = recovered.get(&recovery_use.edge_id).cloned() else {
                    trace_recovery(edge_id, "hint_fail=second_missing_edge");
                    return None;
                };
                curve
            } else {
                let Ok(curve) = recovery_use.edge.parse_curve3d() else {
                    trace_recovery(edge_id, "hint_fail=known_edge_parse");
                    return None;
                };
                curve
            };
            let Some(points) = sample_oriented_curve(curve, recovery_use.orientation) else {
                trace_recovery(edge_id, "hint_fail=known_edge_division");
                return None;
            };
            for (sample_index, point) in points.into_iter().enumerate() {
                if !projected.is_empty() && sample_index == 0 {
                    // Adjacent STEP edges share this endpoint; projecting it twice
                    // only adds numerical noise to the periodic unwrapping.
                    continue;
                }
                let Some(uv) = project_on_surface(surface, point, previous) else {
                    trace_recovery(edge_id, "hint_fail=known_point_projection");
                    return None;
                };
                previous = Some(uv);
                projected.push(uv);
            }
        }

        let (Some(&rest_start), Some(&rest_end)) = (projected.first(), projected.last()) else {
            trace_recovery(edge_id, "hint_fail=empty_complement");
            return None;
        };
        if !point_close(surface.evaluate(rest_start.x, rest_start.y), face_end)
            || !point_close(surface.evaluate(rest_end.x, rest_end.y), face_start)
        {
            trace_recovery(edge_id, "hint_fail=complement_endpoint_mismatch");
            return None;
        }

        let period_u = surface.period_u();
        let period_v = surface.period_v();
        if !span_is_local(&projected, period_u, |uv| uv.x)
            || !span_is_local(&projected, period_v, |uv| uv.y)
        {
            // A complementary chain that itself occupies a full period does not
            // prove which winding the omitted edge used.  Do not guess.
            trace_recovery(edge_id, "hint_fail=full_period_complement");
            return None;
        }

        // Closing in the same unwrapped chart gives the only local zero-wrap
        // edge compatible with the known complementary chain.
        let (use_start_uv, use_end_uv) = (rest_end, rest_start);
        let (start_uv, end_uv) = if missing.orientation {
            (use_start_uv, use_end_uv)
        } else {
            (use_end_uv, use_start_uv)
        };

        Some(SurfaceArcHint {
            surface: surface.clone(),
            evidence: SurfaceHintEvidence::ComplementPath,
            start_uv,
            end_uv,
            period_u,
            period_v,
        })
    }

    fn support_only_hint(surface: &Surface, missing: &RecoveryUse) -> Option<SurfaceArcHint> {
        let front = Point3::from(&missing.edge.edge_start.vertex_geometry);
        let back = Point3::from(&missing.edge.edge_end.vertex_geometry);
        let start_uv = project_on_surface(surface, front, None)?;
        let mut end_uv = project_on_surface(surface, back, Some(start_uv))?;
        end_uv = unwrap_uv_near(surface, end_uv, start_uv);
        Some(SurfaceArcHint {
            surface: surface.clone(),
            evidence: SurfaceHintEvidence::SupportOnly,
            start_uv,
            end_uv,
            period_u: surface.period_u(),
            period_v: surface.period_v(),
        })
    }

    fn missing_edge_hints(
        &self,
        shell: &ShellHolder,
        recovered: &HashMap<u64, Curve3D>,
    ) -> HashMap<u64, Vec<SurfaceArcHint>> {
        let mut hints = HashMap::<u64, Vec<SurfaceArcHint>>::new();

        for face in shell
            .cfs_faces_holder(self)
            .filter_map(|face| self.face_any_to_orientation_and_face(face))
        {
            let (orientation, face) = face;
            let Ok(step_surface) = face.face_geometry.clone().into_owned(self) else {
                continue;
            };
            let Ok(surface) = Surface::try_from(&step_surface) else {
                continue;
            };
            let face_sense = orientation == face.same_sense;

            for (_, bound) in face.bounds_holder(self) {
                let Some(bound) = bound else {
                    continue;
                };
                let Some(uses) = self.ordered_recovery_uses(&bound, face_sense) else {
                    continue;
                };
                let is_unresolved = |use_: &&RecoveryUse| {
                    use_.edge.edge_geometry.is_none() && !recovered.contains_key(&use_.edge_id)
                };
                let missing_count = uses.iter().filter(is_unresolved).count();
                if missing_count == 0 {
                    continue;
                }

                for (missing_index, missing_use) in uses
                    .iter()
                    .enumerate()
                    .filter(|(_, use_)| is_unresolved(use_))
                {
                    // A loop with exactly one absent curve proves the missing
                    // UV path from its complete complementary chain.  With two
                    // or more absent curves that path is underconstrained, but
                    // the ADVANCED_FACE still independently proves the native
                    // support and its incident edge endpoints.  Keep that weaker
                    // evidence rather than throwing the whole face away.
                    let path_hint = (missing_count == 1)
                        .then(|| {
                            Self::arc_hint_from_complement(
                                &surface,
                                &uses,
                                missing_index,
                                recovered,
                            )
                        })
                        .flatten();
                    let hint = path_hint.or_else(|| Self::support_only_hint(&surface, missing_use));
                    if let Some(hint) = hint {
                        trace_recovery(
                            missing_use.edge_id,
                            match hint.evidence {
                                SurfaceHintEvidence::ComplementPath => "hint_ok=path",
                                SurfaceHintEvidence::SupportOnly => "hint_ok=support",
                            },
                        );
                        hints.entry(missing_use.edge_id).or_default().push(hint);
                    } else {
                        trace_recovery(missing_use.edge_id, "hint_fail");
                    }
                }
            }
        }

        hints
    }

    fn support_pair_sample(
        edge_id: u64,
        first: &SurfaceArcHint,
        second: &SurfaceArcHint,
        intersection: &RecoveryIntersection,
        fraction: f64,
        first_sample: bool,
        previous: MarchState,
    ) -> Option<MarchState> {
        let expected0 = first.start_uv + (first.end_uv - first.start_uv) * fraction;
        let expected1 = second.start_uv + (second.end_uv - second.start_uv) * fraction;
        let point0 = first.surface.evaluate(expected0.x, expected0.y);
        let point1 = second.surface.evaluate(expected1.x, expected1.y);
        let guide_point = point0.midpoint(point1);

        let (point, raw_uv0, raw_uv1) = if first_sample {
            (previous.point, previous.uv0, previous.uv1)
        } else {
            // Both march directions solve the same transverse plane equation;
            // the previous UVs only seed Newton and select the local root.
            let predictor_normal0 = first.surface.normal(expected0.x, expected0.y);
            let predictor_normal1 = second.surface.normal(expected1.x, expected1.y);
            let tangent = predictor_normal0.cross(predictor_normal1);
            let normal_scale = predictor_normal0.magnitude() * predictor_normal1.magnitude();
            if !normal_scale.is_finite()
                || normal_scale <= TOLERANCE
                || tangent.magnitude() <= 1.0e-8 * normal_scale
            {
                trace_recovery(edge_id, "pair_fail=predictor_tangent");
                return None;
            }
            let Some(solution) = intersection.search_intersection_in_plane(
                guide_point,
                tangent,
                Some((previous.uv0.x, previous.uv0.y)),
                Some((previous.uv1.x, previous.uv1.y)),
                PROJECTION_TRIALS,
            ) else {
                trace_recovery(edge_id, "pair_fail=plane_ssi_search");
                return None;
            };
            solution
        };

        let uv0 = unwrap_uv_near(&first.surface, raw_uv0, previous.uv0);
        let uv1 = unwrap_uv_near(&second.surface, raw_uv1, previous.uv1);
        if !point_close(first.surface.evaluate(uv0.x, uv0.y), point)
            || !point_close(second.surface.evaluate(uv1.x, uv1.y), point)
        {
            trace_recovery(edge_id, "pair_fail=support_residual");
            return None;
        }

        let normal0 = first.surface.normal(uv0.x, uv0.y);
        let normal1 = second.surface.normal(uv1.x, uv1.y);
        let normal_scale = normal0.magnitude() * normal1.magnitude();
        if !normal_scale.is_finite()
            || normal_scale <= TOLERANCE
            || normal0.cross(normal1).magnitude() <= 1.0e-8 * normal_scale
        {
            trace_recovery(edge_id, "pair_fail=tangent");
            return None;
        }

        if !first_sample {
            let guide_step = previous.guide.distance(guide_point);
            let actual_step = previous.point.distance(point);
            if actual_step > 8.0 * guide_step.max(GEOMETRY_TOLERANCE) {
                trace_recovery(edge_id, "pair_fail=continuation_jump");
                return None;
            }
        }

        for (value, expected, period, evidence) in [
            (uv0.x, expected0.x, first.period_u, first.evidence),
            (uv0.y, expected0.y, first.period_v, first.evidence),
            (uv1.x, expected1.x, second.period_u, second.evidence),
            (uv1.y, expected1.y, second.period_v, second.evidence),
        ] {
            if evidence == SurfaceHintEvidence::ComplementPath
                && let Some(period) =
                    period.filter(|period| period.is_finite() && *period > TOLERANCE)
                && (unwrap_near(value, expected, Some(period)) - expected).abs() >= period * 0.5
            {
                trace_recovery(edge_id, "pair_fail=periodic_branch_drift");
                return None;
            }
        }

        Some(MarchState {
            point,
            uv0,
            uv1,
            guide: guide_point,
        })
    }

    fn march_support_pair(
        edge_id: u64,
        front: Point3,
        back: Point3,
        first: &SurfaceArcHint,
        second: &SurfaceArcHint,
        intersection: &RecoveryIntersection,
        forward: bool,
    ) -> Option<Vec<Point3>> {
        let mut points = Vec::with_capacity(RECOVERY_SAMPLES);
        let mut state = MarchState {
            point: if forward { front } else { back },
            uv0: if forward {
                first.start_uv
            } else {
                first.end_uv
            },
            uv1: if forward {
                second.start_uv
            } else {
                second.end_uv
            },
            guide: if forward { front } else { back },
        };
        let intervals = u32::try_from(RECOVERY_SAMPLES.checked_sub(1)?).ok()?;

        for sample_index in 0..RECOVERY_SAMPLES {
            let sample_index = u32::try_from(sample_index).ok()?;
            let unit = f64::from(sample_index) / f64::from(intervals);
            let fraction = if forward { unit } else { 1.0 - unit };
            state = Self::support_pair_sample(
                edge_id,
                first,
                second,
                intersection,
                fraction,
                sample_index == 0,
                state,
            )?;
            points.push(state.point);
        }

        let target_point = if forward { back } else { front };
        let target_uv0 = if forward {
            first.end_uv
        } else {
            first.start_uv
        };
        let target_uv1 = if forward {
            second.end_uv
        } else {
            second.start_uv
        };
        if !point_close(state.point, target_point) {
            trace_recovery(edge_id, "pair_fail=endpoint_point");
            return None;
        }
        for (value, target, period) in [
            (state.uv0.x, target_uv0.x, first.period_u),
            (state.uv0.y, target_uv0.y, first.period_v),
            (state.uv1.x, target_uv1.x, second.period_u),
            (state.uv1.y, target_uv1.y, second.period_v),
        ] {
            let value = unwrap_near(value, target, period);
            let tolerance = period
                .map_or(GEOMETRY_TOLERANCE, |period| period.abs() * 1.0e-6)
                .max(GEOMETRY_TOLERANCE);
            if (value - target).abs() > tolerance {
                trace_recovery(edge_id, "pair_fail=endpoint_uv");
                return None;
            }
        }
        Some(points)
    }

    fn validate_recovery_candidate(edge_id: u64, candidate: &Curve3D) -> Option<()> {
        let Curve3D::IntersectionCurve(intersection) = candidate else {
            return None;
        };
        let (t0, t1) = intersection.range_tuple();
        let intervals = u32::try_from(RECOVERY_SAMPLES.checked_sub(1)?).ok()?;
        for index in 0..RECOVERY_SAMPLES - 1 {
            let index = u32::try_from(index).ok()?;
            let fraction = (f64::from(index) + 0.5) / f64::from(intervals);
            let t = f64::mul_add(t1 - t0, fraction, t0);
            let Some((point, uv0, uv1)) = intersection.search_triple(t, PROJECTION_TRIALS) else {
                trace_recovery(edge_id, "pair_fail=final_search_triple");
                return None;
            };
            if !point_close(intersection.surface0().evaluate(uv0.x, uv0.y), point)
                || !point_close(intersection.surface1().evaluate(uv1.x, uv1.y), point)
            {
                trace_recovery(edge_id, "pair_fail=final_support_residual");
                return None;
            }
        }
        Some(())
    }

    fn recover_from_support_pair(
        edge_id: u64,
        front: Point3,
        back: Point3,
        first: &SurfaceArcHint,
        second: &SurfaceArcHint,
    ) -> Option<Curve3D> {
        if first.evidence != SurfaceHintEvidence::ComplementPath
            && second.evidence != SurfaceHintEvidence::ComplementPath
        {
            trace_recovery(edge_id, "pair_fail=insufficient_path_evidence");
            return None;
        }
        if step_geometry::same_surface_domain(
            &first.surface,
            &second.surface,
            GEOMETRY_TOLERANCE,
            1.0e-8,
        )
        .is_some()
        {
            trace_recovery(edge_id, "pair_fail=coincident_supports");
            return None;
        }
        for hint in [first, second] {
            if !point_close(
                hint.surface.evaluate(hint.start_uv.x, hint.start_uv.y),
                front,
            ) || !point_close(hint.surface.evaluate(hint.end_uv.x, hint.end_uv.y), back)
            {
                trace_recovery(edge_id, "pair_fail=hint_endpoint_orientation");
                return None;
            }
        }

        let intersection = RecoveryIntersection::new(
            Box::new(first.surface.clone()),
            Box::new(second.surface.clone()),
            Box::new(Curve3D::Polyline(PolylineCurve(vec![front, back]))),
        );
        let Some(mut forward) =
            Self::march_support_pair(edge_id, front, back, first, second, &intersection, true)
        else {
            trace_recovery(edge_id, "pair_fail=forward_march");
            return None;
        };
        let Some(mut backward) =
            Self::march_support_pair(edge_id, front, back, first, second, &intersection, false)
        else {
            trace_recovery(edge_id, "pair_fail=backward_march");
            return None;
        };
        backward.reverse();
        if forward
            .iter()
            .zip(&backward)
            .any(|(left, right)| !point_close(*left, *right))
        {
            trace_recovery(edge_id, "pair_fail=bidirectional_disagreement");
            return None;
        }

        forward[0] = front;
        *forward.last_mut()? = back;
        let candidate = Curve3D::IntersectionCurve(RecoveryIntersection::new(
            Box::new(first.surface.clone()),
            Box::new(second.surface.clone()),
            Box::new(Curve3D::Polyline(PolylineCurve(forward))),
        ));
        Self::validate_recovery_candidate(edge_id, &candidate)?;
        trace_recovery(edge_id, "pair_proof_ok");
        Some(candidate)
    }

    fn shell_has_missing_edge_geometry(&self, shell: &ShellHolder) -> bool {
        shell
            .cfs_faces_holder(self)
            .filter_map(|face| self.face_any_to_orientation_and_face(face))
            .flat_map(|(_, face)| face.bounds_holder(self))
            .filter_map(|(_, bound)| bound)
            .filter_map(|bound| self.ordered_recovery_uses(&bound, true))
            .flatten()
            .any(|use_| use_.edge.edge_geometry.is_none())
    }

    pub(super) fn recover_missing_edge_curves(&self, shell: &ShellHolder) -> HashMap<u64, Curve3D> {
        if !self.shell_has_missing_edge_geometry(shell) {
            return HashMap::new();
        }

        let mut recovered = HashMap::<u64, Curve3D>::new();

        loop {
            let hints = self.missing_edge_hints(shell, &recovered);
            let mut newly_recovered = Vec::<(u64, Curve3D)>::new();

            for (edge_id, hints) in hints {
                if std::env::var_os("MT_STEP_MISSING_EDGE_TRACE").is_some() {
                    eprintln!(
                        "trace step missing_edge edge=#{edge_id} hints={} pass_recovered={}",
                        hints.len(),
                        recovered.len(),
                    );
                }
                let Some(holder) = self.edge_curve.get(&edge_id) else {
                    continue;
                };
                let Ok(edge) = holder.clone().into_owned(self) else {
                    continue;
                };
                if edge.edge_geometry.is_some() || recovered.contains_key(&edge_id) {
                    continue;
                }
                let front = Point3::from(&edge.edge_start.vertex_geometry);
                let back = Point3::from(&edge.edge_end.vertex_geometry);

                let mut accepted = None::<Curve3D>;
                let mut conflicting = false;
                'pairs: for first in 0..hints.len() {
                    for second in first + 1..hints.len() {
                        let Some(candidate) = Self::recover_from_support_pair(
                            edge_id,
                            front,
                            back,
                            &hints[first],
                            &hints[second],
                        ) else {
                            trace_recovery(edge_id, "pair_fail");
                            // A support pair can legitimately be unhelpful, for
                            // example when it is coincident, tangent, or carries
                            // only weak support evidence. That is absence of
                            // evidence, not contradiction.
                            continue;
                        };
                        if let Some(existing) = &accepted {
                            let (e0, e1) = existing.range_tuple();
                            let (c0, c1) = candidate.range_tuple();
                            let agrees = [0.0, 0.25, 0.5, 0.75, 1.0].into_iter().all(|fraction| {
                                point_close(
                                    existing.evaluate(f64::mul_add(e1 - e0, fraction, e0)),
                                    candidate.evaluate(f64::mul_add(c1 - c0, fraction, c0)),
                                )
                            });
                            if !agrees {
                                conflicting = true;
                                break 'pairs;
                            }
                        } else {
                            accepted = Some(candidate);
                        }
                        trace_recovery(edge_id, "pair_ok");
                    }
                }

                if !conflicting && let Some(candidate) = accepted {
                    trace_recovery(edge_id, "recovered");
                    newly_recovered.push((edge_id, candidate));
                } else {
                    trace_recovery(edge_id, "not_recovered");
                }
            }

            if newly_recovered.is_empty() {
                break;
            }
            for (edge_id, curve) in newly_recovered {
                recovered.insert(edge_id, curve);
            }
        }

        recovered
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::step::load::step_geometry::ElementarySurface;
    use monstertruck_geometry::prelude::{Line, Plane, Processor, RevolutionSurface, Vector3};

    fn cylinder() -> Surface {
        let axis = Vector3::unit_z();
        let origin = Point3::origin();
        let point = Point3::new(1.0, 0.0, 0.0);
        Surface::Elementary(ElementarySurface::CylindricalSurface(Processor::new(
            RevolutionSurface::by_revolution(Line(point, point + axis), origin, axis),
        )))
    }

    fn xy_plane() -> Surface { Surface::Elementary(ElementarySurface::Plane(Plane::xy())) }

    #[test]
    fn coincident_supports_do_not_define_a_recoverable_edge() {
        let surface = cylinder();
        let front = Point3::new(1.0, 0.0, 0.0);
        let back = Point3::new(0.0, 1.0, 0.0);
        let start = surface
            .search_parameter(front, None, PROJECTION_TRIALS)
            .expect("front on cylinder");
        let end = surface
            .search_parameter(back, Some(start), PROJECTION_TRIALS)
            .expect("back on cylinder");
        let hint = SurfaceArcHint {
            surface: surface.clone(),
            evidence: SurfaceHintEvidence::ComplementPath,
            start_uv: start.into(),
            end_uv: end.into(),
            period_u: surface.period_u(),
            period_v: surface.period_v(),
        };
        assert!(Table::recover_from_support_pair(0, front, back, &hint, &hint).is_none());
    }

    #[test]
    fn a_full_period_complement_is_not_local_evidence() {
        let values = [
            Point2::new(0.0, 0.0),
            Point2::new(std::f64::consts::PI, 0.0),
            Point2::new(std::f64::consts::TAU, 0.0),
        ];
        assert!(!span_is_local(&values, Some(std::f64::consts::TAU), |uv| {
            uv.x
        },));
    }

    #[test]
    fn support_pair_recovers_the_proven_quarter_circle_branch() {
        let cylinder = cylinder();
        let plane = xy_plane();
        let front = Point3::new(1.0, 0.0, 0.0);
        let back = Point3::new(0.0, 1.0, 0.0);

        let c0 = cylinder
            .search_parameter(front, None, PROJECTION_TRIALS)
            .expect("front on cylinder");
        let c1 = cylinder
            .search_parameter(back, Some(c0), PROJECTION_TRIALS)
            .expect("back on cylinder");
        let p0 = plane
            .search_parameter(front, None, PROJECTION_TRIALS)
            .expect("front on plane");
        let p1 = plane
            .search_parameter(back, Some(p0), PROJECTION_TRIALS)
            .expect("back on plane");

        let cylinder_hint = SurfaceArcHint {
            surface: cylinder.clone(),
            evidence: SurfaceHintEvidence::ComplementPath,
            start_uv: c0.into(),
            end_uv: c1.into(),
            period_u: cylinder.period_u(),
            period_v: cylinder.period_v(),
        };
        let plane_hint = SurfaceArcHint {
            surface: plane.clone(),
            evidence: SurfaceHintEvidence::SupportOnly,
            start_uv: p0.into(),
            end_uv: p1.into(),
            period_u: plane.period_u(),
            period_v: plane.period_v(),
        };

        let curve = Table::recover_from_support_pair(0, front, back, &cylinder_hint, &plane_hint)
            .expect("one proven path plus one native support should recover the edge");
        let (t0, t1) = curve.range_tuple();
        let midpoint = curve.evaluate((t0 + t1) * 0.5);
        assert!(midpoint.x > 0.5);
        assert!(midpoint.y > 0.5);
        assert!(midpoint.z.abs() <= GEOMETRY_TOLERANCE);
        assert!((midpoint.to_vec().magnitude() - 1.0).abs() <= GEOMETRY_TOLERANCE);
    }
    #[test]
    fn support_only_pair_does_not_invent_a_path() {
        let cylinder = cylinder();
        let plane = xy_plane();
        let front = Point3::new(1.0, 0.0, 0.0);
        let back = Point3::new(0.0, 1.0, 0.0);
        let c0 = cylinder
            .search_parameter(front, None, PROJECTION_TRIALS)
            .expect("front lies on support");
        let c1 = cylinder
            .search_parameter(back, Some(c0), PROJECTION_TRIALS)
            .expect("back lies on cylinder");
        let p0 = plane
            .search_parameter(front, None, PROJECTION_TRIALS)
            .expect("front lies on support");
        let p1 = plane
            .search_parameter(back, Some(p0), PROJECTION_TRIALS)
            .expect("back lies on plane");
        let cylinder_hint = SurfaceArcHint {
            surface: cylinder.clone(),
            evidence: SurfaceHintEvidence::SupportOnly,
            start_uv: c0.into(),
            end_uv: c1.into(),
            period_u: cylinder.period_u(),
            period_v: cylinder.period_v(),
        };
        let plane_hint = SurfaceArcHint {
            surface: plane.clone(),
            evidence: SurfaceHintEvidence::SupportOnly,
            start_uv: p0.into(),
            end_uv: p1.into(),
            period_u: plane.period_u(),
            period_v: plane.period_v(),
        };
        assert!(
            Table::recover_from_support_pair(0, front, back, &cylinder_hint, &plane_hint).is_none(),
            "native supports plus endpoints alone do not prove which intersection arc was intended",
        );
    }

    #[test]
    fn support_pair_preserves_short_branch_across_periodic_seam() {
        let cylinder = cylinder();
        let plane = xy_plane();
        let angle = 10.0_f64.to_radians();
        let front = Point3::new(angle.cos(), -angle.sin(), 0.0);
        let back = Point3::new(angle.cos(), angle.sin(), 0.0);

        let c0 = Point2::from(
            cylinder
                .search_parameter(front, None, PROJECTION_TRIALS)
                .expect("front on cylinder"),
        );
        let c1_raw = Point2::from(
            cylinder
                .search_parameter(back, Some((c0.x, c0.y)), PROJECTION_TRIALS)
                .expect("back on cylinder"),
        );
        let c1 = unwrap_uv_near(&cylinder, c1_raw, c0);
        let p0 = Point2::from(
            plane
                .search_parameter(front, None, PROJECTION_TRIALS)
                .expect("front on plane"),
        );
        let p1 = Point2::from(
            plane
                .search_parameter(back, Some((p0.x, p0.y)), PROJECTION_TRIALS)
                .expect("back on plane"),
        );

        let cylinder_hint = SurfaceArcHint {
            surface: cylinder.clone(),
            evidence: SurfaceHintEvidence::ComplementPath,
            start_uv: c0,
            end_uv: c1,
            period_u: cylinder.period_u(),
            period_v: cylinder.period_v(),
        };
        let plane_hint = SurfaceArcHint {
            surface: plane.clone(),
            evidence: SurfaceHintEvidence::SupportOnly,
            start_uv: p0,
            end_uv: p1,
            period_u: plane.period_u(),
            period_v: plane.period_v(),
        };

        let curve = Table::recover_from_support_pair(0, front, back, &cylinder_hint, &plane_hint)
            .expect("the short seam-crossing branch is uniquely proven");
        let (t0, t1) = curve.range_tuple();
        let midpoint = curve.evaluate((t0 + t1) * 0.5);
        assert!(
            midpoint.x > 0.99,
            "recovery took the long periodic branch: {midpoint:?}"
        );
        assert!(midpoint.y.abs() <= 1.0e-5);
        assert!(midpoint.z.abs() <= GEOMETRY_TOLERANCE);
    }

    #[test]
    fn support_pair_refuses_a_misoriented_path_hint() {
        let cylinder = cylinder();
        let plane = xy_plane();
        let front = Point3::new(1.0, 0.0, 0.0);
        let back = Point3::new(0.0, 1.0, 0.0);
        let c0 = Point2::from(
            cylinder
                .search_parameter(front, None, PROJECTION_TRIALS)
                .expect("front on cylinder"),
        );
        let c1 = Point2::from(
            cylinder
                .search_parameter(back, Some((c0.x, c0.y)), PROJECTION_TRIALS)
                .expect("back on cylinder"),
        );
        let p0 = Point2::from(
            plane
                .search_parameter(front, None, PROJECTION_TRIALS)
                .expect("front on plane"),
        );
        let p1 = Point2::from(
            plane
                .search_parameter(back, Some((p0.x, p0.y)), PROJECTION_TRIALS)
                .expect("back on plane"),
        );
        let reversed_path = SurfaceArcHint {
            surface: cylinder.clone(),
            evidence: SurfaceHintEvidence::ComplementPath,
            start_uv: c1,
            end_uv: c0,
            period_u: cylinder.period_u(),
            period_v: cylinder.period_v(),
        };
        let plane_hint = SurfaceArcHint {
            surface: plane.clone(),
            evidence: SurfaceHintEvidence::SupportOnly,
            start_uv: p0,
            end_uv: p1,
            period_u: plane.period_u(),
            period_v: plane.period_v(),
        };

        assert!(
            Table::recover_from_support_pair(0, front, back, &reversed_path, &plane_hint).is_none()
        );
    }
}

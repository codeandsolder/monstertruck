//! Boundary signature for regular intersections of two cylinders.
//!
//! Carrier identity, the transverse root sign, and the increasing angular
//! interval define the physical locus. The polynomial curve must be proven
//! close to that locus over its full domain, not just at endpoint samples.

use super::{
    Arc, BoundaryChainKey, BoundedEdge, DIRECTION_Q, EdgeGeometry, OrientedSupportKey,
    SupportLocusKey, TOL, decode, dequantize_scalar, dequantize_vector, normalize_angular_chains,
    point, support_locus_key, vector,
};
use crate::math3::{distance, dot, mul, normalize, sub};
use monstertruck_geometry::{
    analytic_cylinder_intersection::{
        BoundedSourceIntersection, CylinderCarrier, IntersectionBranch,
    },
    prelude::{NurbsCurve, Point3, Vector4},
};

/// Chart selection and root are part of the canonical chain identity.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub(in crate::boundary_proof) struct CylinderIntersectionKey {
    pub primary_is_second: bool,
    pub positive_branch: bool,
}

#[derive(Clone, Copy)]
struct CylinderChart {
    origin: [f64; 3],
    axis: [f64; 3],
    radius: f64,
    carrier: CylinderCarrier,
}

impl CylinderChart {
    fn new(locus: &SupportLocusKey) -> Option<Self> {
        let SupportLocusKey::Cylinder {
            axis,
            closest_mm,
            radius_mm,
        } = locus
        else {
            return None;
        };
        let origin = decode(*closest_mm)?;
        let axis = normalize(dequantize_vector(*axis, DIRECTION_Q)?, 0.0)?;
        let radius = dequantize_scalar(*radius_mm)?;
        if radius <= TOL {
            return None;
        }
        let carrier = CylinderCarrier::new(point(origin), vector(axis), radius)?;
        Some(Self {
            origin,
            axis,
            radius,
            carrier,
        })
    }

    fn radial(self, p: [f64; 3]) -> [f64; 3] {
        let delta = sub(p, self.origin);
        sub(delta, mul(self.axis, dot(delta, self.axis)))
    }
}

struct IntersectionChart {
    primary: CylinderChart,
    secondary: CylinderChart,
    primary_is_second: bool,
}

impl IntersectionChart {
    fn new(loci: &[SupportLocusKey; 2]) -> Option<Self> {
        let first = CylinderChart::new(&loci[0])?;
        let second = CylinderChart::new(&loci[1])?;
        // A deterministic primary chart is needed across edge subdivision.
        // Smaller radius usually gives a better separated transverse root.
        let primary_is_second = second.radius < first.radius;
        let (primary, secondary) = if primary_is_second {
            (second, first)
        } else {
            (first, second)
        };
        Some(Self {
            primary,
            secondary,
            primary_is_second,
        })
    }

    fn certify(&self, edge: &EdgeGeometry, curve: &NurbsCurve<Vector4>) -> Option<BoundedEdge> {
        let bound = BoundedSourceIntersection::from_spline(
            self.primary.carrier,
            self.secondary.carrier,
            curve,
            TOL,
            TOL / 64.0,
        )?;
        // The bidirectional separation includes the amplification of both
        // carrier deviations by the transverse intersection conditioning.
        // Independent incidence within TOL alone is not enough near tangency.
        if bound.separation_upper > TOL {
            return None;
        }
        let arc = bound.proposal.arc;
        let first = arc.evaluate(0.0)?;
        let last = arc.evaluate(arc.sweep())?;
        let start = decode(edge.start)?;
        let end = decode(edge.end)?;
        let near = |p: Point3, q: [f64; 3]| distance([p.x, p.y, p.z], q) <= TOL;
        let mut ends = if near(first, start) && near(last, end) {
            [edge.start, edge.end]
        } else if near(first, end) && near(last, start) {
            [edge.end, edge.start]
        } else {
            return None;
        };
        let tolerance = TOL / self.primary.radius;
        let sweep = arc.sweep().abs();
        if sweep <= tolerance
            || (ends[0] == ends[1]) != ((sweep - std::f64::consts::TAU).abs() <= tolerance)
            || (ends[0] == ends[1] && edge.start_vertex != edge.end_vertex)
        {
            return None;
        }
        if arc.sweep() < 0.0 {
            ends.reverse();
        }
        Some(BoundedEdge::CylinderIntersection(
            CylinderIntersectionKey {
                primary_is_second: self.primary_is_second,
                positive_branch: arc.branch() == IntersectionBranch::Positive,
            },
            Arc {
                start: ends[0],
                end: ends[1],
                sweep,
            },
        ))
    }
}

pub(super) fn certify_crossing_cylinders(
    edge: &EdgeGeometry,
    curve: &NurbsCurve<Vector4>,
    loci: &[SupportLocusKey; 2],
) -> Option<BoundedEdge> {
    IntersectionChart::new(loci)?.certify(edge, curve)
}

pub(in crate::boundary_proof) fn normalize_intersections(
    supports: &[OrientedSupportKey; 2],
    key: CylinderIntersectionKey,
    arcs: &[Arc],
) -> Result<Vec<BoundaryChainKey>, String> {
    let chart = IntersectionChart::new(&supports.each_ref().map(support_locus_key))
        .ok_or("invalid cylinder intersection chart")?;
    if chart.primary_is_second != key.primary_is_second {
        return Err("inconsistent primary cylinder chart".into());
    }
    normalize_angular_chains(arcs, TOL / chart.primary.radius, chart.primary.axis, |p| {
        chart.primary.radial(p)
    })
    .map(|chains| {
        chains
            .into_iter()
            .map(|endpoints_mm| BoundaryChainKey::CylinderIntersection {
                supports: supports.clone(),
                primary_is_second: key.primary_is_second,
                positive_branch: key.positive_branch,
                endpoints_mm,
            })
            .collect()
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::boundary_proof::point_key;
    use crate::brep::CurveSupport;
    use anyhow::{Context, Result, bail};
    use monstertruck_geometry::{
        analytic_cylinder_intersection::CylinderIntersectionArc,
        prelude::{InnerProduct, Vector3},
    };
    use std::f64::consts::PI;

    fn supports() -> [OrientedSupportKey; 2] {
        [
            OrientedSupportKey::Cylinder {
                axis: [10_000_000_000, 0, 0],
                closest_mm: [0; 3],
                radius_mm: 5_000_000,
                radial_outward: true,
            },
            OrientedSupportKey::Cylinder {
                axis: [0, 0, 10_000_000_000],
                closest_mm: [0; 3],
                radius_mm: 2_000_000,
                radial_outward: false,
            },
        ]
    }

    fn certify(
        arc: CylinderIntersectionArc,
        loci: &[SupportLocusKey; 2],
        id: u64,
    ) -> Result<(CylinderIntersectionKey, Arc)> {
        let leader = arc.bounded_leader(1e-9, 8192).context("bounded leader")?;
        let key =
            |point: Point3| point_key([point.x, point.y, point.z]).context("invalid point key");
        let e = EdgeGeometry {
            start: key(arc.evaluate(0.0).context("start")?)?,
            end: key(arc.evaluate(arc.sweep()).context("end")?)?,
            start_vertex: id,
            end_vertex: id + 1,
            parameter_forward: true,
            curve_id: id,
            support: CurveSupport::Other { entity_id: id },
        };
        let rational = NurbsCurve::from(leader.curve);
        let Some(BoundedEdge::CylinderIntersection(branch, arc)) =
            certify_crossing_cylinders(&e, &rational, loci)
        else {
            bail!("analytical crossing-cylinder locus was not certified");
        };
        Ok((branch, arc))
    }

    #[test]
    fn partition_and_reversal_preserve_crossing_cylinder_locus() -> Result<()> {
        let supports = supports();
        let loci = supports.each_ref().map(support_locus_key);
        let chart = IntersectionChart::new(&loci).context("chart")?;
        let new_arc = |start, sweep, branch| {
            CylinderIntersectionArc::new(
                chart.primary.carrier,
                chart.secondary.carrier,
                start,
                sweep,
                branch,
                1e-8,
            )
            .context("arc")
        };
        let origin = Point3::new(0.2, 0.0, 0.5);
        let whole = new_arc(origin, PI, IntersectionBranch::Positive)?;
        let mid = whole.evaluate(PI / 2.0).context("mid")?;
        let end = whole.evaluate(PI).context("end")?;
        let half_a = new_arc(origin, PI / 2.0, IntersectionBranch::Positive)?;
        let half_b = new_arc(mid, PI / 2.0, IntersectionBranch::Positive)?;
        let reversed = new_arc(end, -PI, IntersectionBranch::Positive)?;
        let complement = new_arc(origin, -PI, IntersectionBranch::Positive)?;
        let (key, arc) = certify(whole, &loci, 1)?;
        let (key_a, a) = certify(half_a, &loci, 2)?;
        let (key_b, b) = certify(half_b, &loci, 3)?;
        let (key_r, r) = certify(reversed, &loci, 4)?;
        let (key_c, c) = certify(complement, &loci, 5)?;
        assert_eq!(key, key_a);
        assert_eq!(key, key_b);
        assert_eq!(key, key_r);
        assert_eq!(key, key_c);
        let expected =
            normalize_intersections(&supports, key, &[arc]).map_err(anyhow::Error::msg)?;
        assert_eq!(
            normalize_intersections(&supports, key, &[b, a]).map_err(anyhow::Error::msg)?,
            expected
        );
        assert_eq!(
            normalize_intersections(&supports, key, &[r]).map_err(anyhow::Error::msg)?,
            expected
        );
        assert_ne!(
            normalize_intersections(&supports, key, &[c]).map_err(anyhow::Error::msg)?,
            expected
        );
        assert!(normalize_intersections(&supports, key, &[arc, a, b]).is_err());
        // Same angular domain, different transverse root.
        let negative = new_arc(
            Point3::new(0.2, 0.0, -0.5),
            PI,
            IntersectionBranch::Negative,
        )?;
        let (negative_key, negative_arc) = certify(negative, &loci, 6)?;
        assert_ne!(key, negative_key);
        assert_ne!(
            normalize_intersections(&supports, negative_key, &[negative_arc])
                .map_err(anyhow::Error::msg)?,
            expected
        );
        // Endpoints match, but an off-carrier interior must not qualify.
        let mut leader = whole.bounded_leader(1e-9, 8192).context("leader")?.curve;
        let mid_control = leader.control_points().len() / 2;
        *leader.control_point_mut(mid_control) += Vector3::z() * 0.001;
        let e = EdgeGeometry {
            start: arc.start,
            end: arc.end,
            start_vertex: 1,
            end_vertex: 2,
            parameter_forward: true,
            curve_id: 1,
            support: CurveSupport::Other { entity_id: 1 },
        };
        assert!(certify_crossing_cylinders(&e, &NurbsCurve::from(leader), &loci).is_none());
        assert!((end - origin).norm() > 0.3);
        Ok(())
    }
}

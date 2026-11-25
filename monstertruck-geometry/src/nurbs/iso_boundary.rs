//! Algebraic evidence for actual clamped isoparametric boundaries.
//!
//! Inverse surface search is singular at collapsed edges and poles. A
//! regular boundary curve still has a perfectly well-defined locus there.
//! Use the exact boundary control net, never a sampled inverse search,
//! to establish these particular curve-on-surface incidences.

use crate::{
    analytic_line_segment::{SegmentLocus, certify_nurbs_segment_locus},
    prelude::*,
    proof_bezier_spans,
};

struct IsoBoundary {
    curve: NurbsCurve<Vector4>,
    uv: Line<Point2>,
}

fn same_knot(a: f64, b: f64) -> bool {
    a.total_cmp(&b).is_eq() || (a.abs().to_bits() == 0 && b.abs().to_bits() == 0)
}

// Tolerance-based knot multiplicity is inappropriate for an exact
// boundary proof: a near-clamped control row need not be a surface edge.
fn clamped_ends(knots: &KnotVector, count: usize) -> Option<(f64, f64, bool, bool)> {
    let degree = knots.len().checked_sub(count + 1)?;
    if degree == 0
        || degree >= count
        || !knots.iter().all(|k| k.is_finite())
        || !knots.iter().zip(knots.iter().skip(1)).all(|(a, b)| a <= b)
    {
        return None;
    }
    let lo = knots[degree];
    let hi = knots[count];
    if lo >= hi {
        return None;
    }
    let front = knots.iter().take(degree + 1).all(|&k| same_knot(k, lo));
    let back = knots.iter().skip(count).all(|&k| same_knot(k, hi));
    Some((lo, hi, front, back))
}

// The same complete clamped control row/column and the same transverse
// basis prove pointwise equality along the entire seam. Exact equality is
// required; a near-coincident boundary is not a periodic chart.
pub(super) fn closed_chart_periods<P: PartialEq>(
    surface: &BsplineSurface<P>,
    valid: impl Fn(&P) -> bool,
) -> (Option<f64>, Option<f64>) {
    let net = surface.control_points();
    let Some(first) = net.first() else {
        return (None, None);
    };
    let rows = net.len();
    let columns = first.len();
    if rows < 2
        || columns < 2
        || net.iter().any(|row| row.len() != columns)
        || !net.iter().flatten().all(valid)
    {
        return (None, None);
    }
    let Some((u0, u1, uf, ub)) = clamped_ends(surface.knot_vector_u(), rows) else {
        return (None, None);
    };
    let Some((v0, v1, vf, vb)) = clamped_ends(surface.knot_vector_v(), columns) else {
        return (None, None);
    };
    let u = u1 - u0;
    let v = v1 - v0;
    (
        (uf && ub && u.is_finite() && net[0] == net[rows - 1]).then_some(u),
        (vf && vb && v.is_finite() && net.iter().all(|row| row[0] == row[columns - 1]))
            .then_some(v),
    )
}

fn boundaries(surface: &NurbsSurface<Vector4>) -> Option<Vec<IsoBoundary>> {
    let net = surface.control_points();
    let rows = net.len();
    let columns = net.first()?.len();
    if rows < 2
        || columns < 2
        || net.iter().any(|r| r.len() != columns)
        || net
            .iter()
            .flatten()
            .any(|p| p.w <= 0.0 || ![p.x, p.y, p.z, p.w].into_iter().all(f64::is_finite))
    {
        return None;
    }
    let (u0, u1, uf, ub) = clamped_ends(surface.knot_vector_u(), rows)?;
    let (v0, v1, vf, vb) = clamped_ends(surface.knot_vector_v(), columns)?;
    let mut result = Vec::with_capacity(4);
    for (enabled, index, u) in [(uf, 0, u0), (ub, rows - 1, u1)] {
        if enabled {
            result.push(IsoBoundary {
                curve: surface.curve_v(index),
                uv: Line(Point2::new(u, v0), Point2::new(u, v1)),
            });
        }
    }
    for (enabled, index, v) in [(vf, 0, v0), (vb, columns - 1, v1)] {
        if enabled {
            result.push(IsoBoundary {
                curve: surface.curve_u(index),
                uv: Line(Point2::new(u0, v), Point2::new(u1, v)),
            });
        }
    }
    Some(result)
}

fn segment_fraction(point: Point3, start: Point3, end: Point3, tolerance: f64) -> Option<f64> {
    let edge = end - start;
    let length2 = edge.magnitude2();
    if !length2.is_finite() || length2 <= 0.0 {
        return None;
    }
    let t = ((point - start).dot(edge) / length2).clamp(0.0, 1.0);
    ((point - (start + edge * t)).magnitude() <= tolerance).then_some(t)
}

impl IsoBoundary {
    fn includes(&self, curve: &NurbsCurve<Vector4>, tolerance: f64) -> bool {
        if proof_bezier_spans(self.curve.non_rationalized()).is_none() {
            return false;
        }
        if curve == &self.curve || curve == &self.curve.inverse() {
            return true;
        }
        // A straight rational boundary need not have the same degree,
        // knots or projective weights as a topological line edge. Prove
        // both loci independently. Candidate endpoints may be subsegments
        // of this boundary, so compare them to the complete finite line.
        let [a, b] = [self.curve.front(), self.curve.back()];
        let [c, d] = [curve.front(), curve.back()];
        let budget = tolerance * 0.25;
        segment_fraction(c, a, b, budget).is_some()
            && segment_fraction(d, a, b, budget).is_some()
            && certify_nurbs_segment_locus(&self.curve, a, b, budget)
                == SegmentLocus::WithinTolerance
            && certify_nurbs_segment_locus(curve, c, d, budget) == SegmentLocus::WithinTolerance
    }

    // An exact UV line for a 3-D line requires an affine boundary
    // parameterization, not merely a straight locus. Test every span's
    // homogeneous controls against the affine parameter law. This also
    // recognizes degree-elevated and multiply-knotted straight sides.
    fn affine_line(&self, line: &Line<Point3>) -> Option<Line<Point2>> {
        let spans = proof_bezier_spans(self.curve.non_rationalized())?;
        let start = spans.first()?.start;
        let end = spans.last()?.end;
        let first = self.curve.front();
        let last = self.curve.back();
        let controls = self.curve.non_rationalized().control_points();
        let weight = controls.first()?.w;
        if !controls.iter().all(|p| same_knot(p.w, weight)) {
            return None;
        }
        let tolerance = KERNEL_EPSILON * 0.25;
        for span in spans {
            let degree = crate::usize_to_f64(span.controls.len() - 1);
            for (i, [x, y, z, w]) in span.controls.into_iter().enumerate() {
                let t =
                    (span.end - span.start).mul_add(crate::usize_to_f64(i) / degree, span.start);
                let expected = first + (last - first) * ((t - start) / (end - start));
                let actual = Point3::new(x / w, y / w, z / w);
                let scale = actual
                    .to_vec()
                    .magnitude()
                    .max(expected.to_vec().magnitude())
                    .max(1.0);
                let floor = 256.0 * f64::EPSILON * (degree + 1.0).powi(2) * scale;
                if (actual - expected).magnitude() + floor > tolerance {
                    return None;
                }
            }
        }
        let front = segment_fraction(line.0, first, last, tolerance)?;
        let back = segment_fraction(line.1, first, last, tolerance)?;
        Some(Line(self.uv.evaluate(front), self.uv.evaluate(back)))
    }
}

pub fn includes_clamped_boundary(
    surface: &NurbsSurface<Vector4>,
    curve: &NurbsCurve<Vector4>,
    tolerance: f64,
) -> bool {
    tolerance.is_finite()
        && tolerance > 0.0
        && proof_bezier_spans(curve.non_rationalized()).is_some()
        && boundaries(surface).is_some_and(|boundaries| {
            boundaries
                .iter()
                .any(|boundary| boundary.includes(curve, tolerance))
        })
}

pub fn affine_clamped_boundary(
    surface: &NurbsSurface<Vector4>,
    line: &Line<Point3>,
) -> Option<Line<Point2>> {
    boundaries(surface)?
        .iter()
        .find_map(|boundary| boundary.affine_line(line))
}

/// Prove collapse using every positive-weight homogeneous boundary control.
pub fn collapsed_boundary_parameter(
    surface: &NurbsSurface<Vector4>,
    point: Point3,
    hint: Point2,
    tolerance: f64,
) -> Option<Point2> {
    if !tolerance.is_finite() || tolerance <= 0.0 {
        return None;
    }
    let net = surface.control_points();
    let near = |control: &Vector4| {
        control.w > 0.0
            && (Point3::new(
                control.x / control.w,
                control.y / control.w,
                control.z / control.w,
            ) - point)
                .magnitude()
                <= tolerance
    };
    // A collapsed clamped side has both of its corner controls at the
    // pole. This cheap necessary condition avoids constructing or
    // extracting any boundary curve for regular interior samples.
    let first = net.first()?;
    let last = net.last()?;
    if ![
        (first.first()?, first.last()?),
        (last.first()?, last.last()?),
        (first.first()?, last.first()?),
        (first.last()?, last.last()?),
    ]
    .into_iter()
    .any(|(a, b)| near(a) && near(b))
    {
        return None;
    }
    boundaries(surface)?.iter().find_map(|boundary| {
        let spans = proof_bezier_spans(boundary.curve.non_rationalized())?;
        let certified = spans
            .iter()
            .flat_map(|span| &span.controls)
            .all(|&[x, y, z, w]| {
                let control = Point3::new(x / w, y / w, z / w);
                let margin = 256.0 * f64::EPSILON * control.to_vec().magnitude().max(1.0);
                (control - point).magnitude() + margin <= tolerance
            });
        if !certified {
            return None;
        }
        let delta = boundary.uv.1 - boundary.uv.0;
        let t = ((hint - boundary.uv.0).dot(delta) / delta.magnitude2()).clamp(0.0, 1.0);
        Some(boundary.uv.evaluate(t))
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn collapsed_surface() -> Result<NurbsSurface<Vector4>> {
        Ok(NurbsSurface::new(BsplineSurface::new(
            (KnotVector::bezier_knot(2), KnotVector::bezier_knot(1)),
            vec![
                vec![
                    Vector4::new(0.0, 0.0, 0.0, 1.0),
                    Vector4::new(0.5, 1.0, 0.0, 1.0),
                ],
                vec![
                    Vector4::new(0.6, 0.0, 0.3, 0.6),
                    Vector4::new(0.3, 0.6, 0.0, 0.6),
                ],
                vec![
                    Vector4::new(2.0, 0.0, 0.0, 1.0),
                    Vector4::new(0.5, 1.0, 0.0, 1.0),
                ],
            ],
        )?))
    }

    #[test]
    fn collapsed_side_and_subsegment_have_structural_incidence() -> Result<()> {
        let surface = collapsed_surface()?;
        let side = Line(Point3::new(2.0, 0.0, 0.0), Point3::new(0.5, 1.0, 0.0));
        for line in [
            side,
            Line(side.evaluate(0.25), side.evaluate(0.75)),
            side.inverse(),
        ] {
            let curve = NurbsCurve::from(BsplineCurve::from(line));
            assert!(includes_clamped_boundary(&surface, &curve, KERNEL_EPSILON));
            let uv = affine_clamped_boundary(&surface, &line).expect("affine boundary");
            for i in 0..=10 {
                let t = f64::from(i) / 10.0;
                let p = uv.evaluate(t);
                assert!((surface.evaluate(p.x, p.y) - line.evaluate(t)).magnitude() < 1.0e-12);
            }
        }
        let off = Line(side.0 + Vector3::unit_z() * 0.01, side.1);
        assert!(!includes_clamped_boundary(
            &surface,
            &NurbsCurve::from(BsplineCurve::from(off)),
            KERNEL_EPSILON
        ));
        assert!(affine_clamped_boundary(&surface, &off).is_none());
        Ok(())
    }
}

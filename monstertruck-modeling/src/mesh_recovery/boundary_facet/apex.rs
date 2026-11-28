//! Deterministic proposals for non-centroid facet fans.
//!
//! The continuous rational-curve certificate remains authoritative.

use super::{
    BoundaryFacetApproximation, BoundaryFacetError, FacetFan, Point3, Vector3, Vector4, centroid,
    cross, directed_pairs, finite, proof_bezier_spans, rational_point,
    subdivide_homogeneous_bezier,
};

// Clip an apex proposal polygon to the left of a directed projected line.
// This is candidate generation only; it cannot waive the final certificate.
fn clip_apex_polygon(points: &mut Vec<Point3>, normal: Vector3, origin: Point3, tangent: Vector3) {
    if tangent.norm_squared() == 0.0 || points.is_empty() {
        return;
    }
    let signed = |p: Point3| cross(normal, tangent, p - origin);
    let old = std::mem::take(points);
    for (&a, &b) in directed_pairs(&old) {
        let (da, db) = (signed(a), signed(b));
        if da >= 0.0 {
            points.push(a);
        }
        if (da < 0.0 && db > 0.0) || (da > 0.0 && db < 0.0) {
            points.push(a + (b - a) * (da / (da - db)));
        }
    }
}

// Endpoint control-polygon directions are the true rational tangent
// directions for positive weights. Subdivision gives additional samples
// without relying on the curve's parameter speed.
fn clip_span_tangents(
    points: &mut Vec<Point3>,
    normal: Vector3,
    controls: &[[f64; 4]],
    depth: usize,
    remaining: &mut usize,
) -> Option<()> {
    *remaining = remaining.checked_sub(1)?;
    if points.is_empty() {
        return Some(());
    }
    if depth > 0 {
        let (left, right) = subdivide_homogeneous_bezier(controls)?;
        clip_span_tangents(points, normal, &left, depth - 1, remaining)?;
        return clip_span_tangents(points, normal, &right, depth - 1, remaining);
    }
    let point = |[x, y, z, w]: [f64; 4]| rational_point(Vector4::new(x, y, z, w));
    let [a, b, ..] = controls else {
        return None;
    };
    let [.., c, d] = controls else {
        return None;
    };
    let (a, b, c, d) = (point(*a), point(*b), point(*c), point(*d));
    clip_apex_polygon(points, normal, a, b - a);
    clip_apex_polygon(points, normal, d, d - c);
    Some(())
}

impl FacetFan<'_> {
    pub(super) fn recover(
        &self,
        reference: Point3,
        original_error: BoundaryFacetError,
    ) -> Result<BoundaryFacetApproximation, BoundaryFacetError> {
        // Restrict the proposal to the intersection of the projected source
        // and native chord polygons, then progressively add exact curve
        // tangent samples. A sampled proposal never substitutes for proof.
        let mut polygon = self.source.to_vec();
        for (&a, &b) in directed_pairs(self.native) {
            clip_apex_polygon(&mut polygon, self.normal, a, b - a);
        }
        let spans = self
            .curves
            .iter()
            .map(proof_bezier_spans)
            .collect::<Option<Vec<_>>>()
            .ok_or(original_error)?;
        let mut remaining = 32768;
        for depth in [0, 2, 4, 6] {
            for span in spans.iter().flatten() {
                clip_span_tangents(
                    &mut polygon,
                    self.normal,
                    &span.controls,
                    depth,
                    &mut remaining,
                )
                .ok_or(original_error)?;
                if polygon.len() < 3 {
                    return Err(original_error);
                }
            }
            let candidate = centroid(&polygon);
            if !finite(candidate) {
                return Err(original_error);
            }
            let penalty = self.source_fan_displacement(reference, candidate);
            if let Ok(result) = self.construct(candidate, penalty) {
                return Ok(result);
            }
        }
        Err(original_error)
    }

    // Both projected source fans cover the same convex polygon.
    // Their heights lie in the interval spanned by the original vertices
    // and their two apices. Charge the entire width (with a numerical
    // margin) to the budget; planar triangles incur only roundoff.
    fn source_fan_displacement(&self, reference: Point3, candidate: Point3) -> f64 {
        let (low, high, scale) = self
            .source
            .iter()
            .copied()
            .chain([candidate])
            .map(|p| p - reference)
            .fold((0.0_f64, 0.0_f64, 0.0_f64), |(lo, hi, scale), p| {
                let h = self.normal.dot(&p);
                (lo.min(h), hi.max(h), scale.max(p.norm()))
            });
        (64.0 * f64::EPSILON).mul_add(scale, high - low)
    }
}

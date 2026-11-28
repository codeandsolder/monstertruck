//! Oblique projection proposals for facets whose curved boundaries are not
//! star-shaped in the original facet plane. A proposal only changes the
//! projection used for the injectivity proof, not the constructed geometry.

use super::tangent::visit_span_tangents;
use super::{
    BoundaryFacetApproximation, BoundaryFacetError, FacetFan, InnerSpace, Point3, Vector3,
    directed_pairs, proof_bezier_spans,
};

// A direction in R^3 is represented by its intersection with one of the six
// faces of a unit cube. Positive scaling preserves every homogeneous
// orientation inequality. Clip each cube face to the sampled normal cone.
fn clip(polygon: &mut Vec<Vector3>, constraint: Vector3) {
    let old = std::mem::take(polygon);
    for (&a, &b) in directed_pairs(&old) {
        let (da, db) = (constraint.dot(a), constraint.dot(b));
        if da >= 0.0 {
            polygon.push(a);
        }
        if (da < 0.0 && db > 0.0) || (da > 0.0 && db < 0.0) {
            polygon.push(a + (b - a) * (da / (da - db)));
        }
    }
}

struct NormalCone {
    faces: [Vec<Vector3>; 6],
}

impl NormalCone {
    fn new() -> Self {
        Self {
            faces: std::array::from_fn(|i| {
                let axis = i / 2;
                let side = if i % 2 == 0 { -1.0 } else { 1.0 };
                [(-1.0, -1.0), (1.0, -1.0), (1.0, 1.0), (-1.0, 1.0)]
                    .into_iter()
                    .map(|(u, v)| {
                        let mut p = Vector3::new(0.0, 0.0, 0.0);
                        p[axis] = side;
                        p[(axis + 1) % 3] = u;
                        p[(axis + 2) % 3] = v;
                        p
                    })
                    .collect()
            }),
        }
    }

    fn constrain(&mut self, constraint: Vector3) {
        let length = constraint.magnitude();
        // Degenerate samples are not a proof of orientation or its absence.
        // The authoritative continuous certificate handles them.
        if length > 0.0 && length.is_finite() {
            let normal = constraint / length;
            for polygon in &mut self.faces {
                clip(polygon, normal);
            }
        }
    }

    fn candidates(&self) -> impl Iterator<Item = Vector3> + '_ {
        self.faces.iter().filter_map(|polygon| {
            if polygon.len() < 3 {
                return None;
            }
            // An average lies in the convex feasible face, avoiding the
            // sampled inequalities' extreme rays. It is only a proposal.
            let sum = polygon.iter().copied().sum::<Vector3>();
            (sum.magnitude2() > 0.0).then(|| sum.normalize())
        })
    }
}

impl FacetFan<'_> {
    pub(super) fn recover_projection(
        &self,
        apex: Point3,
        original_error: BoundaryFacetError,
    ) -> Result<BoundaryFacetApproximation, BoundaryFacetError> {
        let mut cone = NormalCone::new();
        // Keep the source winding and both projected chord polygons convex.
        cone.constrain(self.normal);
        for polygon in [self.source, self.native] {
            for (&a, &b) in directed_pairs(polygon) {
                cone.constrain((b - a).cross(apex - a));
                for &p in polygon {
                    cone.constrain((b - a).cross(p - a));
                }
            }
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
                let viable = visit_span_tangents(
                    &span.controls,
                    depth,
                    &mut remaining,
                    &mut |p, tangent| {
                        cone.constrain((p - apex).cross(tangent));
                        cone.faces.iter().any(|polygon| polygon.len() >= 3)
                    },
                )
                .ok_or(original_error)?;
                if !viable {
                    return Err(original_error);
                }
            }
            for normal in cone.candidates() {
                let fan = FacetFan { normal, ..*self };
                // The apex, exact curves, and reference source fans are
                // unchanged. There is no source displacement penalty.
                if let Ok(result) = fan.construct(apex, 0.0) {
                    return Ok(result);
                }
            }
        }
        Err(original_error)
    }
}

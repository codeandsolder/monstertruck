//! Interval Bernstein controls of the actual homogeneous spline. Use the
//! original coefficients directly for C0 Bezier chains; otherwise evaluate
//! the exact positive-knot-span blossom with interval coefficients.
use super::interval::Interval;
use crate::{BsplineCurve, Vector4};
use monstertruck_geometry::proof_bezier_spans;

pub(super) type Control = [Interval; 4];

pub(super) struct Controls {
    pub(super) values: Vec<Control>,
    pub(super) parameter: [f64; 2],
}

impl Controls {
    pub(super) fn split(&self) -> [Self; 2] {
        let mut work = self.values.clone();
        let mut first = vec![work[0]];
        let mut last = vec![work[work.len() - 1]];
        for count in (1..work.len()).rev() {
            for i in 0..count {
                work[i] = std::array::from_fn(|k| work[i][k].add(work[i + 1][k]).scale(0.5));
            }
            first.push(work[0]);
            last.push(work[count - 1]);
        }
        last.reverse();
        let middle = self.parameter[0].midpoint(self.parameter[1]);
        [
            Self {
                values: first,
                parameter: [self.parameter[0], middle],
            },
            Self {
                values: last,
                parameter: [middle, self.parameter[1]],
            },
        ]
    }
}

fn point(value: Vector4) -> Control { [value.x, value.y, value.z, value.w].map(Interval::point) }

fn blend(a: Control, b: Control, alpha: Interval) -> Control {
    // The true coefficient is a convex combination. Intersect the interval
    // evaluation with that exact hull to keep enclosures tight.
    std::array::from_fn(|i| {
        let hull = a[i].union(b[i]);
        a[i].add(b[i].sub(a[i]).mul(alpha))
            .intersection(hull)
            .unwrap_or(hull)
    })
}

fn same_knot(a: f64, b: f64) -> bool {
    a.total_cmp(&b).is_eq() || (a.abs().to_bits() == 0 && b.abs().to_bits() == 0)
}

pub(super) fn spans(curve: &BsplineCurve<Vector4>) -> Option<Vec<Controls>> {
    // Validate the active layout, positive weights, degree and proof work
    // budget using the existing canonical knot partition. No knot is
    // clustered using the geometric tolerance.
    let partition = proof_bezier_spans(curve)?;
    let knots = curve.knot_vector();
    let points = curve.control_points();
    let degree = curve.degree();
    let mut result = Vec::with_capacity(partition.len());
    let active = (degree..points.len()).filter(|&id| knots[id] < knots[id + 1]);
    for (span, index) in partition.into_iter().zip(active) {
        if !same_knot(knots[index], span.start) || !same_knot(knots[index + 1], span.end) {
            return None;
        }
        let source = &points[index - degree..=index];
        // Both span ends are clamped to at least the degree: the source
        // coefficients already are the exact Bernstein coefficients.
        let direct = (0..degree).all(|i| {
            same_knot(knots[index - i], span.start) && same_knot(knots[index + 1 + i], span.end)
        });
        let values = if direct {
            source.iter().copied().map(point).collect()
        } else {
            let mut values = Vec::with_capacity(degree + 1);
            for repeats in 0..=degree {
                let mut work = source.iter().copied().map(point).collect::<Vec<_>>();
                for round in 1..=degree {
                    let parameter = if round <= degree - repeats {
                        span.start
                    } else {
                        span.end
                    };
                    for local in (round..=degree).rev() {
                        let left = Interval::point(knots[index - degree + local]);
                        let right = Interval::point(knots[index + 1 + local - round]);
                        let alpha = Interval::point(parameter)
                            .sub(left)
                            .div(right.sub(left))?
                            .intersection(Interval::range(0.0, 1.0))?;
                        work[local] = blend(work[local - 1], work[local], alpha);
                    }
                }
                values.push(work[degree]);
            }
            values
        };
        result.push(Controls {
            values,
            parameter: [span.start, span.end],
        });
    }
    Some(result)
}

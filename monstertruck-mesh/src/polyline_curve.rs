use crate::{
    BoundedCurve, CurveParameter, Cut, EuclideanSpace, HashGen, InnerSpace, Invertible, Origin,
    ParameterDivision1D, ParameterRange, ParametricCurve, Point2, Point3, PolylineCurve,
    SearchNearestParameter, SearchParameter, SearchParameterHint1D, SnapCurveEndpoints, Tolerance,
    Transform, Transformed, Vector2, Zero,
};
use itertools::Itertools;
use monstertruck_core::cgmath64::control_point::ControlPoint;
use num_traits::ToPrimitive;
use std::ops::{Bound, Deref, DerefMut};

fn usize_to_f64(value: usize) -> f64 { value.to_f64().unwrap_or(f64::INFINITY) }

fn interior_vertex_indices(
    range: (f64, f64),
    len: usize,
) -> Option<std::ops::RangeInclusive<usize>> {
    if len == 0 || !range.0.is_finite() || !range.1.is_finite() || range.0 >= range.1 {
        return None;
    }

    let first = if range.0 < 0.0 {
        0
    } else {
        range.0.floor().to_usize()?.saturating_add(1)
    };
    if range.1 <= 0.0 {
        return None;
    }
    let last = range.1.ceil().to_usize()?.checked_sub(1)?.min(len - 1);
    (first <= last).then_some(first..=last)
}

impl PolylineCurve<Point2> {
    /// Signed area of area enclosed when endpoints are connected
    /// # Example
    /// ```
    /// use monstertruck_mesh::*;
    /// use std::f64::consts::PI;
    /// let mut hexagon = (0..6).map(|i| {
    ///         let t = PI / 3.0 * i as f64;
    ///         Point2::new(f64::cos(t), f64::sin(t))
    ///     })
    ///     .collect::<PolylineCurve<_>>();
    /// assert_near!(hexagon.area(), 1.5 * f64::sqrt(3.0));
    /// hexagon.invert();
    /// assert_near!(hexagon.area(), -1.5 * f64::sqrt(3.0));
    /// ```
    #[inline]
    #[must_use]
    pub fn area(&self) -> f64 {
        let sum = |sum, (p, q): (&Point2, &Point2)| (q.x + p.x).mul_add(q.y - p.y, sum);
        self.iter().circular_tuple_windows().fold(0.0, sum) / 2.0
    }

    /// whether `c` is included in enclosed domain when endpoints are connected
    /// # Example
    /// ```
    /// use monstertruck_mesh::*;
    /// use std::f64::consts::PI;
    /// let hexagon = (0..6).map(|i| {
    ///         let t = PI / 3.0 * i as f64;
    ///         Point2::new(f64::cos(t), f64::sin(t))
    ///     })
    ///     .collect::<PolylineCurve<_>>();
    /// let p0 = Point2::new(0.9, 0.0);
    /// assert!(hexagon.include(p0));
    /// let p1 = Point2::new(0.0, 1.0);
    /// assert!(!hexagon.include(p1));
    /// ```
    #[allow(
        clippy::many_single_char_names,
        reason = "2D winding calculation uses conventional point, vector, and cross-product notation"
    )]
    #[must_use]
    pub fn include(&self, c: Point2) -> bool {
        if self.iter().any(|p| (*p - c).so_small()) {
            return true;
        }
        let t = 2.0 * std::f64::consts::PI * HashGen::hash1(c);
        let r = Vector2::new(f64::cos(t), f64::sin(t));
        self.iter()
            .circular_tuple_windows()
            .try_fold(0_i32, move |counter, (p0, p1)| {
                if (*p0 - c).so_small() {
                    return None;
                }
                let a = p0 - c;
                let b = p1 - c;
                let s0 = r.x.mul_add(a.y, -(r.y * a.x)); // v times a
                let s1 = r.x.mul_add(b.y, -(r.y * b.x)); // v times b
                let s2 = a.x.mul_add(b.y, -(a.y * b.x)); // a times b
                let x = s2 / (s1 - s0);
                if x.so_small() && s0 * s1 < 0.0 {
                    None
                } else if x > 0.0 && s0 <= 0.0 && s1 > 0.0 {
                    Some(counter + 1)
                } else if x > 0.0 && s0 >= 0.0 && s1 < 0.0 {
                    Some(counter - 1)
                } else {
                    Some(counter)
                }
            })
            .is_none_or(|counter| counter > 0)
    }
}

/// Calculate the area of a region bounded by multiple polylines
/// # Example
/// ```
/// use monstertruck_mesh::*;
/// use std::f64::consts::PI;
/// // the outer boundary is counter-clockwise
/// let large_hexagon = (0..6).map(|i| {
///         let t = PI / 3.0 * i as f64;
///         Point2::new(2.0 * f64::cos(t), 2.0 * f64::sin(t))
///     })
///     .collect::<PolylineCurve<_>>();
/// // the inner boundary is clockwise
/// let small_hexagon = (0..6).map(|i| {
///         let t = PI / 3.0 * i as f64;
///         Point2::new(f64::sin(t), f64::cos(t))
///     })
///     .collect::<PolylineCurve<_>>();
/// let boundaries = [large_hexagon, small_hexagon];
/// assert_near!(polyline_curve::area(&boundaries), 4.5 * f64::sqrt(3.0));
/// ```
pub fn area<'a>(boundaries: impl IntoIterator<Item = &'a PolylineCurve<Point2>>) -> f64 {
    boundaries
        .into_iter()
        .map(super::PolylineCurve::area)
        .sum::<f64>()
}

/// whether `c` is included in enclosed multiple domains when endpoints of each polyline are connected
/// # Example
/// ```
/// use monstertruck_mesh::*;
/// use std::f64::consts::PI;
/// // the outer boundary is counter-clockwise
/// let large_hexagon = (0..6).map(|i| {
///         let t = PI / 3.0 * i as f64;
///         Point2::new(2.0 * f64::cos(t), 2.0 * f64::sin(t))
///     })
///     .collect::<PolylineCurve<_>>();
/// // the inner boundary is clockwise
/// let small_hexagon = (0..6).map(|i| {
///         let t = PI / 3.0 * i as f64;
///         Point2::new(f64::sin(t), f64::cos(t))
///     })
///     .collect::<PolylineCurve<_>>();
/// let boundaries = [large_hexagon, small_hexagon];
/// let p0 = Point2::new(1.5, 0.0);
/// assert!(polyline_curve::include(&boundaries, p0));
/// let p1 = Point2::new(0.0, 0.0);
/// assert!(!polyline_curve::include(&boundaries, p1));
/// ```
#[allow(
    clippy::many_single_char_names,
    reason = "2D winding calculation uses conventional point, vector, and cross-product notation"
)]
pub fn include<'a>(
    boundaries: impl IntoIterator<Item = &'a PolylineCurve<Point2>>,
    c: Point2,
) -> bool {
    let t = 2.0 * std::f64::consts::PI * HashGen::hash1(c);
    let r = Vector2::new(f64::cos(t), f64::sin(t));
    boundaries
        .into_iter()
        .flat_map(|boundary| boundary.iter().circular_tuple_windows())
        .try_fold(0_i32, move |counter, (p0, p1)| {
            if (*p0 - c).so_small() {
                return None;
            }
            let a = p0 - c;
            let b = p1 - c;
            let s0 = r.x.mul_add(a.y, -(r.y * a.x)); // v times a
            let s1 = r.x.mul_add(b.y, -(r.y * b.x)); // v times b
            let s2 = a.x.mul_add(b.y, -(a.y * b.x)); // a times b
            let x = s2 / (s1 - s0);
            if x.so_small() && s0 * s1 < 0.0 {
                None
            } else if x > 0.0 && s0 <= 0.0 && s1 > 0.0 {
                Some(counter + 1)
            } else if x > 0.0 && s0 >= 0.0 && s1 < 0.0 {
                Some(counter - 1)
            } else {
                Some(counter)
            }
        })
        .is_none_or(|counter| counter > 0)
}

impl<P> AsRef<Vec<P>> for PolylineCurve<P> {
    #[inline]
    fn as_ref(&self) -> &Vec<P> { &self.0 }
}

impl<P> AsMut<Vec<P>> for PolylineCurve<P> {
    #[inline]
    fn as_mut(&mut self) -> &mut Vec<P> { &mut self.0 }
}

impl<P> AsRef<[P]> for PolylineCurve<P> {
    #[inline]
    fn as_ref(&self) -> &[P] { &self.0 }
}

impl<P> AsMut<[P]> for PolylineCurve<P> {
    #[inline]
    fn as_mut(&mut self) -> &mut [P] { &mut self.0 }
}

impl<P> Deref for PolylineCurve<P> {
    type Target = Vec<P>;
    #[inline]
    fn deref(&self) -> &Self::Target { &self.0 }
}

impl<P> DerefMut for PolylineCurve<P> {
    #[inline]
    fn deref_mut(&mut self) -> &mut Self::Target { &mut self.0 }
}

impl<P> From<Vec<P>> for PolylineCurve<P> {
    #[inline]
    fn from(v: Vec<P>) -> Self { Self(v) }
}

impl<P> From<PolylineCurve<P>> for Vec<P> {
    #[inline]
    fn from(v: PolylineCurve<P>) -> Self { v.0 }
}

impl<P> FromIterator<P> for PolylineCurve<P> {
    #[inline]
    fn from_iter<I: IntoIterator<Item = P>>(iter: I) -> Self { Self(Vec::from_iter(iter)) }
}

impl<P> IntoIterator for PolylineCurve<P> {
    type Item = P;
    type IntoIter = std::vec::IntoIter<P>;
    #[inline]
    fn into_iter(self) -> Self::IntoIter { self.0.into_iter() }
}

impl<'a, P> IntoIterator for &'a PolylineCurve<P> {
    type Item = &'a P;
    type IntoIter = std::slice::Iter<'a, P>;
    #[inline]
    fn into_iter(self) -> Self::IntoIter { self.0.iter() }
}

impl<P: ControlPoint<f64>> ParametricCurve for PolylineCurve<P> {
    type Point = P;
    type Vector = P::Diff;
    #[inline]
    fn derivative_n(&self, n: usize, t: f64) -> Self::Vector {
        match n {
            0 => self.evaluate(t).to_vec(),
            1 => self.derivative(t),
            _ => Self::Vector::zero(),
        }
    }
    #[inline]
    fn evaluate(&self, t: f64) -> P {
        match self.len() {
            0 => P::origin(),
            1 => self[0],
            _ if t.is_nan() || t <= 0.0 => self[0],
            _ if t + 1.0 >= usize_to_f64(self.len()) => self[self.len() - 1],
            _ => {
                let n = t.to_usize().unwrap_or(0).min(self.len() - 2);
                let local_t = t - usize_to_f64(n);
                self[n] + (self[n + 1] - self[n]) * local_t
            }
        }
    }
    #[inline]
    fn derivative(&self, t: f64) -> P::Diff {
        if self.len() <= 1 || t < 0.0 || usize_to_f64(self.len()) < t + 1.0 {
            P::Diff::zero()
        } else {
            let n = t.to_usize().unwrap_or(0).min(self.len() - 1);
            if n + 1 == self.len() {
                self[n] - self[n - 1]
            } else {
                self[n + 1] - self[n]
            }
        }
    }
    #[inline]
    fn derivative_2(&self, _: f64) -> P::Diff { P::Diff::zero() }
    #[inline]
    fn parameter_range(&self) -> ParameterRange {
        (
            Bound::Included(0.0),
            Bound::Included(usize_to_f64(self.len()) - 1.0),
        )
    }
}

impl<P: ControlPoint<f64>> BoundedCurve for PolylineCurve<P> {
    fn range_tuple(&self) -> (f64, f64) { (0.0, usize_to_f64(self.len()) - 1.0) }
}

impl<P: Clone> Invertible for PolylineCurve<P> {
    #[inline]
    fn invert(&mut self) { self.reverse(); }
    #[inline]
    fn inverse(&self) -> Self { Self(self.iter().rev().cloned().collect()) }
}

impl SnapCurveEndpoints for PolylineCurve<Point3> {
    fn snap_endpoints(&mut self, front: Point3, back: Point3) {
        if let Some(point) = self.first_mut() {
            *point = front;
        }
        if let Some(point) = self.last_mut() {
            *point = back;
        }
    }
}

impl<P: ControlPoint<f64>> Cut for PolylineCurve<P> {
    fn cut(&mut self, t: f64) -> Self {
        if t < 0.0 {
            Self(Vec::new())
        } else if t + 1.0 > usize_to_f64(self.len()) {
            let mut v = Vec::new();
            v.append(&mut self.0);
            Self(v)
        } else {
            let n = t.to_usize().unwrap_or(0).min(self.len().saturating_sub(1));
            if t.near(&usize_to_f64(n)) {
                let mut v = Vec::new();
                v.extend(&self[n..]);
                self.truncate(n + 1);
                Self(v)
            } else {
                let p = self.evaluate(t);
                let mut v = vec![p];
                v.extend(&self[(n + 1)..]);
                self.truncate(n + 1);
                self.push(p);
                Self(v)
            }
        }
    }
}

impl<P> SearchParameter<CurveParameter> for PolylineCurve<P>
where
    P: ControlPoint<f64>,
    P::Diff: InnerSpace<Scalar = f64> + Tolerance,
{
    type Point = P;
    fn search_parameter<H: Into<SearchParameterHint1D>>(
        &self,
        point: P,
        _: H,
        _: usize,
    ) -> Option<f64> {
        for (i, p) in self.0.windows(2).enumerate() {
            let a = point - p[0];
            let b = p[1] - p[0];
            let t = f64::clamp(a.dot(b) / b.dot(b), 0.0, 1.0);
            let h = a - b * t;
            if h.so_small() {
                return Some(t + usize_to_f64(i));
            }
        }
        None
    }
}

impl<P> SearchNearestParameter<CurveParameter> for PolylineCurve<P>
where
    P: ControlPoint<f64>,
    P::Diff: InnerSpace<Scalar = f64>,
{
    type Point = P;
    fn search_nearest_parameter<H: Into<SearchParameterHint1D>>(
        &self,
        point: P,
        _: H,
        _: usize,
    ) -> Option<f64> {
        let (mut t0, mut dist2) = (0.0, f64::INFINITY);
        for (i, p) in self.0.windows(2).enumerate() {
            let a = point - p[0];
            let b = p[1] - p[0];
            let t = f64::clamp(a.dot(b) / b.dot(b), 0.0, 1.0);
            let h = a - b * t;
            if h.dot(h) < dist2 {
                t0 = t + usize_to_f64(i);
                dist2 = h.dot(h);
            }
        }
        Some(t0)
    }
}

impl<P: ControlPoint<f64>> ParameterDivision1D for PolylineCurve<P> {
    type Point = P;
    #[inline]
    fn parameter_division(&self, range: (f64, f64), _: f64) -> (Vec<f64>, Vec<P>) {
        let mut res = (vec![range.0], vec![self.evaluate(range.0)]);
        if let Some(indices) = interior_vertex_indices(range, self.len()) {
            for index in indices {
                res.0.push(usize_to_f64(index));
                res.1.push(self[index]);
            }
        }
        res.0.push(range.1);
        res.1.push(self.evaluate(range.1));
        res
    }
}

impl<P, T> Transformed<T> for PolylineCurve<P>
where
    P: EuclideanSpace,
    T: Transform<P>,
{
    fn transform_by(&mut self, trans: T) {
        self.0
            .iter_mut()
            .for_each(|p| *p = trans.transform_point(*p));
    }
}

#[test]
fn polyline_test() {
    let vec = vec![
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(1.0, 0.0, 0.0),
        Point3::new(0.0, 1.0, 0.0),
        Point3::new(0.0, 0.0, 1.0),
        Point3::new(0.0, 1.0, 1.0),
        Point3::new(1.0, 0.0, 1.0),
        Point3::new(1.0, 1.0, 0.0),
        Point3::new(1.0, 1.0, 1.0),
    ];
    let polyline = PolylineCurve(vec);
    monstertruck_core::assert_near!(polyline.evaluate(2.5), Point3::new(0.0, 0.5, 0.5));

    let mut part0 = polyline.clone();
    let part1 = part0.cut(2.5);
    assert_eq!(part0.len(), 4);
    assert_eq!(part1.len(), 6);
    monstertruck_core::assert_near!(part0.evaluate(2.5), Point3::new(0.0, 0.75, 0.25));
    monstertruck_core::assert_near!(part1.evaluate(0.5), Point3::new(0.0, 0.25, 0.75));
    let mut part0 = polyline.clone();
    let part1 = part0.cut(4.0);
    assert_eq!(part0.len(), 5);
    assert_eq!(part1.len(), 4);
    monstertruck_core::assert_near!(part0.evaluate(2.5), Point3::new(0.0, 0.5, 0.5));
    monstertruck_core::assert_near!(part1.evaluate(0.5), Point3::new(0.5, 0.5, 1.0));

    let pt = polyline.evaluate(2.13);
    let t = polyline.search_parameter(pt, None, 1).unwrap();
    monstertruck_core::assert_near!(t, 2.13);

    let pt = Point3::new(2.0, 0.0, 0.0);
    assert!(polyline.search_parameter(pt, None, 1).is_none());
    monstertruck_core::assert_near!(polyline.search_nearest_parameter(pt, None, 1).unwrap(), 1.0);
    let pt = Point3::new(0.5, 0.5, 0.51);
    assert!(polyline.search_parameter(pt, None, 1).is_none());
    let t = polyline.search_nearest_parameter(pt, None, 1).unwrap();
    assert!(
        polyline
            .derivative(t)
            .dot(pt - polyline.evaluate(t))
            .so_small()
    );

    let div = polyline.parameter_division((1.5, 6.2), 0.0);
    assert_eq!(div.0, vec![1.5, 2.0, 3.0, 4.0, 5.0, 6.0, 6.2]);
    assert_eq!(div.0.len(), div.1.len());
    monstertruck_core::assert_near!(div.1[0], Point3::new(0.5, 0.5, 0.0));
    monstertruck_core::assert_near!(div.1[1], Point3::new(0.0, 1.0, 0.0));
    monstertruck_core::assert_near!(div.1[2], Point3::new(0.0, 0.0, 1.0));
    monstertruck_core::assert_near!(div.1[3], Point3::new(0.0, 1.0, 1.0));
    monstertruck_core::assert_near!(div.1[4], Point3::new(1.0, 0.0, 1.0));
    monstertruck_core::assert_near!(div.1[5], Point3::new(1.0, 1.0, 0.0));
    monstertruck_core::assert_near!(div.1[6], Point3::new(1.0, 1.0, 0.2));
}

#[cfg(test)]
mod parameter_totality_tests {
    use super::*;

    fn line() -> PolylineCurve<Point3> {
        PolylineCurve(vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(2.0, 0.0, 0.0),
            Point3::new(3.0, 0.0, 0.0),
        ])
    }

    #[test]
    fn parameter_division_clamps_interior_vertices_to_curve_domain() {
        let curve = line();
        let (parameters, points) = curve.parameter_division((-2.5, 2.5), 0.0);
        assert_eq!(parameters, vec![-2.5, 0.0, 1.0, 2.0, 2.5]);
        assert_eq!(points.len(), parameters.len());
        assert_eq!(points[0], curve[0]);
        assert_eq!(points[1], curve[0]);
        assert_eq!(points[4], curve.evaluate(2.5));
    }

    #[test]
    fn parameter_division_does_not_duplicate_integer_endpoints() {
        let curve = line();
        let (parameters, points) = curve.parameter_division((1.0, 3.0), 0.0);
        assert_eq!(parameters, vec![1.0, 2.0, 3.0]);
        assert_eq!(points, vec![curve[1], curve[2], curve[3]]);
    }

    #[test]
    fn evaluate_single_point_and_nan_are_total() {
        let single = PolylineCurve(vec![Point3::new(1.0, 2.0, 3.0)]);
        assert_eq!(single.evaluate(f64::NAN), single[0]);
        assert_eq!(single.evaluate(f64::INFINITY), single[0]);
        assert_eq!(single.evaluate(f64::NEG_INFINITY), single[0]);

        let curve = line();
        assert_eq!(curve.evaluate(f64::NAN), curve[0]);
    }

    #[test]
    fn parameter_division_nonfinite_range_does_not_expand_integer_indices() {
        let curve = line();
        let (parameters, points) =
            curve.parameter_division((f64::NEG_INFINITY, f64::INFINITY), 0.0);
        assert_eq!(parameters, vec![f64::NEG_INFINITY, f64::INFINITY]);
        assert_eq!(points, vec![curve[0], curve[3]]);
    }
}

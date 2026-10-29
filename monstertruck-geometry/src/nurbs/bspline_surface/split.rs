#![allow(
    clippy::many_single_char_names,
    reason = "compact symbols in surface splitting follow the standard B-spline recurrence notation"
)]

use super::{
    BoundingBox, BsplineCurve, BsplineSurface, Concat, ControlPoint, CurveCollector, Invertible,
    KnotVector, Tolerance, Vector2,
};

impl<P: ControlPoint<f64> + Tolerance> BsplineSurface<P> {
    /// Cuts the surface into two surfaces at the parameter `u`
    /// # Examples
    /// ```
    /// use monstertruck_geometry::prelude::*;
    ///
    /// let knot_vec0 = KnotVector::uniform_knot(2, 2);
    /// let knot_vec1 = KnotVector::uniform_knot(2, 2);
    /// let control_points0 = vec![
    ///     Vector2::new(0.0, 0.0), Vector2::new(0.5, 0.0), Vector2::new(2.0, 0.0), Vector2::new(2.5, 0.0),
    /// ];
    /// let control_points1 = vec![
    ///     Vector2::new(0.0, 1.0), Vector2::new(0.5, 1.0), Vector2::new(2.0, 1.0), Vector2::new(2.5, 1.0),
    /// ];
    /// let control_points2 = vec![
    ///     Vector2::new(0.0, 1.5), Vector2::new(0.5, 1.5), Vector2::new(2.0, 1.5), Vector2::new(2.5, 1.5),
    /// ];
    /// let control_points3 = vec![
    ///     Vector2::new(0.0, 2.5), Vector2::new(0.5, 2.5), Vector2::new(2.0, 2.5), Vector2::new(2.5, 2.5),
    /// ];
    /// let control_points = vec![control_points0, control_points1, control_points2, control_points3];
    /// let bspsurface = BsplineSurface::new_unchecked((knot_vec0, knot_vec1), control_points);
    ///
    /// let mut part0 = bspsurface.clone();
    /// let part1 = part0.cut_u(0.68);
    /// const N: usize = 100;
    /// for i in 0..=N {
    ///     for j in 0..=N {
    ///         let u = 0.68 * (i as f64) / (N as f64);
    ///         let v = 1.0 * (j as f64) / (N as f64);
    ///         assert_near2!(bspsurface.subs(u, v), part0.subs(u, v));
    ///     }
    /// }
    /// for i in 0..=N {
    ///     for j in 0..=N {
    ///         let u = 0.68 + 0.32 * (i as f64) / (N as f64);
    ///         let v = 1.0 * (j as f64) / (N as f64);
    ///         assert_near2!(bspsurface.subs(u, v), part1.subs(u, v));
    ///     }
    /// }
    /// ```
    #[must_use]
    pub fn cut_u(&mut self, mut u: f64) -> Self {
        let degree = self.udegree();
        let u_start = self.knot_vector_u()[0];
        let u_end = self.knot_vector_u()[self.knot_vector_u().len() - 1];

        // Match curve cutting: parameters within tolerance of a domain endpoint
        // represent that endpoint. If we let such values fall through to the
        // generic knot-snap path, a full-multiplicity endpoint can yield an
        // empty control-net half.
        if u <= u_start && u.near(&u_start) {
            let right = self.clone();
            let knot_vector_v = self.knot_vector_v().clone();
            let control_points = vec![self.control_points[0].clone()];
            *self = Self::new_unchecked(
                (KnotVector::from(vec![u_start, u_start]), knot_vector_v),
                control_points,
            );
            return right;
        } else if u >= u_end && u.near(&u_end) {
            let knot_vector_v = self.knot_vector_v().clone();
            let row = self.control_points[self.control_points.len() - 1].clone();
            return Self::new_unchecked(
                (KnotVector::from(vec![u_end, u_end]), knot_vector_v),
                vec![row],
            );
        }

        let Some(idx) = self.knot_vector_u().floor(u) else {
            let bspline = self.clone();
            let knot_vector_u = KnotVector::from(vec![u, self.knot_vector_u()[0]]);
            let knot_vector_v = self.knot_vector_v().clone();
            let control_points = vec![vec![P::origin(); knot_vector_v.len()]];
            *self = Self::new_unchecked((knot_vector_u, knot_vector_v), control_points);
            return bspline;
        };
        let candidate_knot = self.knot_vector_u()[idx];
        let candidate_is_start = candidate_knot.total_cmp(&u_start).is_eq();
        let candidate_is_end = candidate_knot.total_cmp(&u_end).is_eq();
        let interior_near_endpoint =
            (candidate_is_start && u > u_start) || (candidate_is_end && u < u_end);
        let s = if u.near(&candidate_knot) && !interior_near_endpoint {
            u = candidate_knot;
            self.knot_vector_u().multiplicity(idx)
        } else {
            0
        };

        let mut k = idx;
        for _ in s..=degree {
            k = self.add_knot_u_with_index(u);
        }

        let knot_vector_v = self.knot_vector_v().clone();
        let m = self.knot_vector_u().len();
        let n = self.control_points.len();
        let knot_vec0 = self.knot_vector_u().sub_vec(0..=k);
        let knot_vec1 = self.knot_vector_u().sub_vec((k - degree)..m);
        let control_points0 = Vec::from(&self.control_points[0..(k - degree)]);
        let control_points1 = Vec::from(&self.control_points[(k - degree)..n]);
        *self = Self::new_unchecked((knot_vec0, knot_vector_v.clone()), control_points0);
        Self::new_unchecked((knot_vec1, knot_vector_v), control_points1)
    }
    /// Cuts the curve to two curves at the parameter `t`
    /// # Examples
    /// ```
    /// use monstertruck_geometry::prelude::*;
    ///
    /// let knot_vec0 = KnotVector::uniform_knot(2, 2);
    /// let knot_vec1 = KnotVector::uniform_knot(2, 2);
    /// let control_points0 = vec![
    ///     Vector2::new(0.0, 0.0), Vector2::new(0.5, 0.0), Vector2::new(2.0, 0.0), Vector2::new(2.5, 0.0),
    /// ];
    /// let control_points1 = vec![
    ///     Vector2::new(0.0, 1.0), Vector2::new(0.5, 1.0), Vector2::new(2.0, 1.0), Vector2::new(2.5, 1.0),
    /// ];
    /// let control_points2 = vec![
    ///     Vector2::new(0.0, 1.5), Vector2::new(0.5, 1.5), Vector2::new(2.0, 1.5), Vector2::new(2.5, 1.5),
    /// ];
    /// let control_points3 = vec![
    ///     Vector2::new(0.0, 2.5), Vector2::new(0.5, 2.5), Vector2::new(2.0, 2.5), Vector2::new(2.5, 2.5),
    /// ];
    /// let control_points = vec![control_points0, control_points1, control_points2, control_points3];
    /// let bspsurface = BsplineSurface::new_unchecked((knot_vec0, knot_vec1), control_points);
    ///
    /// let mut part0 = bspsurface.clone();
    /// let part1 = part0.cut_v(0.68);
    /// const N: usize = 100;
    /// for i in 0..=N {
    ///     for j in 0..=N {
    ///         let u = 1.0 * (i as f64) / (N as f64);
    ///         let v = 0.68 * (j as f64) / (N as f64);
    ///         assert_near2!(bspsurface.subs(u, v), part0.subs(u, v));
    ///     }
    /// }
    /// for i in 0..=N {
    ///     for j in 0..=N {
    ///         let u = 1.0 * (i as f64) / (N as f64);
    ///         let v = 0.68 + 0.32 * (j as f64) / (N as f64);
    ///         assert_near2!(bspsurface.subs(u, v), part1.subs(u, v));
    ///     }
    /// }
    /// ```
    #[must_use]
    pub fn cut_v(&mut self, v: f64) -> Self {
        self.swap_axes();
        let mut res = self.cut_u(v);
        self.swap_axes();
        res.swap_axes();
        res
    }

    /// Creates a sectional curve with normalized knot vector from the parameter `p` to the parameter `q`.
    /// # Examples
    /// ```
    /// use monstertruck_geometry::prelude::*;
    ///
    /// // a parabola surface: x = 2u - 1, y = 2v - 1, z = x^2 + y^z
    /// let knot_vecs = (KnotVector::bezier_knot(2), KnotVector::bezier_knot(2));
    /// let control_points = vec![
    ///     vec![Vector3::new(-1.0, -1.0, 2.0), Vector3::new(-1.0, 0.0, 0.0), Vector3::new(-1.0, 1.0, 2.0)],
    ///     vec![Vector3::new(0.0, -1.0, 0.0), Vector3::new(0.0, 0.0, -2.0), Vector3::new(0.0, 1.0, 0.0)],
    ///     vec![Vector3::new(1.0, -1.0, 2.0), Vector3::new(1.0, 0.0, 0.0), Vector3::new(1.0, 1.0, 2.0)],
    /// ];
    /// let mut bspsurface = BsplineSurface::new_unchecked(knot_vecs, control_points);
    ///
    /// // add some knots for the test!
    /// bspsurface.add_knot_u(0.26);
    /// # bspsurface.add_knot_u(0.64);
    /// bspsurface.add_knot_v(0.23);
    /// # bspsurface.add_knot_v(0.82);
    ///
    /// let bnd_box = BoundingBox::from_iter(&[Vector2::new(0.2, 0.3), Vector2::new(0.8, 0.6)]);
    /// let curve = bspsurface.sectional_curve(bnd_box).expect("section is valid");
    /// const N: usize = 100;
    /// assert_near2!(curve.subs(0.0), bspsurface.subs(0.2, 0.3));
    /// assert_near2!(curve.subs(1.0), bspsurface.subs(0.8, 0.6));
    /// for i in 0..=N {
    ///     let t = i as f64 / N as f64;
    ///     let pt = curve.subs(t);
    ///     assert_near2!(pt[1], pt[0] * 0.5 - 0.1);
    ///     assert_near2!(pt[2], pt[0] * pt[0] + pt[1] * pt[1]);
    /// }
    /// ```
    #[must_use]
    pub fn sectional_curve(&self, bnd_box: BoundingBox<Vector2>) -> Option<BsplineCurve<P>> {
        let p = bnd_box.min();
        let q = bnd_box.max();
        let mut bspsurface = self.clone();
        if !p[0].near(&bspsurface.knot_u(0)) {
            bspsurface = bspsurface.cut_u(p[0]);
        }
        if !q[0].near(&bspsurface.knot_u(bspsurface.knot_vector_u().len() - 1)) {
            let _ = bspsurface.cut_u(q[0]);
        }
        if !p[1].near(&bspsurface.knot_v(0)) {
            bspsurface = bspsurface.cut_v(p[1]);
        }
        if !q[1].near(&bspsurface.knot_v(bspsurface.knot_vector_v().len() - 1)) {
            let _ = bspsurface.cut_v(q[1]);
        }
        bspsurface.syncro_uvdegrees();
        bspsurface.syncro_uvknots();
        let degree = bspsurface.udegree();
        let log_comb = log_combinatorial(degree);
        let log_comb2 = log_combinatorial(degree * 2);
        let (knots, _) = bspsurface.knot_vector_u().to_single_multi();
        let mut cc = CurveCollector::Singleton;
        for p in 1..knots.len() {
            let mut backup = None;
            if p + 1 != knots.len() {
                backup = Some(bspsurface.cut_u(knots[p]));
                let _ = bspsurface.cut_v(knots[p]);
            }
            let mut knot_vec = KnotVector::bezier_knot(degree * 2);
            knot_vec.translate(crate::usize_to_f64(p) - 1.0);
            let control_points: Vec<_> = (0..=degree * 2)
                .map(|k| {
                    (0..=k).fold(P::origin(), |sum, i| {
                        if i <= degree && k - i <= degree {
                            let coefficient = (log_comb[i] + log_comb[k - i] - log_comb2[k]).exp();
                            sum + bspsurface.control_points[i][k - i].to_vec() * coefficient
                        } else {
                            sum
                        }
                    })
                })
                .collect();
            cc.try_concat(&BsplineCurve::new_unchecked(knot_vec, control_points))
                .ok()?;
            if let Some(mut backup) = backup {
                bspsurface = backup.cut_v(knots[p]);
            }
        }
        let mut curve: BsplineCurve<P> = cc.into_curve()?;
        curve.knot_normalize();
        Some(curve)
    }

    /// Gets the boundary by four splitted curves.
    /// # Examples
    /// ```
    /// use monstertruck_geometry::prelude::*;
    /// let knot_vecs = (KnotVector::bezier_knot(3), KnotVector::bezier_knot(2));
    /// let control_points = vec![
    ///     vec![Vector2::new(0.0, 0.0), Vector2::new(0.5, -1.0), Vector2::new(1.0, 0.0)],
    ///     vec![Vector2::new(0.0, 1.0), Vector2::new(0.5, 1.0), Vector2::new(1.0, 1.0)],
    ///     vec![Vector2::new(0.0, 2.0), Vector2::new(0.5, 2.0), Vector2::new(1.0, 2.0)],
    ///     vec![Vector2::new(0.0, 3.0), Vector2::new(0.5, 3.5), Vector2::new(1.0, 3.0)],
    /// ];
    /// let bspsurface = BsplineSurface::new_unchecked(knot_vecs, control_points);
    /// let curves = bspsurface.splitted_boundary();
    /// assert_eq!(
    ///     curves[0].control_points(),
    ///     &vec![Vector2::new(0.0, 0.0), Vector2::new(0.0, 1.0), Vector2::new(0.0, 2.0), Vector2::new(0.0, 3.0)],
    /// );
    /// assert_eq!(
    ///     curves[1].control_points(),
    ///     &vec![Vector2::new(0.0, 3.0), Vector2::new(0.5, 3.5), Vector2::new(1.0, 3.0)],
    /// );
    /// assert_eq!(
    ///     curves[2].control_points(),
    ///     &vec![Vector2::new(1.0, 3.0), Vector2::new(1.0, 2.0), Vector2::new(1.0, 1.0), Vector2::new(1.0, 0.0)],
    /// );
    /// assert_eq!(
    ///     curves[3].control_points(),
    ///     &vec![Vector2::new(1.0, 0.0), Vector2::new(0.5, -1.0), Vector2::new(0.0, 0.0)],
    /// );
    /// ```
    #[must_use]
    pub fn splitted_boundary(&self) -> [BsplineCurve<P>; 4] {
        let (knot_vector_u, knot_vector_v) = self.knot_vecs.clone();
        let control_points0 = self.control_points.iter().map(|x| x[0]).collect();
        let control_points1 = self.control_points[self.control_points.len() - 1].clone();
        let control_points2 = self
            .control_points
            .iter()
            .map(|row| row[row.len() - 1])
            .collect();
        let control_points3 = self.control_points[0].clone();
        let curve0 = BsplineCurve::new_unchecked(knot_vector_u.clone(), control_points0);
        let curve1 = BsplineCurve::new_unchecked(knot_vector_v.clone(), control_points1);
        let mut curve2 = BsplineCurve::new_unchecked(knot_vector_u, control_points2);
        let mut curve3 = BsplineCurve::new_unchecked(knot_vector_v, control_points3);
        curve2.invert();
        curve3.invert();
        [curve0, curve1, curve2, curve3]
    }

    /// Extracts the boundary of surface
    #[must_use]
    pub fn boundary(&self) -> Option<BsplineCurve<P>> {
        let (knot_vector_u, knot_vector_v) = self.knot_vecs.clone();
        let (range0, range1) = (knot_vector_u.range_length(), knot_vector_v.range_length());
        let [bspline0, mut bspline1, mut bspline2, mut bspline3] = self.splitted_boundary();
        bspline2.invert();
        bspline3.invert();
        let curve = bspline0.try_concat(bspline1.knot_translate(range0)).ok()?;
        let curve = curve
            .try_concat(bspline2.knot_translate(range0 + range1))
            .ok()?;
        curve
            .try_concat(bspline3.knot_translate(f64::mul_add(range0, 2.0, range1)))
            .ok()
    }
}

fn log_combinatorial(n: usize) -> Vec<f64> {
    let mut result = Vec::with_capacity(n + 1);
    result.push(0.0);
    for i in 1..=n {
        let next =
            result[i - 1] + crate::usize_to_f64(n - i + 1).ln() - crate::usize_to_f64(i).ln();
        result.push(next);
    }
    result
}

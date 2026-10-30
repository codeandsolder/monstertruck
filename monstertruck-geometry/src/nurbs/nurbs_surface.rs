#![allow(
    clippy::similar_names,
    reason = "paired U/V range and basis names mirror standard NURBS surface notation"
)]

use super::{
    BasisWindow, Bounded, BoundedCurve, BoundedSurface, BoundingBox, BsplineCurve, BsplineSurface,
    ControlPoint, D2, Error, EuclideanSpace, HashGen, Homogeneous, INCLUDE_CURVE_TRIALS,
    IncludeCurve, InnerSpace, Invertible, KnotVector, MetricSpace, NurbsCurve, NurbsSurface,
    PRESEARCH_DIVISION, ParameterDivision2D, ParameterRange, ParametricCurve, ParametricSurface,
    ParametricSurface3D, Point3, Result, SearchNearestParameter, SearchParameter,
    SearchParameterHint2D, SurfaceParameter, SurfaceParameterRange, TOLERANCE, TOLERANCE2,
    Tolerance, Transformed, Vector3, Vector4, Zero, algo, multi_rational_derivative,
    rational_derivative,
};
#[cfg(test)]
use crate::prelude::Vector2;
use algo::surface::{SearchNearestParameterVector, SearchParameterVector};

impl<V> NurbsSurface<V> {
    /// constructor
    #[inline]
    #[must_use]
    pub const fn new(bspsurface: BsplineSurface<V>) -> Self { Self(bspsurface) }

    /// Returns the nurbs surface before rationalized
    #[inline]
    #[must_use]
    pub const fn non_rationalized(&self) -> &BsplineSurface<V> { &self.0 }
    /// Returns the nurbs surface before rationalized
    #[inline]
    pub const fn non_rationalized_mut(&mut self) -> &mut BsplineSurface<V> { &mut self.0 }

    /// Returns the reference of the knot vectors.
    #[inline]
    #[must_use]
    pub const fn knot_vectors(&self) -> &(KnotVector, KnotVector) { &self.0.knot_vecs }

    /// Returns the u knot vector.
    #[inline]
    #[must_use]
    pub const fn knot_vector_u(&self) -> &KnotVector { &self.0.knot_vecs.0 }
    /// Returns the v knot vector.
    #[inline]
    #[must_use]
    pub const fn knot_vector_v(&self) -> &KnotVector { &self.0.knot_vecs.1 }

    /// Returns the `idx`th u knot.
    #[inline]
    #[must_use]
    pub fn knot_u(&self, idx: usize) -> f64 { self.0.knot_vecs.0[idx] }
    /// returns the `idx`th v knot.
    #[inline]
    #[must_use]
    pub fn knot_v(&self, idx: usize) -> f64 { self.0.knot_vecs.1[idx] }

    /// Returns the reference of the vector of the control points
    #[inline]
    #[must_use]
    pub const fn control_points(&self) -> &Vec<Vec<V>> { &self.0.control_points }

    /// Returns the reference of the control point corresponding to the index `(idx0, idx1)`.
    #[inline]
    #[must_use]
    pub fn control_point(&self, idx0: usize, idx1: usize) -> &V {
        &self.0.control_points[idx0][idx1]
    }

    /// Apply the given transformation to all control points.
    #[inline]
    pub fn transform_control_points<F: FnMut(&mut V)>(&mut self, f: F) {
        self.0.transform_control_points(f);
    }

    /// Returns the iterator over the control points in the `column_idx`th row.
    /// # Examples
    /// ```
    /// use monstertruck_geometry::prelude::*;
    /// let knot_vector_u = KnotVector::bezier_knot(1);
    /// let knot_vector_v = KnotVector::bezier_knot(2);
    /// let knot_vecs = (knot_vector_u, knot_vector_v);
    /// let control_points = vec![
    ///     vec![Vector3::new(0.0, 0.0, 0.0), Vector3::new(1.0, 0.0, 1.0), Vector3::new(2.0, 0.0, 2.0)],
    ///     vec![Vector3::new(0.0, 1.0, 0.0), Vector3::new(1.0, 1.0, 1.0), Vector3::new(2.0, 1.0, 2.0)],
    /// ];
    /// let bspsurface = BsplineSurface::new_unchecked(knot_vecs, control_points);
    /// let mut iter = bspsurface.control_points_row_iter(1);
    /// assert_eq!(iter.next(), Some(&Vector3::new(1.0, 0.0, 1.0)));
    /// assert_eq!(iter.next(), Some(&Vector3::new(1.0, 1.0, 1.0)));
    /// assert_eq!(iter.next(), None);
    /// ```
    #[inline]
    pub fn control_points_row_iter(
        &self,
        column_idx: usize,
    ) -> impl ExactSizeIterator<Item = &V> + std::iter::FusedIterator<Item = &V> {
        self.0.control_points_row_iter(column_idx)
    }

    /// Returns the iterator over the control points in the `row_idx`th row.
    /// # Examples
    /// ```
    /// use monstertruck_geometry::prelude::*;
    /// let knot_vector_u = KnotVector::bezier_knot(1);
    /// let knot_vector_v = KnotVector::bezier_knot(2);
    /// let knot_vecs = (knot_vector_u, knot_vector_v);
    /// let control_points = vec![
    ///     vec![Vector3::new(0.0, 0.0, 0.0), Vector3::new(1.0, 0.0, 1.0), Vector3::new(2.0, 0.0, 2.0)],
    ///     vec![Vector3::new(0.0, 1.0, 0.0), Vector3::new(1.0, 1.0, 1.0), Vector3::new(2.0, 1.0, 2.0)],
    /// ];
    /// let bspsurface = BsplineSurface::new_unchecked(knot_vecs, control_points);
    /// let mut iter = bspsurface.control_points_column_iter(1);
    /// assert_eq!(iter.next(), Some(&Vector3::new(0.0, 1.0, 0.0)));
    /// assert_eq!(iter.next(), Some(&Vector3::new(1.0, 1.0, 1.0)));
    /// assert_eq!(iter.next(), Some(&Vector3::new(2.0, 1.0, 2.0)));
    /// assert_eq!(iter.next(), None);
    /// ```
    #[inline]
    pub fn control_points_column_iter(&self, row_idx: usize) -> std::slice::Iter<'_, V> {
        self.0.control_points[row_idx].iter()
    }

    /// Returns the mutable reference of the control point corresponding to index `(idx0, idx1)`.
    #[inline]
    pub fn control_point_mut(&mut self, idx0: usize, idx1: usize) -> &mut V {
        &mut self.0.control_points[idx0][idx1]
    }

    /// Returns the iterator on all control points
    #[inline]
    pub fn control_points_mut(&mut self) -> impl Iterator<Item = &mut V> {
        self.0.control_points.iter_mut().flatten()
    }
    /// Returns the degrees of B-spline surface
    /// # Examples
    /// ```
    /// use monstertruck_geometry::prelude::*;
    /// let knot_vector_u = KnotVector::from(vec![0.0, 0.0, 1.0, 1.0]);
    /// let knot_vector_v = KnotVector::from(vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0]);
    /// let knot_vecs = (knot_vector_u, knot_vector_v);
    /// let control_points = vec![
    ///     vec![Vector3::new(0.0, 0.0, 0.0), Vector3::new(1.0, 0.0, 1.0), Vector3::new(2.0, 0.0, 2.0)],
    ///     vec![Vector3::new(0.0, 1.0, 0.0), Vector3::new(1.0, 1.0, 1.0), Vector3::new(2.0, 1.0, 2.0)],
    /// ];
    /// let bspsurface = BsplineSurface::new_unchecked(knot_vecs, control_points);
    /// assert_eq!(bspsurface.udegree(), 1);
    /// ```
    #[inline]
    #[must_use]
    pub fn udegree(&self) -> usize { self.0.udegree() }

    /// Returns the degrees of B-spline surface
    /// # Examples
    /// ```
    /// use monstertruck_geometry::prelude::*;
    /// let knot_vector_u = KnotVector::from(vec![0.0, 0.0, 1.0, 1.0]);
    /// let knot_vector_v = KnotVector::from(vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0]);
    /// let knot_vecs = (knot_vector_u, knot_vector_v);
    /// let control_points = vec![
    ///     vec![Vector3::new(0.0, 0.0, 0.0), Vector3::new(1.0, 0.0, 1.0), Vector3::new(2.0, 0.0, 2.0)],
    ///     vec![Vector3::new(0.0, 1.0, 0.0), Vector3::new(1.0, 1.0, 1.0), Vector3::new(2.0, 1.0, 2.0)],
    /// ];
    /// let bspsurface = BsplineSurface::new_unchecked(knot_vecs, control_points);
    /// assert_eq!(bspsurface.vdegree(), 2);
    /// ```
    #[inline]
    #[must_use]
    pub fn vdegree(&self) -> usize { self.0.vdegree() }

    /// Returns the degrees of B-spline surface
    /// # Examples
    /// ```
    /// use monstertruck_geometry::prelude::*;
    /// let knot_vector_u = KnotVector::from(vec![0.0, 0.0, 1.0, 1.0]);
    /// let knot_vector_v = KnotVector::from(vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0]);
    /// let knot_vecs = (knot_vector_u, knot_vector_v);
    /// let control_points = vec![
    ///     vec![Vector3::new(0.0, 0.0, 0.0), Vector3::new(1.0, 0.0, 1.0), Vector3::new(2.0, 0.0, 2.0)],
    ///     vec![Vector3::new(0.0, 1.0, 0.0), Vector3::new(1.0, 1.0, 1.0), Vector3::new(2.0, 1.0, 2.0)],
    /// ];
    /// let bspsurface = BsplineSurface::new_unchecked(knot_vecs, control_points);
    /// assert_eq!(bspsurface.degrees(), (1, 2));
    /// ```
    #[inline]
    #[must_use]
    pub fn degrees(&self) -> (usize, usize) { (self.udegree(), self.vdegree()) }
    /// Returns whether the knot vectors are clamped or not.
    #[inline]
    #[must_use]
    pub fn is_clamped(&self) -> bool { self.0.is_clamped() }
    /// Swaps two parameters.
    pub fn swap_axes(&mut self) -> &mut Self
    where V: Clone {
        self.0.swap_axes();
        self
    }
    /// The range of the parameter of the surface.
    #[inline]
    #[must_use]
    pub fn parameter_range(&self) -> (ParameterRange, ParameterRange) { self.0.parameter_range() }
    /// Creates the sectional curve along v, at the given u control-point index.
    #[inline]
    #[must_use]
    pub fn curve_v(&self, index_u: usize) -> NurbsCurve<V>
    where V: Clone {
        NurbsCurve(self.0.curve_v(index_u))
    }
    /// Creates the sectional curve along u, at the given v control-point index.
    #[inline]
    #[must_use]
    pub fn curve_u(&self, index_v: usize) -> NurbsCurve<V>
    where V: Clone {
        NurbsCurve(self.0.curve_u(index_v))
    }
}

impl<V: Homogeneous<Scalar = f64>> NurbsSurface<V> {
    /// Constructs a rationalization surface from the non-rationalized surface and weights.
    /// # Errors
    /// the length of `surface.control_points()` and `weights` must be the same.
    #[inline]
    pub fn try_from_bspline_and_weights(
        surface: BsplineSurface<V::Point>,
        weights: Vec<Vec<f64>>,
    ) -> Result<Self> {
        let BsplineSurface {
            knot_vecs,
            control_points,
        } = surface;
        if control_points.len() != weights.len() {
            return Err(Error::DifferentLength);
        }
        let control_points = control_points
            .into_iter()
            .zip(weights)
            .map(|(control_points, weights)| {
                if control_points.len() != weights.len() {
                    return Err(Error::DifferentLength);
                }
                Ok(control_points
                    .into_iter()
                    .zip(weights)
                    .map(|(pt, w)| V::from_point_weight(pt, w))
                    .collect::<Vec<_>>())
            })
            .collect::<Result<Vec<_>>>()?;
        Ok(Self(BsplineSurface::new_unchecked(
            knot_vecs,
            control_points,
        )))
    }
}

impl<V: Homogeneous<Scalar = f64> + ControlPoint<f64, Diff = V>> NurbsSurface<V> {
    /// Returns the closure of substitution.
    #[inline]
    pub fn closure(&self) -> impl Fn(f64, f64) -> V::Point + '_ { move |u, v| self.subs(u, v) }
}

impl<V: Homogeneous<Scalar = f64> + ControlPoint<f64, Diff = V>> NurbsSurface<V>
where V::Point: Tolerance
{
    /// Returns whether constant curve or not, i.e. all control points are same or not.
    /// # Examples
    /// ```
    /// use monstertruck_geometry::prelude::*;
    /// let knot_vector_u = KnotVector::bezier_knot(1);
    /// let knot_vector_v = KnotVector::bezier_knot(2);
    /// let pt = Vector3::new(1.0, 2.0, 1.0);
    /// // allows differences upto scalars
    /// let control_points = vec![
    ///     vec![pt.clone(), pt.clone() * 2.0, pt.clone() * 3.0],
    ///     vec![pt.clone() * 0.5, pt.clone() * 0.25, pt.clone() * 0.125],
    /// ];
    /// let mut surface = NurbsSurface::new(BsplineSurface::new_unchecked((knot_vector_u, knot_vector_v), control_points));
    /// assert!(surface.is_const());
    ///
    /// *surface.control_point_mut(1, 2) = Vector3::new(2.0, 3.0, 1.0);
    /// assert!(!surface.is_const());
    /// ```
    #[inline]
    #[must_use]
    pub fn is_const(&self) -> bool {
        let pt = self.0.control_points[0][0].to_point();
        for vec in self.0.control_points.iter().flat_map(|pts| pts.iter()) {
            if !vec.to_point().near(&pt) {
                return false;
            }
        }
        true
    }
    /// Determines whether `self` and `other` is near as the B-spline rational surfaces or not.
    ///
    /// Divides each knot domain into the number of degree equal parts,
    /// and check `|self(u, v) - other(u, v)| < TOLERANCE` for each end points `(u, v)`.
    /// # Examples
    /// ```
    /// use monstertruck_geometry::prelude::*;
    /// let knot_vecs = (KnotVector::bezier_knot(3), KnotVector::bezier_knot(2));
    /// let control_points = vec![
    ///     vec![Vector3::new(0.0, 0.0, 1.0), Vector3::new(0.5, -1.0, 2.0), Vector3::new(1.0, 0.0, 1.0)],
    ///     vec![Vector3::new(0.0, 1.0, 1.0), Vector3::new(0.5, 1.0, 1.0), Vector3::new(1.0, 1.0, 1.0)],
    ///     vec![Vector3::new(0.0, 2.0, 1.0), Vector3::new(0.5, 2.0, 3.0), Vector3::new(1.0, 2.0, 1.0)],
    ///     vec![Vector3::new(0.0, 3.0, 1.0), Vector3::new(0.5, 3.5, 2.0), Vector3::new(1.0, 3.0, 1.0)],
    /// ];
    /// let surface0 = NurbsSurface::new(BsplineSurface::new_unchecked(knot_vecs, control_points));
    /// let mut surface1 = surface0.clone();
    /// assert!(surface0.near_as_surface(&surface1));
    ///
    /// *surface1.control_point_mut(1, 1) = Vector3::new(0.5, 1.0, 0.9);
    /// assert!(!surface0.near_as_surface(&surface1));
    /// ```
    #[inline]
    #[must_use]
    pub fn near_as_surface(&self, other: &Self) -> bool {
        self.0
            .sub_near_as_surface(&other.0, 2, move |x, y| x.to_point().near(&y.to_point()))
    }
    /// Determines whether `self` and `other` is near in square order as the B-spline rational
    /// surfaces or not.
    ///
    /// Divides each knot domain into the number of degree equal parts,
    /// and check `|self(u, v) - other(u, v)| < TOLERANCE` for each end points `(u, v)`.
    /// # Examples
    /// ```
    /// use monstertruck_geometry::prelude::*;
    /// let eps = TOLERANCE;
    /// let knot_vecs = (KnotVector::bezier_knot(3), KnotVector::bezier_knot(2));
    /// let control_points = vec![
    ///     vec![Vector3::new(0.0, 0.0, 1.0), Vector3::new(0.5, -1.0, 2.0), Vector3::new(1.0, 0.0, 1.0)],
    ///     vec![Vector3::new(0.0, 1.0, 1.0), Vector3::new(0.5, 1.0, 1.0), Vector3::new(1.0, 1.0, 1.0)],
    ///     vec![Vector3::new(0.0, 2.0, 1.0), Vector3::new(0.5, 2.0, 3.0), Vector3::new(1.0, 2.0, 1.0)],
    ///     vec![Vector3::new(0.0, 3.0, 1.0), Vector3::new(0.5, 3.5, 2.0), Vector3::new(1.0, 3.0, 1.0)],
    /// ];
    /// let surface0 = NurbsSurface::new(BsplineSurface::new_unchecked(knot_vecs, control_points));
    /// let mut surface1 = surface0.clone();
    /// assert!(surface0.near_as_surface(&surface1));
    ///
    /// *surface1.control_point_mut(1, 1) = Vector3::new(0.5, 1.0, 1.0 - eps);
    /// assert!(surface0.near_as_surface(&surface1));
    /// assert!(!surface0.near2_as_surface(&surface1));
    /// ```
    #[inline]
    #[must_use]
    pub fn near2_as_surface(&self, other: &Self) -> bool {
        self.0
            .sub_near_as_surface(&other.0, 2, move |x, y| x.to_point().near2(&y.to_point()))
    }
}

impl<V: Homogeneous<Scalar = f64> + ControlPoint<f64, Diff = V> + Tolerance> NurbsSurface<V> {
    /// Adds a knot `x` of the first parameter `u`, and do not change `self` as a surface.
    #[inline]
    pub fn add_knot_u(&mut self, x: f64) -> &mut Self {
        self.0.add_knot_u(x);
        self
    }
    /// Adds a knot `x` of the first parameter `u`, and do not change `self` as a surface.
    #[inline]
    pub fn add_knot_v(&mut self, x: f64) -> &mut Self {
        self.0.add_knot_v(x);
        self
    }
    /// Removes the `knot_u` corresponding to the indice `idx`, and do not change `self` as a curve.
    /// If the knot cannot be removed, returns [`Error::CannotRemoveKnot`].
    /// # Errors
    /// Returns [`Error::CannotRemoveKnot`] when the underlying B-spline U knot cannot be removed.
    ///
    #[inline]
    pub fn try_remove_knot_u(&mut self, idx: usize) -> Result<&mut Self> {
        match self.0.try_remove_knot_u(idx) {
            Ok(_) => Ok(self),
            Err(error) => Err(error),
        }
    }
    /// Removes the `knot_u` corresponding to the indices `idx`, and do not change `self` as a curve.
    /// If cannot remove the knot, do not change `self` and return `self`.
    #[inline]
    pub fn remove_knot_u(&mut self, idx: usize) -> &mut Self {
        self.0.remove_knot_u(idx);
        self
    }
    /// Removes the `knot_u` corresponding to the indice `idx`, and do not change `self` as a curve.
    /// If the knot cannot be removed, returns [`Error::CannotRemoveKnot`].
    ///
    /// # Errors
    /// Returns [`Error::CannotRemoveKnot`] when the underlying B-spline V knot cannot be removed.
    #[inline]
    pub fn try_remove_knot_v(&mut self, idx: usize) -> Result<&mut Self> {
        match self.0.try_remove_knot_v(idx) {
            Ok(_) => Ok(self),
            Err(error) => Err(error),
        }
    }
    /// Removes the `knot_u` corresponding to the indices `idx`, and do not change `self` as a curve.
    /// If cannot remove the knot, do not change `self` and return `self`.
    #[inline]
    pub fn remove_knot_v(&mut self, idx: usize) -> &mut Self {
        self.0.remove_knot_v(idx);
        self
    }
    /// Elevates the udegree.
    #[inline]
    pub fn elevate_udegree(&mut self) -> &mut Self {
        self.0.elevate_udegree();
        self
    }
    /// Elevates the vdegree.
    #[inline]
    pub fn elevate_vdegree(&mut self) -> &mut Self {
        self.0.elevate_vdegree();
        self
    }
    /// Aligns the udegree with the same degrees.
    #[inline]
    pub fn syncro_uvdegrees(&mut self) -> &mut Self {
        self.0.syncro_uvdegrees();
        self
    }
    /// Makes the `knot_u` vector and the `knot_v` vector the same knot vector.
    #[inline]
    pub fn syncro_uvknots(&mut self) -> &mut Self {
        self.0.syncro_uvknots();
        self
    }

    /// Cuts the surface into two surfaces at the parameter `u`
    #[inline]
    #[must_use]
    pub fn cut_u(&mut self, u: f64) -> Self { Self::new(self.0.cut_u(u)) }

    /// Cuts the surface into two surfaces at the parameter `v`
    #[inline]
    #[must_use]
    pub fn cut_v(&mut self, v: f64) -> Self { Self::new(self.0.cut_v(v)) }

    /// Normalizes the knot vectors
    #[inline]
    pub fn knot_normalize(&mut self) -> &mut Self {
        self.0.knot_normalize();
        self
    }
    /// Translates the knot vectors.
    #[inline]
    pub fn knot_translate(&mut self, x: f64, y: f64) -> &mut Self {
        self.0.knot_translate(x, y);
        self
    }

    /// Removes knots in order from the back
    #[inline]
    pub fn optimize(&mut self) -> &mut Self {
        self.0.optimize();
        self
    }

    /// Get the boundary by four splitted curves.
    /// # Example
    /// ```
    /// use monstertruck_geometry::prelude::*;
    /// let knot_vecs = (KnotVector::bezier_knot(3), KnotVector::bezier_knot(2));
    /// let control_points = vec![
    ///     vec![Vector3::new(0.0, 0.0, 1.0), Vector3::new(0.5, -1.0, 2.0), Vector3::new(1.0, 0.0, 1.0)],
    ///     vec![Vector3::new(0.0, 1.0, 2.0), Vector3::new(0.5, 1.0, 3.0), Vector3::new(1.0, 1.0, 2.0)],
    ///     vec![Vector3::new(0.0, 2.0, 2.0), Vector3::new(0.5, 2.0, 3.0), Vector3::new(1.0, 2.0, 2.0)],
    ///     vec![Vector3::new(0.0, 3.0, 1.0), Vector3::new(0.5, 3.5, 2.0), Vector3::new(1.0, 3.0, 1.0)],
    /// ];
    /// let bspsurface = NurbsSurface::new(BsplineSurface::new_unchecked(knot_vecs, control_points));
    /// let curves = bspsurface.splitted_boundary();
    /// assert_eq!(
    ///     curves[0].control_points(),
    ///     &vec![
    ///         Vector3::new(0.0, 0.0, 1.0),
    ///         Vector3::new(0.0, 1.0, 2.0),
    ///         Vector3::new(0.0, 2.0, 2.0),
    ///         Vector3::new(0.0, 3.0, 1.0),
    ///     ],
    /// );
    /// assert_eq!(
    ///     curves[1].control_points(),
    ///     &vec![
    ///         Vector3::new(0.0, 3.0, 1.0),
    ///         Vector3::new(0.5, 3.5, 2.0),
    ///         Vector3::new(1.0, 3.0, 1.0),
    ///     ],
    /// );
    /// assert_eq!(
    ///     curves[2].control_points(),
    ///     &vec![
    ///         Vector3::new(1.0, 3.0, 1.0),
    ///         Vector3::new(1.0, 2.0, 2.0),
    ///         Vector3::new(1.0, 1.0, 2.0),
    ///         Vector3::new(1.0, 0.0, 1.0),
    ///     ],
    /// );
    /// assert_eq!(
    ///     curves[3].control_points(),
    ///     &vec![
    ///         Vector3::new(1.0, 0.0, 1.0),
    ///         Vector3::new(0.5, -1.0, 2.0),
    ///         Vector3::new(0.0, 0.0, 1.0),
    ///     ],
    /// );
    /// ```
    #[inline]
    pub fn splitted_boundary(&self) -> [NurbsCurve<V>; 4] {
        self.0.splitted_boundary().map(NurbsCurve::new)
    }

    /// Extracts the boundary of surface
    #[inline]
    pub fn boundary(&self) -> Option<NurbsCurve<V>> { self.0.boundary().map(NurbsCurve::new) }
}

impl<V: Homogeneous<Scalar = f64> + ControlPoint<f64, Diff = V>> NurbsSurface<V>
where V::Point: MetricSpace<Metric = f64> + Copy
{
    /// Separable-basis presearch: bit-identical to [`algo::surface::presearch`]
    /// for this surface, but computes each grid row's `u`-direction basis window
    /// and each grid column's `v`-direction basis window exactly once and reuses
    /// them across the whole `(division + 1) x (division + 1)` grid. The
    /// tensor-product B-spline basis is separable, so the generic scan's
    /// per-point re-evaluation of both bases repeats each window `division + 1`
    /// times; this hoists that redundant work out of the inner loop.
    ///
    /// Everything that affects the result is unchanged: the same grid nodes
    /// (`u`, `v` computed by the identical expressions), the same evaluation
    /// kernel ([`BsplineSurface::combine_basis`] + [`Homogeneous::to_point`]),
    /// the same `distance2` metric, and the same strict-`<` first-minimum
    /// tie-break. The returned `(u, v)` -- and therefore every downstream Newton
    /// seed -- is byte-for-byte identical to the generic presearch.
    fn presearch_separable(
        &self,
        point: V::Point,
        (urange, vrange): SurfaceParameterRange,
        division: usize,
    ) -> (f64, f64) {
        let bsp = &self.0;
        let (degree0, degree1) = bsp.degrees();
        let (knot_u, knot_v) = &bsp.knot_vecs;
        // The separable kernel reproduces only the non-degenerate branch of
        // `derivative_mn`; defer any zero-range surface to the generic scan so
        // the degenerate corner cases stay byte-for-byte unchanged.
        if knot_u.range_length().is_zero() || knot_v.range_length().is_zero() {
            return algo::surface::presearch(self, point, (urange, vrange), division);
        }
        // Spec 014 W3: the separable fast path scans the SAME grid as
        // `algo::surface::presearch`, so it charges the same nodes -- otherwise
        // the unit would read lower purely because a surface took the fast path.
        // (The degenerate early return above delegates and is charged there.)
        algo::surface::charge_presearch_nodes(algo::surface::presearch_nodes(division));
        let ((u0, u1), (v0, v1)) = (urange, vrange);
        // Every grid column's v-basis window, computed once.
        let vbasis: Vec<BasisWindow> = (0..=division)
            .map(|j| {
                let q = crate::usize_to_f64(j) / crate::usize_to_f64(division);
                let v = f64::mul_add(v1, q, v0 * (1.0 - q));
                knot_v.bspline_basis_functions_unchecked(degree1, 0, v)
            })
            .collect();
        let mut res = (0.0, 0.0);
        let mut min = f64::INFINITY;
        for i in 0..=division {
            let p = crate::usize_to_f64(i) / crate::usize_to_f64(division);
            let u = f64::mul_add(u1, p, u0 * (1.0 - p));
            // This grid row's u-basis window, computed once and reused for
            // every column below.
            let ubasis = knot_u.bspline_basis_functions_unchecked(degree0, 0, u);
            for (j, vb) in vbasis.iter().enumerate() {
                let q = crate::usize_to_f64(j) / crate::usize_to_f64(division);
                let v = f64::mul_add(v1, q, v0 * (1.0 - q));
                let dist = V::from_vec(bsp.combine_basis(&ubasis, vb))
                    .to_point()
                    .distance2(point);
                if dist < min {
                    min = dist;
                    res = (u, v);
                }
            }
        }
        res
    }
}

impl<V: Homogeneous<Scalar = f64> + ControlPoint<f64, Diff = V>>
    SearchNearestParameter<SurfaceParameter> for NurbsSurface<V>
where
    Self: ParametricSurface<Point = V::Point, Vector = <V::Point as EuclideanSpace>::Diff>,
    V::Point: EuclideanSpace<Scalar = f64> + MetricSpace<Metric = f64>,
    <V::Point as EuclideanSpace>::Diff: SearchNearestParameterVector<Point = V::Point>,
{
    type Point = V::Point;
    /// Searches the parameter `(u, v)` which minimize `|self(u, v) - point|` by Newton's method
    /// with initial guess `(u0, v0)`. If the repeated trial does not converge, then returns `None`.
    /// # Examples
    /// ```
    /// use monstertruck_geometry::prelude::*;
    /// let knot_vecs = (KnotVector::bezier_knot(3), KnotVector::bezier_knot(2));
    /// let control_points = vec![
    ///     vec![Vector3::new(0.0, 0.0, 1.0), Vector3::new(1.0, -2.0, 2.0), Vector3::new(1.0, 0.0, 1.0)],
    ///     vec![Vector3::new(0.0, 2.0, 2.0), Vector3::new(2.0, 4.0, 4.0), Vector3::new(2.0, 2.0, 2.0)],
    ///     vec![Vector3::new(0.0, 4.0, 2.0), Vector3::new(2.0, 8.0, 4.0), Vector3::new(2.0, 4.0, 2.0)],
    ///     vec![Vector3::new(0.0, 3.0, 1.0), Vector3::new(1.0, 7.0, 2.0), Vector3::new(1.0, 3.0, 1.0)],
    /// ];
    /// let surface = NurbsSurface::new(BsplineSurface::new_unchecked(knot_vecs, control_points));
    /// let pt = surface.subs(0.3, 0.7);
    /// let (u, v) = surface.search_nearest_parameter(pt, Some((0.5, 0.5)), 100).unwrap();
    /// assert!(u.near(&0.3) && v.near(&0.7));
    /// ```
    /// # Remarks
    /// It may converge to a local solution depending on the hint.
    /// cf. [`BsplineCurve::search_rational_nearest_parameter`](struct.BsplineCurve.html#method.search_rational_nearest_parameter)
    #[inline]
    fn search_nearest_parameter<H: Into<SearchParameterHint2D>>(
        &self,
        point: V::Point,
        hint: H,
        trials: usize,
    ) -> Option<(f64, f64)> {
        let hint = match hint.into() {
            SearchParameterHint2D::Parameter(x, y) => (x, y),
            SearchParameterHint2D::Range(range0, range1) => {
                self.presearch_separable(point, (range0, range1), PRESEARCH_DIVISION)
            }
            SearchParameterHint2D::None => {
                self.presearch_separable(point, self.range_tuple(), PRESEARCH_DIVISION)
            }
        };
        algo::surface::search_nearest_parameter(self, point, hint, trials)
    }
}

impl<V: Homogeneous<Scalar = f64>> NurbsSurface<V>
where V::Point: Bounded<Scalar = f64>
{
    /// Returns the bounding box including all control points.
    #[inline]
    #[must_use]
    pub fn roughly_bounding_box(&self) -> BoundingBox<V::Point> {
        self.0
            .control_points
            .iter()
            .flatten()
            .map(|pt| pt.to_point())
            .collect()
    }
}

impl<V: Clone> Invertible for NurbsSurface<V> {
    #[inline]
    fn invert(&mut self) { self.swap_axes(); }
    #[inline]
    fn inverse(&self) -> Self {
        let mut surface = self.clone();
        surface.swap_axes();
        surface
    }
}

impl<V: Homogeneous<Scalar = f64> + ControlPoint<f64, Diff = V>> ParametricSurface
    for NurbsSurface<V>
{
    type Point = V::Point;
    type Vector = <V::Point as EuclideanSpace>::Diff;
    #[inline]
    fn derivative_mn(&self, m: usize, n: usize, u: f64, v: f64) -> Self::Vector {
        if m < 7 && n < 7 {
            let mut ders = [[V::zero(); 8]; 8];
            (0..=m)
                .for_each(|i| (0..=n).for_each(|j| ders[i][j] = self.0.derivative_mn(i, j, u, v)));
            let ders = std::array::from_fn::<_, 8, _>(|i| &ders[i][..=n]);
            multi_rational_derivative(&ders[..=m])
                .unwrap_or_else(|_| Self::Vector::zero() * f64::NAN)
        } else {
            let ders = (0..=m)
                .map(|i| (0..=n).map(|j| self.0.derivative_mn(i, j, u, v)).collect())
                .collect::<Vec<Vec<_>>>();
            multi_rational_derivative(&ders).unwrap_or_else(|_| Self::Vector::zero() * f64::NAN)
        }
    }
    #[inline]
    fn evaluate(&self, u: f64, v: f64) -> V::Point { self.0.evaluate(u, v).to_point() }
    #[inline]
    fn derivative_u(&self, u: f64, v: f64) -> Self::Vector {
        rational_derivative(&[self.0.evaluate(u, v), self.0.derivative_u(u, v)])
    }
    #[inline]
    fn derivative_v(&self, u: f64, v: f64) -> <V::Point as EuclideanSpace>::Diff {
        rational_derivative(&[self.0.evaluate(u, v), self.0.derivative_v(u, v)])
    }
    #[inline]
    fn derivative_uu(&self, u: f64, v: f64) -> <V::Point as EuclideanSpace>::Diff {
        rational_derivative(&[
            self.0.evaluate(u, v),
            self.0.derivative_u(u, v),
            self.0.derivative_uu(u, v),
        ])
    }
    #[inline]
    fn derivative_uv(&self, u: f64, v: f64) -> <V::Point as EuclideanSpace>::Diff {
        multi_rational_derivative(&[
            [self.0.evaluate(u, v), self.0.derivative_v(u, v)],
            [self.0.derivative_u(u, v), self.0.derivative_uv(u, v)],
        ])
        .unwrap_or_else(|_| Self::Vector::zero() * f64::NAN)
    }
    #[inline]
    fn derivative_vv(&self, u: f64, v: f64) -> <V::Point as EuclideanSpace>::Diff {
        rational_derivative(&[
            self.0.evaluate(u, v),
            self.0.derivative_v(u, v),
            self.0.derivative_vv(u, v),
        ])
    }
    #[inline]
    fn parameter_range(&self) -> (ParameterRange, ParameterRange) { self.parameter_range() }
}

impl ParametricSurface3D for NurbsSurface<Vector4> {
    fn normal(&self, u: f64, v: f64) -> Vector3 {
        let pt = self.0.evaluate(u, v);
        let ud = self.0.derivative_u(u, v);
        let vd = self.0.derivative_v(u, v);
        let n = rational_derivative(&[pt, ud]).cross(rational_derivative(&[pt, vd]));
        if n.magnitude2() > TOLERANCE2 {
            return n.normalize();
        }
        // Degenerate point (e.g. pole of a revolution surface).
        // Sample nearby non-degenerate points and average their normals.
        let eps = 1.0e-4;
        let kvu = self.0.knot_vector_u();
        let kvv = self.0.knot_vector_v();
        let u_lo = kvu[0];
        let u_hi = kvu[kvu.len() - 1];
        let v_lo = kvv[0];
        let v_hi = kvv[kvv.len() - 1];
        let mut sum = Vector3::zero();
        let mut count = 0u32;
        for &(du, dv) in &[(eps, 0.0), (-eps, 0.0), (0.0, eps), (0.0, -eps)] {
            let u2 = (u + du).clamp(u_lo, u_hi);
            let v2 = (v + dv).clamp(v_lo, v_hi);
            if (u2 - u).abs() < eps * 0.5 && (v2 - v).abs() < eps * 0.5 {
                continue; // Clamped back to the original point.
            }
            let pt2 = self.0.evaluate(u2, v2);
            let ud2 = self.0.derivative_u(u2, v2);
            let vd2 = self.0.derivative_v(u2, v2);
            let n2 = rational_derivative(&[pt2, ud2]).cross(rational_derivative(&[pt2, vd2]));
            if n2.magnitude2() > TOLERANCE2 {
                sum += n2.normalize();
                count += 1;
            }
        }
        if count > 0 {
            (sum / f64::from(count)).normalize()
        } else {
            Vector3::unit_z()
        }
    }
}

impl<V: Homogeneous<Scalar = f64> + ControlPoint<f64, Diff = V>> ParameterDivision2D
    for NurbsSurface<V>
where V::Point: MetricSpace<Metric = f64> + HashGen<f64>
{
    #[inline]
    fn parameter_division(
        &self,
        range: ((f64, f64), (f64, f64)),
        tol: f64,
    ) -> (Vec<f64>, Vec<f64>) {
        algo::surface::parameter_division(self, range, tol)
    }
}

impl<V> BoundedSurface for NurbsSurface<V>
where Self: ParametricSurface
{
    fn range_tuple(&self) -> SurfaceParameterRange {
        let u_degree = self.0.knot_vecs.0.len() - self.0.control_points.len() - 1;
        let v_degree = self.0.knot_vecs.1.len() - self.0.control_points[0].len() - 1;
        (
            (
                self.0.knot_vecs.0[u_degree],
                self.0.knot_vecs.0[self.0.control_points.len()],
            ),
            (
                self.0.knot_vecs.1[v_degree],
                self.0.knot_vecs.1[self.0.control_points[0].len()],
            ),
        )
    }
}

// -- v2 scalar-generic impls ------------------------------------------------

use monstertruck_core::scalar::HasScalar;
use monstertruck_traits::v2;

impl<V> v2::ParametricSurface for NurbsSurface<V>
where V: HasScalar<Scalar = f64> + Homogeneous<Scalar = f64> + ControlPoint<f64, Diff = V>
{
    type Scalar = f64;
    type Point = V::Point;
    type Vector = <V::Point as EuclideanSpace>::Diff;

    #[inline]
    fn evaluate(&self, u: Self::Scalar, v: Self::Scalar) -> Self::Point {
        ParametricSurface::evaluate(self, u, v)
    }
    #[inline]
    fn derivative_u(&self, u: Self::Scalar, v: Self::Scalar) -> Self::Vector {
        ParametricSurface::derivative_u(self, u, v)
    }
    #[inline]
    fn derivative_v(&self, u: Self::Scalar, v: Self::Scalar) -> Self::Vector {
        ParametricSurface::derivative_v(self, u, v)
    }
    #[inline]
    fn derivative_uu(&self, u: Self::Scalar, v: Self::Scalar) -> Self::Vector {
        ParametricSurface::derivative_uu(self, u, v)
    }
    #[inline]
    fn derivative_uv(&self, u: Self::Scalar, v: Self::Scalar) -> Self::Vector {
        ParametricSurface::derivative_uv(self, u, v)
    }
    #[inline]
    fn derivative_vv(&self, u: Self::Scalar, v: Self::Scalar) -> Self::Vector {
        ParametricSurface::derivative_vv(self, u, v)
    }
    #[inline]
    fn period_u(&self) -> Option<Self::Scalar> { ParametricSurface::period_u(self) }
    #[inline]
    fn period_v(&self) -> Option<Self::Scalar> { ParametricSurface::period_v(self) }
}

impl<V> v2::BoundedSurface for NurbsSurface<V>
where
    V: HasScalar<Scalar = f64> + Homogeneous<Scalar = f64> + ControlPoint<f64, Diff = V>,
    Self: ParametricSurface + BoundedSurface,
{
    #[inline]
    fn range_tuple(&self) -> ((f64, f64), (f64, f64)) { BoundedSurface::range_tuple(self) }
}

impl v2::ParametricSurface3D for NurbsSurface<Vector4> {}

impl<V> v2::SearchNearestParameter<v2::SurfaceParameter<f64>> for NurbsSurface<V>
where
    V: HasScalar<Scalar = f64> + Homogeneous<Scalar = f64> + ControlPoint<f64, Diff = V>,
    Self: ParametricSurface<Point = V::Point, Vector = <V::Point as EuclideanSpace>::Diff>,
    V::Point: EuclideanSpace<Scalar = f64> + MetricSpace<Metric = f64>,
    <V::Point as EuclideanSpace>::Diff: SearchNearestParameterVector<Point = V::Point>,
{
    type Point = V::Point;
    #[inline]
    fn search_nearest_parameter<H: Into<v2::SearchParameterHint2D<f64>>>(
        &self,
        pt: V::Point,
        _: H,
        trials: usize,
    ) -> Option<(f64, f64)> {
        SearchNearestParameter::<D2>::search_nearest_parameter(self, pt, None, trials)
    }
}

impl<V> v2::SearchParameter<v2::SurfaceParameter<f64>> for NurbsSurface<V>
where
    V: HasScalar<Scalar = f64> + Homogeneous<Scalar = f64> + ControlPoint<f64, Diff = V>,
    V::Point: ControlPoint<f64, Diff = <V::Point as EuclideanSpace>::Diff>
        + MetricSpace<Metric = f64>
        + Tolerance,
    <V::Point as EuclideanSpace>::Diff: SearchParameterVector<Point = V::Point>,
{
    type Point = V::Point;
    #[inline]
    fn search_parameter<H: Into<v2::SearchParameterHint2D<f64>>>(
        &self,
        pt: V::Point,
        _: H,
        trials: usize,
    ) -> Option<(f64, f64)> {
        SearchParameter::<D2>::search_parameter(self, pt, None, trials)
    }
}

impl IncludeCurve<NurbsCurve<Vector3>> for NurbsSurface<Vector3> {
    #[inline]
    fn include(&self, curve: &NurbsCurve<Vector3>) -> bool {
        let pt = curve.subs(curve.knot_vector()[0]);
        let Some(mut hint) = self.search_parameter(pt, None, INCLUDE_CURVE_TRIALS) else {
            return false;
        };
        let knot_vector_u = self.knot_vector_u();
        let knot_vector_v = self.knot_vector_v();
        let degree = curve.degree() * 6;
        let (knots, _) = curve.knot_vector().to_single_multi();
        for i in 1..knots.len() {
            for j in 1..=degree {
                let p = crate::usize_to_f64(j) / crate::usize_to_f64(degree);
                let t = f64::mul_add(knots[i], p, knots[i - 1] * (1.0 - p));
                let pt = curve.subs(t);
                hint = match self.search_parameter(pt, Some(hint), INCLUDE_CURVE_TRIALS) {
                    Some(got) => got,
                    None => return false,
                };
                if !self.subs(hint.0, hint.1).near(&pt)
                    || hint.0 < knot_vector_u[0] - TOLERANCE
                    || hint.0 - knot_vector_u[0] > knot_vector_u.range_length() + TOLERANCE
                    || hint.1 < knot_vector_v[0] - TOLERANCE
                    || hint.1 - knot_vector_v[0] > knot_vector_v.range_length() + TOLERANCE
                {
                    return false;
                }
            }
        }
        true
    }
}

impl IncludeCurve<BsplineCurve<Point3>> for NurbsSurface<Vector4> {
    #[inline]
    fn include(&self, curve: &BsplineCurve<Point3>) -> bool {
        let pt = curve.front();
        let Some(mut hint) = self.search_parameter(pt, None, INCLUDE_CURVE_TRIALS) else {
            return false;
        };
        let knot_vector_u = self.knot_vector_u();
        let knot_vector_v = self.knot_vector_v();
        let degree = curve.degree() * 6;
        let (knots, _) = curve.knot_vector().to_single_multi();
        for i in 1..knots.len() {
            for j in 1..=degree {
                let p = crate::usize_to_f64(j) / crate::usize_to_f64(degree);
                let t = f64::mul_add(knots[i], p, knots[i - 1] * (1.0 - p));
                let pt = curve.subs(t);
                hint = match self.search_parameter(pt, Some(hint), INCLUDE_CURVE_TRIALS) {
                    Some(got) => got,
                    None => return false,
                };
                if !self.subs(hint.0, hint.1).near(&pt)
                    || hint.0 < knot_vector_u[0] - TOLERANCE
                    || hint.0 - knot_vector_u[0] > knot_vector_u.range_length() + TOLERANCE
                    || hint.1 < knot_vector_v[0] - TOLERANCE
                    || hint.1 - knot_vector_v[0] > knot_vector_v.range_length() + TOLERANCE
                {
                    return false;
                }
            }
        }
        true
    }
}

impl IncludeCurve<NurbsCurve<Vector4>> for NurbsSurface<Vector4> {
    #[inline]
    fn include(&self, curve: &NurbsCurve<Vector4>) -> bool {
        let pt = curve.front();
        let Some(mut hint) = self.search_parameter(pt, None, INCLUDE_CURVE_TRIALS) else {
            return false;
        };
        let knot_vector_u = self.knot_vector_u();
        let knot_vector_v = self.knot_vector_v();
        let degree = curve.degree() * 6;
        let (knots, _) = curve.knot_vector().to_single_multi();
        for i in 1..knots.len() {
            for j in 1..=degree {
                let p = crate::usize_to_f64(j) / crate::usize_to_f64(degree);
                let t = f64::mul_add(knots[i], p, knots[i - 1] * (1.0 - p));
                let pt = curve.subs(t);
                hint = match self.search_parameter(pt, Some(hint), INCLUDE_CURVE_TRIALS) {
                    Some(got) => got,
                    None => return false,
                };
                if !self.subs(hint.0, hint.1).near(&pt)
                    || hint.0 < knot_vector_u[0] - TOLERANCE
                    || hint.0 - knot_vector_u[0] > knot_vector_u.range_length() + TOLERANCE
                    || hint.1 < knot_vector_v[0] - TOLERANCE
                    || hint.1 - knot_vector_v[0] > knot_vector_v.range_length() + TOLERANCE
                {
                    return false;
                }
            }
        }
        true
    }
}

impl<M, V: Copy> Transformed<M> for NurbsSurface<V>
where M: Copy + std::ops::Mul<V, Output = V>
{
    #[inline]
    fn transform_by(&mut self, trans: M) {
        self.0
            .control_points
            .iter_mut()
            .flatten()
            .for_each(move |v| *v = trans * *v);
    }
}

impl<V> SearchParameter<SurfaceParameter> for NurbsSurface<V>
where
    V: Homogeneous<Scalar = f64> + ControlPoint<f64, Diff = V>,
    V::Point: ControlPoint<f64, Diff = <V::Point as EuclideanSpace>::Diff>
        + MetricSpace<Metric = f64>
        + Tolerance,
    <V::Point as EuclideanSpace>::Diff: SearchParameterVector<Point = V::Point>,
{
    type Point = V::Point;
    /// Search the parameter `(u, v)` such that `self.subs(u, v).rational_projection()` is near `pt`.
    /// If cannot find, then return `None`.
    /// # Examples
    /// ```
    /// use monstertruck_geometry::prelude::*;
    /// let knot_vec = KnotVector::uniform_knot(2, 2);
    /// let control_points = vec![
    ///     vec![Vector4::new(0.0, 0.0, 0.0, 1.0), Vector4::new(0.1, 0.0, 0.5, 0.4), Vector4::new(1.0, 0.0, 0.6, 2.0), Vector4::new(0.4, 0.0, 0.4, 0.4)],
    ///     vec![Vector4::new(0.0, 0.2, 0.2, 2.0), Vector4::new(0.24, 0.24, 0.1, 1.2), Vector4::new(2.4, 1.8, 2.4, 0.6), Vector4::new(1.4, 0.42, 0.98, 1.4)],
    ///     vec![Vector4::new(0.0, 1.5, 1.2, 3.0), Vector4::new(1.02, 2.04, 1.7, 3.4), Vector4::new(0.42, 0.28, 0.7, 0.7), Vector4::new(0.6, 0.3, 0.0, 0.6)],
    ///     vec![Vector4::new(0.0, 1.0, 1.0, 1.0), Vector4::new(0.2, 2.0, 2.0, 2.0), Vector4::new(0.85, 1.7, 0.85, 1.7), Vector4::new(1.0, 1.0, 0.3, 1.0)],
    /// ];
    /// let bspsurface = BsplineSurface::new_unchecked((knot_vec.clone(), knot_vec), control_points);
    /// let surface = NurbsSurface::new(bspsurface);
    ///
    /// let pt = surface.subs(0.3, 0.7);
    /// let (u, v) = surface.search_parameter(pt, Some((0.5, 0.5)), 100).unwrap();
    /// assert_near!(surface.subs(u, v), pt);
    /// ```
    fn search_parameter<H: Into<SearchParameterHint2D>>(
        &self,
        point: V::Point,
        hint: H,
        trials: usize,
    ) -> Option<(f64, f64)> {
        let hint = match hint.into() {
            SearchParameterHint2D::Parameter(x, y) => (x, y),
            SearchParameterHint2D::Range(range0, range1) => {
                self.presearch_separable(point, (range0, range1), PRESEARCH_DIVISION)
            }
            SearchParameterHint2D::None => {
                self.presearch_separable(point, self.range_tuple(), PRESEARCH_DIVISION)
            }
        };
        algo::surface::search_parameter(self, point, hint, trials)
    }
}

impl<V: Homogeneous<Scalar = f64>> From<BsplineSurface<V::Point>> for NurbsSurface<V> {
    fn from(bsp: BsplineSurface<V::Point>) -> Self {
        let control_points = bsp
            .control_points
            .into_iter()
            .map(|vec| vec.into_iter().map(|p| V::from_point(p)).collect())
            .collect();
        Self(BsplineSurface {
            knot_vecs: bsp.knot_vecs,
            control_points,
        })
    }
}

impl<V> From<NurbsSurface<V>> for BsplineSurface<V> {
    #[inline]
    fn from(nurbs: NurbsSurface<V>) -> Self { nurbs.0 }
}

#[test]
fn test_include2d() {
    let knot_vec = KnotVector::uniform_knot(2, 3);
    let control_points = vec![
        vec![
            Vector3::new(0.0, 0.0, 1.0),
            Vector3::new(0.05, 0.0, 0.5),
            Vector3::new(0.15, 0.0, 0.3),
            Vector3::new(1.0, 0.0, 1.0),
        ],
        vec![
            Vector3::new(0.0, 0.01, 0.1),
            Vector3::new(0.02, 0.02, 0.1),
            Vector3::new(0.16, 0.12, 0.4),
            Vector3::new(0.7, 0.21, 0.7),
        ],
        vec![
            Vector3::new(0.0, 0.02, 0.4),
            Vector3::new(0.15, 0.3, 0.5),
            Vector3::new(0.6, 0.4, 1.0),
            Vector3::new(0.4, 0.2, 0.4),
        ],
        vec![
            Vector3::new(0.0, 1.0, 1.0),
            Vector3::new(0.1, 1.0, 1.0),
            Vector3::new(0.25, 0.5, 0.5),
            Vector3::new(0.3, 0.3, 0.3),
        ],
    ];
    let surface = BsplineSurface::new_unchecked((knot_vec.clone(), knot_vec), control_points);
    let bnd_box = BoundingBox::from_iter(&[Vector2::new(0.2, 0.3), Vector2::new(0.8, 0.6)]);
    let mut curve = (surface.sectional_curve(bnd_box)).expect("section is valid");
    curve.control_points_mut().for_each(|pt| *pt *= 3.0);
    let surface = NurbsSurface::new(surface);
    let curve = NurbsCurve::new(curve);
    assert!(surface.include(&curve));
}

#[test]
fn test_include3d() {
    let knot_vec = KnotVector::bezier_knot(2);
    let control_points = vec![
        vec![
            Vector4::new(-1.0, -1.0, 2.0, 1.0),
            Vector4::new(-1.0, 0.0, 0.0, 1.0),
            Vector4::new(-1.0, 1.0, 2.0, 1.0),
        ],
        vec![
            Vector4::new(0.0, -1.0, 0.0, 1.0),
            Vector4::new(0.0, 0.0, -2.0, 1.0),
            Vector4::new(0.0, 1.0, 0.0, 1.0),
        ],
        vec![
            Vector4::new(1.0, -1.0, 2.0, 1.0),
            Vector4::new(1.0, 0.0, 0.0, 1.0),
            Vector4::new(1.0, 1.0, 2.0, 1.0),
        ],
    ];
    let surface = NurbsSurface::new(BsplineSurface::new_unchecked(
        (knot_vec.clone(), knot_vec),
        control_points,
    ));

    let knot_vec = KnotVector::from(vec![
        0.0, 0.0, 0.0, 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1.0, 1.0, 1.0,
    ]);
    let control_points = vec![
        // the vector of the indices of control points
        Vector4::new(0.0, -2.0, 2.0, 2.0),
        Vector4::new(1.0, -1.0, 1.0, 1.0),
        Vector4::new(1.0, 0.0, 1.0, 1.0),
        Vector4::new(1.0, 1.0, 1.0, 1.0),
        Vector4::new(0.0, 2.0, 2.0, 2.0),
        Vector4::new(-1.0, 1.0, 1.0, 1.0),
        Vector4::new(-1.0, 0.0, 1.0, 1.0),
        Vector4::new(-1.0, -1.0, 1.0, 1.0),
        Vector4::new(0.0, -2.0, 2.0, 2.0),
    ];
    let mut curve = NurbsCurve::new(BsplineCurve::new_unchecked(knot_vec, control_points));
    assert!(surface.include(&curve));
    *curve.control_point_mut(1) += Vector4::new(0.0, 0.0, 0.00001, 0.0);
    assert!(!surface.include(&curve));
}

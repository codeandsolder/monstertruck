#![expect(
    clippy::many_single_char_names,
    reason = "compact U/V derivative-order symbols mirror the tensor-product surface equations"
)]

use super::basis::{DIFF_EPS, basis_function_der};
use super::{
    Arc, BoundedSurface, BsplineCurve, BsplineSurface, Error, InnerSpace, Invertible, KnotVector,
    Matrix4, Origin, ParameterDivision2D, ParameterRange, ParametricSurface, ParametricSurface3D,
    Point3, Result, RwLock, SearchNearestParameter, SearchParameter, SearchParameterHint2D,
    SurfaceParameter, SurfaceParameterRange, Tmesh, TmeshConnectionType, TmeshControlPoint,
    TmeshDirection, Transformed, Vector3, algo,
};
use serde::{Deserialize, Serialize};
use std::panic::{AssertUnwindSafe, catch_unwind};

impl Tmesh<Point3> {
    /// Evaluates the analytical derivative d^(m+n)S / du^m dv^n at `(u, v)` using the quotient
    /// rule on the rational surface `S = N / W` where `N = sum(B_i * P_i)` and `W = sum(B_i)`.
    ///
    /// Supports analytical derivatives up to 2nd order in each parameter direction.
    /// Falls back to finite differences for higher orders.
    pub(super) fn try_analytical_der_mn(
        &self,
        m: usize,
        n: usize,
        u: f64,
        v: f64,
    ) -> Result<Vector3> {
        // Orders are at most 2 in the analytical path, so max binomial index is 2.
        const BINOM: [[f64; 3]; 3] = [[1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [1.0, 2.0, 1.0]];

        if self.knot_vectors.read().is_none() {
            self.generate_knot_vectors()?;
        }

        let num_points = self.control_points.len();
        let n_cols = n + 1;

        // Precompute all needed partial basis derivatives: B^(p,q)_i for p in 0..=m, q in 0..=n.
        // Flat layout with strided indexing [p][q][i] populated in [i][p][q] order -- the
        // transposed iteration prevents a clean iterator chain, so imperative indexing is used.
        let bd_stride = n_cols * num_points;
        let basis_derivs = {
            let knot_vectors = self.knot_vectors.read();
            let all_kvs = knot_vectors.as_ref().ok_or(Error::TmeshMalformedMesh)?;
            let mut basis_derivs = vec![0.0f64; (m + 1) * bd_stride];
            for (i, (s_kv, t_kv)) in all_kvs.iter().enumerate().take(num_points) {
                let s_slice = s_kv.as_slice();
                let t_slice = t_kv.as_slice();
                for p in 0..=m {
                    let s_val = basis_function_der(u, s_slice, p);
                    for q in 0..=n {
                        let t_val = basis_function_der(v, t_slice, q);
                        basis_derivs[p * bd_stride + q * num_points + i] = s_val * t_val;
                    }
                }
            }
            drop(knot_vectors);
            basis_derivs
        };

        // Compute partial derivatives of the numerator N and denominator W.
        // Flat layout: index [p * n_cols + q].
        let pq_size = (m + 1) * n_cols;
        let mut n_derivs = vec![Vector3::new(0.0, 0.0, 0.0); pq_size];
        let mut w_derivs = vec![0.0f64; pq_size];

        for p in 0..=m {
            for q in 0..=n {
                let bd_base = p * bd_stride + q * num_points;
                let mut nx = 0.0;
                let mut ny = 0.0;
                let mut nz = 0.0;
                let mut w = 0.0;
                for (i, cp) in self.control_points.iter().enumerate() {
                    let b = basis_derivs[bd_base + i];
                    let pt = *cp.read().point();
                    nx = f64::mul_add(b, pt.x, nx);
                    ny = f64::mul_add(b, pt.y, ny);
                    nz = f64::mul_add(b, pt.z, nz);
                    w += b;
                }
                let idx = p * n_cols + q;
                n_derivs[idx] = Vector3::new(nx, ny, nz);
                w_derivs[idx] = w;
            }
        }

        // Apply the general Leibniz rule for the derivative of a quotient S = N/W:
        // S^(m,n) = (1/W) * (N^(m,n) - sum_{(j,k) != (0,0)} C(m,j)*C(n,k) * W^(j,k) * S^(m-j,n-k)).
        let mut s_derivs = vec![Vector3::new(0.0, 0.0, 0.0); pq_size];
        let w0 = w_derivs[0];

        #[expect(
            clippy::needless_range_loop,
            reason = "tensor derivative assembly uses p and q as mathematical derivative orders"
        )]
        for p in 0..=m {
            for q in 0..=n {
                let idx = p * n_cols + q;
                let mut val = n_derivs[idx];
                for j in 0..=p {
                    for k in 0..=q {
                        if j == 0 && k == 0 {
                            continue;
                        }
                        let s_idx = (p - j) * n_cols + (q - k);
                        val -= s_derivs[s_idx]
                            * (BINOM[p][j] * BINOM[q][k] * w_derivs[j * n_cols + k]);
                    }
                }
                s_derivs[idx] = val / w0;
            }
        }

        Ok(s_derivs[m * n_cols + n])
    }

    pub(super) fn analytical_der_mn(&self, m: usize, n: usize, u: f64, v: f64) -> Vector3 {
        self.try_analytical_der_mn(m, n, u, v)
            .unwrap_or_else(|_| Vector3::new(f64::NAN, f64::NAN, f64::NAN))
    }

    fn try_gaussian_curvature(&self, u: f64, v: f64) -> Result<f64> {
        let su = self.try_analytical_der_mn(1, 0, u, v)?;
        let sv = self.try_analytical_der_mn(0, 1, u, v)?;
        let suu = self.try_analytical_der_mn(2, 0, u, v)?;
        let suv = self.try_analytical_der_mn(1, 1, u, v)?;
        let svv = self.try_analytical_der_mn(0, 2, u, v)?;

        let normal = su.cross(sv);
        let normal_len = normal.magnitude();
        let curvature = if normal_len.so_small() {
            0.0
        } else {
            let n = normal / normal_len;

            // First fundamental form coefficients.
            let cap_e = su.dot(su);
            let cap_f = su.dot(sv);
            let cap_g = sv.dot(sv);

            // Second fundamental form coefficients.
            let e = suu.dot(n);
            let f = suv.dot(n);
            let g = svv.dot(n);

            let cap_f_squared = cap_f * cap_f;
            let denom = f64::mul_add(cap_e, cap_g, -cap_f_squared);
            if denom.abs().so_small() {
                0.0
            } else {
                let f_squared = f * f;
                f64::mul_add(e, g, -f_squared) / denom
            }
        };
        Ok(curvature)
    }

    /// Computes the Gaussian curvature K at `(u, v)`.
    #[cfg(test)]
    pub(super) fn gaussian_curvature(&self, u: f64, v: f64) -> f64 {
        self.try_gaussian_curvature(u, v).unwrap_or(f64::NAN)
    }

    /// Adaptively refines the T-mesh by inserting knots where Gaussian curvature exceeds the threshold.
    ///
    /// Knot insertion uses `try_absolute_local_knot_insertion`, which requires the target
    /// coordinate to lie on an existing edge. For each high-curvature cell, the method
    /// inserts a midpoint on the nearest straddling edge in both the u and v directions.
    ///
    /// # Arguments
    /// * `curvature_threshold` - Minimum absolute Gaussian curvature to trigger refinement.
    /// * `max_iterations` - Maximum number of refinement passes.
    /// * `initial_samples` - Grid density in each direction for the first pass (doubles each iteration).
    ///
    /// # Errors
    /// Total number of control points inserted, or an error if refinement fails.
    pub fn adaptive_refine(
        &mut self,
        curvature_threshold: f64,
        max_iterations: usize,
        initial_samples: usize,
    ) -> Result<usize> {
        let mut total_insertions = 0usize;
        let mut samples = initial_samples;

        for _ in 0..max_iterations {
            // Collect all unique knot lines (s and t) from existing control points.
            let mut s_lines: Vec<f64> = self
                .control_points
                .iter()
                .map(|cp| cp.read().knot_coordinates().0)
                .collect();
            let mut t_lines: Vec<f64> = self
                .control_points
                .iter()
                .map(|cp| cp.read().knot_coordinates().1)
                .collect();
            s_lines.sort_by(f64::total_cmp);
            s_lines.dedup_by(|a, b| (*a - *b).so_small());
            t_lines.sort_by(f64::total_cmp);
            t_lines.dedup_by(|a, b| (*a - *b).so_small());

            // Sample curvature on a grid and collect cells that exceed the threshold.
            let step = 1.0 / crate::usize_to_f64(samples);
            let mut high_curvature_cells = Vec::new();
            for i in 0..samples {
                let u = (crate::usize_to_f64(i) + 0.5) * step;
                for j in 0..samples {
                    let v = (crate::usize_to_f64(j) + 0.5) * step;
                    if self.try_gaussian_curvature(u, v)?.abs() > curvature_threshold {
                        high_curvature_cells.push((u, v));
                    }
                }
            }

            if high_curvature_cells.is_empty() {
                break;
            }

            // For each high-curvature cell, find the straddling edge intervals and
            // insert midpoints on the existing knot lines.
            let mut targets: Vec<(f64, f64)> = high_curvature_cells
                .iter()
                .flat_map(|&(u, v)| {
                    // Nearest existing t-line -> insert at (u, t_val).
                    let on_t = t_lines
                        .iter()
                        .copied()
                        .min_by(|a, b| (a - v).abs().total_cmp(&(b - v).abs()))
                        .map(|t_val| (u, t_val));
                    // Nearest existing s-line -> insert at (s_val, v).
                    let on_s = s_lines
                        .iter()
                        .copied()
                        .min_by(|a, b| (a - u).abs().total_cmp(&(b - u).abs()))
                        .map(|s_val| (s_val, v));
                    on_t.into_iter().chain(on_s)
                })
                .collect();
            targets.sort_by(|a, b| a.0.total_cmp(&b.0).then(a.1.total_cmp(&b.1)));
            targets.dedup_by(|a, b| (a.0 - b.0).so_small() && (a.1 - b.1).so_small());

            let mut insertions = 0usize;
            for (u, v) in targets {
                // Clone before attempting insertion so a panic in LKI
                // doesn't corrupt the mesh.
                let backup = self.clone();
                let result = catch_unwind(AssertUnwindSafe(|| {
                    self.try_absolute_local_knot_insertion((u, v))
                }));
                match result {
                    Ok(Ok(_)) => insertions += 1,
                    Ok(Err(
                        Error::TmeshExistingControlPoint
                        | Error::TmeshConnectionNotFound
                        | Error::TmeshControlPointNotFound
                        | Error::TmeshKnotVectorsNotEqual,
                    ))
                    | Err(_) => {
                        // Restore from backup on any structural error or panic.
                        *self = backup;
                    }
                    Ok(Err(e)) => return Err(e),
                }
            }

            if insertions == 0 {
                break;
            }

            total_insertions += insertions;
            samples *= 2;
        }

        Ok(total_insertions)
    }
}

impl ParametricSurface for Tmesh<Point3> {
    type Point = Point3;
    type Vector = Vector3;

    fn evaluate(&self, u: f64, v: f64) -> Point3 {
        Self::try_evaluate(self, u, v).unwrap_or_else(|_| Point3::new(f64::NAN, f64::NAN, f64::NAN))
    }

    fn derivative_u(&self, u: f64, v: f64) -> Vector3 { self.derivative_mn(1, 0, u, v) }
    fn derivative_v(&self, u: f64, v: f64) -> Vector3 { self.derivative_mn(0, 1, u, v) }
    fn derivative_uu(&self, u: f64, v: f64) -> Vector3 { self.derivative_mn(2, 0, u, v) }
    fn derivative_uv(&self, u: f64, v: f64) -> Vector3 { self.derivative_mn(1, 1, u, v) }
    fn derivative_vv(&self, u: f64, v: f64) -> Vector3 { self.derivative_mn(0, 2, u, v) }

    fn derivative_mn(&self, m: usize, n: usize, u: f64, v: f64) -> Vector3 {
        if m == 0 && n == 0 {
            let p = <Self as ParametricSurface>::evaluate(self, u, v);
            return Vector3::new(p.x, p.y, p.z);
        }
        // Use analytical derivatives for orders up to 2.
        if m <= 2 && n <= 2 {
            return self.analytical_der_mn(m, n, u, v);
        }
        // Fall back to finite differences for higher orders.
        let h = DIFF_EPS;
        if m > 0 {
            let forward = self.derivative_mn(m - 1, n, u + h, v);
            let backward = self.derivative_mn(m - 1, n, u - h, v);
            (forward - backward) / (2.0 * h)
        } else {
            let forward = self.derivative_mn(m, n - 1, u, v + h);
            let backward = self.derivative_mn(m, n - 1, u, v - h);
            (forward - backward) / (2.0 * h)
        }
    }

    fn parameter_range(&self) -> (ParameterRange, ParameterRange) {
        use std::ops::Bound::Included;
        (
            (Included(0.0), Included(1.0)),
            (Included(0.0), Included(1.0)),
        )
    }
}

impl ParametricSurface3D for Tmesh<Point3> {}

impl BoundedSurface for Tmesh<Point3> {
    fn range_tuple(&self) -> SurfaceParameterRange { ((0.0, 1.0), (0.0, 1.0)) }
}

impl ParameterDivision2D for Tmesh<Point3> {
    fn parameter_division(
        &self,
        range: ((f64, f64), (f64, f64)),
        tol: f64,
    ) -> (Vec<f64>, Vec<f64>) {
        algo::surface::parameter_division(self, range, tol)
    }
}

impl Invertible for Tmesh<Point3> {
    fn invert(&mut self) {
        // Swap u and v by swapping Right<->Up and Left<->Down connections for every control point,
        // and swapping the (s, t) knot coordinates.
        for cp in &self.control_points {
            let mut w = cp.write();
            w.connections
                .swap(TmeshDirection::Up as usize, TmeshDirection::Right as usize);
            w.connections
                .swap(TmeshDirection::Down as usize, TmeshDirection::Left as usize);
            w.knot_coordinates = (w.knot_coordinates.1, w.knot_coordinates.0);
        }
        // Invalidate cached knot vectors.
        *self.knot_vectors.write() = None;
    }
}

impl Transformed<Matrix4> for Tmesh<Point3> {
    fn transform_by(&mut self, trans: Matrix4) {
        use monstertruck_core::cgmath64::Transform;
        for cp in &self.control_points {
            let mut w = cp.write();
            let p = *w.point();
            w.set_point(trans.transform_point(p));
        }
        // Invalidate cached knot vectors.
        *self.knot_vectors.write() = None;
    }
}

impl SearchParameter<SurfaceParameter> for Tmesh<Point3> {
    type Point = Point3;
    fn search_parameter<H: Into<SearchParameterHint2D>>(
        &self,
        point: Point3,
        hint: H,
        trials: usize,
    ) -> Option<(f64, f64)> {
        let hint = match hint.into() {
            SearchParameterHint2D::Parameter(u, v) => (u, v),
            SearchParameterHint2D::Range(x, y) => {
                algo::surface::presearch(self, point, (x, y), 100)
            }
            SearchParameterHint2D::None => {
                algo::surface::presearch(self, point, self.range_tuple(), 100)
            }
        };
        algo::surface::search_parameter(self, point, hint, trials)
    }
}

impl SearchNearestParameter<SurfaceParameter> for Tmesh<Point3> {
    type Point = Point3;
    fn search_nearest_parameter<H: Into<SearchParameterHint2D>>(
        &self,
        point: Point3,
        hint: H,
        trials: usize,
    ) -> Option<(f64, f64)> {
        let hint = match hint.into() {
            SearchParameterHint2D::Parameter(u, v) => (u, v),
            SearchParameterHint2D::Range(x, y) => {
                algo::surface::presearch(self, point, (x, y), 100)
            }
            SearchParameterHint2D::None => {
                algo::surface::presearch(self, point, self.range_tuple(), 100)
            }
        };
        algo::surface::search_nearest_parameter(self, point, hint, trials)
    }
}

/// Serializable representation of a single T-mesh connection.
type TmeshSerdeConnection = Option<(Option<usize>, f64)>;

/// Flat serialization helper for `Tmesh<P>`.
#[derive(Serialize, Deserialize)]
struct TmeshSerde<P> {
    /// Control point positions and knot coordinates.
    points: Vec<(P, (f64, f64))>,
    /// For each point, 4 connections (Up, Right, Down, Left).
    /// `None` = T-junction, `Some((None, ki))` = edge, `Some((Some(idx), ki))` = point connection.
    connections: Vec<[TmeshSerdeConnection; 4]>,
}

impl Serialize for Tmesh<Point3> {
    fn serialize<S: serde::Serializer>(
        &self,
        serializer: S,
    ) -> std::result::Result<S::Ok, S::Error> {
        let mut points = Vec::with_capacity(self.control_points.len());
        let mut connections = Vec::with_capacity(self.control_points.len());

        for cp in &self.control_points {
            let (position, knot_coordinates, point_connections) = {
                let control_point = cp.read();
                (
                    *control_point.point(),
                    control_point.knot_coordinates(),
                    [
                        control_point.get(TmeshDirection::Up).cloned(),
                        control_point.get(TmeshDirection::Right).cloned(),
                        control_point.get(TmeshDirection::Down).cloned(),
                        control_point.get(TmeshDirection::Left).cloned(),
                    ],
                )
            };
            points.push((position, knot_coordinates));

            let mut cons = [None; 4];
            for (index, connection) in point_connections.iter().enumerate() {
                cons[index] = match connection.as_ref() {
                    None => None,
                    Some((None, interval)) => Some((None, *interval)),
                    Some((Some(connected), interval)) => {
                        let connected_index = self
                            .control_points
                            .iter()
                            .position(|point| Arc::ptr_eq(point, connected))
                            .ok_or_else(|| {
                                serde::ser::Error::custom(
                                    "T-mesh connection points outside the serialized mesh",
                                )
                            })?;
                        Some((Some(connected_index), *interval))
                    }
                };
            }
            connections.push(cons);
        }

        TmeshSerde {
            points,
            connections,
        }
        .serialize(serializer)
    }
}

impl<'de> Deserialize<'de> for Tmesh<Point3> {
    fn deserialize<D: serde::Deserializer<'de>>(
        deserializer: D,
    ) -> std::result::Result<Self, D::Error> {
        let data = TmeshSerde::<Point3>::deserialize(deserializer)?;

        // Create control points with dummy edge conditions.
        let points: Vec<Arc<RwLock<TmeshControlPoint<Point3>>>> = data
            .points
            .iter()
            .map(|(p, _kc)| Arc::new(RwLock::new(TmeshControlPoint::new(*p, 1.0))))
            .collect();

        // Establish connections (same logic as Clone impl).
        for (point_index, cons) in data.connections.iter().enumerate() {
            for dir in TmeshDirection::iter() {
                let con = &cons[dir as usize];
                if let Some((maybe_idx, ki)) = con {
                    if let Some(con_index) = maybe_idx {
                        // Point connection -- skip if already established from the other side.
                        if points[point_index].read().con_type(dir) == TmeshConnectionType::Point {
                            continue;
                        }
                        points[point_index].write().remove_connection(dir).ok();
                        points[*con_index]
                            .write()
                            .remove_connection(dir.flip())
                            .ok();
                        TmeshControlPoint::connect(
                            &points[point_index],
                            &points[*con_index],
                            dir,
                            *ki,
                        )
                        .map_err(serde::de::Error::custom)?;
                    } else {
                        // Edge condition.
                        points[point_index]
                            .write()
                            .set_edge_con_weight(dir, *ki)
                            .ok();
                    }
                } else {
                    // T-junction.
                    points[point_index].write().remove_connection(dir).ok();
                }
            }
        }

        // Set knot coordinates.
        for (i, (_, kc)) in data.points.iter().enumerate() {
            points[i].write().knot_coordinates = *kc;
        }

        Ok(Self {
            control_points: points,
            knot_vectors: RwLock::new(None),
        })
    }
}

/// Computes the Greville abscissae for a knot vector of given degree.
/// These are the optimal parameter values for B-spline interpolation.
fn greville_abscissae(knots: &KnotVector, degree: usize) -> Vec<f64> {
    let n = knots.len() - degree - 1;
    (0..n)
        .map(|i| (1..=degree).map(|j| knots[i + j]).sum::<f64>() / crate::usize_to_f64(degree))
        .collect()
}

impl Tmesh<Point3> {
    /// Converts this T-spline surface to an approximate `BsplineSurface`.
    ///
    /// STEP (ISO 10303) has no T-spline entity, so T-spline surfaces must be
    /// decomposed into B-spline patches for export. This method evaluates the
    /// T-spline at Greville abscissae and uses tensor-product interpolation
    /// to find the B-spline control points.
    ///
    /// `division` controls the number of spans in each parametric direction.
    /// Higher values give better approximation at the cost of more control
    /// points: `division + 3` control points per direction.
    ///
    /// # Errors
    /// Propagates T-spline evaluation and B-spline interpolation errors.
    pub fn to_bspline_surface(&self, division: usize) -> Result<BsplineSurface<Point3>> {
        let u_knots = KnotVector::uniform_knot(3, division);
        let v_knots = KnotVector::uniform_knot(3, division);
        let n = division + 3;

        let u_grev = greville_abscissae(&u_knots, 3);
        let v_grev = greville_abscissae(&v_knots, 3);

        // Evaluate T-spline at the grid of Greville abscissae.
        let surface_points: Vec<Vec<Point3>> = u_grev
            .iter()
            .map(|&u| {
                v_grev
                    .iter()
                    .map(|&v| ParametricSurface::evaluate(self, u, v))
                    .collect()
            })
            .collect();

        // Tensor-product interpolation: first interpolate each row (v-direction).
        let row_curves: Vec<BsplineCurve<Point3>> = surface_points
            .iter()
            .map(|row| {
                let params: Vec<(f64, Point3)> =
                    v_grev.iter().copied().zip(row.iter().copied()).collect();
                BsplineCurve::try_interpolate(v_knots.clone(), params)
            })
            .collect::<Result<_>>()?;

        // Collect intermediate control points (one row per u-sample).
        let intermediate: Vec<Vec<Point3>> = row_curves
            .iter()
            .map(|c| c.control_points().clone())
            .collect();

        // Interpolate each column (u-direction) through the intermediate control points.
        // col_cps[j] contains the U-direction control points for V-index j.
        let col_cps: Vec<Vec<Point3>> = (0..n)
            .map(|j| {
                let params: Vec<(f64, Point3)> = u_grev
                    .iter()
                    .copied()
                    .zip(intermediate.iter().map(|row| row[j]))
                    .collect();
                BsplineCurve::try_interpolate(u_knots.clone(), params)
                    .map(|curve| curve.control_points().clone())
            })
            .collect::<Result<_>>()?;

        // Transpose from [V][U] to [U][V] for BsplineSurface.
        let control_points: Vec<Vec<Point3>> = (0..n)
            .map(|i| (0..n).map(|j| col_cps[j][i]).collect())
            .collect();

        BsplineSurface::new((u_knots, v_knots), control_points)
    }
}

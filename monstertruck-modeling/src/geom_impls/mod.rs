use crate::errors::Error;
use crate::{
    AffinePoint, BoundingBox, Connector, ExtrusionSurface, GeometricMapping, GeometryTransform,
    InnerProduct, Line, Matrix3, Matrix4, Origin, Plane, Point3, Processor, RevolutionSurface,
    SquareLinearMap, ToSameGeometry, Transformed, TrimmedCurve, UnitCircle, Vector3, Zero,
};
use itertools::Itertools;
use monstertruck_core::circumcenter;
use std::f64::consts::PI;

pub fn circle_arc_by_three_points(
    point0: Point3,
    point1: Point3,
    transit: Point3,
) -> std::result::Result<Processor<TrimmedCurve<UnitCircle<Point3>>, Matrix4>, Error> {
    let finite = |point: Point3| point.x.is_finite() && point.y.is_finite() && point.z.is_finite();
    if !finite(point0) || !finite(point1) || !finite(transit) {
        return Err(Error::DegenerateCircularArcPoints);
    }

    let (vec0, vec1) = (point0 - transit, point1 - transit);
    let axis_raw = vec1.cross(&vec0);
    let lengths_squared = vec0.norm_squared() * vec1.norm_squared();
    // Absolute cross-product thresholds scale as radius^4 and reject
    // perfectly valid small circles. Bound the dimensionless sin^2 of
    // the angle between the defining chords instead.
    if !lengths_squared.is_finite()
        || lengths_squared <= 0.0
        || !axis_raw.norm_squared().is_finite()
        || axis_raw.norm_squared() <= f64::EPSILON * lengths_squared
    {
        return Err(Error::DegenerateCircularArcPoints);
    }
    let origin = circumcenter(point0, point1, transit);
    if !finite(origin) {
        return Err(Error::DegenerateCircularArcPoints);
    }
    let axis = axis_raw.normalize();
    let angle = PI - monstertruck_core::angle_between(vec0, vec1);
    if !angle.is_finite() {
        return Err(Error::DegenerateCircularArcPoints);
    }
    Ok(circle_arc(point0, origin, axis, angle * 2.0))
}

/// Constructs the unique circular arc that starts at `point0` with the
/// given start `tangent`, ends at `point1`, and turns less than a full
/// revolution.
///
/// Returns:
/// - [`Error::DegenerateCircularArcTangent`] if `tangent` is zero
///   or near-zero.
/// - [`Error::CircularArcTangentParallelToChord`] if `tangent` is
///   parallel (or anti-parallel) to `point1 - point0` -- the plane of
///   the arc would be under-determined and the radius infinite.
pub fn try_circle_arc_by_start_tangent(
    point0: Point3,
    point1: Point3,
    tangent: Vector3,
) -> std::result::Result<Processor<TrimmedCurve<UnitCircle<Point3>>, Matrix4>, Error> {
    let chord = point1 - point0;
    if tangent.norm_squared().so_small() {
        return Err(Error::DegenerateCircularArcTangent);
    }
    let tangent = tangent.normalize();
    let axis_raw = tangent.cross(&chord);
    if axis_raw.norm_squared().so_small() {
        return Err(Error::CircularArcTangentParallelToChord);
    }
    let axis = axis_raw.normalize();
    let to_origin = axis.cross(&tangent);
    let radius = chord.dot(&chord) / (2.0 * chord.dot(&to_origin));
    let origin = point0 + radius * to_origin;
    let (vec0, vec1) = (point0 - origin, point1 - origin);
    let mut angle = f64::atan2(axis.dot(&vec0.cross(&vec1)), vec0.dot(&vec1));
    if angle <= 0.0 {
        angle = 2.0f64.mul_add(PI, angle);
    }
    Ok(circle_arc(point0, origin, axis, angle))
}

pub fn circle_arc(
    point: Point3,
    origin: Point3,
    axis: Vector3,
    angle: f64,
) -> Processor<TrimmedCurve<UnitCircle<Point3>>, Matrix4> {
    let origin = origin + (axis.dot(&(point - origin))) * axis;
    let diag = point - origin;
    let axis_trsf = Matrix4::from_columns(&[
        diag.push(0.0),
        axis.cross(&diag).push(0.0),
        axis.push(0.0),
        origin.to_homogeneous(),
    ]);
    let unit_arc = TrimmedCurve::new(UnitCircle::new(), (0.0, angle));
    Processor::with_transform(unit_arc, axis_trsf)
}

fn closed_polyline_orientation<'a>(pts: impl IntoIterator<Item = &'a Vec<Point3>>) -> bool {
    pts.into_iter()
        .flat_map(|vec| vec.iter().circular_tuple_windows())
        .map(|(p0, p1)| (p1[0] + p0[0]) * (p1[1] - p0[1]))
        .sum::<f64>()
        >= 0.0
}

fn take_one_axis_by_normal(n: Vector3) -> Vector3 {
    let a = n.map(f64::abs);
    if a.x > a.z || a.y > a.z {
        Vector3::new(-n.y, n.x, 0.0).normalize()
    } else {
        Vector3::new(-n.z, 0.0, n.x).normalize()
    }
}

pub fn attach_plane(mut pts: Vec<Vec<Point3>>) -> Option<Plane> {
    let point_count = pts.iter().map(Vec::len).sum::<usize>();
    if point_count == 0 {
        return None;
    }
    let center = pts
        .iter()
        .flatten()
        .fold(Point3::origin(), |sum, pt| sum + pt.position_vector())
        / crate::usize_to_f64(point_count);
    let normal = pts
        .iter()
        .flat_map(|vec| vec.iter().circular_tuple_windows())
        .fold(Vector3::zero(), |sum, (p0, p1)| {
            sum + (p0 - center).cross(&(p1 - center))
        });
    let n = if normal.so_small() {
        return None;
    } else {
        normal.normalize()
    };
    let a = take_one_axis_by_normal(n);
    let mat: Matrix4 = Matrix3::from_columns(&[a, n.cross(&a), n]).into();
    let inverse = mat.try_inverse()?;
    pts.iter_mut()
        .flatten()
        .for_each(|pt| *pt = inverse.apply_point(*pt));
    let bnd_box: BoundingBox<Point3> = pts.iter().flatten().collect();
    let diag = bnd_box.diagonal();
    if !diag[2].so_small() {
        return None;
    }
    let (max, min) = if closed_polyline_orientation(&pts) {
        (bnd_box.max(), bnd_box.min())
    } else {
        (bnd_box.min(), bnd_box.max())
    };
    let plane = Plane::new(
        Point3::new(min[0], min[1], min[2]),
        Point3::new(max[0], min[1], min[2]),
        Point3::new(min[0], max[1], min[2]),
    )
    .transformed(mat);
    Some(plane)
}

#[cfg(test)]
mod test_geom_impl;

impl<T: Clone> GeometricMapping<T> for () {
    #[inline]
    fn mapping(self) -> impl Fn(&T) -> T { Clone::clone }
}
impl<T: Transformed<Self>> GeometricMapping<T> for Matrix4 {
    #[inline]
    fn mapping(self) -> impl Fn(&T) -> T { move |t| t.transformed(self) }
}
impl<T> GeometricMapping<T> for fn(&T) -> T {
    #[inline]
    fn mapping(self) -> impl Fn(&T) -> T { self }
}

impl<T, H> Connector<T, H> for fn(&T, &T) -> H {
    #[inline]
    fn connector(self) -> impl Fn(&T, &T) -> H { self }
}

#[derive(Debug, Clone, Copy)]
pub struct LineConnector;

impl<C> Connector<Point3, C> for LineConnector
where Line<Point3>: ToSameGeometry<C>
{
    fn connector(self) -> impl Fn(&Point3, &Point3) -> C { |p, q| Line(*p, *q).to_same_geometry() }
}

#[derive(Debug, Clone, Copy)]
pub struct ExtrudeConnector {
    pub vector: Vector3,
}

impl<C, S> Connector<C, S> for ExtrudeConnector
where
    C: Clone,
    ExtrusionSurface<C, Vector3>: ToSameGeometry<S>,
{
    fn connector(self) -> impl Fn(&C, &C) -> S {
        move |curve0, _| {
            ExtrusionSurface::by_extrusion(curve0.clone(), self.vector).to_same_geometry()
        }
    }
}

#[derive(Clone, Copy, Debug)]
pub struct ArcConnector {
    pub origin: Point3,
    pub axis: Vector3,
    pub angle: f64,
}

impl<C> Connector<Point3, C> for ArcConnector
where Processor<TrimmedCurve<UnitCircle<Point3>>, Matrix4>: ToSameGeometry<C>
{
    fn connector(self) -> impl Fn(&Point3, &Point3) -> C {
        let Self {
            origin,
            axis,
            angle,
        } = self;
        move |p, _| circle_arc(*p, origin, axis, angle).to_same_geometry()
    }
}

#[derive(Clone, Copy, Debug)]
pub struct RevoluteConnector {
    pub origin: Point3,
    pub axis: Vector3,
}

impl<C, S> Connector<C, S> for RevoluteConnector
where
    C: Clone,
    RevolutionSurface<C>: ToSameGeometry<S>,
{
    fn connector(self) -> impl Fn(&C, &C) -> S {
        let Self { origin, axis } = self;
        move |curve, _| {
            RevolutionSurface::by_revolution(curve.clone(), origin, axis).to_same_geometry()
        }
    }
}

//! Points, directions, vectors and axis placements.
//!
//! The coordinate-carrying entities every other family refers to: a
//! `cartesian_point` is a position, a `direction` a unit vector, and the
//! `axis*_placement` entities the local frames that curves, surfaces and
//! transformations are expressed in.

use super::{
    Angle, Array, Deserialize, EuclideanSpace, Holder, InnerSpace, Matrix, Matrix3, Matrix4,
    Origin, Point2, Point3, Result, Serialize, StepConvertingError, Table, Tolerance, Vector2,
    Vector3, Zero,
};

/// `cartesian_point`
#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize, Holder)]
#[holder(table = Table)]
#[holder(field = cartesian_point)]
#[holder(generate_deserialize)]
pub struct CartesianPoint {
    pub label: String,
    pub coordinates: Vec<f64>,
}
impl From<&CartesianPoint> for Point2 {
    fn from(pt: &CartesianPoint) -> Self {
        let pt = &pt.coordinates;
        match pt.len() {
            0 => Self::origin(),
            1 => Self::new(pt[0], 0.0),
            _ => Self::new(pt[0], pt[1]),
        }
    }
}
impl From<&CartesianPoint> for Point3 {
    fn from(pt: &CartesianPoint) -> Self {
        let pt = &pt.coordinates;
        match pt.len() {
            0 => Self::origin(),
            1 => Self::new(pt[0], 0.0, 0.0),
            2 => Self::new(pt[0], pt[1], 0.0),
            _ => Self::new(pt[0], pt[1], pt[2]),
        }
    }
}

/// `direction`
#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize, Holder)]
#[holder(table = Table)]
#[holder(field = direction)]
#[holder(generate_deserialize)]
pub struct Direction {
    pub label: String,
    pub direction_ratios: Vec<f64>,
}
impl From<&Direction> for Vector2 {
    fn from(dir: &Direction) -> Self {
        let dir = &dir.direction_ratios;
        match dir.len() {
            0 => Self::zero(),
            1 => Self::new(dir[0], 0.0),
            _ => Self::new(dir[0], dir[1]),
        }
    }
}
impl From<&Direction> for Vector3 {
    fn from(dir: &Direction) -> Self {
        let dir = &dir.direction_ratios;
        match dir.len() {
            0 => Self::zero(),
            1 => Self::new(dir[0], 0.0, 0.0),
            2 => Self::new(dir[0], dir[1], 0.0),
            _ => Self::new(dir[0], dir[1], dir[2]),
        }
    }
}

/// `vector`
#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize, Holder)]
#[holder(table = Table)]
#[holder(field = vector)]
#[holder(generate_deserialize)]
pub struct Vector {
    pub label: String,
    #[holder(use_place_holder)]
    pub orientation: Direction,
    pub magnitude: f64,
}
impl From<&Vector> for Vector2 {
    fn from(vec: &Vector) -> Self { Self::from(&vec.orientation) * vec.magnitude }
}
impl From<&Vector> for Vector3 {
    fn from(vec: &Vector) -> Self { Self::from(&vec.orientation) * vec.magnitude }
}

/// `placement`
#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize, Holder)]
#[holder(table = Table)]
#[holder(field = placement)]
#[holder(generate_deserialize)]
pub struct Placement {
    pub label: String,
    #[holder(use_place_holder)]
    pub location: CartesianPoint,
}
impl From<&Placement> for Point2 {
    fn from(p: &Placement) -> Self { Self::from(&p.location) }
}
impl From<&Placement> for Point3 {
    fn from(p: &Placement) -> Self { Self::from(&p.location) }
}

/// `axis1_placement`
#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize, Holder)]
#[holder(table = Table)]
#[holder(field = axis1_placement)]
#[holder(generate_deserialize)]
pub struct Axis1Placement {
    pub label: String,
    #[holder(use_place_holder)]
    pub location: CartesianPoint,
    #[holder(use_place_holder)]
    pub direction: Option<Direction>,
}

impl Axis1Placement {
    pub fn direction(&self) -> Vector3 {
        self.direction
            .as_ref()
            .map_or_else(Vector3::unit_z, Vector3::from)
    }
}

/// `axis2_placement`
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize, Holder)]
#[holder(table = Table)]
#[holder(generate_deserialize)]
pub enum Axis2Placement {
    #[holder(use_place_holder)]
    Axis2Placement2d(Axis2Placement2d),
    #[holder(use_place_holder)]
    Axis2Placement3d(Axis2Placement3d),
}

impl TryFrom<&Axis2Placement> for Matrix3 {
    type Error = StepConvertingError;
    fn try_from(axis: &Axis2Placement) -> Result<Self, StepConvertingError> {
        use Axis2Placement::{Axis2Placement2d, Axis2Placement3d};
        match axis {
            Axis2Placement2d(axis) => Ok(Self::from(axis)),
            Axis2Placement3d(_) => Err("This is not a 2D axis placement.".into()),
        }
    }
}
impl TryFrom<&Axis2Placement> for Matrix4 {
    type Error = StepConvertingError;
    fn try_from(axis: &Axis2Placement) -> Result<Self, StepConvertingError> {
        use Axis2Placement::{Axis2Placement2d, Axis2Placement3d};
        match axis {
            Axis2Placement2d(_) => Err("This is not a 3D axis placement.".into()),
            Axis2Placement3d(axis) => Ok(Self::from(axis)),
        }
    }
}

/// `axis2_placement_2d`
#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize, Holder)]
#[holder(table = Table)]
#[holder(field = axis2_placement_2d)]
#[holder(generate_deserialize)]
pub struct Axis2Placement2d {
    pub label: String,
    #[holder(use_place_holder)]
    pub location: CartesianPoint,
    #[holder(use_place_holder)]
    pub ref_direction: Option<Direction>,
}

impl From<&Axis2Placement2d> for Matrix3 {
    fn from(axis: &Axis2Placement2d) -> Self {
        let z = Point2::from(&axis.location);
        let x = match &axis.ref_direction {
            Some(axis) => Vector2::from(axis),
            None => Vector2::unit_x(),
        };
        let y = Vector2::new(-x.y, x.x);
        Self::from_cols(x.extend(0.0), y.extend(0.0), z.to_vec().extend(1.0))
    }
}

/// `axis2_placement_3d`
#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize, Holder)]
#[holder(table = Table)]
#[holder(field = axis2_placement_3d)]
#[holder(generate_deserialize)]
pub struct Axis2Placement3d {
    pub label: String,
    #[holder(use_place_holder)]
    pub location: CartesianPoint,
    #[holder(use_place_holder)]
    pub axis: Option<Direction>,
    #[holder(use_place_holder)]
    pub ref_direction: Option<Direction>,
}

impl From<&Axis2Placement3d> for Matrix4 {
    fn from(axis: &Axis2Placement3d) -> Self {
        let w = Point3::from(&axis.location);
        let z = match &axis.axis {
            Some(axis) => Vector3::from(axis),
            None => Vector3::unit_z(),
        };
        // Pick a fallback reference direction that is not parallel to `z`.
        let fallback = if z.near(&Vector3::unit_x()) {
            Vector3::unit_y()
        } else {
            Vector3::unit_x()
        };
        let x = match &axis.ref_direction {
            Some(axis) => Vector3::from(axis),
            None => fallback,
        };
        // Gram-Schmidt: remove the `z` component of `x`. ISO 10303 permits
        // `ref_direction` to be parallel to `axis`; in that case the projected
        // vector is zero and `normalize()` would produce `NaN`, which later
        // panics in `nonpositive_tolerance!` during meshing. Fall back to an
        // arbitrary direction orthogonal to `z` when the projection degenerates.
        let projected = x - x.dot(z) * z;
        let x = if projected.magnitude2().so_small() {
            (fallback - fallback.dot(z) * z).normalize()
        } else {
            projected.normalize()
        };
        let y = z.cross(x);
        Self::from_cols(
            x.extend(0.0),
            y.extend(0.0),
            z.extend(0.0),
            w.to_vec().extend(1.0),
        )
    }
}

//! Shared three-dimensional vector operations for the independent oracle.

pub type Point = [f64; 3];

pub fn dot(a: Point, b: Point) -> f64 { a.into_iter().zip(b).map(|(x, y)| x * y).sum() }

pub fn cross(a: Point, b: Point) -> Point {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

pub fn sub(a: Point, b: Point) -> Point { std::array::from_fn(|i| a[i] - b[i]) }
pub fn scale(a: Point, factor: f64) -> Point { a.map(|x| x * factor) }
pub fn norm(a: Point) -> f64 { dot(a, a).sqrt() }

pub fn unit(a: Point, minimum_length: f64) -> Option<Point> {
    let length = norm(a);
    (length.is_finite() && length > minimum_length).then(|| a.map(|x| x / length))
}

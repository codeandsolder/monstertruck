use crate::{AbstractShape, Edge, Face, IntoWasm, Vertex, Wire, wasm_bindgen};
use monstertruck_modeling::*;
use wasm_bindgen::JsValue;

fn js_error(message: impl AsRef<str>) -> JsValue { JsValue::from_str(message.as_ref()) }

fn point3(values: &[f64], name: &str) -> std::result::Result<Point3, JsValue> {
    let [x, y, z] = values else {
        return Err(js_error(format!("{name} must contain exactly 3 numbers")));
    };
    Ok(Point3::new(*x, *y, *z))
}

fn vector3(values: &[f64], name: &str) -> std::result::Result<Vector3, JsValue> {
    let [x, y, z] = values else {
        return Err(js_error(format!("{name} must contain exactly 3 numbers")));
    };
    Ok(Vector3::new(*x, *y, *z))
}

/// Creates and returns a vertex by a three dimensional point.
#[wasm_bindgen]
pub fn vertex(x: f64, y: f64, z: f64) -> Vertex { builder::vertex(Point3::new(x, y, z)).into() }

/// Returns a line from `vertex0` to `vertex1`.
#[wasm_bindgen]
pub fn line(vertex0: &Vertex, vertex1: &Vertex) -> Edge { builder::line(vertex0, vertex1).into() }

/// Returns a circle arc from `vertex0` to `vertex1` via `transit`.
///
/// # Errors
/// Returns a JavaScript error for malformed coordinates or degenerate arc geometry.
#[wasm_bindgen]
pub fn circle_arc(
    vertex0: &Vertex,
    vertex1: &Vertex,
    transit: &[f64],
) -> std::result::Result<Edge, JsValue> {
    let transit = point3(transit, "transit")?;
    builder::circle_arc(vertex0, vertex1, transit)
        .map(Into::into)
        .map_err(|error| js_error(error.to_string()))
}

/// Returns a Bezier curve from `vertex0` to `vertex1` with inter control points `inter_points`.
///
/// # Errors
/// Returns a JavaScript error unless `inter_points` contains complete XYZ triples.
#[wasm_bindgen]
pub fn bezier(
    vertex0: &Vertex,
    vertex1: &Vertex,
    inter_points: &[f64],
) -> std::result::Result<Edge, JsValue> {
    let (points, remainder) = inter_points.as_chunks::<3>();
    if !remainder.is_empty() {
        return Err(js_error(
            "inter_points length must be a multiple of 3 coordinates",
        ));
    }
    let inter_points = points
        .iter()
        .map(|[x, y, z]| Point3::new(*x, *y, *z))
        .collect();
    Ok(builder::bezier(vertex0, vertex1, inter_points).into())
}

/// Returns a homotopic face from `edge0` to `edge1`.
#[wasm_bindgen]
pub fn homotopy(edge0: &Edge, edge1: &Edge) -> Face { builder::homotopy(edge0, edge1).into() }

/// Try attaching a plane whose boundary is `wire`.
#[wasm_bindgen]
pub fn try_attach_plane(wire: &Wire) -> Option<Face> {
    builder::try_attach_plane(&[wire.as_ref().clone()])
        .map(Into::into)
        .map_err(|error| {
            gloo::console::error!(format!("{error}"));
        })
        .ok()
}

macro_rules! transform_if_chain {
    ($shape: expr, $function: expr, ($($arg: expr),*), $exception: expr, $member: ident) => {
        if let Some(entity) = AbstractShape::$member($shape) {
            Ok($function(entity.as_ref(), $($arg),*).into_wasm().upcast())
        } else {
            $exception
        }
    };
    ($shape: expr, $function: expr, ($($arg: expr),*), $exception: expr, $member: ident, $($a: ident),*) => {
        if let Some(entity) = AbstractShape::$member($shape) {
            Ok($function(entity.as_ref(), $($arg),*).into_wasm().upcast())
        } else {
            transform_if_chain!($shape, $function, ($($arg),*), $exception, $($a),*)
        }
    };
}

macro_rules! try_transform_if_chain {
    ($shape: expr, $function: expr, ($($arg: expr),*), $exception: expr, $member: ident) => {
        if let Some(entity) = AbstractShape::$member($shape) {
            $function(entity.as_ref(), $($arg),*)
                .map(|value| value.into_wasm().upcast())
                .map_err(|error| js_error(error.to_string()))
        } else {
            $exception
        }
    };
    ($shape: expr, $function: expr, ($($arg: expr),*), $exception: expr, $member: ident, $($a: ident),*) => {
        if let Some(entity) = AbstractShape::$member($shape) {
            $function(entity.as_ref(), $($arg),*)
                .map(|value| value.into_wasm().upcast())
                .map_err(|error| js_error(error.to_string()))
        } else {
            try_transform_if_chain!($shape, $function, ($($arg),*), $exception, $($a),*)
        }
    };
}

macro_rules! derive_all_shape {
    ($shape: expr, $function: expr, ($($arg: expr),*)) => {
        transform_if_chain!(
            $shape,
            $function,
            ($($arg),*),
            Err(js_error("unsupported abstract shape")),
            as_vertex,
            as_edge,
            as_wire,
            as_face,
            as_shell,
            as_solid
        )
    };
}

macro_rules! derive_all_sweepable {
    ($shape: expr, $function: expr, ($($arg: expr),*)) => {
        transform_if_chain!(
            $shape,
            $function,
            ($($arg),*),
            Err(js_error("sweep is only supported for Vertex, Edge, Wire, and Face")),
            as_vertex,
            as_edge,
            as_wire,
            as_face
        )
    };
}

macro_rules! derive_all_fallible_sweepable {
    ($shape: expr, $function: expr, ($($arg: expr),*)) => {
        try_transform_if_chain!(
            $shape,
            $function,
            ($($arg),*),
            Err(js_error("sweep is only supported for Vertex, Edge, Wire, and Face")),
            as_vertex,
            as_edge,
            as_wire,
            as_face
        )
    };
}

/// Returns a translated vertex, edge, wire, face, shell or solid.
///
/// # Errors
/// Returns a JavaScript error for malformed vectors or unsupported shape variants.
#[wasm_bindgen]
pub fn translated(
    shape: &AbstractShape,
    vector: &[f64],
) -> std::result::Result<AbstractShape, JsValue> {
    let vector = vector3(vector, "vector")?;
    derive_all_shape!(shape, builder::translated, (vector))
}

/// Returns a rotated vertex, edge, wire, face, shell or solid.
///
/// # Errors
/// Returns a JavaScript error for malformed transform inputs or unsupported shape variants.
#[wasm_bindgen]
pub fn rotated(
    shape: &AbstractShape,
    origin: &[f64],
    axis: &[f64],
    angle: f64,
) -> std::result::Result<AbstractShape, JsValue> {
    let origin = point3(origin, "origin")?;
    let axis = vector3(axis, "axis")?;
    derive_all_shape!(shape, builder::rotated, (origin, axis, Rad(angle)))
}

/// Returns a scaled vertex, edge, wire, face, shell or solid.
///
/// # Errors
/// Returns a JavaScript error for malformed scale/origin inputs or unsupported shape variants.
#[wasm_bindgen]
pub fn scaled(
    shape: &AbstractShape,
    origin: &[f64],
    scalars: &[f64],
) -> std::result::Result<AbstractShape, JsValue> {
    let origin = point3(origin, "origin")?;
    let scale = match scalars {
        [scalar] => Vector3::new(*scalar, *scalar, *scalar),
        [x, y, z] => Vector3::new(*x, *y, *z),
        _ => return Err(js_error("scalars must contain either 1 or 3 numbers")),
    };
    derive_all_shape!(shape, builder::scaled, (origin, scale))
}

/// Sweeps a vertex, an edge, a wire, or a face by a vector.
///
/// # Errors
/// Returns a JavaScript error for malformed vectors or unsupported shape variants.
#[wasm_bindgen]
pub fn extrude(
    shape: &AbstractShape,
    vector: &[f64],
) -> std::result::Result<AbstractShape, JsValue> {
    let vector = vector3(vector, "vector")?;
    derive_all_sweepable!(shape, builder::extrude, (vector))
}

/// Sweeps a vertex, an edge, a wire, or a face by a rotation.
///
/// # Errors
/// Returns a JavaScript error for malformed/degenerate revolution parameters or unsupported shapes.
#[wasm_bindgen]
pub fn revolve(
    shape: &AbstractShape,
    origin: &[f64],
    axis: &[f64],
    angle: f64,
    division: usize,
) -> std::result::Result<AbstractShape, JsValue> {
    let origin = point3(origin, "origin")?;
    let axis = vector3(axis, "axis")?;
    derive_all_fallible_sweepable!(
        shape,
        builder::revolve,
        (
            origin,
            axis,
            builder::SweepAngle::Partial(Rad(angle)),
            division
        )
    )
}

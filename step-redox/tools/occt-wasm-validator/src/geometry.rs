use anyhow::{Result, ensure};
use occt_wasm::{BoundingBox, OcctKernel, ShapeHandle, StepImport};
use serde::Serialize;
use std::path::Path;

#[derive(Debug, Serialize)]
struct Bounds {
    bbox: [f64; 6],
    bbox_center: [f64; 3],
    bbox_extent: [f64; 3],
}

impl Bounds {
    fn new(bounds: &BoundingBox) -> Result<Self> {
        let low = [bounds.min.x, bounds.min.y, bounds.min.z];
        let high = [bounds.max.x, bounds.max.y, bounds.max.z];
        ensure!(
            low.into_iter()
                .zip(high)
                .all(|(lo, hi)| lo.is_finite() && hi.is_finite() && lo <= hi),
            "invalid BRep bounds"
        );
        Ok(Self {
            bbox: [low[0], low[1], low[2], high[0], high[1], high[2]],
            bbox_center: std::array::from_fn(|i| (low[i] + high[i]) * 0.5),
            bbox_extent: std::array::from_fn(|i| high[i] - low[i]),
        })
    }
}

#[derive(Debug, Serialize)]
struct Topology {
    faces: u32,
    edges: u32,
    vertices: u32,
}

impl Topology {
    fn measure(kernel: &mut OcctKernel, shape: ShapeHandle) -> Result<Self> {
        Ok(Self {
            faces: u32::try_from(kernel.sub_shape_occurrence_count(shape, "face")?)?,
            edges: u32::try_from(kernel.sub_shape_occurrence_count(shape, "edge")?)?,
            vertices: u32::try_from(kernel.sub_shape_occurrence_count(shape, "vertex")?)?,
        })
    }
}

#[derive(Debug, Serialize)]
struct Solid {
    volume: f64,
    area: f64,
    bbox_center: [f64; 3],
    bbox_extent: [f64; 3],
    #[serde(flatten)]
    topology: Topology,
}

#[derive(Debug, Serialize)]
pub struct GeometryReport<'a> {
    path: &'a Path,
    bytes: usize,
    roots: u32,
    transferred: u32,
    brep_valid: bool,
    volume: f64,
    area: f64,
    edge_length: f64,
    center_of_mass: [f64; 3],
    inertia: [[f64; 3]; 3],
    #[serde(flatten)]
    bounds: Bounds,
    #[serde(flatten)]
    topology: Topology,
    solids: u32,
    shells: u32,
    solid_volume_sum: f64,
    solid_area_sum: f64,
    solid_inventory: Vec<Solid>,
}

pub fn analyze<'a>(
    kernel: &mut OcctKernel,
    path: &'a Path,
    bytes: usize,
    import: StepImport,
) -> Result<GeometryReport<'a>> {
    let shape = import.shape;
    let props = kernel.get_b_rep_properties(shape)?;
    let bounds = Bounds::new(&kernel.get_bounding_box_loose(shape, true)?)?;
    let topology = Topology::measure(kernel, shape)?;
    let solids = kernel.get_sub_shape_occurrences(shape, "solid")?;
    let mut inventory = Vec::with_capacity(solids.len());
    for solid in solids {
        let props = kernel.get_b_rep_properties(solid)?;
        let bounds = Bounds::new(&kernel.get_bounding_box_loose(solid, true)?)?;
        inventory.push(Solid {
            volume: props.volume,
            area: props.area,
            bbox_center: bounds.bbox_center,
            bbox_extent: bounds.bbox_extent,
            topology: Topology::measure(kernel, solid)?,
        });
        kernel.release(solid)?;
    }
    let solid_count = u32::try_from(kernel.sub_shape_occurrence_count(shape, "solid")?)?;
    ensure!(
        usize::try_from(solid_count)? == inventory.len(),
        "solid occurrence mismatch"
    );
    Ok(GeometryReport {
        path,
        bytes,
        roots: import.roots,
        transferred: import.transferred,
        brep_valid: kernel.is_valid(shape)?,
        volume: props.volume,
        area: props.area,
        edge_length: props.edge_length,
        center_of_mass: props.center_of_mass,
        inertia: props.inertia,
        bounds,
        topology,
        solids: solid_count,
        shells: u32::try_from(kernel.sub_shape_occurrence_count(shape, "shell")?)?,
        solid_volume_sum: inventory.iter().map(|solid| solid.volume).sum(),
        solid_area_sum: inventory.iter().map(|solid| solid.area).sum(),
        solid_inventory: inventory,
    })
}

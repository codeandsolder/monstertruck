use crate::PolygonMesh;
use derive_more::{Deref, DerefMut, From, Into};
use gloo::console;
use monstertruck_io::step::load::Table as StepTable;
use monstertruck_io::step::load::step_geometry::{Curve3D, Point3, Surface};
use monstertruck_meshing::tessellation::{MeshedShape, RobustMeshableShape};
use monstertruck_topology::compress::{CompressedShell, CompressedSolid};
use wasm_bindgen::prelude::wasm_bindgen;

/// STEP parse table.
#[derive(Clone, Debug, Deref, DerefMut, From, Into)]
#[wasm_bindgen]
pub struct Table(StepTable);

#[derive(Clone, Debug)]
enum SubShapeFromStep {
    Shell(CompressedShell<Point3, Curve3D, Surface>),
    #[expect(
        dead_code,
        reason = "WASM STEP input enum preserves solid support even when the current conversion path does not construct it"
    )]
    Solid(CompressedSolid<Point3, Curve3D, Surface>),
}

/// Shell and solid parsed from STEP.
#[derive(Clone, Debug, From, Into)]
#[wasm_bindgen]
pub struct ShapeFromStep(SubShapeFromStep);

#[wasm_bindgen]
impl ShapeFromStep {
    /// Meshes a shape from STEP.
    #[must_use]
    pub fn to_polygon(&self, tol: f64) -> PolygonMesh {
        match &self.0 {
            SubShapeFromStep::Shell(x) => x.robust_triangulation(tol).to_polygon().into(),
            SubShapeFromStep::Solid(x) => x.robust_triangulation(tol).to_polygon().into(),
        }
    }
}

#[wasm_bindgen]
impl Table {
    /// Reads a STEP file.
    #[must_use]
    pub fn from_step(step_str: &str) -> Option<Self> {
        StepTable::from_step(step_str)
            .map(Table)
            .map_err(|e| {
                console::error!(format!("{e}"));
            })
            .ok()
    }
    /// Gets shell indices.
    #[must_use]
    pub fn shell_indices(&self) -> Vec<u64> { self.0.shell.keys().copied().collect() }
    /// Gets a shape from an entity index.
    #[must_use]
    pub fn shape(&self, idx: u64) -> Option<ShapeFromStep> {
        let stepshell = self.shell.get(&idx)?;
        let shell = self
            .to_compressed_shell(stepshell)
            .map_err(|e| {
                console::error!(format!("{e}"));
            })
            .ok()?;
        Some(SubShapeFromStep::Shell(shell).into())
    }
}

//! Tessellate a shape and output an obj file.
//!
//! ```bash
//! usage: tessellate-shape <input json file> <output obj file>
//! ```
//!
//! The default `<output file>` is output.obj.

use anyhow::{Context as _, bail};
use monstertruck_mesh::{KERNEL_EPSILON, obj};
use monstertruck_meshing::{
    analyzers::Topology,
    filters::OptimizingFilter,
    tessellation::{MeshableShape, MeshedShape},
};
use monstertruck_modeling::{
    Point3,
    geometry::{Curve, Surface},
};
use monstertruck_topology::compress::{CompressedShell, CompressedSolid};

type CShell = CompressedShell<Point3, Curve, Surface>;
type CSolid = CompressedSolid<Point3, Curve, Surface>;

fn main() -> anyhow::Result<()> {
    let mut args = std::env::args();
    let _program = args.next();
    let Some(input) = args.next() else {
        bail!(
            "usage: tessellate-shape <input json file> <output obj file>\nThe default <output file> is output.obj."
        );
    };
    let output = args.next().unwrap_or_else(|| "output.obj".to_owned());
    let json =
        std::fs::read_to_string(&input).with_context(|| format!("failed to read {input}"))?;
    let mut poly = if let Ok(solid) = serde_json::from_str::<CSolid>(&json) {
        solid.triangulation(0.005).to_polygon()
    } else if let Ok(shell) = serde_json::from_str::<CShell>(&json) {
        shell.triangulation(0.005).to_polygon()
    } else {
        bail!("input is neither a compressed solid nor a compressed shell");
    };
    poly.put_together_same_attrs(KERNEL_EPSILON * 10.0)
        .remove_degenerate_faces()
        .remove_unused_attrs();
    println!("polygon shell condition: {:?}", poly.shell_condition());
    let file =
        std::fs::File::create(&output).with_context(|| format!("failed to create {output}"))?;
    obj::write(&poly, file)?;
    Ok(())
}

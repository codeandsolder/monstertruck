//! Export selected imported faces and their shared tessellation samples.
//! Retains the complete shell's edge sampling, rather than tessellating
//! each selected face in isolation.

use anyhow::{Context, Result};
use clap::Parser;
use monstertruck_io::step::load::Table;
use monstertruck_meshing::tessellation::{MeshableShape, stitch_compressed_shell};
use monstertruck_modeling::KERNEL_EPSILON;
use serde_json::json;
use std::{
    collections::BTreeSet,
    fs::File,
    io::{BufWriter, Write},
    path::PathBuf,
};

#[derive(Parser)]
struct Args {
    input: PathBuf,
    output: PathBuf,
    #[arg(long, default_value_t = 0.005)]
    tolerance: f64,
    /// STEP shell entity ID; default is the shell with the most faces.
    #[arg(long)]
    shell: Option<u64>,
    /// Zero-based face index within the imported shell.
    #[arg(long, required = true)]
    face: Vec<usize>,
}

fn main() -> Result<()> {
    let args = Args::parse();
    anyhow::ensure!(
        args.tolerance.is_finite() && args.tolerance >= KERNEL_EPSILON,
        "invalid tessellation tolerance"
    );
    let table = Table::from_step(&std::fs::read_to_string(&args.input)?)?;
    let mut selected = None;
    for (&id, shell) in &table.shell {
        if args.shell.is_some_and(|wanted| wanted != id) {
            continue;
        }
        let shell = table.to_compressed_trimmed_shell(shell)?;
        if selected.as_ref().is_none_or(
            |(_, previous): &(
                u64,
                monstertruck_io::step::load::convert::StepCompressedTrimmedShell,
            )| shell.faces.len() > previous.faces.len(),
        ) {
            selected = Some((id, shell));
        }
    }
    let (id, shell) = selected.context("no matching STEP shell")?;
    anyhow::ensure!(
        args.face.iter().all(|&i| i < shell.faces.len()),
        "face index out of range"
    );
    let meshed = shell.triangulation(args.tolerance);
    let mut faces = Vec::new();
    for &i in &args.face {
        let face = &shell.faces[i];
        let native_edges = face
            .boundaries
            .iter()
            .flatten()
            .map(|e| e.index)
            .collect::<BTreeSet<_>>();
        let edges = native_edges.into_iter().map(|j|{
            json!({"index":j,"native":shell.edges[j],
                "samples":meshed.edges[j].curve.0,
                "vertices":[shell.vertices[shell.edges[j].vertices.0],shell.vertices[shell.edges[j].vertices.1]]})
        }).collect::<Vec<_>>();
        let polygon = meshed.faces[i].surface.as_ref().map(|mesh|{
            let polygons=mesh.face_iter().map(|p|p.iter().map(|v|v.pos).collect::<Vec<_>>()).collect::<Vec<_>>();
            json!({"positions":mesh.positions(),"uvs":mesh.uv_coords(),"polygons":polygons,"orientation":meshed.faces[i].orientation})
        });
        faces.push(json!({"index":i,"native":face,"edges":edges,"mesh":polygon}));
    }
    let error = stitch_compressed_shell(&meshed, KERNEL_EPSILON)
        .err()
        .map(|e| e.to_string());
    let mut writer = BufWriter::new(File::create_new(&args.output)?);
    serde_json::to_writer(
        &mut writer,
        &json!({"shell":id,"face_count":shell.faces.len(),"error":error,"faces":faces}),
    )?;
    writer.flush()?;
    Ok(())
}

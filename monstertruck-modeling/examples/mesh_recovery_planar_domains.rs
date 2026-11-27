//! Certify selected planar native trim domains from a compressed shell.
use anyhow::{Context, Result};
use monstertruck_modeling::{
    Curve, Point3, Shell, Surface,
    mesh_recovery::{PlanarFanOptions, certify_projected_planar_domain, simplify_planar_fan_shell},
};
use monstertruck_topology::compress::CompressedShell;
use serde_json::json;
use std::{
    fs::File,
    io::{BufReader, BufWriter, Write},
    time::Instant,
};

fn main() -> Result<()> {
    let mut args = std::env::args_os().skip(1);
    let path = args
        .next()
        .context("expected compressed shell JSON or --simplify")?;
    if path == "--simplify" {
        let input = args.next().context("expected input shell JSON")?;
        let output = args.next().context("expected output shell JSON")?;
        let limit = args
            .next()
            .context("expected additional deviation budget")?
            .to_str()
            .context("non-UTF8 budget")?
            .parse::<f64>()?;
        anyhow::ensure!(args.next().is_none(), "--simplify INPUT OUTPUT BUDGET");
        let compressed: CompressedShell<Point3, Curve, Surface> =
            serde_json::from_reader(BufReader::new(File::open(input)?))?;
        let shell = Shell::extract(compressed)?;
        let start = Instant::now();
        let result = simplify_planar_fan_shell(&shell, PlanarFanOptions::new(limit))
            .context("shell simplification failed validation")?;
        let mut writer = BufWriter::new(File::create_new(output)?);
        serde_json::to_writer(&mut writer, &result.shell.compress())?;
        writer.flush()?;
        println!(
            "{}",
            json!({"faces_before":shell.len(),"faces_after":result.shell.len(),"merged_groups":result.merged_groups,"rejected_groups":result.rejected_groups,"deviation_bound":result.deviation_upper_bound,"seam_bound":result.seam_deviation_upper_bound,"original_to_simplified":result.original_to_simplified,"seconds":start.elapsed().as_secs_f64()})
        );
        return Ok(());
    }
    let input: CompressedShell<Point3, Curve, Surface> =
        serde_json::from_reader(BufReader::new(File::open(path)?))?;
    let shell = Shell::extract(input)?;
    let indices = args
        .map(|arg| {
            arg.to_str()
                .context("non-UTF8 face index")?
                .parse::<usize>()
                .context("face index")
        })
        .collect::<Result<Vec<_>>>()?;
    let start = Instant::now();
    let mut results = Vec::new();
    for (index, face) in shell.iter().enumerate() {
        if !indices.is_empty() && !indices.contains(&index) {
            continue;
        }
        let Surface::Plane(plane) = face.surface() else {
            continue;
        };
        let normal = if face.orientation() {
            plane.normal()
        } else {
            -plane.normal()
        };
        let result = certify_projected_planar_domain(
            &face.boundaries(),
            plane.origin(),
            normal,
            monstertruck_modeling::KERNEL_EPSILON,
        );
        results.push(json!({"face":index,"wires":face.boundaries().len(),"accepted":result.is_some(),"certificate":result.map(|c|json!({"outer_wire":c.outer_wire,"pieces":c.pieces,"seam_bound":c.seam_displacement_upper_bound}))}));
    }
    println!(
        "{}",
        json!({"seconds":start.elapsed().as_secs_f64(),"attempted":results.len(),"accepted":results.iter().filter(|r|r["accepted"]==true).count(),"results":results})
    );
    Ok(())
}

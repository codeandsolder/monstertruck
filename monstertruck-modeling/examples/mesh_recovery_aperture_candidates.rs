//! Enumerate explicit source-only topology hypotheses for collapsed apertures.
use anyhow::{Context, Result, ensure};
use monstertruck_mesh::obj;
use monstertruck_modeling::mesh_recovery::{
    AnalyticRecoveryOptions, AperturePolygonOrigin, NativeRecoveryOptions,
    collapsed_aperture_candidates,
};
use serde_json::json;
use std::{
    fs::{self, File},
    io::{BufReader, BufWriter},
    path::PathBuf,
};

fn main() -> Result<()> {
    let mut args = std::env::args_os().skip(1);
    let source = PathBuf::from(args.next().context("source OBJ")?);
    let budget = args
        .next()
        .context("source chord budget")?
        .into_string()
        .map_err(|_| anyhow::anyhow!("invalid budget"))?
        .parse::<f64>()?;
    let out = PathBuf::from(args.next().context("output directory")?);
    let recover = match args.next() {
        None => false,
        Some(flag) if flag == "--recover" => true,
        Some(_) => anyhow::bail!("unknown option"),
    };
    ensure!(
        args.next().is_none(),
        "usage: mesh_recovery_aperture_candidates input.obj budget output_directory [--recover]"
    );
    let mesh = obj::read(BufReader::new(File::open(&source)?))?;
    let candidates = collapsed_aperture_candidates(&mesh, budget)?;
    fs::create_dir_all(&out)?;
    let mut reports = Vec::new();
    for (i, candidate) in candidates.into_iter().enumerate() {
        let path = out.join(format!("candidate-{i}.obj"));
        obj::write(candidate.mesh(), BufWriter::new(File::create(&path)?))?;
        let e = candidate.evidence();
        let provenance = candidate
            .polygon_origins()
            .iter()
            .map(|origin| match origin {
                AperturePolygonOrigin::Source(id) => json!({"source":id}),
                AperturePolygonOrigin::ReconstructedCap => json!("cap"),
                AperturePolygonOrigin::ReconstructedWall => json!("wall"),
            })
            .collect::<Vec<_>>();
        let mut report = json!({"path":path,"depth":e.depth,"deviation_upper_bound":e.deviation_upper_bound,
            "cap_polygons":e.cap_polygons,"wall_polygons":e.wall_polygons,"floor_polygons":e.floor_polygons,
            "floor_ring":e.floor_ring,"contact_edges":e.contact_edges,"polygon_origins":provenance});
        if recover {
            let options = NativeRecoveryOptions {
                analytic: AnalyticRecoveryOptions {
                    source_chord_budget: Some(budget),
                    ..Default::default()
                },
                hybrid: None,
            };
            let remaining = candidate.remaining_recovery_options(&options)?;
            report["remaining_chord_budget"] = json!(remaining.analytic.source_chord_budget);
            let start = std::time::Instant::now();
            let native = candidate.recover_native(&options, |_| {})?;
            report["native_seconds"] = json!(start.elapsed().as_secs_f64());
            report["native_faces"] = json!(native.shell.as_ref().map(|s| s.len()));
            report["native_components"] = json!(native.component_count);
            if let Some(shell) = native.shell {
                let shell_path = out.join(format!("candidate-{i}.shell.json"));
                serde_json::to_writer(
                    BufWriter::new(File::create_new(&shell_path)?),
                    &shell.compress(),
                )?;
                report["native_shell"] = json!(shell_path);
            }
        }
        reports.push(report);
    }
    println!(
        "{}",
        json!({"source":source,"source_chord_budget":budget,"candidates":reports})
    );
    Ok(())
}

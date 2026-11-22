use super::{
    IdentityKind, ManifestArgs, ScanOutcome, Source, load_manifest, parent, scan_source, write_line,
};
use crate::{detect_solid_surface_signatures_bytes, solid_revolutions::SolidSurfaceSignature};
use anyhow::{Result, ensure};
use clap::Args;
use rayon::prelude::*;
use serde::Serialize;
use std::{
    fs::{self, File},
    io::{BufWriter, Write},
    path::PathBuf,
};

#[derive(Debug, Args)]
pub struct CensusArgs {
    #[command(flatten)]
    manifest: ManifestArgs,
    #[arg(long)]
    out: PathBuf,
    #[arg(long)]
    face_count: Option<usize>,
    #[arg(long = "contains-support")]
    supports: Vec<String>,
}

#[derive(Serialize)]
struct CensusRecord<'a> {
    source: &'a Source,
    solids: Vec<SolidSurfaceSignature>,
    #[serde(skip_serializing_if = "Option::is_none")]
    error: Option<String>,
}

fn scan_supports<'a>(args: &CensusArgs, source: &'a Source) -> CensusRecord<'a> {
    let signatures = fs::read(&source.path)
        .map_err(anyhow::Error::from)
        .and_then(|bytes| detect_solid_surface_signatures_bytes(&bytes));
    let (solids, error) = match signatures {
        Ok(mut solids) => {
            solids.retain(|solid| {
                args.face_count
                    .is_none_or(|count| solid.face_count == count)
                    && args.supports.iter().all(|name| {
                        solid
                            .support_counts
                            .get(name)
                            .is_some_and(|count| *count > 0)
                    })
            });
            (solids, None)
        }
        Err(error) => (Vec::new(), Some(format!("{error:#}"))),
    };
    CensusRecord {
        source,
        solids,
        error,
    }
}

#[derive(Serialize)]
struct CensusSummary {
    sources: usize,
    matching_solids: usize,
    failures: usize,
}

pub(super) fn census(args: &CensusArgs) -> Result<bool> {
    let manifest = load_manifest(&args.manifest.manifest)?;
    let pool = args.manifest.parallel.pool()?;
    fs::create_dir_all(parent(&args.out))?;
    let mut writer = BufWriter::new(File::create_new(&args.out)?);
    let mut summary = CensusSummary {
        sources: manifest.sources.len(),
        matching_solids: 0,
        failures: 0,
    };
    for chunk in manifest.sources.chunks(args.manifest.parallel.jobs.get()) {
        let records = pool.install(|| {
            chunk
                .par_iter()
                .map(|source| scan_supports(args, source))
                .collect::<Vec<_>>()
        });
        for record in records {
            summary.failures += usize::from(record.error.is_some());
            summary.matching_solids += record.solids.len();
            if !record.solids.is_empty() || record.error.is_some() {
                write_line(&mut writer, &record)?;
            }
        }
    }
    writer.flush()?;
    println!("{}", serde_json::to_string_pretty(&summary)?);
    Ok(summary.failures == 0)
}

#[derive(Debug, Args)]
pub struct GapArgs {
    #[command(flatten)]
    manifest: ManifestArgs,
    /// Start at this zero-based manifest row.
    #[arg(long, default_value_t = 0)]
    start: usize,
}

pub(super) fn find_gap(args: &GapArgs) -> Result<bool> {
    let manifest = load_manifest(&args.manifest.manifest)?;
    ensure!(
        args.start <= manifest.sources.len(),
        "start outside manifest"
    );
    let pool = args.manifest.parallel.pool()?;
    let mut index = args.start;
    for chunk in manifest.sources[args.start..].chunks(args.manifest.parallel.jobs.get()) {
        let reports = pool.install(|| {
            chunk
                .par_iter()
                .map(|source| scan_source(source, None))
                .collect::<Vec<_>>()
        });
        for record in reports {
            let failure = match &record.outcome {
                ScanOutcome::Error { .. } => true,
                ScanOutcome::Scanned { report, .. } => {
                    report.skipped_solids > 0
                        || report
                            .solids
                            .iter()
                            .any(|solid| solid.identity_kind != IdentityKind::Canonical)
                }
            };
            if failure {
                println!(
                    "{}",
                    serde_json::to_string_pretty(
                        &serde_json::json!({"index":index,"record":record})
                    )?
                );
                return Ok(false);
            }
            index += 1;
        }
    }
    println!(
        "{}",
        serde_json::json!({"all_canonical":true,"models":manifest.sources.len(),"checked_from":args.start})
    );
    Ok(true)
}

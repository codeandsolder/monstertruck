use super::{
    Body, ExtractionOutcome, ExtractionReport, Manifest, ManifestArgs, Occurrence, Origin, SCHEMA,
    ScanCache, Summary, Totals, digest, executable_digest, load_manifest, parent, scan_source,
    write_line,
};
use crate::body_corpus::{ExtractedBody, Selection, extract};
use anyhow::{Context, Result, ensure};
use clap::Args;
use rayon::{ThreadPool, prelude::*};
use std::{
    collections::{BTreeMap, btree_map::Entry},
    fs::{self, File},
    io::{BufWriter, Write},
    path::{Path, PathBuf},
};

#[derive(Debug, Args)]
pub struct BuildArgs {
    #[command(flatten)]
    manifest: ManifestArgs,
    /// New output directory; existing corpora are never overwritten.
    #[arg(long)]
    out: PathBuf,
    /// Shared content-addressed scan cache outside the published corpus.
    #[arg(long)]
    cache: PathBuf,
    /// Publish manifests without materializing representative STEP objects.
    #[arg(long)]
    no_extract: bool,
    /// Fail on any unsupported body identity or skipped solid.
    #[arg(long)]
    require_canonical: bool,
}

fn preference(origin: &Origin) -> (bool, usize, &str, usize, u64) {
    let solid = &origin.solid;
    (
        !solid.has_shape_representation,
        solid.closure_entities,
        &origin.source.uuid,
        solid.data_section,
        solid.solid_id,
    )
}

pub(super) fn insert_body(
    groups: &mut BTreeMap<String, Body>,
    occurrence: Occurrence,
) -> Result<()> {
    let Occurrence {
        fingerprint,
        origin,
    } = occurrence;
    match groups.entry(fingerprint.clone()) {
        Entry::Vacant(entry) => {
            entry.insert(Body {
                fingerprint,
                occurrences: 1,
                representative: origin,
            });
        }
        Entry::Occupied(mut entry) => {
            let body = entry.get_mut();
            let a = &body.representative.solid;
            let b = &origin.solid;
            ensure!(
                a.topology() == b.topology()
                    && a.geometry_fingerprint == b.geometry_fingerprint
                    && a.identity_kind == b.identity_kind
                    && a.unit_scale_mm.map(f64::to_bits) == b.unit_scale_mm.map(f64::to_bits),
                "fingerprint collision or inconsistent body geometry {fingerprint}"
            );
            body.occurrences += 1;
            if preference(&origin) < preference(&body.representative) {
                body.representative = origin;
            }
        }
    }
    Ok(())
}

fn extract_group(
    source: &Path,
    source_sha256: &str,
    selections: &[Selection],
    directory: &Path,
) -> Result<Vec<ExtractedBody>> {
    let bytes = fs::read(source)?;
    ensure!(
        digest(&bytes) == source_sha256,
        "source changed since scan: {}",
        source.display()
    );
    extract(&bytes, &directory.join("objects"), selections)
}

fn materialize<'a>(
    bodies: impl Iterator<Item = &'a Body>,
    directory: &Path,
    pool: &ThreadPool,
) -> Result<usize> {
    fs::create_dir(directory.join("objects"))?;
    let mut groups = BTreeMap::<(&Path, &str), Vec<Selection>>::new();
    for body in bodies.filter(|body| body.extractable()) {
        let origin = &body.representative;
        groups
            .entry((&origin.source.path, &origin.source_sha256))
            .or_default()
            .push(Selection {
                section: origin.solid.data_section,
                solid_id: origin.solid.solid_id,
                fingerprint: body.fingerprint.clone(),
            });
    }
    let groups: Vec<_> = groups.into_iter().collect();
    let reports = pool.install(|| {
        groups
            .par_iter()
            .map(|((source, hash), selections)| {
                let outcome = match extract_group(source, hash, selections, directory) {
                    Ok(objects) => ExtractionOutcome::Extracted { objects },
                    Err(error) => ExtractionOutcome::Error {
                        error: format!("{error:#}"),
                    },
                };
                ExtractionReport {
                    source: source.to_path_buf(),
                    source_sha256: (*hash).into(),
                    requested: selections.len(),
                    outcome,
                }
            })
            .collect::<Vec<_>>()
    });
    serde_json::to_writer_pretty(
        File::create(directory.join("extraction-reports.json"))?,
        &reports,
    )?;
    Ok(reports
        .iter()
        .filter(|report| matches!(report.outcome, ExtractionOutcome::Error { .. }))
        .count())
}

fn collect_inventory(
    manifest: &Manifest,
    cache: &ScanCache,
    pool: &ThreadPool,
    directory: &Path,
) -> Result<(Totals, BTreeMap<String, Body>)> {
    let mut scans = BufWriter::new(File::create(directory.join("model-scans.jsonl"))?);
    let mut occurrences = BufWriter::new(File::create(directory.join("occurrences.jsonl"))?);
    let mut groups = BTreeMap::new();
    let mut totals = Totals::default();
    let jobs = pool.current_num_threads();
    // Keep at most one wave of full scan results in memory; JSONL and
    // representative-only state scale independently of corpus occurrence count.
    for chunk in manifest.sources.chunks(jobs) {
        let records = pool.install(|| {
            chunk
                .par_iter()
                .map(|source| scan_source(source, Some(cache)))
                .collect::<Vec<_>>()
        });
        for record in records {
            totals.record(&record);
            write_line(&mut scans, &record)?;
            for occurrence in record.into_occurrences() {
                write_line(&mut occurrences, &occurrence)?;
                insert_body(&mut groups, occurrence)?;
            }
        }
    }
    scans.flush()?;
    occurrences.flush()?;
    drop(scans);
    drop(occurrences);
    let mut bodies = BufWriter::new(File::create(directory.join("bodies.jsonl"))?);
    for body in groups.values() {
        write_line(&mut bodies, body)?;
    }
    bodies.flush()?;
    drop(bodies);
    Ok((totals, groups))
}

pub(super) fn run(args: &BuildArgs) -> Result<bool> {
    ensure!(
        !args.out.try_exists()?,
        "output already exists: {}",
        args.out.display()
    );
    fs::create_dir_all(parent(&args.out))?;
    let directory = tempfile::Builder::new()
        .prefix(".body-corpus-")
        .tempdir_in(parent(&args.out))?;
    let manifest = load_manifest(&args.manifest.manifest)?;
    let executable = executable_digest()?;
    let cache = ScanCache::new(&args.cache, &executable)?;
    let pool = args.manifest.parallel.pool()?;

    let (totals, groups) = collect_inventory(&manifest, &cache, &pool, directory.path())?;
    let extraction_failures = if args.no_extract {
        0
    } else {
        materialize(groups.values(), directory.path(), &pool)?
    };
    let pass = totals.succeeded(args.require_canonical, extraction_failures);
    let summary = Summary {
        schema: SCHEMA.into(),
        source_manifest_sha256: manifest.sha256,
        executable_sha256: executable,
        totals,
        unique_bodies: groups.len(),
        extractable_bodies: groups.values().filter(|body| body.extractable()).count(),
        objects_materialized: !args.no_extract,
        extraction_failures,
        require_canonical: args.require_canonical,
        pass,
    };
    serde_json::to_writer_pretty(
        File::create(directory.path().join("summary.json"))?,
        &summary,
    )?;
    // Unexpected I/O errors discard staging. Observed scan/extraction failures
    // still publish inspectable reports, marked FAILED, never COMPLETE.
    File::create(
        directory
            .path()
            .join(if pass { "COMPLETE" } else { "FAILED" }),
    )?;
    ensure!(
        !args.out.try_exists()?,
        "output created concurrently: {}",
        args.out.display()
    );
    fs::rename(directory.path(), &args.out)
        .with_context(|| format!("publish {}", args.out.display()))?;
    println!("{}", serde_json::to_string_pretty(&summary)?);
    Ok(pass)
}

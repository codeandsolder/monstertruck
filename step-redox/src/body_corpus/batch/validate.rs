use super::{
    Body, ExtractionOutcome, ExtractionReport, IdentityKind, Occurrence, Parallel, SCHEMA,
    ScanOutcome, ScanRecord, Summary, Totals, build::insert_body, digest, fingerprint, read_jsonl,
    valid_digest, valid_scan,
};
use crate::body_corpus::{ExtractedBody, scan_isolated_object};
use anyhow::{Context, Result, ensure};
use clap::Args;
use rayon::prelude::*;
use serde::Serialize;
use std::{
    collections::{BTreeMap, BTreeSet, HashSet},
    fs::{self, File},
    path::{Path, PathBuf},
};

#[derive(Debug, Args)]
pub struct ValidateArgs {
    #[arg(long)]
    corpus: PathBuf,
    #[command(flatten)]
    parallel: Parallel,
}

#[derive(Default, Serialize)]
struct Errors {
    count: usize,
    examples: Vec<String>,
}

impl Errors {
    fn record(&mut self, error: impl Into<String>) {
        self.count += 1;
        if self.examples.len() < 20 {
            self.examples.push(error.into());
        }
    }

    fn check(&mut self, condition: bool, error: impl Into<String>) {
        if !condition {
            self.record(error);
        }
    }
}

fn reconstruct(directory: &Path, errors: &mut Errors) -> Result<(Totals, BTreeMap<String, Body>)> {
    let mut totals = Totals::default();
    let mut groups = BTreeMap::new();
    let mut occurrences = read_jsonl::<Occurrence>(&directory.join("occurrences.jsonl"))?;
    let mut uuids = HashSet::new();
    for record in read_jsonl::<ScanRecord>(&directory.join("model-scans.jsonl"))? {
        let record = record?;
        ensure!(
            !record.source.uuid.is_empty() && uuids.insert(record.source.uuid.clone()),
            "missing or duplicate model UUID"
        );
        totals.record(&record);
        match &record.outcome {
            ScanOutcome::Error { error } => {
                errors.record(format!("{} scan failed: {error}", record.source.uuid));
            }
            ScanOutcome::Scanned {
                source_sha256,
                report,
                ..
            } => {
                ensure!(
                    valid_digest(source_sha256) && valid_scan(report),
                    "invalid saved scan {}",
                    record.source.uuid
                );
            }
        }
        for expected in record.into_occurrences() {
            let actual = occurrences.next().transpose()?;
            errors.check(
                actual.is_some_and(|actual| {
                    actual.fingerprint == expected.fingerprint && actual.origin == expected.origin
                }),
                format!(
                    "occurrence differs from scan {} #{}",
                    expected.origin.source.uuid, expected.origin.solid.solid_id
                ),
            );
            insert_body(&mut groups, expected)?;
        }
    }
    for extra in occurrences {
        let extra = extra?;
        errors.record(format!("unexpected occurrence {}", extra.fingerprint));
    }
    Ok((totals, groups))
}

fn check_body_manifest(
    directory: &Path,
    expected: &BTreeMap<String, Body>,
    errors: &mut Errors,
) -> Result<()> {
    let mut seen = HashSet::new();
    for body in read_jsonl::<Body>(&directory.join("bodies.jsonl"))? {
        let body = body?;
        errors.check(
            seen.insert(body.fingerprint.clone()),
            format!("duplicate body {}", body.fingerprint),
        );
        errors.check(
            expected.get(&body.fingerprint) == Some(&body),
            format!(
                "body record differs from reconstructed representative/count {}",
                body.fingerprint
            ),
        );
    }
    for key in expected.keys().filter(|key| !seen.contains(*key)) {
        errors.record(format!("missing body {key}"));
    }
    Ok(())
}

fn actual_files(directory: &Path, relative: &Path, files: &mut BTreeSet<PathBuf>) -> Result<()> {
    for entry in fs::read_dir(directory.join(relative))? {
        let entry = entry?;
        let path = relative.join(entry.file_name());
        let kind = entry.file_type()?;
        if kind.is_dir() {
            actual_files(directory, &path, files)?;
        } else {
            ensure!(
                kind.is_file(),
                "object is not a regular file: {}",
                path.display()
            );
            files.insert(path);
        }
    }
    Ok(())
}

fn validate_object(directory: &Path, body: &Body, object: &ExtractedBody) -> Result<()> {
    let path = directory.join(body.object_path());
    let bytes = fs::read(&path)?;
    ensure!(
        bytes.len() == object.bytes && digest(&bytes) == object.sha256,
        "object byte digest changed"
    );
    let report = scan_isolated_object(&bytes)?;
    ensure!(valid_scan(&report), "invalid extracted scan");
    ensure!(
        report.manifold_solids == 1 && report.analyzed_solids == 1,
        "expected exactly one analyzable solid"
    );
    let solid = report
        .solids
        .first()
        .context("extracted file has no solid")?;
    let reference = &body.representative.solid;
    ensure!(
        solid.identity_kind == reference.identity_kind,
        "body identity kind changed"
    );
    ensure!(
        solid.unit_scale_mm.map(f64::to_bits) == reference.unit_scale_mm.map(f64::to_bits),
        "physical unit scale changed"
    );
    if reference.identity_kind == IdentityKind::Canonical {
        ensure!(
            fingerprint(&body.representative.source, solid) == body.fingerprint,
            "canonical body identity changed"
        );
    }
    ensure!(
        solid.topology() == reference.topology(),
        "body topology counts changed"
    );
    Ok(())
}

fn extraction_records(
    directory: &Path,
    bodies: &BTreeMap<String, Body>,
    expected_failures: usize,
    errors: &mut Errors,
) -> Result<BTreeMap<String, ExtractedBody>> {
    let reports: Vec<ExtractionReport> =
        serde_json::from_reader(File::open(directory.join("extraction-reports.json"))?)?;
    let mut expected = BTreeMap::<(PathBuf, String), BTreeSet<&str>>::new();
    for body in bodies.values().filter(|body| body.extractable()) {
        let origin = &body.representative;
        expected
            .entry((origin.source.path.clone(), origin.source_sha256.clone()))
            .or_default()
            .insert(&body.fingerprint);
    }
    let mut records = BTreeMap::new();
    let mut failures = 0;
    for report in reports {
        let Some(fingerprints) =
            expected.remove(&(report.source.clone(), report.source_sha256.clone()))
        else {
            errors.record(format!(
                "unexpected or duplicate extraction group {}",
                report.source.display()
            ));
            continue;
        };
        errors.check(
            report.requested == fingerprints.len(),
            format!(
                "wrong extraction request count for {}",
                report.source.display()
            ),
        );
        match report.outcome {
            ExtractionOutcome::Error { error } => {
                failures += 1;
                errors.record(format!(
                    "{} extraction failed: {error}",
                    report.source.display()
                ));
            }
            ExtractionOutcome::Extracted { objects } => {
                let mut missing = fingerprints;
                for object in objects {
                    errors.check(
                        missing.remove(object.fingerprint.as_str()),
                        format!(
                            "unexpected or repeated object digest {}",
                            object.fingerprint
                        ),
                    );
                    errors.check(
                        valid_digest(&object.sha256) && object.bytes > 0,
                        format!("invalid object digest {}", object.fingerprint),
                    );
                    let key = object.fingerprint.clone();
                    errors.check(
                        records.insert(key.clone(), object).is_none(),
                        format!("duplicate object digest {key}"),
                    );
                }
                for key in missing {
                    errors.record(format!("missing object digest {key}"));
                }
            }
        }
    }
    errors.check(expected.is_empty(), "missing extraction group reports");
    errors.check(
        failures == expected_failures,
        "summary extraction failure count differs",
    );
    Ok(records)
}

fn validate_objects(
    args: &ValidateArgs,
    bodies: &BTreeMap<String, Body>,
    failures: usize,
    errors: &mut Errors,
) -> Result<usize> {
    let directory = &args.corpus;
    let digests = extraction_records(directory, bodies, failures, errors)?;
    let expected: BTreeSet<_> = bodies
        .values()
        .filter(|body| body.extractable())
        .map(Body::object_path)
        .collect();
    let mut actual = BTreeSet::new();
    let object_dir = Path::new("objects");
    if fs::symlink_metadata(directory.join(object_dir)).is_ok_and(|metadata| metadata.is_dir()) {
        actual_files(directory, object_dir, &mut actual)?;
    } else {
        errors.record("missing regular object directory");
    }
    for missing in expected.difference(&actual) {
        errors.record(format!("missing object {}", missing.display()));
    }
    for extra in actual.difference(&expected) {
        errors.record(format!("unexpected object {}", extra.display()));
    }
    let jobs: Vec<_> = bodies
        .values()
        .filter(|body| body.extractable() && actual.contains(&body.object_path()))
        .filter_map(|body| digests.get(&body.fingerprint).map(|digest| (body, digest)))
        .collect();
    let pool = args.parallel.pool()?;
    for chunk in jobs.chunks(args.parallel.jobs.get()) {
        let results = pool.install(|| {
            chunk
                .par_iter()
                .map(|(body, digest)| {
                    validate_object(directory, body, digest)
                        .err()
                        .map(|error| format!("{}: {error:#}", body.fingerprint))
                })
                .collect::<Vec<_>>()
        });
        for error in results.into_iter().flatten() {
            errors.record(error);
        }
    }
    Ok(jobs.len())
}

#[derive(Serialize)]
struct Validation {
    bodies: usize,
    occurrences: usize,
    extractable_bodies: usize,
    objects_checked: usize,
    pass: bool,
    errors: Errors,
}

pub(super) fn run(args: &ValidateArgs) -> Result<bool> {
    let directory = &args.corpus;
    let summary: Summary = serde_json::from_reader(File::open(directory.join("summary.json"))?)?;
    ensure!(
        summary.schema == SCHEMA,
        "unsupported body corpus schema {}",
        summary.schema
    );
    ensure!(
        valid_digest(&summary.source_manifest_sha256) && valid_digest(&summary.executable_sha256),
        "invalid corpus provenance digest"
    );
    let mut errors = Errors::default();
    let (totals, bodies) = reconstruct(directory, &mut errors)?;
    check_body_manifest(directory, &bodies, &mut errors)?;
    let extractable = bodies.values().filter(|body| body.extractable()).count();
    errors.check(summary.totals == totals, "summary scan totals differ");
    errors.check(
        summary.unique_bodies == bodies.len(),
        "summary unique body count differs",
    );
    errors.check(
        summary.extractable_bodies == extractable,
        "summary extractable body count differs",
    );
    let pass = totals.succeeded(summary.require_canonical, summary.extraction_failures);
    errors.check(summary.pass == pass, "summary pass status is inconsistent");
    errors.check(summary.pass, "corpus build reported failure");
    errors.check(
        directory.join("COMPLETE").is_file() == summary.pass
            && directory.join("FAILED").is_file() != summary.pass,
        "incorrect corpus completion markers",
    );
    let checked = if summary.objects_materialized {
        validate_objects(args, &bodies, summary.extraction_failures, &mut errors)?
    } else {
        errors.check(
            !directory.join("objects").try_exists()?,
            "metadata-only corpus contains materialized objects",
        );
        errors.check(
            !directory.join("extraction-reports.json").try_exists()?,
            "metadata-only corpus contains extraction reports",
        );
        errors.check(
            summary.extraction_failures == 0,
            "metadata-only corpus reports extraction failures",
        );
        0
    };
    let validation = Validation {
        bodies: bodies.len(),
        occurrences: totals.analyzed_solids,
        extractable_bodies: extractable,
        objects_checked: checked,
        pass: errors.count == 0,
        errors,
    };
    println!("{}", serde_json::to_string_pretty(&validation)?);
    Ok(validation.pass)
}

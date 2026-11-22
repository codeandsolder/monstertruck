//! Native, content-addressed body corpus orchestration.
mod build;
mod diagnostics;
mod validate;

use super::{ExtractedBody, IdentityKind, ScanReport, SolidScan, digest, scan_bytes};
use anyhow::{Context, Result, ensure};
use clap::{Args, Subcommand};
use rayon::{ThreadPool, ThreadPoolBuilder};
use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha256};
use std::{
    collections::HashSet,
    fs::{self, File},
    io::{BufRead, BufReader, Write},
    num::NonZeroUsize,
    path::{Path, PathBuf},
    process::ExitCode,
    sync::Mutex,
};
use tempfile::NamedTempFile;

const SCHEMA: &str = "step-redox.body-corpus.v3";

#[derive(Debug, Subcommand)]
pub enum CorpusCommand {
    /// Build one representative object per physical body identity.
    Build(build::BuildArgs),
    /// Validate published metadata and all materialized body objects.
    Validate(validate::ValidateArgs),
    /// Filter analytic support signatures in a source manifest.
    Census(diagnostics::CensusArgs),
    /// Find the first source with a noncanonical or missing body identity.
    FindGap(diagnostics::GapArgs),
}

/// Execute the selected native corpus operation.
///
/// # Errors
/// Returns an error for invalid inputs, I/O failures or unexpected internal errors.
pub fn run(command: CorpusCommand) -> Result<ExitCode> {
    let pass = match command {
        CorpusCommand::Build(args) => build::run(&args)?,
        CorpusCommand::Validate(args) => validate::run(&args)?,
        CorpusCommand::Census(args) => diagnostics::census(&args)?,
        CorpusCommand::FindGap(args) => diagnostics::find_gap(&args)?,
    };
    Ok(if pass {
        ExitCode::SUCCESS
    } else {
        ExitCode::FAILURE
    })
}

fn default_jobs() -> NonZeroUsize {
    std::thread::available_parallelism().unwrap_or(NonZeroUsize::MIN)
}

#[derive(Debug, Args)]
struct Parallel {
    /// Bound the number of simultaneous input parses.
    #[arg(long, default_value_t = default_jobs())]
    jobs: NonZeroUsize,
}

impl Parallel {
    fn pool(&self) -> Result<ThreadPool> {
        Ok(ThreadPoolBuilder::new()
            .num_threads(self.jobs.get())
            .build()?)
    }
}

#[derive(Debug, Args)]
struct ManifestArgs {
    /// UTF-8 TSV with mandatory uuid and path columns.
    #[arg(long)]
    manifest: PathBuf,
    #[command(flatten)]
    parallel: Parallel,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
struct Source {
    uuid: String,
    path: PathBuf,
}

struct Manifest {
    sources: Vec<Source>,
    sha256: String,
}

fn file_digest(path: &Path) -> Result<String> {
    let mut reader = BufReader::new(File::open(path)?);
    let mut hash = Sha256::new();
    loop {
        let bytes = reader.fill_buf()?;
        if bytes.is_empty() {
            break;
        }
        hash.update(bytes);
        let length = bytes.len();
        reader.consume(length);
    }
    Ok(hex::encode(hash.finalize()))
}

fn executable_digest() -> Result<String> { file_digest(&std::env::current_exe()?) }

fn parent(path: &Path) -> &Path {
    path.parent()
        .filter(|p| !p.as_os_str().is_empty())
        .unwrap_or_else(|| Path::new("."))
}

fn load_manifest(path: &Path) -> Result<Manifest> {
    let bytes = fs::read(path).with_context(|| format!("read {}", path.display()))?;
    let base = parent(path).canonicalize()?;
    let mut reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_reader(bytes.as_slice());
    let mut seen = HashSet::new();
    let mut sources = Vec::new();
    for row in reader.deserialize::<Source>() {
        let mut source =
            row.with_context(|| format!("malformed TSV manifest {}", path.display()))?;
        ensure!(
            !source.uuid.is_empty() && !source.path.as_os_str().is_empty(),
            "empty uuid or source path"
        );
        ensure!(
            seen.insert(source.uuid.clone()),
            "duplicate source UUID {}",
            source.uuid
        );
        if source.path.is_relative() {
            source.path = base.join(source.path);
        }
        sources.push(source);
    }
    ensure!(!sources.is_empty(), "manifest has no source rows");
    Ok(Manifest {
        sources,
        sha256: digest(&bytes),
    })
}

fn valid_digest(value: &str) -> bool {
    value.len() == 64
        && value
            .bytes()
            .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
}

fn valid_scan(report: &ScanReport) -> bool {
    report.analyzed_solids == report.solids.len()
        && report.analyzed_solids.checked_add(report.skipped_solids) == Some(report.manifold_solids)
        && report.fallback_identity_solids
            == report
                .solids
                .iter()
                .filter(|s| s.identity_kind == IdentityKind::FileLocalFallback)
                .count()
        && report.solids_without_shape_representation
            == report
                .solids
                .iter()
                .filter(|s| !s.has_shape_representation)
                .count()
        && report.solids.windows(2).all(|pair| {
            (pair[0].data_section, pair[0].solid_id) < (pair[1].data_section, pair[1].solid_id)
        })
        && report.solids.iter().all(|s| {
            s.unit_scale_mm.is_none_or(|v| v.is_finite() && v > 0.0)
                && match s.identity_kind {
                    IdentityKind::Canonical => {
                        s.geometry_fingerprint.as_ref().is_some_and(|v| {
                            v.len() == 32 && v.bytes().all(|b| b.is_ascii_hexdigit())
                        }) && s
                            .center_source_units
                            .is_some_and(|v| v.iter().all(|x| x.is_finite()))
                            && s.canonical_quarter_turn.is_some_and(|v| v < 4)
                    }
                    IdentityKind::FileLocalFallback => {
                        s.geometry_fingerprint.is_none()
                            && s.center_source_units.is_none()
                            && s.canonical_quarter_turn.is_none()
                    }
                }
        })
}

#[derive(Debug, Serialize, Deserialize)]
struct CachedScan {
    source_sha256: String,
    report: ScanReport,
}

struct ScanCache {
    namespace: PathBuf,
    // Same-content concurrent inputs share one parse within this process.
    shards: [Mutex<()>; 256],
}

impl ScanCache {
    fn new(root: &Path, executable: &str) -> Result<Self> {
        ensure!(valid_digest(executable), "invalid scanner SHA-256");
        let namespace = root.join(SCHEMA).join(executable);
        fs::create_dir_all(&namespace)?;
        Ok(Self {
            namespace,
            shards: std::array::from_fn(|_| Mutex::new(())),
        })
    }

    fn scan(&self, bytes: &[u8], source_sha256: &str) -> Result<(ScanReport, bool)> {
        let prefix = &source_sha256[..2];
        let index = usize::from(u8::from_str_radix(prefix, 16)?);
        let _guard = self.shards[index]
            .lock()
            .map_err(|_| anyhow::anyhow!("scan cache shard lock poisoned"))?;
        let shard = self.namespace.join(prefix);
        let path = shard.join(format!("{source_sha256}.json"));
        if let Ok(file) = File::open(&path)
            && let Ok(cached) = serde_json::from_reader::<_, CachedScan>(BufReader::new(file))
            && cached.source_sha256 == source_sha256
            && valid_scan(&cached.report)
        {
            return Ok((cached.report, true));
        }
        let report = scan_bytes(bytes)?;
        ensure!(valid_scan(&report), "invalid scan report invariants");
        fs::create_dir_all(&shard)?;
        let mut temp = NamedTempFile::new_in(&shard)?;
        let cached = CachedScan {
            source_sha256: source_sha256.into(),
            report,
        };
        serde_json::to_writer(&mut temp, &cached)?;
        temp.persist(path)?;
        Ok((cached.report, false))
    }
}

#[derive(Debug, Serialize, Deserialize)]
#[serde(tag = "status", rename_all = "snake_case")]
enum ScanOutcome {
    Scanned {
        source_sha256: String,
        cache_hit: bool,
        report: ScanReport,
    },
    Error {
        error: String,
    },
}

#[derive(Debug, Serialize, Deserialize)]
struct ScanRecord {
    source: Source,
    #[serde(flatten)]
    outcome: ScanOutcome,
}

impl ScanRecord {
    fn into_occurrences(self) -> impl Iterator<Item = Occurrence> {
        let data = match self.outcome {
            ScanOutcome::Scanned {
                source_sha256,
                report,
                ..
            } => Some((source_sha256, report.solids)),
            ScanOutcome::Error { .. } => None,
        };
        data.into_iter().flat_map(move |(source_sha256, solids)| {
            let source = self.source.clone();
            solids.into_iter().map(move |solid| Occurrence {
                fingerprint: fingerprint(&source, &solid),
                origin: Origin {
                    source: source.clone(),
                    source_sha256: source_sha256.clone(),
                    solid,
                },
            })
        })
    }
}

fn scan_source(source: &Source, cache: Option<&ScanCache>) -> ScanRecord {
    let outcome = (|| -> Result<ScanOutcome> {
        let bytes =
            fs::read(&source.path).with_context(|| format!("read {}", source.path.display()))?;
        let source_sha256 = digest(&bytes);
        let (report, cache_hit) = match cache {
            Some(cache) => cache.scan(&bytes, &source_sha256)?,
            None => (scan_bytes(&bytes)?, false),
        };
        ensure!(valid_scan(&report), "invalid scan report invariants");
        Ok(ScanOutcome::Scanned {
            source_sha256,
            cache_hit,
            report,
        })
    })()
    .unwrap_or_else(|error| ScanOutcome::Error {
        error: format!("{error:#}"),
    });
    ScanRecord {
        source: source.clone(),
        outcome,
    }
}

// New v3 identities use domain-separated binary unit bits. Unknown scales
// remain file-local; unsupported identities are never merged across sources.
fn fingerprint(source: &Source, solid: &SolidScan) -> String {
    let mut hash = Sha256::new();
    hash.update(b"step-redox/body/v3\0");
    match (&solid.geometry_fingerprint, solid.unit_scale_mm) {
        (Some(geometry), Some(scale)) => {
            hash.update(b"physical\0");
            hash.update(scale.to_bits().to_le_bytes());
            hash.update(geometry.as_bytes());
        }
        (Some(geometry), None) => {
            hash.update(b"unknown-unit\0");
            hash.update(source.uuid.len().to_string().as_bytes());
            hash.update([0]);
            hash.update(source.uuid.as_bytes());
            hash.update(geometry.as_bytes());
        }
        (None, _) => {
            hash.update(b"fallback\0");
            hash.update(source.uuid.as_bytes());
            hash.update([0]);
            hash.update(solid.data_section.to_string().as_bytes());
            hash.update([0]);
            hash.update(solid.solid_id.to_le_bytes());
        }
    }
    hex::encode(hash.finalize())
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
struct Origin {
    source: Source,
    source_sha256: String,
    solid: SolidScan,
}

#[derive(Debug, Serialize, Deserialize)]
struct Occurrence {
    fingerprint: String,
    origin: Origin,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
struct Body {
    fingerprint: String,
    occurrences: usize,
    representative: Origin,
}

impl Body {
    fn object_path(&self) -> PathBuf {
        Path::new("objects")
            .join(&self.fingerprint[..2])
            .join(format!("{}.step", self.fingerprint))
    }

    const fn extractable(&self) -> bool { self.representative.solid.has_shape_representation }
}

#[derive(Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
struct Totals {
    source_models: usize,
    scan_failures: usize,
    cache_hits: usize,
    manifold_solids: usize,
    analyzed_solids: usize,
    skipped_solids: usize,
    fallback_identity_solids: usize,
    solids_without_shape_representation: usize,
    unknown_unit_solids: usize,
}

impl Totals {
    const fn succeeded(&self, require_canonical: bool, extraction_failures: usize) -> bool {
        self.scan_failures == 0
            && extraction_failures == 0
            && (!require_canonical
                || (self.fallback_identity_solids == 0 && self.skipped_solids == 0))
    }

    fn record(&mut self, record: &ScanRecord) {
        self.source_models += 1;
        match &record.outcome {
            ScanOutcome::Error { .. } => self.scan_failures += 1,
            ScanOutcome::Scanned {
                cache_hit, report, ..
            } => {
                self.cache_hits += usize::from(*cache_hit);
                self.manifold_solids += report.manifold_solids;
                self.analyzed_solids += report.analyzed_solids;
                self.skipped_solids += report.skipped_solids;
                self.fallback_identity_solids += report.fallback_identity_solids;
                self.solids_without_shape_representation +=
                    report.solids_without_shape_representation;
                self.unknown_unit_solids += report
                    .solids
                    .iter()
                    .filter(|solid| solid.unit_scale_mm.is_none())
                    .count();
            }
        }
    }
}

#[derive(Debug, Serialize, Deserialize)]
struct ExtractionReport {
    source: PathBuf,
    source_sha256: String,
    requested: usize,
    #[serde(flatten)]
    outcome: ExtractionOutcome,
}

#[derive(Debug, Serialize, Deserialize)]
#[serde(tag = "status", rename_all = "snake_case")]
enum ExtractionOutcome {
    Extracted { objects: Vec<ExtractedBody> },
    Error { error: String },
}

#[derive(Debug, Serialize, Deserialize)]
struct Summary {
    schema: String,
    source_manifest_sha256: String,
    executable_sha256: String,
    totals: Totals,
    unique_bodies: usize,
    extractable_bodies: usize,
    objects_materialized: bool,
    extraction_failures: usize,
    require_canonical: bool,
    pass: bool,
}

fn write_line(writer: &mut impl Write, value: &impl Serialize) -> Result<()> {
    serde_json::to_writer(&mut *writer, value)?;
    writer.write_all(b"\n")?;
    Ok(())
}

fn read_jsonl<T: serde::de::DeserializeOwned>(
    path: &Path,
) -> Result<impl Iterator<Item = Result<T>> + use<T>> {
    let reader = BufReader::new(File::open(path)?);
    let path = path.to_path_buf();
    Ok(reader.lines().enumerate().map(move |(index, line)| {
        let line = line?;
        serde_json::from_str(&line)
            .with_context(|| format!("{}:{}: invalid JSONL record", path.display(), index + 1))
    }))
}

#[cfg(test)]
mod tests {
    use super::*;
    use clap::Parser;

    const FRUSTUM: &[u8] = include_bytes!("../../validation/fixtures/native_conical_frustum.step");
    const BOX: &[u8] = include_bytes!("../../validation/fixtures/native_body_corpus_box.step");
    const TORUS: &[u8] = include_bytes!("../../validation/fixtures/native_torus.step");

    #[derive(Parser)]
    struct TestCli {
        #[command(subcommand)]
        command: CorpusCommand,
    }

    fn execute(args: &[&str]) -> Result<bool> {
        let command =
            TestCli::try_parse_from(std::iter::once("corpus").chain(args.iter().copied()))?.command;
        Ok(run(command)? == ExitCode::SUCCESS)
    }

    fn source(path: &Path, uuid: &str) -> Source {
        Source {
            uuid: uuid.into(),
            path: path.to_path_buf(),
        }
    }

    #[test]
    fn scan_cache_tracks_content_and_executable_identity() -> Result<()> {
        let temp = tempfile::tempdir()?;
        let path = temp.path().join("part.step");
        fs::write(&path, FRUSTUM)?;
        let source = source(&path, "A");
        let executable = digest(b"scanner revision one");
        let cache = ScanCache::new(temp.path(), &executable)?;
        let first = scan_source(&source, Some(&cache));
        let ScanOutcome::Scanned {
            source_sha256,
            cache_hit,
            report,
        } = first.outcome
        else {
            anyhow::bail!("frustum scan failed")
        };
        assert!(!cache_hit);
        assert_eq!(report.manifold_solids, 1);
        assert_eq!(source_sha256, digest(FRUSTUM));
        assert!(matches!(
            scan_source(&source, Some(&cache)).outcome,
            ScanOutcome::Scanned {
                cache_hit: true,
                ..
            }
        ));
        fs::write(&path, TORUS)?;
        assert!(matches!(
            scan_source(&source, Some(&cache)).outcome,
            ScanOutcome::Scanned {
                cache_hit: false,
                ..
            }
        ));
        let other = ScanCache::new(temp.path(), &digest(b"scanner revision two"))?;
        assert!(matches!(
            scan_source(&source, Some(&other)).outcome,
            ScanOutcome::Scanned {
                cache_hit: false,
                ..
            }
        ));
        let torus_hash = digest(TORUS);
        let entry = cache
            .namespace
            .join(&torus_hash[..2])
            .join(format!("{torus_hash}.json"));
        fs::write(entry, b"not JSON")?;
        assert!(matches!(
            scan_source(&source, Some(&cache)).outcome,
            ScanOutcome::Scanned {
                cache_hit: false,
                ..
            }
        ));
        Ok(())
    }

    #[test]
    fn physical_fingerprint_distinguishes_units_and_local_fallbacks() -> Result<()> {
        let a = source(Path::new("a.step"), "a");
        let b = source(Path::new("b.step"), "b");
        let report = scan_bytes(FRUSTUM)?;
        let mut solid = report
            .solids
            .first()
            .context("fixture has no solid")?
            .clone();
        // Define a synthetic canonical identity independently of the fixture's
        // geometry grammar, to exercise every identity domain explicitly.
        solid.geometry_fingerprint = Some("0123456789abcdef0123456789abcdef".into());
        solid.unit_scale_mm = Some(1.0);
        let mm = fingerprint(&a, &solid);
        assert_eq!(mm, fingerprint(&b, &solid));
        solid.unit_scale_mm = Some(25.4);
        assert_ne!(mm, fingerprint(&a, &solid));
        solid.unit_scale_mm = None;
        assert_ne!(fingerprint(&a, &solid), fingerprint(&b, &solid));
        solid.geometry_fingerprint = None;
        let fallback = fingerprint(&a, &solid);
        assert_ne!(fallback, fingerprint(&b, &solid));
        solid.solid_id += 1;
        assert_ne!(fallback, fingerprint(&a, &solid));
        Ok(())
    }

    #[test]
    fn corpus_build_validate_and_detect_tampering() -> Result<()> {
        let temp = tempfile::tempdir()?;
        let root = temp.path();
        fs::write(root.join("part A.step"), BOX)?;
        fs::write(root.join("part B.step"), BOX)?;
        let manifest = root.join("sources.tsv");
        fs::write(
            &manifest,
            "uuid\tpath\textra\nB\t\"part B.step\"\t1\nA\t\"part A.step\"\t2\n",
        )?;
        let corpus = root.join("corpus");
        let cache = root.join("cache");
        let m = manifest.to_str().context("manifest path not UTF-8")?;
        let c = corpus.to_str().context("corpus path not UTF-8")?;
        let k = cache.to_str().context("cache path not UTF-8")?;
        assert!(execute(&[
            "build",
            "--manifest",
            m,
            "--out",
            c,
            "--cache",
            k,
            "--jobs",
            "1"
        ])?);
        assert!(execute(&["validate", "--corpus", c, "--jobs", "2"])?);
        let summary: Summary = serde_json::from_reader(File::open(corpus.join("summary.json"))?)?;
        assert_eq!(summary.totals.analyzed_solids, 2);
        assert_eq!(summary.totals.cache_hits, 1);
        assert_eq!(summary.extractable_bodies, summary.unique_bodies);
        assert_eq!(summary.totals.fallback_identity_solids, 0);
        assert_eq!(summary.totals.unknown_unit_solids, 0);
        assert_eq!(summary.unique_bodies, 1);
        let body = read_jsonl::<Body>(&corpus.join("bodies.jsonl"))?
            .next()
            .context("missing body")??;
        assert_eq!(body.representative.source.uuid, "A");
        assert_eq!(body.occurrences, 2);
        assert!(execute(&["build", "--manifest", m, "--out", c, "--cache", k]).is_err());
        let bodies = corpus.join("bodies.jsonl");
        let original = fs::read(&bodies)?;
        fs::write(&bodies, b"")?;
        assert!(!execute(&["validate", "--corpus", c, "--jobs", "1"])?);
        fs::write(&bodies, original)?;
        let body = read_jsonl::<Body>(&bodies)?
            .next()
            .context("missing body")??;
        let object = corpus.join(body.object_path());
        let original_bytes = fs::read(&object)?;
        let mut modified = original_bytes.clone();
        modified.push(b'\n');
        fs::write(&object, modified)?;
        assert!(!execute(&["validate", "--corpus", c, "--jobs", "1"])?);
        fs::write(&object, original_bytes)?;
        assert!(execute(&["validate", "--corpus", c, "--jobs", "1"])?);
        fs::remove_file(object)?;
        assert!(!execute(&["validate", "--corpus", c, "--jobs", "1"])?);
        let metadata = root.join("metadata");
        let md = metadata.to_str().context("metadata path not UTF-8")?;
        assert!(execute(&[
            "build",
            "--manifest",
            m,
            "--out",
            md,
            "--cache",
            k,
            "--no-extract",
            "--jobs",
            "2"
        ])?);
        assert!(execute(&["validate", "--corpus", md])?);
        assert!(!metadata.join("objects").exists());
        Ok(())
    }

    #[test]
    fn fallback_bodies_remain_source_local() -> Result<()> {
        let temp = tempfile::tempdir()?;
        let root = temp.path();
        fs::write(root.join("a.step"), FRUSTUM)?;
        fs::write(root.join("b.step"), FRUSTUM)?;
        let manifest = root.join("sources.tsv");
        fs::write(&manifest, "uuid\tpath\na\ta.step\nb\tb.step\n")?;
        let corpus = root.join("corpus");
        let cache = root.join("cache");
        let m = manifest.to_str().context("manifest path not UTF-8")?;
        let c = corpus.to_str().context("corpus path not UTF-8")?;
        let k = cache.to_str().context("cache path not UTF-8")?;
        assert!(execute(&[
            "build",
            "--manifest",
            m,
            "--out",
            c,
            "--cache",
            k,
            "--jobs",
            "2"
        ])?);
        assert!(execute(&["validate", "--corpus", c, "--jobs", "2"])?);
        let summary: Summary = serde_json::from_reader(File::open(corpus.join("summary.json"))?)?;
        assert_eq!(summary.totals.fallback_identity_solids, 2);
        assert_eq!(summary.unique_bodies, 2);
        assert_eq!(summary.totals.cache_hits, 1);
        let strict = root.join("strict");
        let t = strict.to_str().context("strict path not UTF-8")?;
        assert!(!execute(&[
            "build",
            "--manifest",
            m,
            "--out",
            t,
            "--cache",
            k,
            "--require-canonical",
            "--no-extract"
        ])?);
        assert!(strict.join("FAILED").is_file());
        assert!(!strict.join("COMPLETE").exists());
        Ok(())
    }

    #[test]
    fn diagnostics_report_first_failed_source_in_manifest_order() -> Result<()> {
        let temp = tempfile::tempdir()?;
        let root = temp.path();
        fs::write(root.join("valid.step"), FRUSTUM)?;
        let manifest = root.join("sources.tsv");
        fs::write(
            &manifest,
            "uuid\tpath\nfirst\tmissing.step\nsecond\tvalid.step\n",
        )?;
        let m = manifest.to_str().context("manifest path not UTF-8")?;
        let out = root.join("supports.jsonl");
        let o = out.to_str().context("census path not UTF-8")?;
        assert!(!execute(&["find-gap", "--manifest", m, "--jobs", "2"])?);
        assert!(!execute(&[
            "census",
            "--manifest",
            m,
            "--out",
            o,
            "--jobs",
            "2"
        ])?);
        let lines = fs::read_to_string(out)?;
        assert!(lines.contains("\"uuid\":\"first\""));
        assert!(lines.contains("\"uuid\":\"second\""));
        let duplicate = root.join("duplicate.tsv");
        fs::write(&duplicate, "uuid\tpath\na\tvalid.step\na\tvalid.step\n")?;
        assert!(load_manifest(&duplicate).is_err());
        Ok(())
    }
}

use crate::{
    FullComparison, compare_steps, fixtures,
    geometry::{
        compare::{Alignment, Thresholds},
        render::View,
    },
};
use anyhow::{Context, Result, ensure};
use clap::Args;
use occt_wasm::{OcctKernel, wasmtime::Engine};
use serde::{Deserialize, Serialize};
use std::{
    collections::{BTreeMap, BTreeSet},
    fmt::Write as _,
    fs,
    path::{Path, PathBuf},
    process::Command,
};

#[derive(Debug, Deserialize)]
struct Manifest {
    schema: u32,
    fixtures: BTreeMap<String, fixtures::Fixture>,
    cases: Vec<Case>,
}

fn enabled() -> bool { true }

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
struct Case {
    name: String,
    #[serde(default = "enabled")]
    enabled: bool,
    input_fixture: String,
    reference_fixture: Option<String>,
    command: Vec<String>,
    #[serde(default)]
    alignment: Alignment,
    #[serde(default)]
    thresholds: Thresholds,
    #[serde(default)]
    views: Vec<View>,
}

#[derive(Debug, Args)]
pub struct FixtureArgs {
    #[arg(long)]
    manifest: PathBuf,
    #[arg(long)]
    cache: PathBuf,
    /// Explicit subset of cases; omitted selects all enabled cases.
    #[arg(long = "case")]
    cases: Vec<String>,
}

#[derive(Debug, Args)]
pub struct RunArgs {
    #[command(flatten)]
    fixtures: FixtureArgs,
    #[arg(long)]
    out: PathBuf,
    #[arg(long)]
    step_redox: PathBuf,
    #[arg(long)]
    step_count_resize: Option<PathBuf>,
    #[arg(long)]
    step_chain_resize: Option<PathBuf>,
}

#[derive(Serialize)]
struct CaseReport<'a> {
    name: &'a str,
    command: Vec<String>,
    thresholds: &'a Thresholds,
    returncode: Option<i32>,
    stdout: String,
    stderr: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    comparison: Option<FullComparison>,
    #[serde(skip_serializing_if = "Option::is_none")]
    error: Option<String>,
    pass: bool,
}

fn expand(argument: &str, values: &BTreeMap<&str, Option<&Path>>) -> Result<String> {
    let mut result = String::new();
    let mut remaining = argument;
    while let Some(index) = remaining.find(['{', '}']) {
        result.push_str(&remaining[..index]);
        remaining = &remaining[index..];
        if let Some(tail) = remaining.strip_prefix("{{") {
            result.push('{');
            remaining = tail;
        } else if let Some(tail) = remaining.strip_prefix("}}") {
            result.push('}');
            remaining = tail;
        } else {
            let tail = remaining
                .strip_prefix('{')
                .context("unmatched closing template brace")?;
            let (key, rest) = tail
                .split_once('}')
                .context("unclosed command placeholder")?;
            let path = values
                .get(key)
                .with_context(|| format!("unknown command placeholder {key}"))?
                .with_context(|| format!("command requires {key} but no path was supplied"))?;
            result.push_str(path.to_str().context("non-UTF-8 command path")?);
            remaining = rest;
        }
    }
    result.push_str(remaining);
    Ok(result)
}

fn tail_text(bytes: &[u8]) -> String {
    // Retain the final diagnostic text even for unusually verbose generators.
    let start = bytes.len().saturating_sub(10_000);
    String::from_utf8_lossy(&bytes[start..]).into_owned()
}

impl FixtureArgs {
    fn load(&self) -> Result<Manifest> {
        let manifest: Manifest = serde_json::from_reader(fs::File::open(&self.manifest)?)?;
        ensure!(
            manifest.schema == 1,
            "unsupported STEP manifest schema {}",
            manifest.schema
        );
        Ok(manifest)
    }

    fn select<'a>(&self, manifest: &'a Manifest) -> Result<Vec<&'a Case>> {
        let available: BTreeSet<_> = manifest
            .cases
            .iter()
            .map(|case| case.name.as_str())
            .collect();
        ensure!(
            available.len() == manifest.cases.len(),
            "duplicate manifest case names"
        );
        for name in &self.cases {
            ensure!(
                manifest
                    .cases
                    .iter()
                    .any(|case| case.enabled && case.name == *name),
                "unknown or disabled validation case {name}"
            );
        }
        let selected: Vec<_> = manifest
            .cases
            .iter()
            .filter(|case| {
                case.enabled && (self.cases.is_empty() || self.cases.contains(&case.name))
            })
            .collect();
        ensure!(!selected.is_empty(), "no enabled validation cases selected");
        Ok(selected)
    }

    fn acquire<'a>(
        &self,
        manifest: &'a Manifest,
        cases: &[&'a Case],
    ) -> Result<BTreeMap<&'a str, PathBuf>> {
        let base = self
            .manifest
            .parent()
            .context("manifest has no containing directory")?;
        let mut needed = BTreeSet::new();
        for case in cases {
            fixtures::safe_filename(&case.name)?;
            needed.insert(case.input_fixture.as_str());
            needed.insert(
                case.reference_fixture
                    .as_deref()
                    .unwrap_or(&case.input_fixture),
            );
        }
        needed
            .into_iter()
            .map(|name| {
                let spec = manifest
                    .fixtures
                    .get(name)
                    .with_context(|| format!("unknown fixture {name}"))?;
                let path = fixtures::fetch(name, spec, base, &self.cache)?;
                Ok((name, path))
            })
            .collect()
    }
}

pub fn fetch(args: &FixtureArgs) -> Result<()> {
    let manifest = args.load()?;
    let cases = args.select(&manifest)?;
    let paths = args.acquire(&manifest, &cases)?;
    println!("{}", serde_json::to_string_pretty(&paths)?);
    Ok(())
}

pub fn run(args: &RunArgs) -> Result<bool> {
    let manifest = args.fixtures.load()?;
    let cases = args.fixtures.select(&manifest)?;
    let paths = args.fixtures.acquire(&manifest, &cases)?;
    fs::create_dir_all(&args.out)?;
    let engine = Engine::new(&OcctKernel::engine_config())?;
    let module = OcctKernel::compile(&engine)?;
    let mut reports = Vec::with_capacity(cases.len());
    for case in cases {
        let directory = args.out.join(&case.name);
        fs::create_dir_all(&directory)?;
        let output = directory.join("candidate.step");
        if output.try_exists()? {
            fs::remove_file(&output)?;
        }
        let input = paths
            .get(case.input_fixture.as_str())
            .context("input fixture was not acquired")?;
        let reference = paths
            .get(
                case.reference_fixture
                    .as_deref()
                    .unwrap_or(&case.input_fixture),
            )
            .context("reference fixture was not acquired")?;
        let values = BTreeMap::from([
            ("step_redox", Some(args.step_redox.as_path())),
            ("step_count_resize", args.step_count_resize.as_deref()),
            ("step_chain_resize", args.step_chain_resize.as_deref()),
            ("input", Some(input.as_path())),
            ("reference", Some(reference.as_path())),
            ("output", Some(output.as_path())),
        ]);
        let command = case
            .command
            .iter()
            .map(|argument| expand(argument, &values))
            .collect::<Result<Vec<_>>>()?;
        let (program, argv) = command.split_first().context("empty generator command")?;
        let mut report = CaseReport {
            name: &case.name,
            thresholds: &case.thresholds,
            command: command.clone(),
            returncode: None,
            stdout: String::new(),
            stderr: String::new(),
            comparison: None,
            error: None,
            pass: false,
        };
        match Command::new(program).args(argv).output() {
            Ok(process) => {
                report.returncode = process.status.code();
                report.stdout = tail_text(&process.stdout);
                report.stderr = tail_text(&process.stderr);
                if process.status.success() && output.is_file() {
                    match compare_steps(
                        &module,
                        reference,
                        &output,
                        &case.thresholds,
                        case.alignment,
                        &case.views,
                        &directory,
                    ) {
                        Ok(comparison) => {
                            report.pass = comparison.pass;
                            report.comparison = Some(comparison);
                        }
                        Err(error) => {
                            report.error = Some(format!("{error:#}"));
                        }
                    }
                } else {
                    report.error = Some(format!(
                        "generator command failed ({}) or produced no STEP file",
                        process.status
                    ));
                }
            }
            Err(error) => report.error = Some(format!("cannot execute {program}: {error}")),
        }
        fs::write(
            directory.join("report.json"),
            serde_json::to_vec_pretty(&report)?,
        )?;
        eprintln!(
            "{}: {}",
            report.name,
            if report.pass { "PASS" } else { "FAIL" }
        );
        reports.push(report);
    }
    fs::write(
        args.out.join("report.json"),
        serde_json::to_vec_pretty(&reports)?,
    )?;
    let passed = reports.iter().filter(|report| report.pass).count();
    let mut summary = format!(
        "# STEP validation\n\n**{passed}/{} cases passed.**\n\n",
        reports.len()
    );
    for report in &reports {
        writeln!(
            summary,
            "- {}: **{}** — [{}/report.json]({}/report.json)",
            report.name,
            if report.pass { "PASS" } else { "FAIL" },
            report.name,
            report.name
        )?;
    }
    fs::write(args.out.join("summary.md"), summary)?;
    Ok(passed == reports.len())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn real_manifest_acquires_gzip_and_zip_fixtures() -> Result<()> {
        let cache = tempfile::tempdir()?;
        let args = FixtureArgs {
            manifest: Path::new(env!("CARGO_MANIFEST_DIR")).join("../../validation/fixtures.json"),
            cache: cache.path().to_path_buf(),
            cases: vec!["compact_header_2x24".into(), "compact_te_3_353293_4".into()],
        };
        let manifest = args.load()?;
        let selected = args.select(&manifest)?;
        let paths = args.acquire(&manifest, &selected)?;
        assert!(paths.contains_key("jlc_header_2x24_2mm"));
        assert!(paths.contains_key("te_3_353293_4"));
        // The second acquisition must come from the verified cache.
        assert_eq!(paths, args.acquire(&manifest, &selected)?);
        Ok(())
    }

    #[test]
    fn template_expansion_preserves_literal_paths_and_rejects_unknown_fields() -> Result<()> {
        let values = BTreeMap::from([
            ("input", Some(Path::new("path with {literal} braces.step"))),
            ("step_count_resize", None),
        ]);
        assert_eq!(
            expand("--input={input} {{literal}}", &values)?,
            "--input=path with {literal} braces.step {literal}"
        );
        for bad in ["{}", "{unknown}", "{input", "x}", "{step_count_resize}"] {
            assert!(expand(bad, &values).is_err(), "{bad}");
        }
        Ok(())
    }
}

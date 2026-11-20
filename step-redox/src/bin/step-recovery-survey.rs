//! Read-only corpus measurement of automatic STEP-to-constructive-IR recovery.
//!
//! Structural recovery and independently source-boundary-verified admission
//! are counted separately, without inferring verification from recovery family.

#[cfg(feature = "cad-kernel-monstertruck")]
mod survey {
    use anyhow::{Context, Result, ensure};
    use clap::Parser;
    use serde::{Deserialize, Serialize};
    use sha2::{Digest, Sha256};
    use std::collections::BTreeMap;
    use std::io::Write;
    use std::path::{Path, PathBuf};
    use std::time::Instant;
    use step_redox::complete_ir::{CompleteIr, KernelRecoveryFailure};
    use step_redox::{SourceEntity, cad_kernel::monstertruck::MonstertruckKernel};

    mod journal;
    const SCHEMA: &str = "step-redox.recovery-survey.v3";

    #[derive(Debug, Parser)]
    #[command(about = "Measure automatic constructive recovery across STEP inputs as JSON")]
    struct Cli {
        #[arg(long, value_name = "JSON")]
        output: Option<PathBuf>,
        /// Resume an interrupted survey with identical binary and ordered inputs.
        #[arg(long, requires = "output")]
        resume: bool,
        #[arg(required = true, value_name = "STEP")]
        inputs: Vec<PathBuf>,
    }

    #[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
    #[serde(rename_all = "snake_case")]
    enum Outcome {
        InputError,
        NoSolids,
        CompleteConstructive,
        PartialConstructive,
        SemanticOnly,
        Unrecovered,
    }

    #[derive(Debug, Serialize, Deserialize)]
    struct Case {
        input: String,
        input_bytes: Option<usize>,
        sha256: Option<String>,
        elapsed_ms: u128,
        outcome: Outcome,
        body_count: Option<usize>,
        constructively_recovered: Option<usize>,
        semantically_recovered: Option<usize>,
        unrecovered: Option<usize>,
        families: BTreeMap<String, usize>,
        source_boundary_verified_sources: Vec<SourceEntity>,
        unrecovered_sources: Vec<SourceEntity>,
        kernel_failures: Vec<KernelRecoveryFailure>,
        error: Option<String>,
    }

    #[derive(Debug, Default, Serialize)]
    struct Totals {
        inputs: usize,
        inputs_with_errors: usize,
        files_with_no_solids: usize,
        fully_constructive_files: usize,
        fully_source_boundary_verified_files: usize,
        partially_constructive_files: usize,
        semantic_only_files: usize,
        unrecovered_files: usize,
        source_solids: usize,
        constructively_recovered_solids: usize,
        source_boundary_verified_solids: usize,
        semantically_recovered_solids: usize,
        unrecovered_solids: usize,
    }

    #[derive(Debug, Serialize)]
    struct Survey {
        schema: &'static str,
        executable_sha256: String,
        /// Only explicitly verified source IDs count as independent boundary evidence.
        claim: &'static str,
        totals: Totals,
        cases: Vec<Case>,
    }

    const fn classify(body_count: usize, constructive: usize, semantic: usize) -> Outcome {
        if body_count == 0 {
            Outcome::NoSolids
        } else if constructive == body_count {
            Outcome::CompleteConstructive
        } else if constructive > 0 {
            Outcome::PartialConstructive
        } else if semantic > 0 {
            Outcome::SemanticOnly
        } else {
            Outcome::Unrecovered
        }
    }

    fn count_families(report: &CompleteIr) -> BTreeMap<String, usize> {
        [
            ("extrusion", report.solid_extrusions.len()),
            ("loft", report.solid_lofts.len()),
            ("revolution", report.solid_revolutions.len()),
            (
                "radial_slot_revolution",
                report.radial_slot_revolutions.len(),
            ),
            ("sweep", report.solid_sweeps.len()),
            ("periodic_chain", report.periodic_chains.len()),
            (
                "boundary_feature_composite",
                report.boundary_feature_composites.len(),
            ),
            (
                "edge_treatment_composite",
                report.edge_treatment_composites.len(),
            ),
            (
                "same_domain_normalized",
                report.same_domain_normalized_solids.len(),
            ),
            ("extrusion_core_union", report.extrusion_core_unions.len()),
        ]
        .into_iter()
        .filter(|(_, count)| *count > 0)
        .map(|(family, count)| (family.to_owned(), count))
        .collect()
    }

    fn recover_case(path: &Path, kernel: &MonstertruckKernel) -> Case {
        let started = Instant::now();
        let mut entry = Case {
            input: path.to_string_lossy().into_owned(),
            input_bytes: None,
            sha256: None,
            elapsed_ms: 0,
            outcome: Outcome::InputError,
            body_count: None,
            constructively_recovered: None,
            semantically_recovered: None,
            unrecovered: None,
            families: BTreeMap::new(),
            source_boundary_verified_sources: Vec::new(),
            unrecovered_sources: Vec::new(),
            kernel_failures: Vec::new(),
            error: None,
        };
        let result = (|| -> Result<CompleteIr> {
            let bytes = std::fs::read(path).with_context(|| format!("read {}", path.display()))?;
            entry.input_bytes = Some(bytes.len());
            entry.sha256 = Some(hex::encode(Sha256::digest(&bytes)));
            step_redox::complete_ir::recover_complete_ir_with_kernel_bytes(&bytes, kernel)
                .with_context(|| format!("recover {}", path.display()))
        })();
        match result {
            Ok(report) => {
                entry.outcome = classify(
                    report.body_count,
                    report.constructively_recovered_solid_count,
                    report.semantically_recovered_solid_count,
                );
                entry.body_count = Some(report.body_count);
                entry.constructively_recovered = Some(report.constructively_recovered_solid_count);
                entry.semantically_recovered = Some(report.semantically_recovered_solid_count);
                entry.unrecovered = Some(report.unrecovered_solid_count);
                entry.families = count_families(&report);
                entry.unrecovered_sources = report
                    .unrecovered_solids
                    .iter()
                    .map(|solid| solid.source)
                    .collect();
                entry.source_boundary_verified_sources = report.source_boundary_verified_solids;
                entry.kernel_failures = report.kernel_recovery_failures;
            }
            Err(error) => entry.error = Some(format!("{error:#}")),
        }
        entry.elapsed_ms = started.elapsed().as_millis();
        entry
    }

    impl Totals {
        fn include(&mut self, entry: &Case) {
            self.inputs += 1;
            match entry.outcome {
                Outcome::InputError => self.inputs_with_errors += 1,
                Outcome::NoSolids => self.files_with_no_solids += 1,
                Outcome::CompleteConstructive => self.fully_constructive_files += 1,
                Outcome::PartialConstructive => self.partially_constructive_files += 1,
                Outcome::SemanticOnly => self.semantic_only_files += 1,
                Outcome::Unrecovered => self.unrecovered_files += 1,
            }
            let verified = entry.source_boundary_verified_sources.len();
            self.source_boundary_verified_solids += verified;
            if verified > 0 && entry.body_count == Some(verified) {
                self.fully_source_boundary_verified_files += 1;
            }
            self.source_solids += entry.body_count.unwrap_or_default();
            self.constructively_recovered_solids +=
                entry.constructively_recovered.unwrap_or_default();
            self.semantically_recovered_solids += entry.semantically_recovered.unwrap_or_default();
            self.unrecovered_solids += entry.unrecovered.unwrap_or_default();
        }
    }

    fn run(cli: &Cli) -> Result<()> {
        let inputs = cli
            .inputs
            .iter()
            .map(std::path::absolute)
            .collect::<std::io::Result<Vec<_>>>()?;
        let output = cli.output.as_ref().map(std::path::absolute).transpose()?;
        let identity = journal::Header::current(&inputs)?;
        let (mut journal, mut cases) = if let Some(output) = &output {
            let checkpoint = output.with_extension("cases.jsonl");
            ensure!(
                output != &checkpoint,
                "report and journal must have distinct paths"
            );
            ensure!(
                !inputs
                    .iter()
                    .any(|input| input == output || input == &checkpoint),
                "report or journal path is also a source input"
            );
            ensure!(
                !output.try_exists()?,
                "report already exists: {}",
                output.display()
            );
            let (journal, cases) = journal::Journal::open(&checkpoint, &identity, cli.resume)?;
            (Some(journal), cases)
        } else {
            (None, Vec::new())
        };
        let kernel = MonstertruckKernel::project_standard();
        let mut totals = Totals::default();
        for entry in &cases {
            totals.include(entry);
        }
        for path in inputs.iter().skip(cases.len()) {
            let entry = recover_case(path, &kernel);
            totals.include(&entry);
            if let Some(journal) = journal.as_mut() {
                journal.append(&entry)?;
            }
            cases.push(entry);
        }
        let report = Survey {
            schema: SCHEMA,
            executable_sha256: identity.executable_sha256,
            claim: "automatic constructive IR selection; source-boundary verified IDs are explicit independent analytic equivalence certificates",
            totals,
            cases,
        };
        let mut bytes = serde_json::to_vec_pretty(&report)?;
        bytes.push(b'\n');
        if let Some(output) = output {
            journal::write_report(&output, &bytes)
        } else {
            std::io::stdout().lock().write_all(&bytes)?;
            Ok(())
        }
    }

    pub fn main() -> Result<()> { run(&Cli::parse()) }

    #[cfg(test)]
    mod tests {
        use super::*;

        #[test]
        fn distinguishes_full_partial_semantic_and_unrecovered_cases() {
            assert_eq!(classify(0, 0, 0), Outcome::NoSolids);
            assert_eq!(classify(2, 2, 2), Outcome::CompleteConstructive);
            assert_eq!(classify(2, 1, 1), Outcome::PartialConstructive);
            assert_eq!(classify(2, 0, 1), Outcome::SemanticOnly);
            assert_eq!(classify(2, 0, 0), Outcome::Unrecovered);
        }

        #[test]
        fn resumed_report_reuses_a_completed_real_source_case() -> Result<()> {
            let dir = tempfile::tempdir()?;
            let first = dir.path().join("loft.step");
            let second = dir.path().join("missing.step");
            let output = dir.path().join("report.json");
            let checkpoint = output.with_extension("cases.jsonl");
            std::fs::write(
                &first,
                include_bytes!("../../validation/fixtures/native_staged_rectangular_loft.step"),
            )?;
            let run_survey = |resume| {
                run(&Cli {
                    output: Some(output.clone()),
                    resume,
                    inputs: vec![first.clone(), second.clone()],
                })
            };
            run_survey(false)?;
            let original: serde_json::Value = serde_json::from_slice(&std::fs::read(&output)?)?;
            assert_eq!(original["totals"]["fully_constructive_files"], 1);
            assert_eq!(
                original["totals"]["fully_source_boundary_verified_files"],
                1
            );
            assert_eq!(original["totals"]["source_boundary_verified_solids"], 1);
            let verified = original["cases"][0]["source_boundary_verified_sources"]
                .as_array()
                .context("missing verified source IDs")?;
            assert_eq!(verified.len(), 1);
            assert_eq!(original["totals"]["inputs_with_errors"], 1);
            assert!(run_survey(false).is_err());
            std::fs::remove_file(&output)?;
            let contents = std::fs::read_to_string(&checkpoint)?;
            let prefix = contents.split_inclusive('\n').take(2).collect::<String>();
            std::fs::write(&checkpoint, format!("{prefix}{{\"input\":"))?;
            run_survey(true)?;
            let resumed: serde_json::Value = serde_json::from_slice(&std::fs::read(&output)?)?;
            assert_eq!(original["cases"][0], resumed["cases"][0]);
            assert_eq!(original["totals"], resumed["totals"]);
            assert_eq!(original["executable_sha256"], resumed["executable_sha256"]);
            Ok(())
        }

        #[test]
        fn missing_input_does_not_abort_the_corpus() {
            let kernel = MonstertruckKernel::project_standard();
            let entry = recover_case(Path::new("/missing/monstertruck-survey.step"), &kernel);
            assert_eq!(entry.outcome, Outcome::InputError);
            assert!(entry.error.is_some());
            assert_eq!(entry.body_count, None);
        }
    }
}

#[cfg(feature = "cad-kernel-monstertruck")]
fn main() -> anyhow::Result<()> { survey::main() }

#[cfg(not(feature = "cad-kernel-monstertruck"))]
fn main() -> anyhow::Result<()> {
    anyhow::bail!("step-recovery-survey requires --features cad-kernel-monstertruck")
}

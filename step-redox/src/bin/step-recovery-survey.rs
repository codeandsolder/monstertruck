//! Read-only, corpus-level measurement of automatic STEP-to-constructive-IR recovery.
//!
//! Constructive recovery and independent source-exact admission are distinct
//! claims: this command reports the former, never silently asserts the latter.

#[cfg(feature = "cad-kernel-monstertruck")]
mod survey {
    use anyhow::{Context, Result};
    use clap::Parser;
    use serde::Serialize;
    use sha2::{Digest, Sha256};
    use std::collections::BTreeMap;
    use std::io::Write;
    use std::path::{Path, PathBuf};
    use std::time::Instant;
    use step_redox::cad_kernel::monstertruck::MonstertruckKernel;
    use step_redox::complete_ir::{CompleteIr, KernelRecoveryFailure};

    #[derive(Debug, Parser)]
    #[command(about = "Measure automatic constructive recovery across STEP inputs as JSON")]
    struct Cli {
        #[arg(long, value_name = "JSON")]
        output: Option<PathBuf>,
        #[arg(required = true, value_name = "STEP")]
        inputs: Vec<PathBuf>,
    }

    #[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
    #[serde(rename_all = "snake_case")]
    enum Outcome {
        InputError,
        NoSolids,
        CompleteConstructive,
        PartialConstructive,
        SemanticOnly,
        Unrecovered,
    }

    #[derive(Debug, Serialize)]
    struct SourceIdentity {
        data_section: usize,
        ordinal: usize,
        id: u64,
    }

    #[derive(Debug, Serialize)]
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
        families: BTreeMap<&'static str, usize>,
        unrecovered_sources: Vec<SourceIdentity>,
        kernel_failures: Vec<KernelRecoveryFailure>,
        error: Option<String>,
    }

    #[derive(Debug, Default, Serialize)]
    struct Totals {
        inputs: usize,
        inputs_with_errors: usize,
        files_with_no_solids: usize,
        fully_constructive_files: usize,
        partially_constructive_files: usize,
        semantic_only_files: usize,
        unrecovered_files: usize,
        source_solids: usize,
        constructively_recovered_solids: usize,
        semantically_recovered_solids: usize,
        unrecovered_solids: usize,
    }

    #[derive(Debug, Serialize)]
    struct Survey {
        schema: &'static str,
        /// Recovery counts are not independent proof of source-equivalent geometry.
        claim: &'static str,
        totals: Totals,
        cases: Vec<Case>,
    }

    fn classify(body_count: usize, constructive: usize, semantic: usize) -> Outcome {
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

    fn count_families(report: &CompleteIr) -> BTreeMap<&'static str, usize> {
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
                    .map(|solid| SourceIdentity {
                        data_section: solid.source.data_section,
                        ordinal: solid.source.ordinal,
                        id: solid.source.id,
                    })
                    .collect();
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
            self.source_solids += entry.body_count.unwrap_or_default();
            self.constructively_recovered_solids +=
                entry.constructively_recovered.unwrap_or_default();
            self.semantically_recovered_solids += entry.semantically_recovered.unwrap_or_default();
            self.unrecovered_solids += entry.unrecovered.unwrap_or_default();
        }
    }

    fn run(cli: Cli) -> Result<()> {
        let kernel = MonstertruckKernel::project_standard();
        let mut totals = Totals::default();
        let cases = cli
            .inputs
            .iter()
            .map(|path| {
                let entry = recover_case(path, &kernel);
                totals.include(&entry);
                entry
            })
            .collect();
        let report = Survey {
            schema: "step-redox.recovery-survey.v1",
            claim: "automatic constructive IR selection; independent source-boundary equivalence not measured",
            totals,
            cases,
        };
        let mut bytes = serde_json::to_vec_pretty(&report)?;
        bytes.push(b'\n');
        if let Some(output) = cli.output {
            std::fs::write(&output, bytes).with_context(|| format!("write {}", output.display()))
        } else {
            std::io::stdout().lock().write_all(&bytes)?;
            Ok(())
        }
    }

    pub fn main() -> Result<()> { run(Cli::parse()) }

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

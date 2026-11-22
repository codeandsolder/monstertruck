use anyhow::{Context, Result, ensure};
use clap::Parser;
use serde::{Deserialize, de::DeserializeOwned};
use serde_json::{Value, json};
use sha2::{Digest, Sha256};
use std::{
    collections::BTreeMap,
    fs,
    path::{Path, PathBuf},
    process::Command,
};

#[derive(Parser)]
#[command(about = "Validate independently specified STEP recovery expectations against CLI output")]
struct Cli {
    #[arg(long)]
    manifest: PathBuf,
    #[arg(long)]
    cache: PathBuf,
    #[arg(long)]
    out: PathBuf,
    #[arg(long)]
    step_redox: PathBuf,
}

#[derive(Deserialize)]
struct Manifest {
    fixtures: BTreeMap<String, Fixture>,
}

#[derive(Deserialize)]
struct Fixture {
    filename: String,
    sha256: String,
    periodic_chain_expectation: Option<ChainExpectation>,
    instance_pattern_expectation: Option<PatternExpectation>,
    formed_sheet_expectation: Option<SheetExpectation>,
}

#[derive(Debug, Deserialize, PartialEq, Eq)]
struct ChainFacts {
    read_only_proven: bool,
    complete_partition: bool,
    faces_without_geometry: usize,
    nonmanifold_edges: usize,
    cross_site_edges: usize,
    interior_site_face_count: usize,
    interior_gap_face_count: usize,
}

#[derive(Deserialize)]
struct ChainExpectation {
    sites: usize,
    pitch_mm: f64,
    stretch_faces: usize,
    #[serde(flatten)]
    facts: ChainFacts,
}

#[derive(Deserialize)]
struct Chain {
    sites: usize,
    pitch_mm: f64,
    #[serde(flatten)]
    facts: ChainFacts,
    lattice_face_family_votes: usize,
    stretch_face_ids: Vec<u64>,
    fixed_negative_face_ids: Vec<u64>,
    fixed_positive_face_ids: Vec<u64>,
    repeat_coverage_ratio: f64,
}

#[derive(Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case")]
enum FragmentSource {
    PeriodicChain {
        solid_id: u64,
        sites: usize,
    },
    #[serde(other)]
    Other,
}

#[derive(Deserialize)]
struct Fragment {
    source: FragmentSource,
    model: FragmentModel,
    root: usize,
}

#[derive(Deserialize)]
struct FragmentModel {
    nodes: Vec<Value>,
    provenance: BTreeMap<usize, Provenance>,
}

#[derive(Deserialize)]
struct Provenance {
    proof: String,
    max_residual_mm: Value,
}

#[derive(Deserialize)]
struct ChainNode {
    source_solid_id: u64,
    sites: usize,
    pitch_mm: f64,
    interior_site_face_count: usize,
    interior_gap_face_count: usize,
    stretch_face_count: usize,
    fixed_negative_face_count: usize,
    fixed_positive_face_count: usize,
    site_face_count_overrides: Value,
    gap_face_count_overrides: Value,
}

#[derive(Deserialize)]
struct PatternExpectation {
    linear_patterns: usize,
    instances_per_pattern: usize,
    step_mm: [f64; 3],
    max_residual_mm: f64,
}

#[derive(Deserialize)]
struct Pattern {
    dimension: u8,
    item_ids: Vec<u64>,
    basis: Vec<[f64; 3]>,
    pitch: Vec<f64>,
    max_residual_mm: f64,
    representation_map: u64,
}

#[derive(Debug, Deserialize, PartialEq, Eq)]
struct SheetCounts {
    cylindrical_faces: u32,
    paired_cylindrical_faces: u32,
    coaxial_radius_pairs: u32,
    parallel_plane_pairs: u32,
}

#[derive(Deserialize)]
struct SheetExpectation {
    candidates: usize,
    thickness_mm: f64,
    #[serde(flatten)]
    counts: SheetCounts,
}

#[derive(Deserialize)]
struct Sheet {
    solid_id: u64,
    thickness_mm: f64,
    paired_cylinder_face_ratio: f64,
    #[serde(flatten)]
    counts: SheetCounts,
}

fn read_json<T: DeserializeOwned>(path: &Path) -> Result<T> {
    let bytes = fs::read(path).with_context(|| format!("read {}", path.display()))?;
    serde_json::from_slice(&bytes).with_context(|| format!("parse {}", path.display()))
}

fn close(a: f64, b: f64, tolerance: f64) -> bool {
    a.is_finite() && b.is_finite() && (a - b).abs() <= tolerance
}

fn verify_fixture(input: &Path, expected: &str) -> Result<()> {
    let bytes =
        fs::read(input).with_context(|| format!("read cached fixture {}", input.display()))?;
    let actual = hex::encode(Sha256::digest(bytes));
    ensure!(
        actual == expected,
        "{}: fixture SHA-256 mismatch: expected {expected}, got {actual}",
        input.display()
    );
    Ok(())
}

fn run_probe(step_redox: &Path, input: &Path, prefix: &Path) -> Result<()> {
    let output = Command::new(step_redox)
        .args(["--profile", "compact"])
        .arg(input)
        .arg(prefix.with_extension("compact.step"))
        .arg("--patterns-json")
        .arg(prefix.with_extension("patterns.json"))
        .arg("--cad-fragments-json")
        .arg(prefix.with_extension("cad-fragments.json"))
        .arg("--formed-sheet-json")
        .arg(prefix.with_extension("formed-sheet.json"))
        .arg("--periodic-bodies-json")
        .arg(prefix.with_extension("bodies.json"))
        .arg("--count-parameters-json")
        .arg(prefix.with_extension("counts.json"))
        .arg("--periodic-chains-json")
        .arg(prefix.with_extension("chains.json"))
        .arg("--json")
        .output()
        .with_context(|| format!("execute {}", step_redox.display()))?;
    fs::write(prefix.with_extension("stdout.txt"), &output.stdout)?;
    fs::write(prefix.with_extension("stats.json"), &output.stderr)?;
    ensure!(
        output.status.success(),
        "{}: semantic probe failed ({}); see {}",
        input.display(),
        output.status,
        prefix.with_extension("stats.json").display()
    );
    Ok(())
}

fn validate_patterns(
    name: &str,
    patterns: &[Pattern],
    expected: &PatternExpectation,
) -> Result<Value> {
    let linear: Vec<_> = patterns
        .iter()
        .filter(|pattern| pattern.dimension == 1)
        .collect();
    ensure!(
        linear.len() == expected.linear_patterns,
        "{name}: expected {} linear patterns, got {} from {} total",
        expected.linear_patterns,
        linear.len(),
        patterns.len()
    );
    let mut observed = Vec::new();
    for (index, pattern) in linear.into_iter().enumerate() {
        ensure!(
            pattern.item_ids.len() == expected.instances_per_pattern,
            "{name}: pattern {index} incorrect instance count"
        );
        let ([step], [pitch]) = (pattern.basis.as_slice(), pattern.pitch.as_slice()) else {
            anyhow::bail!("{name}: pattern {index} inconsistent 1-D lattice metadata");
        };
        ensure!(
            step.iter()
                .zip(expected.step_mm)
                .all(|(&a, b)| close(a, b, 1e-12)),
            "{name}: pattern {index} step {step:?} != {:?}",
            expected.step_mm
        );
        let step_norm = step
            .iter()
            .map(|component| component * component)
            .sum::<f64>()
            .sqrt();
        ensure!(
            close(*pitch, step_norm, 1e-12),
            "{name}: pattern {index} pitch {pitch} != basis norm {step_norm}"
        );
        ensure!(
            pattern.max_residual_mm.is_finite()
                && pattern.max_residual_mm >= 0.0
                && pattern.max_residual_mm <= expected.max_residual_mm,
            "{name}: pattern {index} residual {} outside expected bound {}",
            pattern.max_residual_mm,
            expected.max_residual_mm
        );
        observed.push(json!({
            "count": pattern.item_ids.len(),
            "step_mm": step,
            "max_residual_mm": pattern.max_residual_mm,
            "representation_map": pattern.representation_map,
        }));
    }
    Ok(json!({"fixture": name, "patterns": observed}))
}

fn validate_sheet(name: &str, sheets: &[Sheet], expected: &SheetExpectation) -> Result<Value> {
    ensure!(
        sheets.len() == expected.candidates,
        "{name}: expected {} formed-sheet candidates, got {}",
        expected.candidates,
        sheets.len()
    );
    ensure!(
        expected.counts.cylindrical_faces > 0,
        "{name}: invalid expected cylindrical face count"
    );
    let expected_ratio = f64::from(expected.counts.paired_cylindrical_faces)
        / f64::from(expected.counts.cylindrical_faces);
    let mut observed = Vec::new();
    for (index, candidate) in sheets.iter().enumerate() {
        ensure!(
            candidate.counts == expected.counts,
            "{name}: formed-sheet candidate {index} counts {:?} != {:?}",
            candidate.counts,
            expected.counts
        );
        ensure!(
            close(candidate.thickness_mm, expected.thickness_mm, 1e-9),
            "{name}: candidate {index} thickness {} != {}",
            candidate.thickness_mm,
            expected.thickness_mm
        );
        ensure!(
            close(candidate.paired_cylinder_face_ratio, expected_ratio, 1e-12),
            "{name}: candidate {index} paired ratio {} != {expected_ratio}",
            candidate.paired_cylinder_face_ratio
        );
        observed.push(json!({
            "solid_id": candidate.solid_id,
            "thickness_mm": candidate.thickness_mm,
            "paired_cylinder_face_ratio": candidate.paired_cylinder_face_ratio,
            "coaxial_radius_pairs": candidate.counts.coaxial_radius_pairs,
            "parallel_plane_pairs": candidate.counts.parallel_plane_pairs,
        }));
    }
    Ok(json!({"fixture": name, "candidates": observed}))
}

fn validate_cad_fragment(
    name: &str,
    fragments: &[Fragment],
    expected: &ChainExpectation,
    chain: &Chain,
) -> Result<Value> {
    let mut candidates = fragments
        .iter()
        .filter_map(|fragment| match fragment.source {
            FragmentSource::PeriodicChain { solid_id, sites } if sites == expected.sites => {
                Some((fragment, solid_id))
            }
            _ => None,
        });
    let (fragment, source_solid_id) = candidates.next().with_context(|| {
        format!(
            "{name}: missing {}-site periodic CAD fragment",
            expected.sites
        )
    })?;
    ensure!(
        candidates.next().is_none(),
        "{name}: ambiguous periodic CAD fragment"
    );
    let value = fragment
        .model
        .nodes
        .get(fragment.root)
        .and_then(|node| node.get("PeriodicChain"))
        .with_context(|| format!("{name}: CAD root is not PeriodicChain"))?;
    let node = ChainNode::deserialize(value)?;
    ensure!(
        node.source_solid_id == source_solid_id,
        "{name}: CAD source solid mismatch"
    );
    ensure!(
        node.sites == expected.sites,
        "{name}: CAD site count mismatch"
    );
    ensure!(
        close(node.pitch_mm, expected.pitch_mm, 1e-9),
        "{name}: CAD pitch mismatch"
    );
    ensure!(
        node.interior_site_face_count == expected.facts.interior_site_face_count
            && node.interior_gap_face_count == expected.facts.interior_gap_face_count,
        "{name}: CAD interior face count mismatch"
    );
    ensure!(
        node.stretch_face_count == chain.stretch_face_ids.len()
            && node.fixed_negative_face_count == chain.fixed_negative_face_ids.len()
            && node.fixed_positive_face_count == chain.fixed_positive_face_ids.len(),
        "{name}: CAD terminal face count mismatch"
    );
    let provenance = fragment
        .model
        .provenance
        .get(&fragment.root)
        .with_context(|| format!("{name}: missing CAD root provenance"))?;
    ensure!(
        provenance.proof == "StructurallyProven",
        "{name}: expected structural CAD proof, got {}",
        provenance.proof
    );
    ensure!(
        provenance.max_residual_mm.is_null(),
        "{name}: structural proof must not invent a residual"
    );
    Ok(json!({
        "source_solid_id": node.source_solid_id,
        "sites": node.sites,
        "pitch_mm": node.pitch_mm,
        "site_face_count_overrides": node.site_face_count_overrides,
        "gap_face_count_overrides": node.gap_face_count_overrides,
    }))
}

#[derive(Debug, PartialEq, Eq)]
struct FamilyGrammar {
    slope: usize,
    fixed_negative_faces: usize,
    fixed_positive_faces: usize,
    stretch_faces: usize,
}

struct ChainReport {
    pitch_mm: f64,
    grammar: FamilyGrammar,
    summary: Value,
}

fn validate_chain(
    name: &str,
    chains: &[Chain],
    fragments: &[Fragment],
    expected: &ChainExpectation,
) -> Result<ChainReport> {
    let mut candidates = chains.iter().filter(|chain| {
        chain.sites == expected.sites && close(chain.pitch_mm, expected.pitch_mm, 1e-8)
    });
    let chain = candidates.next().with_context(|| {
        format!(
            "{name}: missing {}-site/{} mm chain",
            expected.sites, expected.pitch_mm
        )
    })?;
    ensure!(
        candidates.next().is_none(),
        "{name}: ambiguous periodic chain"
    );
    ensure!(
        chain.facts == expected.facts,
        "{name}: chain topology/proof {:?} != {:?}",
        chain.facts,
        expected.facts
    );
    ensure!(
        chain.stretch_face_ids.len() == expected.stretch_faces,
        "{name}: incorrect stretch face count"
    );
    let grammar = FamilyGrammar {
        slope: chain.facts.interior_site_face_count + chain.facts.interior_gap_face_count,
        fixed_negative_faces: chain.fixed_negative_face_ids.len(),
        fixed_positive_faces: chain.fixed_positive_face_ids.len(),
        stretch_faces: chain.stretch_face_ids.len(),
    };
    let cad_ir = validate_cad_fragment(name, fragments, expected, chain)?;
    let summary = json!({
        "fixture": name,
        "sites": chain.sites,
        "pitch_mm": chain.pitch_mm,
        "lattice_face_family_votes": chain.lattice_face_family_votes,
        "interior_site_face_count": chain.facts.interior_site_face_count,
        "interior_gap_face_count": chain.facts.interior_gap_face_count,
        "faces_per_added_site": grammar.slope,
        "stretch_faces": grammar.stretch_faces,
        "fixed_negative_faces": grammar.fixed_negative_faces,
        "fixed_positive_faces": grammar.fixed_positive_faces,
        "repeat_coverage_ratio": chain.repeat_coverage_ratio,
        "read_only_proven": chain.facts.read_only_proven,
        "cad_ir": cad_ir,
    });
    Ok(ChainReport {
        pitch_mm: chain.pitch_mm,
        grammar,
        summary,
    })
}

fn validate_family(reports: &[ChainReport]) -> Result<()> {
    if let Some((first, rest)) = reports.split_first() {
        for report in rest {
            // Count variants must share pitch to nine decimal places, per-site
            // topology law and terminal grammar, not absolute fixed-face counts.
            ensure!(
                (report.pitch_mm * 1e9)
                    .round_ties_even()
                    .total_cmp(&(first.pitch_mm * 1e9).round_ties_even())
                    .is_eq()
                    && report.grammar == first.grammar,
                "periodic family mismatch: pitch {} and {:?} vs {} and {:?}",
                report.pitch_mm,
                report.grammar,
                first.pitch_mm,
                first.grammar
            );
        }
    }
    Ok(())
}

fn write_summary(out: &Path, filename: &str, results: &[Value]) -> Result<()> {
    if !results.is_empty() {
        let mut bytes = serde_json::to_vec_pretty(results)?;
        bytes.push(b'\n');
        fs::write(out.join(filename), bytes)?;
    }
    Ok(())
}

fn run(cli: &Cli) -> Result<()> {
    let manifest: Manifest = read_json(&cli.manifest)?;
    fs::create_dir_all(&cli.out)?;
    let mut chains = Vec::new();
    let mut patterns = Vec::new();
    let mut sheets = Vec::new();
    for (name, spec) in manifest.fixtures {
        if spec.periodic_chain_expectation.is_none()
            && spec.instance_pattern_expectation.is_none()
            && spec.formed_sheet_expectation.is_none()
        {
            continue;
        }
        let input = cli.cache.join(spec.filename);
        verify_fixture(&input, &spec.sha256)?;
        let prefix = cli.out.join(&name);
        run_probe(&cli.step_redox, &input, &prefix)?;
        if let Some(expected) = spec.periodic_chain_expectation {
            let detected: Vec<Chain> = read_json(&prefix.with_extension("chains.json"))?;
            let fragments: Vec<Fragment> = read_json(&prefix.with_extension("cad-fragments.json"))?;
            chains.push(validate_chain(&name, &detected, &fragments, &expected)?);
        }
        if let Some(expected) = spec.instance_pattern_expectation {
            let detected: Vec<Pattern> = read_json(&prefix.with_extension("patterns.json"))?;
            patterns.push(validate_patterns(&name, &detected, &expected)?);
        }
        if let Some(expected) = spec.formed_sheet_expectation {
            let detected: Vec<Sheet> = read_json(&prefix.with_extension("formed-sheet.json"))?;
            sheets.push(validate_sheet(&name, &detected, &expected)?);
        }
        println!("{name}: semantic expectations passed");
    }
    ensure!(
        !chains.is_empty() || !patterns.is_empty() || !sheets.is_empty(),
        "manifest contains no semantic expectations"
    );
    validate_family(&chains)?;
    let summaries: Vec<_> = chains.into_iter().map(|report| report.summary).collect();
    write_summary(&cli.out, "periodic-chain-summary.json", &summaries)?;
    write_summary(&cli.out, "cad-fragment-summary.json", &patterns)?;
    write_summary(&cli.out, "formed-sheet-summary.json", &sheets)?;
    Ok(())
}

fn main() -> Result<()> { run(&Cli::parse()) }

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn semantic_report_preserves_binary64_roundtrip() -> Result<()> {
        // Actual TE fixture coverage was decoded one ULP high by the default
        // fast JSON parser, even though the CLI had emitted the right value.
        let decimal = "0.9267461669505963";
        let decoded: f64 = serde_json::from_str(decimal)?;
        let expected = decimal.parse::<f64>()?;
        assert_eq!(decoded.to_bits(), expected.to_bits());
        let encoded = serde_json::to_string(&decoded)?;
        assert_eq!(
            serde_json::from_str::<f64>(&encoded)?.to_bits(),
            expected.to_bits()
        );
        Ok(())
    }

    #[test]
    fn patterns_require_consistent_pitch_and_exact_step_dimensions() -> Result<()> {
        let expected: PatternExpectation = serde_json::from_value(json!({
            "linear_patterns": 1, "instances_per_pattern": 2,
            "step_mm": [2.0, 0.0, 0.0], "max_residual_mm": 1e-7,
        }))?;
        let mut pattern: Pattern = serde_json::from_value(json!({
            "dimension": 1, "item_ids": [1, 2], "basis": [[2.0, 0.0, 0.0]],
            "pitch": [2.0], "max_residual_mm": 0.0, "representation_map": 1,
        }))?;
        assert!(validate_patterns("fixture", std::slice::from_ref(&pattern), &expected).is_ok());
        pattern.pitch[0] = 4.0;
        assert!(validate_patterns("fixture", &[pattern], &expected).is_err());
        assert!(
            serde_json::from_value::<PatternExpectation>(json!({
                "linear_patterns": 1, "instances_per_pattern": 2,
                "step_mm": [2.0], "max_residual_mm": 1e-7,
            }))
            .is_err()
        );
        Ok(())
    }

    #[test]
    fn sheet_checks_both_counts_and_ratio() -> Result<()> {
        let expected: SheetExpectation = serde_json::from_value(json!({
            "candidates": 1, "thickness_mm": 0.2, "cylindrical_faces": 20,
            "paired_cylindrical_faces": 16, "coaxial_radius_pairs": 8, "parallel_plane_pairs": 9,
        }))?;
        let mut sheet: Sheet = serde_json::from_value(json!({
            "solid_id": 1, "thickness_mm": 0.2, "cylindrical_faces": 20,
            "paired_cylindrical_faces": 16, "coaxial_radius_pairs": 8, "parallel_plane_pairs": 9,
            "paired_cylinder_face_ratio": 0.8,
        }))?;
        assert!(validate_sheet("fixture", std::slice::from_ref(&sheet), &expected).is_ok());
        sheet.paired_cylinder_face_ratio = 0.5;
        assert!(validate_sheet("fixture", &[sheet], &expected).is_err());
        Ok(())
    }

    #[test]
    fn family_detects_terminal_and_pitch_drift() {
        let report = |pitch_mm, fixed_positive_faces| ChainReport {
            pitch_mm,
            grammar: FamilyGrammar {
                slope: 73,
                fixed_negative_faces: 3,
                fixed_positive_faces,
                stretch_faces: 6,
            },
            summary: Value::Null,
        };
        assert!(validate_family(&[report(1.5, 3), report(1.5, 3)]).is_ok());
        assert!(validate_family(&[report(1.5, 3), report(1.51, 3)]).is_err());
        assert!(validate_family(&[report(1.5, 3), report(1.5, 4)]).is_err());
    }

    #[test]
    fn structural_cad_proof_cannot_claim_residual_or_ambiguous_root() -> Result<()> {
        let expected: ChainExpectation = serde_json::from_value(json!({
            "sites": 14, "pitch_mm": 1.5, "stretch_faces": 1,
            "interior_site_face_count": 70, "interior_gap_face_count": 3,
            "read_only_proven": true, "complete_partition": true,
            "faces_without_geometry": 0, "nonmanifold_edges": 0, "cross_site_edges": 0,
        }))?;
        let chain: Chain = serde_json::from_value(json!({
            "sites": 14, "pitch_mm": 1.5, "interior_site_face_count": 70, "interior_gap_face_count": 3,
            "read_only_proven": true, "complete_partition": true,
            "faces_without_geometry": 0, "nonmanifold_edges": 0, "cross_site_edges": 0,
            "lattice_face_family_votes": 20, "stretch_face_ids": [1],
            "fixed_negative_face_ids": [2], "fixed_positive_face_ids": [3], "repeat_coverage_ratio": 0.8,
        }))?;
        let value = json!({
            "source": {"kind": "periodic_chain", "solid_id": 42, "sites": 14},
            "root": 0,
            "model": {
                "nodes": [{"PeriodicChain": {
                    "source_solid_id": 42, "sites": 14, "pitch_mm": 1.5,
                    "interior_site_face_count": 70, "interior_gap_face_count": 3,
                    "stretch_face_count": 1, "fixed_negative_face_count": 1, "fixed_positive_face_count": 1,
                    "site_face_count_overrides": [], "gap_face_count_overrides": [],
                }}],
                "provenance": {"0": {"proof": "StructurallyProven", "max_residual_mm": null}},
            },
        });
        let mut fragment: Fragment = serde_json::from_value(value)?;
        assert!(
            validate_chain(
                "fixture",
                std::slice::from_ref(&chain),
                std::slice::from_ref(&fragment),
                &expected
            )
            .is_ok()
        );
        let provenance = fragment
            .model
            .provenance
            .get_mut(&0)
            .context("test provenance")?;
        provenance.max_residual_mm = json!(0.0);
        assert!(validate_chain("fixture", &[chain], &[fragment], &expected).is_err());
        Ok(())
    }
}

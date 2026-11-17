//! Manual Part 21 load profiler; filesystem I/O is outside measured sections.
//!
//! cargo +1.100.0 run -p monstertruck-io --example profile-step-parse --release --offline
//! cargo +1.100.0 run -p monstertruck-io --example profile-step-parse --release --offline -- --fingerprint-corpus
//!
//! Timing depends on host load; no time-based assertions are made.

use anyhow::{Context, Result, ensure};
use monstertruck_io::step::load::{ParsedStepSource, decode_step_input};
use std::hint::black_box;
use std::time::{Duration, Instant};

fn median_duration(runs: usize, mut measure: impl FnMut() -> Result<()>) -> Result<Duration> {
    ensure!(runs > 0, "benchmark needs at least one run");
    for _ in 0..3 {
        measure()?;
    }
    let mut samples = Vec::with_capacity(runs);
    for _ in 0..runs {
        let start = Instant::now();
        measure()?;
        samples.push(start.elapsed());
    }
    samples.sort_unstable();
    Ok(samples[runs / 2])
}

fn profile(name: &str, bytes: &[u8], runs: usize) -> Result<()> {
    let (decoded, encoding) = decode_step_input(bytes);
    let exchange = step_p21::parser::parse(&decoded)?;
    let section_counts = exchange
        .data
        .iter()
        .map(|section| section.entities.len())
        .collect::<Vec<_>>();
    let cold = ParsedStepSource::parse(bytes)?;
    ensure!(
        cold.exchange == exchange,
        "full loader changed the parsed AST"
    );
    ensure!(cold.input_encoding == encoding, "input encoding changed");

    let decode = median_duration(runs, || {
        black_box(decode_step_input(black_box(bytes)));
        Ok(())
    })?;
    let parse = median_duration(runs, || {
        black_box(step_p21::parser::parse(black_box(&decoded))?);
        Ok(())
    })?;
    let wrapper = median_duration(runs, || {
        black_box(ParsedStepSource::parse(black_box(bytes))?);
        Ok(())
    })?;
    let input_bytes = bytes.len();
    let sections = section_counts.len();
    let decode_median_us = decode.as_micros();
    let step_p21_median_us = parse.as_micros();
    let wrapper_median_us = wrapper.as_micros();
    println!(
        "{name}: input_bytes={input_bytes}, sections={sections}, entities={section_counts:?}, encoding={encoding}, runs={runs}, decode_median_us={decode_median_us}, step_p21_median_us={step_p21_median_us}, wrapper_median_us={wrapper_median_us}"
    );
    Ok(())
}

/// Compare complete parsed ASTs across dependency versions. These Debug hashes
/// are comparable only with the same Rust toolchain and fixture content; they
/// are not a stable interchange format or a cryptographic integrity check.
fn fingerprint_corpus() -> Result<()> {
    use std::hash::{Hash, Hasher};

    let directory =
        std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../step-redox/validation/fixtures");
    let mut files = std::fs::read_dir(directory)?
        .map(|entry| entry.map(|item| item.path()))
        .collect::<std::io::Result<Vec<_>>>()?;
    files.retain(|path| path.extension().is_some_and(|ext| ext == "step"));
    files.sort_unstable();
    for path in files {
        let filename = path
            .file_name()
            .context("STEP fixture path has no file name")?
            .to_string_lossy();
        let bytes = std::fs::read(&path)?;
        match ParsedStepSource::parse(&bytes) {
            Ok(document) => {
                let mut hash = std::collections::hash_map::DefaultHasher::new();
                format!("{:?}", document.exchange).hash(&mut hash);
                let count = document
                    .exchange
                    .data
                    .iter()
                    .map(|section| section.entities.len())
                    .sum::<usize>();
                println!(
                    "{filename}\tOK\t{}\t{}\t{count}\t{:016x}",
                    bytes.len(),
                    document.exchange.data.len(),
                    hash.finish(),
                );
            }
            Err(_) => println!("{filename}\tERROR\t{}", bytes.len()),
        }
    }
    Ok(())
}

fn main() -> Result<()> {
    if std::env::args().any(|arg| arg == "--fingerprint-corpus") {
        return fingerprint_corpus();
    }
    let runs = if cfg!(debug_assertions) { 10 } else { 64 };
    profile(
        "production_radial_slot",
        include_bytes!(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../step-redox/validation/fixtures/production_radial_slot.step"
        )),
        runs,
    )?;
    profile(
        "native_split_ring_hollow_cylinder",
        include_bytes!(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../step-redox/validation/fixtures/native_split_ring_hollow_cylinder.step"
        )),
        runs,
    )?;
    profile(
        "nested_boundary_boss",
        include_bytes!(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../step-redox/validation/fixtures/boundary_boss_nested.step"
        )),
        runs,
    )?;
    Ok(())
}

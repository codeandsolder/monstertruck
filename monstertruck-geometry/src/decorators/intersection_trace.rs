//! Opt-in diagnostics for identifying excessive exact SSI re-evaluation.
use std::{cell::RefCell, sync::OnceLock};

#[derive(Default)]
struct Counts {
    evaluate: u64,
    triple: u64,
    closest: u64,
    projection: u64,
    curve_divisions: u64,
    curve_division_samples: u64,
    leader_vertices: u64,
    biggest: (u64, u64, u64),
}

impl Drop for Counts {
    fn drop(&mut self) {
        if enabled() {
            eprintln!(
                "SSI_WORKLOAD evaluate={} search_triple={} nearest={} double_projection={} divisions={} divided_samples={} leader_vertices={} biggest_division=(leader={},out={},tol_bits={})",
                self.evaluate,
                self.triple,
                self.closest,
                self.projection,
                self.curve_divisions,
                self.curve_division_samples,
                self.leader_vertices,
                self.biggest.0,
                self.biggest.1,
                self.biggest.2
            );
        }
    }
}

thread_local! {
    static COUNTS: RefCell<Counts> = RefCell::new(Counts::default());
}

fn enabled() -> bool {
    static ENABLED: OnceLock<bool> = OnceLock::new();
    *ENABLED.get_or_init(|| std::env::var_os("MT_SSI_WORKLOAD").is_some())
}

pub(super) fn mark(kind: u8) {
    if !enabled() {
        return;
    }
    COUNTS.with(|metrics| {
        let mut metrics = metrics.borrow_mut();
        match kind {
            0 => metrics.evaluate += 1,
            1 => metrics.triple += 1,
            2 => metrics.closest += 1,
            _ => metrics.projection += 1,
        }
    });
}

pub(super) fn division(leader: usize, samples: usize, tol: f64) {
    if !enabled() {
        return;
    }
    COUNTS.with(|metrics| {
        let mut metrics = metrics.borrow_mut();
        metrics.curve_divisions += 1;
        metrics.curve_division_samples += samples as u64;
        metrics.leader_vertices += leader as u64;
        if samples as u64 > metrics.biggest.1 {
            metrics.biggest = (leader as u64, samples as u64, tol.to_bits());
        }
        if samples > 2048 {
            eprintln!(
                "SSI_DIVISION_HEAVY leader={} samples={} tol={}",
                leader, samples, tol
            );
        }
    });
}

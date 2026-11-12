//! Opt-in workload diagnostics for real NURBS inverse searches.

use std::{
    cell::RefCell,
    sync::OnceLock,
    time::{Duration, Instant},
};

#[derive(Default)]
struct Metrics {
    queries: u64,
    hits: u64,
    misses: u64,
    indexed: u64,
    kd_builds: u64,
    gpu_queries: u64,
    gpu_mismatches: u64,
    gpu_time: Duration,
    seed_time: Duration,
    refined: u64,
    refinement_time: Duration,
}

impl Metrics {
    fn print(&self, finished: bool) {
        eprintln!(
            "NURBS_PRESEARCH_METRICS final={} queries={} hits={} misses={} indexed={} kd_builds={} seed_ms={:.3} refinements={} refine_ms={:.3} gpu_queries={} gpu_mismatches={} gpu_ms={:.3}",
            finished,
            self.queries,
            self.hits,
            self.misses,
            self.indexed,
            self.kd_builds,
            self.seed_time.as_secs_f64() * 1000.0,
            self.refined,
            self.refinement_time.as_secs_f64() * 1000.0,
            self.gpu_queries,
            self.gpu_mismatches,
            self.gpu_time.as_secs_f64() * 1000.0,
        );
    }
}

impl Drop for Metrics {
    fn drop(&mut self) {
        if enabled() && self.queries > 0 {
            self.print(true);
        }
    }
}

thread_local! {
    static METRICS: RefCell<Metrics> = RefCell::new(Metrics::default());
}

pub(super) fn enabled() -> bool {
    static ENABLED: OnceLock<bool> = OnceLock::new();
    *ENABLED.get_or_init(|| std::env::var_os("MONSTERTRUCK_PRESEARCH_TRACE").is_some())
}

pub(super) fn begin() -> Option<Instant> { enabled().then(Instant::now) }

pub(super) fn seed(hit: bool, indexed: bool, built: bool, started: Option<Instant>) {
    if let Some(started) = started {
        METRICS.with(|state| {
            let mut state = state.borrow_mut();
            state.queries += 1;
            state.hits += u64::from(hit);
            state.misses += u64::from(!hit);
            state.indexed += u64::from(indexed);
            state.kd_builds += u64::from(built);
            state.seed_time += started.elapsed();
            if state.queries.is_multiple_of(25_000) {
                state.print(false);
            }
        });
    }
}

pub(super) fn refine(started: Option<Instant>) {
    if let Some(started) = started {
        METRICS.with(|state| {
            let mut state = state.borrow_mut();
            state.refined += 1;
            state.refinement_time += started.elapsed();
        });
    }
}

pub(super) fn gpu(mismatch: bool, started: Option<Instant>) {
    if let Some(started) = started {
        METRICS.with(|state| {
            let mut state = state.borrow_mut();
            state.gpu_queries += 1;
            state.gpu_mismatches += u64::from(mismatch);
            state.gpu_time += started.elapsed();
        });
    }
}

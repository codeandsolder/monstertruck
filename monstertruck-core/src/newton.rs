//! Implementation of Newton method

use std::ops::{Mul, Sub};

use crate::{
    cgmath64::{One, SquareMatrix, Zero, cgmath},
    tolerance::Tolerance,
};

/// the value and jacobian corresponding to parameter
#[derive(Clone, Debug)]
pub struct CalcOutput<V, M> {
    /// the value of function
    pub value: V,
    /// the jacobian of function
    pub derivation: M,
}

/// jacobian of function
pub trait Jacobian<V>: Mul<V, Output = V> + Sized {
    #[doc(hidden)]
    fn invert(self) -> Option<Self>;
}

/// Scalar-generic Jacobian for 1D (scalar as its own jacobian).
impl<S> Jacobian<S> for S
where S: Mul<S, Output = S> + Zero + One + Copy + PartialEq + std::ops::Div<Output = S>
{
    #[inline]
    fn invert(self) -> Option<Self> {
        if self == S::zero() {
            None
        } else {
            Some(S::one() / self)
        }
    }
}

/// Scalar-generic matrix Jacobians. Requires `S: BaseFloat` for
/// `SquareMatrix::invert`.
macro_rules! impl_jacobian {
    ($matrix: ident, $vector: ident) => {
        impl<S: cgmath::BaseFloat> Jacobian<cgmath::$vector<S>> for cgmath::$matrix<S> {
            #[inline]
            fn invert(self) -> Option<Self> { SquareMatrix::invert(&self) }
        }
    };
}

impl_jacobian!(Matrix2, Vector2);
impl_jacobian!(Matrix3, Vector3);
impl_jacobian!(Matrix4, Vector4);

/// Opt-in iteration histogram for diagnosing expensive CAD projections.
///
/// The solver's success/failure and precision are not changed by collecting
/// the histogram. Set MT_NEWTON_PROFILE=1 when replaying a real CAD operation.
mod profiling {
    use std::{any::type_name, cell::RefCell, collections::BTreeMap, sync::OnceLock};

    #[derive(Default)]
    struct Count {
        calls: u64,
        converged: u64,
        singular: u64,
        capped: u64,
        iterations: u64,
        histogram: [u64; 8],
    }

    #[derive(Default)]
    struct Profile(BTreeMap<&'static str, Count>);

    impl Drop for Profile {
        fn drop(&mut self) {
            if enabled() {
                for (dimension, s) in &self.0 {
                    eprintln!(
                        "NEWTON_PROFILE vector={} calls={} success={} singular={} cap={} iterations={} bins_1_2_4_8_16_32_64_more={:?}",
                        dimension,
                        s.calls,
                        s.converged,
                        s.singular,
                        s.capped,
                        s.iterations,
                        s.histogram
                    );
                }
            }
        }
    }

    thread_local! {
        static STATS: RefCell<Profile> = RefCell::new(Profile::default());
    }

    fn enabled() -> bool {
        static ENABLED: OnceLock<bool> = OnceLock::new();
        *ENABLED.get_or_init(|| std::env::var_os("MT_NEWTON_PROFILE").is_some())
    }

    pub(super) fn record<V>(steps: usize, outcome: u8) {
        if !enabled() {
            return;
        }
        STATS.with(|stats| {
            let mut stats = stats.borrow_mut();
            let count = stats.0.entry(type_name::<V>()).or_default();
            count.calls += 1;
            count.iterations += steps as u64;
            match outcome {
                0 => count.converged += 1,
                1 => count.singular += 1,
                _ => count.capped += 1,
            }
            let bucket = if steps <= 1 {
                0
            } else {
                ((steps - 1).ilog2() + 1).min(7) as usize
            };
            count.histogram[bucket] += 1;
        });
    }
}

/// Solve equation by Newton's method
/// # Errors
/// Returns a [`NewtonLog`] when the Jacobian becomes singular or convergence is not reached within `trials`.
///
/// # Examples
/// ```
/// # fn main() -> anyhow::Result<()> {
/// use monstertruck_core::{newton::*, assert_near2};
///
/// let function = |x: f64| CalcOutput {
///     value: x * x - 2.0,
///     derivation: 2.0 * x,
/// };
/// let sqrt2 = solve(function, 1.0, 10).map_err(|e| anyhow::anyhow!("{e}"))?;
/// assert_near2!(sqrt2, f64::sqrt(2.0));
/// # Ok(())
/// # }
/// ```
pub fn solve<V, M>(
    function: impl Fn(V) -> CalcOutput<V, M>,
    mut hint: V,
    trials: usize,
) -> Result<V, NewtonLog<V>>
where
    V: Sub<Output = V> + Copy + Tolerance,
    M: Jacobian<V>,
{
    let mut log = NewtonLog::new(cfg!(debug_assertions), trials);
    for step in 0..=trials {
        log.push(hint);
        let CalcOutput { value, derivation } = function(hint);
        let Some(inv) = derivation.invert() else {
            log.set_degenerate(true);
            profiling::record::<V>(step + 1, 1);
            return Err(log);
        };
        let next = hint - inv * value;
        if next.near(&hint) {
            profiling::record::<V>(step + 1, 0);
            return Ok(next);
        }
        hint = next;
    }
    profiling::record::<V>(trials.saturating_add(1), 2);
    Err(log)
}

/// Scalar-generic Newton solver.
///
/// Like [`solve`] but uses [`ToleranceV2`](crate::scalar::ToleranceV2) instead
/// of the `f64`-bound [`Tolerance`] trait, allowing convergence checks with any
/// scalar whose epsilon implements [`ToleranceScalar`](crate::scalar::ToleranceScalar).
///
/// # Errors
/// Returns a [`NewtonLog`] when the Jacobian becomes singular or convergence is not reached within `trials`.
pub fn solve_v2<V, M>(
    function: impl Fn(V) -> CalcOutput<V, M>,
    mut hint: V,
    trials: usize,
) -> Result<V, NewtonLog<V>>
where
    V: Sub<Output = V> + Copy + crate::scalar::ToleranceV2,
    V::Epsilon: crate::scalar::ToleranceScalar,
    M: Jacobian<V>,
{
    let mut log = NewtonLog::new(cfg!(debug_assertions), trials);
    for step in 0..=trials {
        log.push(hint);
        let CalcOutput { value, derivation } = function(hint);
        let Some(inv) = derivation.invert() else {
            log.set_degenerate(true);
            profiling::record::<V>(step + 1, 1);
            return Err(log);
        };
        let next = hint - inv * value;
        if next.near_v2(&hint) {
            profiling::record::<V>(step + 1, 0);
            return Ok(next);
        }
        hint = next;
    }
    profiling::record::<V>(trials.saturating_add(1), 2);
    Err(log)
}

mod newtonlog {
    use std::fmt::{Debug, Display, Formatter, Result};
    /// A structure that stores logs for debugging.
    #[derive(Clone, Debug)]
    pub struct NewtonLog<T> {
        log: Option<Vec<T>>,
        degenerate: bool,
    }

    impl<T> NewtonLog<T> {
        /// constructor
        #[inline]
        #[must_use]
        pub fn new(activate: bool, trials: usize) -> Self {
            if activate {
                Self {
                    log: Some(Vec::with_capacity(trials)),
                    degenerate: false,
                }
            } else {
                Self {
                    log: None,
                    degenerate: false,
                }
            }
        }
        /// Returns `true` iff the Newton method terminates due to Jacobian degeneracy.
        #[inline]
        #[must_use]
        pub const fn degenerate(&self) -> bool { self.degenerate }
        #[inline]
        pub(super) fn push(&mut self, log: T) {
            if let Some(vec) = &mut self.log {
                vec.push(log);
            }
        }
        #[inline]
        pub(super) const fn set_degenerate(&mut self, degenerate: bool) {
            self.degenerate = degenerate;
        }
    }

    impl<T: Debug> Display for NewtonLog<T> {
        fn fmt(&self, f: &mut Formatter<'_>) -> Result {
            if self.degenerate {
                f.pad("Jacobian is degenerate. ")?;
            } else {
                f.pad("Newton method did not converge. ")?;
            }
            match &self.log {
                None => f.pad(
                    "If you want to see the Newton log, please re-run it with the debug build.",
                ),
                Some(vec) => {
                    f.pad("Newton Log:\n")?;
                    vec.iter()
                        .try_for_each(|log| f.write_fmt(format_args!("{log:?}\n")))
                }
            }
        }
    }
}
pub use newtonlog::NewtonLog;

//! All input dimensions produce consistent output prefixes. These checks
//! exercise the scalar stage, the recursively hashed vector stages, and
//! their f32/f64 coefficient conversions without pinning platform libm bits.
use monstertruck_core::hash::HashGen;

trait TestHashScalar: Copy + core::fmt::Debug {
    type Bits: Eq + core::fmt::Debug;

    fn bits(self) -> Self::Bits;
    fn in_unit_interval(self) -> bool;
}

impl TestHashScalar for f64 {
    type Bits = u64;
    fn bits(self) -> Self::Bits { self.to_bits() }
    fn in_unit_interval(self) -> bool { (0.0..1.0).contains(&self) }
}

impl TestHashScalar for f32 {
    type Bits = u32;
    fn bits(self) -> Self::Bits { self.to_bits() }
    fn in_unit_interval(self) -> bool { (0.0..1.0).contains(&self) }
}

fn check_prefix<S: TestHashScalar, const N: usize>(seed: [S; N])
where [S; N]: HashGen<S> {
    let h1 = <[S; N] as HashGen<S>>::hash1(seed);
    let h2 = <[S; N] as HashGen<S>>::hash2(seed);
    let h3 = <[S; N] as HashGen<S>>::hash3(seed);
    let h4 = <[S; N] as HashGen<S>>::hash4(seed);
    for (lhs, rhs) in [
        (h1, h4[0]),
        (h2[0], h4[0]),
        (h2[1], h4[1]),
        (h3[0], h4[0]),
        (h3[1], h4[1]),
        (h3[2], h4[2]),
    ] {
        assert_eq!(lhs.bits(), rhs.bits());
    }
    for value in h4 {
        assert!(value.in_unit_interval(), "unexpected hash value {value:?}");
    }
}

#[test]
fn output_width_preserves_every_component_prefix() {
    for (f64_values, f32_values) in [
        ([0.25, -1.5, 3.125, 0.75], [0.25_f32, -1.5, 3.125, 0.75]),
        (
            [-3.25, 0.01, 11.875, 0.125],
            [-3.25_f32, 0.01, 11.875, 0.125],
        ),
        ([1.0, 2.0, 3.0, 4.0], [1.0_f32, 2.0, 3.0, 4.0]),
    ] {
        check_prefix([f64_values[0]]);
        check_prefix([f64_values[0], f64_values[1]]);
        check_prefix([f64_values[0], f64_values[1], f64_values[2]]);
        check_prefix(f64_values);
        check_prefix([f32_values[0]]);
        check_prefix([f32_values[0], f32_values[1]]);
        check_prefix([f32_values[0], f32_values[1], f32_values[2]]);
        check_prefix(f32_values);
    }
}

/// Non-gating local diagnostic. Run with a release build and --ignored to
/// compare the old/new code under identical sccache/toolchain settings.
#[test]
#[ignore = "release-only throughput diagnostic, not a test assertion"]
fn release_hash4_throughput_diagnostic() {
    use std::hint::black_box;
    use std::time::Instant;

    const SAMPLES: usize = 200_000;
    let mut results = [0_u128; 5];
    for duration in &mut results {
        let started = Instant::now();
        let mut sum = 0.0_f64;
        let mut x = 0.0_f64;
        for _ in 0..SAMPLES {
            x += 0.000_173;
            let x = black_box(x);
            let result = <[f64; 4] as HashGen<f64>>::hash4([x, -x, x + 1.3, 2.1 - x]);
            sum += black_box(result[0]);
        }
        black_box(sum);
        *duration = started.elapsed().as_nanos();
    }
    results.sort_unstable();
    println!("hash4 f64 elapsed_ns {results:?}, median {}", results[2]);

    for duration in &mut results {
        let started = Instant::now();
        let mut sum = 0.0_f32;
        let mut x = 0.0_f32;
        for _ in 0..SAMPLES {
            x += 0.000_173;
            let x = black_box(x);
            let result = <[f32; 4] as HashGen<f32>>::hash4([x, -x, x + 1.3, 2.1 - x]);
            sum += black_box(result[0]);
        }
        black_box(sum);
        *duration = started.elapsed().as_nanos();
    }
    results.sort_unstable();
    println!("hash4 f32 elapsed_ns {results:?}, median {}", results[2]);
}

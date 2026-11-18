#![expect(
    clippy::many_single_char_names,
    reason = "hash mixing formulas use conventional coefficient and component symbols"
)]
use cgmath::num_traits::Float;

/// Scalar types for which the deterministic hash constants have an explicit,
/// total representation conversion.
trait HashFloat: Float {
    fn from_hash_constant(value: f64) -> Self;
}

impl HashFloat for f64 {
    #[inline]
    fn from_hash_constant(value: f64) -> Self { value }
}

impl HashFloat for f32 {
    #[inline]
    #[expect(
        clippy::cast_possible_truncation,
        reason = "f32 hashing deliberately uses the f32-rounded form of the fixed f64 mixing constants"
    )]
    fn from_hash_constant(value: f64) -> Self {
        // The f32 hash intentionally runs with f32-rounded constants. Rust's
        // float narrowing is total and preserves NaN/∞ semantics if this
        // private helper is ever reused with such a value.
        value as Self
    }
}

#[inline]
fn f64_to<S: HashFloat>(value: f64) -> S { S::from_hash_constant(value) }

/// Deterministic hash generator
pub trait HashGen<S> {
    /// deterministic hash, output 1-dim
    fn hash1(seed: Self) -> S;
    /// deterministic hash, output 2-dim
    fn hash2(seed: Self) -> [S; 2];
    /// deterministic hash, output 3-dim
    fn hash3(seed: Self) -> [S; 3];
    /// deterministic hash, output 4-dim
    fn hash4(seed: Self) -> [S; 4];
}

// Original hash mixing constants. Stage number is the output component.
const SCALAR_COEFFICIENTS: [[f64; 3]; 4] = [
    [
        61.909_685_033_545_934,
        8.436_303_256_302_796,
        220.678_620_083_637_8,
    ],
    [
        92.388_483_452_867_79,
        2.777_147_670_083_144_3,
        573.304_430_981_608_9,
    ],
    [
        69.611_190_301_379_92,
        0.881_442_274_895_636_3,
        176.561_790_403_821_36,
    ],
    [
        21.423_754_555_191_877,
        0.216_460_113_604_786_9,
        871.064_908_412_064_8,
    ],
];

#[inline]
fn scalar_stage<S: HashFloat, const STAGE: usize>(seed: S) -> S {
    let a = f64_to(SCALAR_COEFFICIENTS[STAGE][0]);
    let b = f64_to(SCALAR_COEFFICIENTS[STAGE][1]);
    let c = f64_to(SCALAR_COEFFICIENTS[STAGE][2]);
    let x = S::sin(seed * a + b) * c;
    x - S::floor(x)
}

// Original hash mixing constants. Stage number is the output component.
const PAIR_COEFFICIENTS: [[f64; 3]; 4] = [
    [
        9.784_225_605_373_198,
        68.948_070_147_109_01,
        81.499_072_897_379_97,
    ],
    [
        80.846_736_527_084_62,
        47.747_983_481_580_206,
        967.624_285_198_662_2,
    ],
    [
        75.650_611_853_748_19,
        0.752_916_243_450_729_7,
        825.239_418_077_631_3,
    ],
    [
        97.270_786_957_904_9,
        85.310_458_882_159_8,
        329.188_368_007_135_47,
    ],
];

#[inline]
fn pair_inputs<S: HashFloat>(seed: [S; 2]) -> [S; 2] {
    [HashGen::hash1(seed[0]), HashGen::hash1(seed[1])]
}

#[inline]
fn pair_stage<S: HashFloat, const STAGE: usize>([a, b]: [S; 2]) -> S {
    let c = f64_to(PAIR_COEFFICIENTS[STAGE][0]);
    let d = f64_to(PAIR_COEFFICIENTS[STAGE][1]);
    let e = f64_to(PAIR_COEFFICIENTS[STAGE][2]);
    let x = S::sin(a * c + b * d) * e;
    x - S::floor(x)
}

// Original hash mixing constants. Stage number is the output component.
const TRIPLE_COEFFICIENTS: [[f64; 4]; 4] = [
    [
        98.540_901_844_300_27,
        67.781_416_841_284_85,
        2.094_806_887_694_622_4,
        280.465_965_479_366_2,
    ],
    [
        31.298_218_789_251_166,
        20.629_869_323_489_913,
        49.454_054_533_406_79,
        65.418_047_936_625_31,
    ],
    [
        24.301_021_506_408_738,
        2.287_647_662_619_474,
        97.795_539_177_359,
        738.377_375_351_762_2,
    ],
    [
        13.335_891_984_218_563,
        77.765_494_753_703_58,
        78.433_323_915_272_21,
        582.218_155_377_012_3,
    ],
];

#[inline]
fn triple_inputs<S: HashFloat>(seed: [S; 3]) -> [S; 3] {
    [
        HashGen::hash1([seed[1], seed[2]]),
        HashGen::hash1([seed[2], seed[0]]),
        HashGen::hash1([seed[0], seed[1]]),
    ]
}

#[inline]
fn triple_stage<S: HashFloat, const STAGE: usize>([a, b, c]: [S; 3]) -> S {
    let d = f64_to(TRIPLE_COEFFICIENTS[STAGE][0]);
    let e = f64_to(TRIPLE_COEFFICIENTS[STAGE][1]);
    let f = f64_to(TRIPLE_COEFFICIENTS[STAGE][2]);
    let g = f64_to(TRIPLE_COEFFICIENTS[STAGE][3]);
    let x = S::sin(a * d + b * e + c * f) * g;
    x - S::floor(x)
}

// Original hash mixing constants. Stage number is the output component.
const QUAD_COEFFICIENTS: [[f64; 5]; 4] = [
    [
        45.487_415_065_002_66,
        61.983_961_049_580_714,
        71.687_107_480_472_83,
        17.018_593_362_218_127,
        233.187_597_250_993_15,
    ],
    [
        94.899_566_224_222_71,
        67.233_540_630_417_12,
        1.703_998_383_970_983_8,
        32.280_513_717_302_526,
        939.884_074_512_732_8,
    ],
    [
        8.536_774_960_497_063,
        88.789_434_660_774_06,
        57.999_707_672_724_874,
        93.952_891_460_388_17,
        310.706_244_154_918_86,
    ],
    [
        46.900_144_129_289_025,
        38.332_917_821_725_6,
        95.109_841_471_202_59,
        38.513_377_313_130_87,
        843.506_386_011_757_2,
    ],
];

#[inline]
fn quad_inputs<S: HashFloat>(seed: [S; 4]) -> [S; 4] {
    [
        HashGen::hash1([seed[1], seed[2], seed[3]]),
        HashGen::hash1([seed[2], seed[3], seed[0]]),
        HashGen::hash1([seed[3], seed[0], seed[1]]),
        HashGen::hash1([seed[0], seed[1], seed[2]]),
    ]
}

#[inline]
fn quad_stage<S: HashFloat, const STAGE: usize>([a, b, c, d]: [S; 4]) -> S {
    let e = f64_to(QUAD_COEFFICIENTS[STAGE][0]);
    let f = f64_to(QUAD_COEFFICIENTS[STAGE][1]);
    let g = f64_to(QUAD_COEFFICIENTS[STAGE][2]);
    let h = f64_to(QUAD_COEFFICIENTS[STAGE][3]);
    let i = f64_to(QUAD_COEFFICIENTS[STAGE][4]);
    let x = S::sin(a * e + b * f + c * g + d * h) * i;
    x - S::floor(x)
}

impl<S: HashFloat> HashGen<S> for S {
    #[inline]
    fn hash1(seed: S) -> S { scalar_stage::<S, 0>(seed) }

    #[inline]
    fn hash2(seed: S) -> [S; 2] { [scalar_stage::<S, 0>(seed), scalar_stage::<S, 1>(seed)] }

    #[inline]
    fn hash3(seed: S) -> [S; 3] {
        [
            scalar_stage::<S, 0>(seed),
            scalar_stage::<S, 1>(seed),
            scalar_stage::<S, 2>(seed),
        ]
    }

    #[inline]
    fn hash4(seed: S) -> [S; 4] {
        [
            scalar_stage::<S, 0>(seed),
            scalar_stage::<S, 1>(seed),
            scalar_stage::<S, 2>(seed),
            scalar_stage::<S, 3>(seed),
        ]
    }
}

/// Every multidimensional seed recursively hashes its input components,
/// then evaluates one through four dimension-specific mixing stages.
/// Macro expansion preserves compile-time stage indices and one input
/// prehash per call, rather than constructing a runtime stage loop.
macro_rules! impl_vector_hashgen {
    ($dimension:literal, $inputs:ident, $stage:ident) => {
        impl<S: HashFloat> HashGen<S> for [S; $dimension] {
            #[inline]
            fn hash1(seed: Self) -> S {
                let components = $inputs(seed);
                $stage::<S, 0>(components)
            }

            #[inline]
            fn hash2(seed: Self) -> [S; 2] {
                let components = $inputs(seed);
                [$stage::<S, 0>(components), $stage::<S, 1>(components)]
            }

            #[inline]
            fn hash3(seed: Self) -> [S; 3] {
                let components = $inputs(seed);
                [
                    $stage::<S, 0>(components),
                    $stage::<S, 1>(components),
                    $stage::<S, 2>(components),
                ]
            }

            #[inline]
            fn hash4(seed: Self) -> [S; 4] {
                let components = $inputs(seed);
                [
                    $stage::<S, 0>(components),
                    $stage::<S, 1>(components),
                    $stage::<S, 2>(components),
                    $stage::<S, 3>(components),
                ]
            }
        }
    };
}
impl_vector_hashgen!(2, pair_inputs, pair_stage);
impl_vector_hashgen!(3, triple_inputs, triple_stage);
impl_vector_hashgen!(4, quad_inputs, quad_stage);

impl<S: HashFloat> HashGen<S> for [S; 1] {
    fn hash1(seed: Self) -> S { S::hash1(seed[0]) }
    fn hash2(seed: Self) -> [S; 2] { S::hash2(seed[0]) }
    fn hash3(seed: Self) -> [S; 3] { S::hash3(seed[0]) }
    fn hash4(seed: Self) -> [S; 4] { S::hash4(seed[0]) }
}

macro_rules! derive_hashgen {
    ($from: ty, $into: ty) => {
        impl<S: HashFloat> HashGen<S> for $from {
            fn hash1(seed: Self) -> S { <$into>::hash1(seed.into()) }
            fn hash2(seed: Self) -> [S; 2] { <$into>::hash2(seed.into()) }
            fn hash3(seed: Self) -> [S; 3] { <$into>::hash3(seed.into()) }
            fn hash4(seed: Self) -> [S; 4] { <$into>::hash4(seed.into()) }
        }
    };
}

use cgmath::{Point1, Point2, Point3, Vector1, Vector2, Vector3, Vector4};
derive_hashgen!(Point1<S>, [S; 1]);
derive_hashgen!(Point2<S>, [S; 2]);
derive_hashgen!(Point3<S>, [S; 3]);
derive_hashgen!(Vector1<S>, [S; 1]);
derive_hashgen!(Vector2<S>, [S; 2]);
derive_hashgen!(Vector3<S>, [S; 3]);
derive_hashgen!(Vector4<S>, [S; 4]);

/// Take one random unit vector
pub fn take_one_unit<G: HashGen<f64>>(seed: G) -> Vector3<f64> {
    let u = HashGen::hash2(seed);
    let theta = 2.0 * std::f64::consts::PI * u[0];
    let z = 2.0f64.mul_add(u[1], -1.0);
    let r = f64::sqrt(f64::mul_add(z, -z, 1.0));
    Vector3::new(r * f64::cos(theta), r * f64::sin(theta), z)
}

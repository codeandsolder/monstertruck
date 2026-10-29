#![allow(
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
    #[allow(
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

impl<S: HashFloat> HashGen<S> for S {
    fn hash1(s: S) -> S {
        let a = f64_to(61.909_685_033_545_934);
        let b = f64_to(8.436_303_256_302_796);
        let c = f64_to(220.678_620_083_637_8);
        let x = S::sin(s * a + b) * c;
        x - S::floor(x)
    }

    fn hash2(s: S) -> [S; 2] {
        let a = f64_to(61.909_685_033_545_934);
        let b = f64_to(8.436_303_256_302_796);
        let c = f64_to(220.678_620_083_637_8);
        let x = S::sin(s * a + b) * c;
        let a = f64_to(92.388_483_452_867_79);
        let b = f64_to(2.777_147_670_083_144_3);
        let c = f64_to(573.304_430_981_608_9);
        let y = S::sin(s * a + b) * c;
        [x - S::floor(x), y - S::floor(y)]
    }

    fn hash3(s: S) -> [S; 3] {
        let a = f64_to(61.909_685_033_545_934);
        let b = f64_to(8.436_303_256_302_796);
        let c = f64_to(220.678_620_083_637_8);
        let x = S::sin(s * a + b) * c;
        let a = f64_to(92.388_483_452_867_79);
        let b = f64_to(2.777_147_670_083_144_3);
        let c = f64_to(573.304_430_981_608_9);
        let y = S::sin(s * a + b) * c;
        let a = f64_to(69.611_190_301_379_92);
        let b = f64_to(0.881_442_274_895_636_3);
        let c = f64_to(176.561_790_403_821_36);
        let z = S::sin(s * a + b) * c;
        [x - S::floor(x), y - S::floor(y), z - S::floor(z)]
    }

    fn hash4(s: S) -> [S; 4] {
        let a = f64_to(61.909_685_033_545_934);
        let b = f64_to(8.436_303_256_302_796);
        let c = f64_to(220.678_620_083_637_8);
        let x = S::sin(s * a + b) * c;
        let a = f64_to(92.388_483_452_867_79);
        let b = f64_to(2.777_147_670_083_144_3);
        let c = f64_to(573.304_430_981_608_9);
        let y = S::sin(s * a + b) * c;
        let a = f64_to(69.611_190_301_379_92);
        let b = f64_to(0.881_442_274_895_636_3);
        let c = f64_to(176.561_790_403_821_36);
        let z = S::sin(s * a + b) * c;
        let a = f64_to(21.423_754_555_191_877);
        let b = f64_to(0.216_460_113_604_786_9);
        let c = f64_to(871.064_908_412_064_8);
        let w = S::sin(s * a + b) * c;
        [
            x - S::floor(x),
            y - S::floor(y),
            z - S::floor(z),
            w - S::floor(w),
        ]
    }
}

impl<S: HashFloat> HashGen<S> for [S; 2] {
    fn hash1(v: [S; 2]) -> S {
        let a = HashGen::hash1(v[0]);
        let b = HashGen::hash1(v[1]);
        let c = f64_to(9.784_225_605_373_198);
        let d = f64_to(68.948_070_147_109_01);
        let e = f64_to(81.499_072_897_379_97);
        let x = S::sin(a * c + b * d) * e;
        x - S::floor(x)
    }

    fn hash2(v: [S; 2]) -> [S; 2] {
        let a = HashGen::hash1(v[0]);
        let b = HashGen::hash1(v[1]);
        let c = f64_to(9.784_225_605_373_198);
        let d = f64_to(68.948_070_147_109_01);
        let e = f64_to(81.499_072_897_379_97);
        let x = S::sin(a * c + b * d) * e;
        let c = f64_to(80.846_736_527_084_62);
        let d = f64_to(47.747_983_481_580_206);
        let e = f64_to(967.624_285_198_662_2);
        let y = S::sin(a * c + b * d) * e;
        [x - S::floor(x), y - S::floor(y)]
    }

    fn hash3(v: [S; 2]) -> [S; 3] {
        let a = HashGen::hash1(v[0]);
        let b = HashGen::hash1(v[1]);
        let c = f64_to(9.784_225_605_373_198);
        let d = f64_to(68.948_070_147_109_01);
        let e = f64_to(81.499_072_897_379_97);
        let x = S::sin(a * c + b * d) * e;
        let c = f64_to(80.846_736_527_084_62);
        let d = f64_to(47.747_983_481_580_206);
        let e = f64_to(967.624_285_198_662_2);
        let y = S::sin(a * c + b * d) * e;
        let c = f64_to(75.650_611_853_748_19);
        let d = f64_to(0.752_916_243_450_729_7);
        let e = f64_to(825.239_418_077_631_3);
        let z = S::sin(a * c + b * d) * e;
        [x - S::floor(x), y - S::floor(y), z - S::floor(z)]
    }

    fn hash4(v: [S; 2]) -> [S; 4] {
        let a = HashGen::hash1(v[0]);
        let b = HashGen::hash1(v[1]);
        let c = f64_to(9.784_225_605_373_198);
        let d = f64_to(68.948_070_147_109_01);
        let e = f64_to(81.499_072_897_379_97);
        let x = S::sin(a * c + b * d) * e;
        let c = f64_to(80.846_736_527_084_62);
        let d = f64_to(47.747_983_481_580_206);
        let e = f64_to(967.624_285_198_662_2);
        let y = S::sin(a * c + b * d) * e;
        let c = f64_to(75.650_611_853_748_19);
        let d = f64_to(0.752_916_243_450_729_7);
        let e = f64_to(825.239_418_077_631_3);
        let z = S::sin(a * c + b * d) * e;
        let c = f64_to(97.270_786_957_904_9);
        let d = f64_to(85.310_458_882_159_8);
        let e = f64_to(329.188_368_007_135_47);
        let w = S::sin(a * c + b * d) * e;
        [
            x - S::floor(x),
            y - S::floor(y),
            z - S::floor(z),
            w - S::floor(w),
        ]
    }
}

impl<S: HashFloat> HashGen<S> for [S; 3] {
    fn hash1(v: [S; 3]) -> S {
        let a = HashGen::hash1([v[1], v[2]]);
        let b = HashGen::hash1([v[2], v[0]]);
        let c = HashGen::hash1([v[0], v[1]]);
        let d = f64_to(98.540_901_844_300_27);
        let e = f64_to(67.781_416_841_284_85);
        let f = f64_to(2.094_806_887_694_622_4);
        let g = f64_to(280.465_965_479_366_2);
        let x = S::sin(a * d + b * e + c * f) * g;
        x - S::floor(x)
    }

    fn hash2(v: [S; 3]) -> [S; 2] {
        let a = HashGen::hash1([v[1], v[2]]);
        let b = HashGen::hash1([v[2], v[0]]);
        let c = HashGen::hash1([v[0], v[1]]);
        let d = f64_to(98.540_901_844_300_27);
        let e = f64_to(67.781_416_841_284_85);
        let f = f64_to(2.094_806_887_694_622_4);
        let g = f64_to(280.465_965_479_366_2);
        let x = S::sin(a * d + b * e + c * f) * g;
        let d = f64_to(31.298_218_789_251_166);
        let e = f64_to(20.629_869_323_489_913);
        let f = f64_to(49.454_054_533_406_79);
        let g = f64_to(65.418_047_936_625_31);
        let y = S::sin(a * d + b * e + c * f) * g;
        [x - S::floor(x), y - S::floor(y)]
    }

    fn hash3(v: [S; 3]) -> [S; 3] {
        let a = HashGen::hash1([v[1], v[2]]);
        let b = HashGen::hash1([v[2], v[0]]);
        let c = HashGen::hash1([v[0], v[1]]);
        let d = f64_to(98.540_901_844_300_27);
        let e = f64_to(67.781_416_841_284_85);
        let f = f64_to(2.094_806_887_694_622_4);
        let g = f64_to(280.465_965_479_366_2);
        let x = S::sin(a * d + b * e + c * f) * g;
        let d = f64_to(31.298_218_789_251_166);
        let e = f64_to(20.629_869_323_489_913);
        let f = f64_to(49.454_054_533_406_79);
        let g = f64_to(65.418_047_936_625_31);
        let y = S::sin(a * d + b * e + c * f) * g;
        let d = f64_to(24.301_021_506_408_738);
        let e = f64_to(2.287_647_662_619_474);
        let f = f64_to(97.795_539_177_359);
        let g = f64_to(738.377_375_351_762_2);
        let z = S::sin(a * d + b * e + c * f) * g;
        [x - S::floor(x), y - S::floor(y), z - S::floor(z)]
    }

    fn hash4(v: [S; 3]) -> [S; 4] {
        let a = HashGen::hash1([v[1], v[2]]);
        let b = HashGen::hash1([v[2], v[0]]);
        let c = HashGen::hash1([v[0], v[1]]);
        let d = f64_to(98.540_901_844_300_27);
        let e = f64_to(67.781_416_841_284_85);
        let f = f64_to(2.094_806_887_694_622_4);
        let g = f64_to(280.465_965_479_366_2);
        let x = S::sin(a * d + b * e + c * f) * g;
        let d = f64_to(31.298_218_789_251_166);
        let e = f64_to(20.629_869_323_489_913);
        let f = f64_to(49.454_054_533_406_79);
        let g = f64_to(65.418_047_936_625_31);
        let y = S::sin(a * d + b * e + c * f) * g;
        let d = f64_to(24.301_021_506_408_738);
        let e = f64_to(2.287_647_662_619_474);
        let f = f64_to(97.795_539_177_359);
        let g = f64_to(738.377_375_351_762_2);
        let z = S::sin(a * d + b * e + c * f) * g;
        let d = f64_to(13.335_891_984_218_563);
        let e = f64_to(77.765_494_753_703_58);
        let f = f64_to(78.433_323_915_272_21);
        let g = f64_to(582.218_155_377_012_3);
        let w = S::sin(a * d + b * e + c * f) * g;
        [
            x - S::floor(x),
            y - S::floor(y),
            z - S::floor(z),
            w - S::floor(w),
        ]
    }
}

impl<S: HashFloat> HashGen<S> for [S; 4] {
    fn hash1(v: [S; 4]) -> S {
        let a = HashGen::hash1([v[1], v[2], v[3]]);
        let b = HashGen::hash1([v[2], v[3], v[0]]);
        let c = HashGen::hash1([v[3], v[0], v[1]]);
        let d = HashGen::hash1([v[0], v[1], v[2]]);
        let e = f64_to(45.487_415_065_002_66);
        let f = f64_to(61.983_961_049_580_714);
        let g = f64_to(71.687_107_480_472_83);
        let h = f64_to(17.018_593_362_218_127);
        let i = f64_to(233.187_597_250_993_15);
        let x = S::sin(a * e + b * f + c * g + d * h) * i;
        x - S::floor(x)
    }

    fn hash2(v: [S; 4]) -> [S; 2] {
        let a = HashGen::hash1([v[1], v[2], v[3]]);
        let b = HashGen::hash1([v[2], v[3], v[0]]);
        let c = HashGen::hash1([v[3], v[0], v[1]]);
        let d = HashGen::hash1([v[0], v[1], v[2]]);
        let e = f64_to(45.487_415_065_002_66);
        let f = f64_to(61.983_961_049_580_714);
        let g = f64_to(71.687_107_480_472_83);
        let h = f64_to(17.018_593_362_218_127);
        let i = f64_to(233.187_597_250_993_15);
        let x = S::sin(a * e + b * f + c * g + d * h) * i;
        let e = f64_to(94.899_566_224_222_71);
        let f = f64_to(67.233_540_630_417_12);
        let g = f64_to(1.703_998_383_970_983_8);
        let h = f64_to(32.280_513_717_302_526);
        let i = f64_to(939.884_074_512_732_8);
        let y = S::sin(a * e + b * f + c * g + d * h) * i;
        [x - S::floor(x), y - S::floor(y)]
    }

    fn hash3(v: [S; 4]) -> [S; 3] {
        let a = HashGen::hash1([v[1], v[2], v[3]]);
        let b = HashGen::hash1([v[2], v[3], v[0]]);
        let c = HashGen::hash1([v[3], v[0], v[1]]);
        let d = HashGen::hash1([v[0], v[1], v[2]]);
        let e = f64_to(45.487_415_065_002_66);
        let f = f64_to(61.983_961_049_580_714);
        let g = f64_to(71.687_107_480_472_83);
        let h = f64_to(17.018_593_362_218_127);
        let i = f64_to(233.187_597_250_993_15);
        let x = S::sin(a * e + b * f + c * g + d * h) * i;
        let e = f64_to(94.899_566_224_222_71);
        let f = f64_to(67.233_540_630_417_12);
        let g = f64_to(1.703_998_383_970_983_8);
        let h = f64_to(32.280_513_717_302_526);
        let i = f64_to(939.884_074_512_732_8);
        let y = S::sin(a * e + b * f + c * g + d * h) * i;
        let e = f64_to(8.536_774_960_497_063);
        let f = f64_to(88.789_434_660_774_06);
        let g = f64_to(57.999_707_672_724_874);
        let h = f64_to(93.952_891_460_388_17);
        let i = f64_to(310.706_244_154_918_86);
        let z = S::sin(a * e + b * f + c * g + d * h) * i;
        [x - S::floor(x), y - S::floor(y), z - S::floor(z)]
    }

    fn hash4(v: [S; 4]) -> [S; 4] {
        let a = HashGen::hash1([v[1], v[2], v[3]]);
        let b = HashGen::hash1([v[2], v[3], v[0]]);
        let c = HashGen::hash1([v[3], v[0], v[1]]);
        let d = HashGen::hash1([v[0], v[1], v[2]]);
        let e = f64_to(45.487_415_065_002_66);
        let f = f64_to(61.983_961_049_580_714);
        let g = f64_to(71.687_107_480_472_83);
        let h = f64_to(17.018_593_362_218_127);
        let i = f64_to(233.187_597_250_993_15);
        let x = S::sin(a * e + b * f + c * g + d * h) * i;
        let e = f64_to(94.899_566_224_222_71);
        let f = f64_to(67.233_540_630_417_12);
        let g = f64_to(1.703_998_383_970_983_8);
        let h = f64_to(32.280_513_717_302_526);
        let i = f64_to(939.884_074_512_732_8);
        let y = S::sin(a * e + b * f + c * g + d * h) * i;
        let e = f64_to(8.536_774_960_497_063);
        let f = f64_to(88.789_434_660_774_06);
        let g = f64_to(57.999_707_672_724_874);
        let h = f64_to(93.952_891_460_388_17);
        let i = f64_to(310.706_244_154_918_86);
        let z = S::sin(a * e + b * f + c * g + d * h) * i;
        let e = f64_to(46.900_144_129_289_025);
        let f = f64_to(38.332_917_821_725_6);
        let g = f64_to(95.109_841_471_202_59);
        let h = f64_to(38.513_377_313_130_87);
        let i = f64_to(843.506_386_011_757_2);
        let w = S::sin(a * e + b * f + c * g + d * h) * i;
        [
            x - S::floor(x),
            y - S::floor(y),
            z - S::floor(z),
            w - S::floor(w),
        ]
    }
}

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

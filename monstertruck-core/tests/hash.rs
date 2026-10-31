use cgmath::*;
use monstertruck_core::hash::*;

fn scalar_from_usize<S: BaseFloat>(mut value: usize) -> S {
    let mut result = S::zero();
    let mut place = S::one();
    while value != 0 {
        if value & 1 != 0 {
            result += place;
        }
        value >>= 1;
        place += place;
    }
    result
}

fn scalar_ratio<S: BaseFloat>(numerator: usize, denominator: usize) -> S {
    scalar_from_usize::<S>(numerator) / scalar_from_usize::<S>(denominator)
}

#[test]
fn hash11_test() {
    fn exec_test<S>()
    where
        S: BaseFloat + HashGen<S>,
        [S; 1]: HashGen<S>,
        [S; 2]: HashGen<S>,
        [S; 3]: HashGen<S>,
        [S; 4]: HashGen<S>,
        Vector2<S>: HashGen<S>,
        Vector3<S>: HashGen<S>,
        Vector4<S>: HashGen<S>,
        Point2<S>: HashGen<S>,
        Point3<S>: HashGen<S>, {
        const N: usize = 10000;
        let mean = (0..N).fold(S::zero(), |sum, i| {
            sum + HashGen::hash1(scalar_from_usize::<S>(i))
        }) / scalar_from_usize::<S>(N);
        let var = (0..N).fold(S::zero(), |sum, i| {
            let x = HashGen::hash1(scalar_from_usize::<S>(i)) - scalar_ratio::<S>(1, 2);
            sum + x * x
        }) / scalar_from_usize::<S>(N);
        assert!(
            S::abs(mean - scalar_ratio::<S>(1, 2)) < scalar_ratio::<S>(1, 40),
            "mean: {mean:?}"
        );
        assert!(
            S::abs(var - scalar_ratio::<S>(1, 12)) < scalar_ratio::<S>(1, 240),
            "var: {var:?}"
        );
    }

    exec_test::<f32>();
    exec_test::<f64>();
}

macro_rules! define_hashx1_test {
    ($func: ident, $vec: ty, $unit: expr, $hash: ident) => {
        #[test]
        fn $func() {
            fn exec_test<S>()
            where
                S: BaseFloat + HashGen<S>,
                [S; 1]: HashGen<S>,
                [S; 2]: HashGen<S>,
                [S; 3]: HashGen<S>,
                [S; 4]: HashGen<S>,
                Vector2<S>: HashGen<S>,
                Vector3<S>: HashGen<S>,
                Vector4<S>: HashGen<S>,
                Point2<S>: HashGen<S>,
                Point3<S>: HashGen<S>, {
                const N: usize = 10000;
                let unit: $vec = $unit;
                let mean = (0..N).fold(<$vec>::zero(), |sum, i| {
                    sum + <$vec>::from(HashGen::$hash(scalar_from_usize::<S>(i)))
                }) / scalar_from_usize::<S>(N);
                let var = (0..N).fold(<$vec>::zero(), |sum, i| {
                    let x = <$vec>::from(HashGen::$hash(scalar_from_usize::<S>(i)))
                        - unit * scalar_ratio::<S>(1, 2);
                    sum + x.mul_element_wise(x)
                }) / scalar_from_usize::<S>(N);
                assert!(
                    (mean - unit * scalar_ratio::<S>(1, 2)).magnitude() < scalar_ratio::<S>(1, 40),
                    "mean: {:?}",
                    mean,
                );
                assert!(
                    (var - unit / scalar_from_usize::<S>(12)).magnitude()
                        < scalar_ratio::<S>(1, 240),
                    "var: {:?}",
                    var,
                );
            }
            exec_test::<f32>();
            exec_test::<f64>();
        }
    };
}

define_hashx1_test!(
    hash21_test,
    Vector2<S>,
    Vector2::new(S::one(), S::one()),
    hash2
);
define_hashx1_test!(
    hash31_test,
    Vector3<S>,
    Vector3::new(S::one(), S::one(), S::one()),
    hash3
);
define_hashx1_test!(
    hash41_test,
    Vector4<S>,
    Vector4::new(S::one(), S::one(), S::one(), S::one()),
    hash4
);

#[test]
fn hash12_test() {
    fn exec_test<S>()
    where
        S: BaseFloat + HashGen<S>,
        [S; 1]: HashGen<S>,
        [S; 2]: HashGen<S>,
        [S; 3]: HashGen<S>,
        [S; 4]: HashGen<S>,
        Vector2<S>: HashGen<S>,
        Vector3<S>: HashGen<S>,
        Vector4<S>: HashGen<S>,
        Point2<S>: HashGen<S>,
        Point3<S>: HashGen<S>, {
        const N: usize = 100;
        let mean = (0..(N * N)).fold(S::zero(), |sum, i| {
            let v = [scalar_from_usize::<S>(i / N), scalar_from_usize::<S>(i % N)];
            sum + HashGen::hash1(v)
        }) / scalar_from_usize::<S>(N * N);
        let var = (0..(N * N)).fold(S::zero(), |sum, i| {
            let v = [scalar_from_usize::<S>(i / N), scalar_from_usize::<S>(i % N)];
            let x = HashGen::hash1(v) - scalar_ratio::<S>(1, 2);
            sum + x * x
        }) / scalar_from_usize::<S>(N * N);
        assert!(
            S::abs(mean - scalar_ratio::<S>(1, 2)) < scalar_ratio::<S>(1, 40),
            "mean: {mean:?}"
        );
        assert!(
            S::abs(var - scalar_ratio::<S>(1, 12)) < scalar_ratio::<S>(1, 240),
            "var: {var:?}"
        );
    }

    exec_test::<f32>();
    exec_test::<f64>();
}

macro_rules! define_hashx2_test {
    ($func: ident, $vec: ty, $unit: expr, $hash: ident) => {
        #[test]
        fn $func() {
            fn exec_test<S>()
            where
                S: BaseFloat + HashGen<S>,
                [S; 1]: HashGen<S>,
                [S; 2]: HashGen<S>,
                [S; 3]: HashGen<S>,
                [S; 4]: HashGen<S>,
                Vector2<S>: HashGen<S>,
                Vector3<S>: HashGen<S>,
                Vector4<S>: HashGen<S>,
                Point2<S>: HashGen<S>,
                Point3<S>: HashGen<S>, {
                const N: usize = 100;
                let unit: $vec = $unit;
                let mean = (0..(N * N)).fold(<$vec>::zero(), |sum, i| {
                    let v = [scalar_from_usize::<S>(i / N), scalar_from_usize::<S>(i % N)];
                    sum + <$vec>::from(HashGen::$hash(v))
                }) / scalar_from_usize::<S>(N * N);
                let var = (0..(N * N)).fold(<$vec>::zero(), |sum, i| {
                    let v = [scalar_from_usize::<S>(i / N), scalar_from_usize::<S>(i % N)];
                    let x = <$vec>::from(HashGen::$hash(v)) - unit * scalar_ratio::<S>(1, 2);
                    sum + x.mul_element_wise(x)
                }) / scalar_from_usize::<S>(N * N);
                assert!(
                    (mean - unit * scalar_ratio::<S>(1, 2)).magnitude() < scalar_ratio::<S>(1, 40),
                    "mean: {:?}",
                    mean,
                );
                assert!(
                    (var - unit / scalar_from_usize::<S>(12)).magnitude()
                        < scalar_ratio::<S>(1, 240),
                    "var: {:?}",
                    var,
                );
            }
            exec_test::<f32>();
            exec_test::<f64>();
        }
    };
}

define_hashx2_test!(
    hash22_test,
    Vector2<S>,
    Vector2::new(S::one(), S::one()),
    hash2
);
define_hashx2_test!(
    hash32_test,
    Vector3<S>,
    Vector3::new(S::one(), S::one(), S::one()),
    hash3
);
define_hashx2_test!(
    hash42_test,
    Vector4<S>,
    Vector4::new(S::one(), S::one(), S::one(), S::one()),
    hash4
);

#[test]
fn hash13_test() {
    fn exec_test<S>()
    where
        S: BaseFloat + HashGen<S>,
        [S; 1]: HashGen<S>,
        [S; 2]: HashGen<S>,
        [S; 3]: HashGen<S>,
        [S; 4]: HashGen<S>,
        Vector2<S>: HashGen<S>,
        Vector3<S>: HashGen<S>,
        Vector4<S>: HashGen<S>,
        Point2<S>: HashGen<S>,
        Point3<S>: HashGen<S>, {
        const N: usize = 21;
        let mean = (0..(N * N * N)).fold(S::zero(), |sum, i| {
            let v = [
                scalar_from_usize::<S>(i / (N * N)),
                scalar_from_usize::<S>((i / N) % N),
                scalar_from_usize::<S>(i % N),
            ];
            sum + HashGen::hash1(v)
        }) / scalar_from_usize::<S>(N * N * N);
        let var = (0..(N * N * N)).fold(S::zero(), |sum, i| {
            let v = [
                scalar_from_usize::<S>(i / (N * N)),
                scalar_from_usize::<S>((i / N) % N),
                scalar_from_usize::<S>(i % N),
            ];
            let x = HashGen::hash1(v) - scalar_ratio::<S>(1, 2);
            sum + x * x
        }) / scalar_from_usize::<S>(N * N * N);
        assert!(
            S::abs(mean - scalar_ratio::<S>(1, 2)) < scalar_ratio::<S>(1, 40),
            "mean: {mean:?}"
        );
        assert!(
            S::abs(var - scalar_ratio::<S>(1, 12)) < scalar_ratio::<S>(1, 240),
            "var: {var:?}"
        );
    }

    exec_test::<f32>();
    exec_test::<f64>();
}

macro_rules! define_hashx3_test {
    ($func: ident, $vec: ty, $unit: expr, $hash: ident) => {
        #[test]
        fn $func() {
            fn exec_test<S>()
            where
                S: BaseFloat + HashGen<S>,
                [S; 1]: HashGen<S>,
                [S; 2]: HashGen<S>,
                [S; 3]: HashGen<S>,
                [S; 4]: HashGen<S>,
                Vector2<S>: HashGen<S>,
                Vector3<S>: HashGen<S>,
                Vector4<S>: HashGen<S>,
                Point2<S>: HashGen<S>,
                Point3<S>: HashGen<S>, {
                const N: usize = 21;
                let unit: $vec = $unit;
                let mean = (0..(N * N * N)).fold(<$vec>::zero(), |sum, i| {
                    let v = [
                        scalar_from_usize::<S>(i / (N * N)),
                        scalar_from_usize::<S>((i / N) % N),
                        scalar_from_usize::<S>(i % N),
                    ];
                    sum + <$vec>::from(HashGen::$hash(v))
                }) / scalar_from_usize::<S>(N * N * N);
                let var = (0..(N * N * N)).fold(<$vec>::zero(), |sum, i| {
                    let v = [
                        scalar_from_usize::<S>(i / (N * N)),
                        scalar_from_usize::<S>((i / N) % N),
                        scalar_from_usize::<S>(i % N),
                    ];
                    let x = <$vec>::from(HashGen::$hash(v)) - unit * scalar_ratio::<S>(1, 2);
                    sum + x.mul_element_wise(x)
                }) / scalar_from_usize::<S>(N * N * N);
                assert!(
                    (mean - unit * scalar_ratio::<S>(1, 2)).magnitude() < scalar_ratio::<S>(1, 40),
                    "mean: {:?}",
                    mean,
                );
                assert!(
                    (var - unit * scalar_ratio::<S>(1, 12)).magnitude() < scalar_ratio::<S>(1, 240),
                    "var: {:?}",
                    var,
                );
            }

            exec_test::<f32>();
            exec_test::<f64>();
        }
    };
}

define_hashx3_test!(
    hash23_test,
    Vector2<S>,
    Vector2::new(S::one(), S::one()),
    hash2
);
define_hashx3_test!(
    hash33_test,
    Vector3<S>,
    Vector3::new(S::one(), S::one(), S::one()),
    hash3
);
define_hashx3_test!(
    hash43_test,
    Vector4<S>,
    Vector4::new(S::one(), S::one(), S::one(), S::one()),
    hash4
);

#[test]
fn hash14_test() {
    fn exec_test<S>()
    where
        S: BaseFloat + HashGen<S>,
        [S; 1]: HashGen<S>,
        [S; 2]: HashGen<S>,
        [S; 3]: HashGen<S>,
        [S; 4]: HashGen<S>,
        Vector2<S>: HashGen<S>,
        Vector3<S>: HashGen<S>,
        Vector4<S>: HashGen<S>,
        Point2<S>: HashGen<S>,
        Point3<S>: HashGen<S>, {
        const N: usize = 10;
        let mean = (0..(N * N * N * N)).fold(S::zero(), |sum, i| {
            let v = [
                scalar_from_usize::<S>(i / (N * N * N)),
                scalar_from_usize::<S>((i / (N * N)) % N),
                scalar_from_usize::<S>((i / N) % N),
                scalar_from_usize::<S>(i % N),
            ];
            sum + HashGen::hash1(v)
        }) / scalar_from_usize::<S>(N * N * N * N);
        let var = (0..(N * N * N * N)).fold(S::zero(), |sum, i| {
            let v = [
                scalar_from_usize::<S>(i / (N * N * N)),
                scalar_from_usize::<S>((i / (N * N)) % N),
                scalar_from_usize::<S>((i / N) % N),
                scalar_from_usize::<S>(i % N),
            ];
            let x = HashGen::hash1(v) - scalar_ratio::<S>(1, 2);
            sum + x * x
        }) / scalar_from_usize::<S>(N * N * N * N);
        assert!(
            S::abs(mean - scalar_ratio::<S>(1, 2)) < scalar_ratio::<S>(1, 40),
            "mean: {mean:?}"
        );
        assert!(
            S::abs(var - scalar_ratio::<S>(1, 12)) < scalar_ratio::<S>(1, 240),
            "var: {var:?}"
        );
    }

    exec_test::<f32>();
    exec_test::<f64>();
}

macro_rules! define_hashx4_test {
    ($func: ident, $vec: ty, $unit: expr, $hash: ident) => {
        #[test]
        fn $func() {
            fn exec_test<S>()
            where
                S: BaseFloat + HashGen<S>,
                [S; 1]: HashGen<S>,
                [S; 2]: HashGen<S>,
                [S; 3]: HashGen<S>,
                [S; 4]: HashGen<S>,
                Vector2<S>: HashGen<S>,
                Vector3<S>: HashGen<S>,
                Vector4<S>: HashGen<S>,
                Point2<S>: HashGen<S>,
                Point3<S>: HashGen<S>, {
                const N: usize = 10;
                let unit: $vec = $unit;
                let mean = (0..(N * N * N * N)).fold(<$vec>::zero(), |sum, i| {
                    let v = [
                        scalar_from_usize::<S>(i / (N * N * N)),
                        scalar_from_usize::<S>((i / (N * N)) % N),
                        scalar_from_usize::<S>((i / N) % N),
                        scalar_from_usize::<S>(i % N),
                    ];
                    sum + <$vec>::from(HashGen::$hash(v))
                }) / scalar_from_usize::<S>(N * N * N * N);
                let var = (0..(N * N * N * N)).fold(<$vec>::zero(), |sum, i| {
                    let v = [
                        scalar_from_usize::<S>(i / (N * N * N)),
                        scalar_from_usize::<S>((i / (N * N)) % N),
                        scalar_from_usize::<S>((i / N) % N),
                        scalar_from_usize::<S>(i % N),
                    ];
                    let x = <$vec>::from(HashGen::$hash(v)) - unit * scalar_ratio::<S>(1, 2);
                    sum + x.mul_element_wise(x)
                }) / scalar_from_usize::<S>(N * N * N * N);
                assert!(
                    (mean - unit * scalar_ratio::<S>(1, 2)).magnitude() < scalar_ratio::<S>(1, 40),
                    "mean: {:?}",
                    mean
                );
                assert!(
                    (var - unit * scalar_ratio::<S>(1, 12)).magnitude() < scalar_ratio::<S>(1, 240),
                    "var: {:?}",
                    var
                );
            }

            exec_test::<f32>();
            exec_test::<f64>();
        }
    };
}

define_hashx4_test!(
    hash24_test,
    Vector2<S>,
    Vector2::new(S::one(), S::one()),
    hash2
);
define_hashx4_test!(
    hash34_test,
    Vector3<S>,
    Vector3::new(S::one(), S::one(), S::one()),
    hash3
);
define_hashx4_test!(
    hash44_test,
    Vector4<S>,
    Vector4::new(S::one(), S::one(), S::one(), S::one()),
    hash4
);

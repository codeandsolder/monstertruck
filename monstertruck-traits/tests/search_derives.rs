//! Runtime coverage for the four generated search-parameter dispatchers.
//! Both enum variant shapes and generic tuple structs must forward the
//! parameter-space dimension, hint and trial count unchanged.
#![cfg(feature = "derive")]

use monstertruck_traits::{
    CurveParameter, SearchNearestParameter, SearchNearestParameterD1, SearchNearestParameterD2,
    SearchParameter, SearchParameterD1, SearchParameterD2, SearchParameterHint1D,
    SearchParameterHint2D, SurfaceParameter,
};

#[derive(Clone, Copy, Debug)]
struct Probe(f64);

impl SearchNearestParameter<CurveParameter> for Probe {
    type Point = f64;

    fn search_nearest_parameter<H: Into<SearchParameterHint1D>>(
        &self,
        pt: f64,
        hint: H,
        trials: usize,
    ) -> Option<f64> {
        let hint = match hint.into() {
            SearchParameterHint1D::Parameter(value) => value,
            SearchParameterHint1D::Range(low, high) => low + high,
            SearchParameterHint1D::None => 0.0,
        };
        Some(self.0 + pt + hint + trials as f64)
    }
}

impl SearchParameter<CurveParameter> for Probe {
    type Point = f64;

    fn search_parameter<H: Into<SearchParameterHint1D>>(
        &self,
        pt: f64,
        hint: H,
        trials: usize,
    ) -> Option<f64> {
        <Self as SearchNearestParameter<CurveParameter>>::search_nearest_parameter(
            self, pt, hint, trials,
        )
        .map(|value| value + 1.0)
    }
}

impl SearchNearestParameter<SurfaceParameter> for Probe {
    type Point = f64;

    fn search_nearest_parameter<H: Into<SearchParameterHint2D>>(
        &self,
        pt: f64,
        hint: H,
        trials: usize,
    ) -> Option<(f64, f64)> {
        let hint = match hint.into() {
            SearchParameterHint2D::Parameter(u, v) => (u, v),
            SearchParameterHint2D::Range((u, _), (v, _)) => (u, v),
            SearchParameterHint2D::None => (0.0, 0.0),
        };
        Some((self.0 + pt + hint.0, self.0 + pt + hint.1 + trials as f64))
    }
}

impl SearchParameter<SurfaceParameter> for Probe {
    type Point = f64;

    fn search_parameter<H: Into<SearchParameterHint2D>>(
        &self,
        pt: f64,
        hint: H,
        trials: usize,
    ) -> Option<(f64, f64)> {
        <Self as SearchNearestParameter<SurfaceParameter>>::search_nearest_parameter(
            self, pt, hint, trials,
        )
        .map(|(u, v)| (u + 1.0, v + 1.0))
    }
}

#[derive(
    SearchNearestParameterD1, SearchNearestParameterD2, SearchParameterD1, SearchParameterD2,
)]
struct Wrapped<T: Copy>(T)
where T: core::fmt::Debug;

#[derive(
    SearchNearestParameterD1, SearchNearestParameterD2, SearchParameterD1, SearchParameterD2,
)]
enum Choice<T: Copy>
where T: core::fmt::Debug {
    Tuple(T),
    Named { value: T },
}

fn assert_matches_probe<T>(derived: &T, source: &Probe)
where T: SearchNearestParameter<CurveParameter, Point = f64>
        + SearchNearestParameter<SurfaceParameter, Point = f64>
        + SearchParameter<CurveParameter, Point = f64>
        + SearchParameter<SurfaceParameter, Point = f64> {
    let point = 1.25;
    let trials = 7;
    assert_eq!(
        <T as SearchNearestParameter<CurveParameter>>::search_nearest_parameter(
            derived,
            point,
            Some(0.5),
            trials
        ),
        <Probe as SearchNearestParameter<CurveParameter>>::search_nearest_parameter(
            source,
            point,
            Some(0.5),
            trials
        )
    );
    assert_eq!(
        <T as SearchParameter<CurveParameter>>::search_parameter(
            derived,
            point,
            (0.25, 0.75),
            trials
        ),
        <Probe as SearchParameter<CurveParameter>>::search_parameter(
            source,
            point,
            (0.25, 0.75),
            trials
        )
    );
    assert_eq!(
        <T as SearchNearestParameter<SurfaceParameter>>::search_nearest_parameter(
            derived,
            point,
            (0.1, 0.2),
            trials
        ),
        <Probe as SearchNearestParameter<SurfaceParameter>>::search_nearest_parameter(
            source,
            point,
            (0.1, 0.2),
            trials
        )
    );
    assert_eq!(
        <T as SearchParameter<SurfaceParameter>>::search_parameter(
            derived,
            point,
            ((0.1, 0.4), (0.2, 0.5)),
            trials
        ),
        <Probe as SearchParameter<SurfaceParameter>>::search_parameter(
            source,
            point,
            ((0.1, 0.4), (0.2, 0.5)),
            trials
        )
    );
}

#[test]
fn generic_tuple_and_both_enum_variant_shapes_delegate_all_searches() {
    let probe = Probe(0.25);
    assert_matches_probe(&Wrapped(probe), &probe);
    assert_matches_probe(&Choice::Tuple(probe), &probe);
    assert_matches_probe(&Choice::Named { value: probe }, &probe);
}

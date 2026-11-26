//! Independent material occupancy witnesses from three STEP operands.
use anyhow::{Context, Result, bail, ensure};
use clap::{Args, ValueEnum};
use num_traits::ToPrimitive;
use occt_wasm::{OcctKernel, ShapeHandle};
use rand_chacha::{
    ChaCha12Rng,
    rand_core::{Rng, SeedableRng},
};
use serde::{Deserialize, Serialize};
use std::{
    fs::{self, File},
    path::{Path, PathBuf},
    process::ExitCode,
};

type Point = [f64; 3];

#[derive(Debug, Clone, Copy, Serialize, ValueEnum)]
#[serde(rename_all = "snake_case")]
enum Operation {
    Difference,
    Intersection,
    Union,
    Xor,
}

impl Operation {
    fn material(self, left: bool, right: bool) -> bool {
        match self {
            Self::Difference => left && !right,
            Self::Intersection => left && right,
            Self::Union => left || right,
            Self::Xor => left != right,
        }
    }
}

#[derive(Debug, Args)]
pub struct OccupancyArgs {
    #[arg(long, value_enum)]
    operation: Operation,
    #[arg(long)]
    left: PathBuf,
    #[arg(long)]
    right: PathBuf,
    #[arg(long)]
    result: PathBuf,
    #[arg(long = "point", num_args = 3, allow_hyphen_values = true, value_names = ["X", "Y", "Z"])]
    points: Vec<f64>,
    #[arg(long)]
    points_json: Option<PathBuf>,
    #[arg(long = "box", num_args = 6, allow_hyphen_values = true, value_names = ["XMIN", "XMAX", "YMIN", "YMAX", "ZMIN", "ZMAX"])]
    bounds: Option<Vec<f64>>,
    #[arg(long, default_value_t = 0)]
    samples: usize,
    #[arg(long, default_value_t = 20261115)]
    seed: u64,
    #[arg(long, default_value_t = 1e-7)]
    tolerance: f64,
    #[arg(long)]
    report: Option<PathBuf>,
}

#[derive(Deserialize)]
#[serde(untagged)]
enum WitnessFile {
    List(Vec<Point>),
    Object { points: Vec<Point> },
}

fn sampled_points(bounds: [f64; 6], count: usize, seed: u64) -> Result<Vec<Point>> {
    let bounds = [
        [bounds[0], bounds[1]],
        [bounds[2], bounds[3]],
        [bounds[4], bounds[5]],
    ];
    ensure!(
        bounds
            .iter()
            .all(|&[lo, hi]| lo.is_finite() && hi.is_finite() && lo < hi && (hi - lo).is_finite()),
        "sample box needs ordered finite bounds and finite spans"
    );
    let mut rng = ChaCha12Rng::seed_from_u64(seed);
    let mut points = Vec::with_capacity(count);
    for _ in 0..count {
        let mut point = [0.0; 3];
        for (coordinate, [lo, hi]) in point.iter_mut().zip(bounds) {
            // Exactly 53 uniformly random mantissa bits, independent of
            // platform-specific floating-point random distribution code.
            let fraction = (rng.next_u64() >> 11)
                .to_f64()
                .context("sample exceeds f64")?
                * 2_f64.powi(-53);
            *coordinate = lo + (hi - lo) * fraction;
        }
        points.push(point);
    }
    Ok(points)
}

impl OccupancyArgs {
    fn witnesses(&self) -> Result<Vec<Point>> {
        ensure!(
            self.tolerance.is_finite() && self.tolerance > 0.0,
            "tolerance must be finite and positive"
        );
        ensure!(
            self.bounds.is_some() == (self.samples > 0),
            "box and positive sample count must appear together"
        );
        let (explicit, trailing) = self.points.as_chunks::<3>();
        ensure!(trailing.is_empty(), "point needs three coordinates");
        let mut points = explicit.to_vec();
        if let Some(path) = &self.points_json {
            let file: WitnessFile = serde_json::from_reader(File::open(path)?)?;
            points.extend(match file {
                WitnessFile::List(points) | WitnessFile::Object { points } => points,
            });
        }
        if let Some(bounds) = &self.bounds {
            points.extend(sampled_points(
                bounds.as_slice().try_into()?,
                self.samples,
                self.seed,
            )?);
        }
        ensure!(!points.is_empty(), "no occupancy witnesses supplied");
        ensure!(
            points.iter().flatten().all(|v| v.is_finite()),
            "points must be finite"
        );
        Ok(points)
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
#[serde(rename_all = "snake_case")]
enum State {
    Outside,
    Inside,
    Boundary,
    Unknown,
}

fn solids(kernel: &mut OcctKernel, shape: ShapeHandle) -> Result<Vec<ShapeHandle>> {
    let solids = kernel.get_sub_shape_occurrences(shape, "solid")?;
    for &solid in &solids {
        ensure!(kernel.is_valid(solid)?, "invalid BRep solid");
    }
    Ok(solids)
}

fn read_solids(kernel: &mut OcctKernel, path: &Path) -> Result<Vec<ShapeHandle>> {
    let (import, _) = crate::load(kernel, path)?;
    ensure!(
        import.transferred > 0,
        "{} transferred no STEP roots",
        path.display()
    );
    solids(kernel, import.shape).with_context(|| format!("classify {}", path.display()))
}

fn classify(
    kernel: &mut OcctKernel,
    solids: &[ShapeHandle],
    [x, y, z]: Point,
    tolerance: f64,
) -> Result<State> {
    let mut boundary = false;
    let mut unknown = false;
    for &solid in solids {
        match kernel.classify_solid_point(solid, x, y, z, tolerance)? {
            0 => {}
            1 => return Ok(State::Inside),
            2 => boundary = true,
            3 => unknown = true,
            value => bail!("invalid OCCT classifier state {value}"),
        }
    }
    Ok(if unknown {
        State::Unknown
    } else if boundary {
        State::Boundary
    } else {
        State::Outside
    })
}

#[derive(Debug, Default, Serialize)]
struct Counts {
    total: usize,
    matched: usize,
    mismatched: usize,
    boundary_ambiguous: usize,
    classifier_unknown: usize,
}

#[derive(Debug, Serialize)]
struct Violation {
    index: usize,
    point: Point,
    left: bool,
    right: bool,
    expected: bool,
    actual: bool,
}

#[derive(Debug, Serialize)]
struct UnknownWitness {
    index: usize,
    point: Point,
    states: [State; 3],
}

#[derive(Debug, Default, Serialize)]
struct Evidence {
    counts: Counts,
    violations: Vec<Violation>,
    unknown_witnesses: Vec<UnknownWitness>,
}

impl Evidence {
    fn record(&mut self, operation: Operation, index: usize, point: Point, states: [State; 3]) {
        self.counts.total += 1;
        if states.contains(&State::Unknown) {
            self.counts.classifier_unknown += 1;
            self.unknown_witnesses.push(UnknownWitness {
                index,
                point,
                states,
            });
            return;
        }
        if states.contains(&State::Boundary) {
            self.counts.boundary_ambiguous += 1;
            return;
        }
        let [left, right, actual] = states.map(|state| state == State::Inside);
        let expected = operation.material(left, right);
        if actual == expected {
            self.counts.matched += 1;
        } else {
            self.counts.mismatched += 1;
            self.violations.push(Violation {
                index,
                point,
                left,
                right,
                expected,
                actual,
            });
        }
    }

    fn exit_code(&self) -> ExitCode {
        if self.counts.mismatched > 0 {
            ExitCode::FAILURE
        } else if self.counts.matched == 0 || self.counts.classifier_unknown > 0 {
            ExitCode::from(2)
        } else {
            ExitCode::SUCCESS
        }
    }
}

fn assess(
    kernel: &mut OcctKernel,
    operands: [&[ShapeHandle]; 3],
    operation: Operation,
    points: &[Point],
    tolerance: f64,
) -> Result<Evidence> {
    let mut evidence = Evidence::default();
    for (index, &point) in points.iter().enumerate() {
        let [left, right, result] = operands;
        let states = [
            classify(kernel, left, point, tolerance)?,
            classify(kernel, right, point, tolerance)?,
            classify(kernel, result, point, tolerance)?,
        ];
        evidence.record(operation, index, point, states);
    }
    Ok(evidence)
}

#[derive(Serialize)]
struct Report<'a> {
    operation: Operation,
    left: &'a Path,
    right: &'a Path,
    result: &'a Path,
    tolerance: f64,
    points: &'a [Point],
    #[serde(flatten)]
    evidence: Evidence,
}

pub fn run(args: &OccupancyArgs) -> Result<ExitCode> {
    let points = args.witnesses()?;
    let mut kernel = OcctKernel::new()?;
    let left = read_solids(&mut kernel, &args.left)?;
    let right = read_solids(&mut kernel, &args.right)?;
    let result = read_solids(&mut kernel, &args.result)?;
    let evidence = assess(
        &mut kernel,
        [&left, &right, &result],
        args.operation,
        &points,
        args.tolerance,
    )?;
    let code = evidence.exit_code();
    let report = Report {
        operation: args.operation,
        left: &args.left,
        right: &args.right,
        result: &args.result,
        tolerance: args.tolerance,
        points: &points,
        evidence,
    };
    let json = serde_json::to_string_pretty(&report)?;
    if let Some(path) = &args.report {
        if let Some(parent) = path.parent().filter(|p| !p.as_os_str().is_empty()) {
            fs::create_dir_all(parent)?;
        }
        fs::write(path, format!("{json}\n"))?;
    }
    println!("{json}");
    Ok(code)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sampled_witnesses_are_reproducible_and_bounded() -> Result<()> {
        let bounds = [-3.0, 7.0, 1.0, 2.0, -8.0, -2.0];
        let a = sampled_points(bounds, 100, 7)?;
        let b = sampled_points(bounds, 100, 7)?;
        assert_eq!(a, b);
        assert_ne!(a, sampled_points(bounds, 100, 8)?);
        for point in a {
            assert!((-3.0..=7.0).contains(&point[0]));
            assert!((1.0..=2.0).contains(&point[1]));
            assert!((-8.0..=-2.0).contains(&point[2]));
        }
        assert!(sampled_points([0.0; 6], 1, 7).is_err());
        Ok(())
    }

    #[test]
    fn ambiguity_and_unknown_are_not_material_evidence() {
        use State::{Boundary, Inside, Outside, Unknown};
        let mut e = Evidence::default();
        let p = [0.0; 3];
        let op = Operation::Difference;
        e.record(op, 0, p, [Inside, Outside, Inside]);
        assert_eq!(e.exit_code(), ExitCode::SUCCESS);
        e.record(op, 1, p, [Boundary, Inside, Outside]);
        assert_eq!(e.counts.matched, 1);
        e.record(op, 2, p, [Inside, Unknown, Outside]);
        assert_eq!(e.exit_code(), ExitCode::from(2));
        e.record(op, 3, p, [Inside, Outside, Outside]);
        assert_eq!(e.exit_code(), ExitCode::FAILURE);
        assert_eq!(e.counts.total, 4);
        assert_eq!(e.counts.boundary_ambiguous, 1);
        assert_eq!(e.counts.classifier_unknown, 1);
        assert_eq!(e.violations.len(), 1);
    }

    #[test]
    fn wasm_boolean_material_classification() -> Result<()> {
        let mut kernel = crate::test_kernel()?;
        let a = kernel.make_box(2.0, 2.0, 2.0)?;
        let b = kernel.translate(a, 1.0, 0.0, 0.0)?;
        let cut = kernel.cut(a, b)?;
        let common = kernel.common(a, b)?;
        let union = kernel.fuse(a, b)?;
        let other_cut = kernel.cut(b, a)?;
        let xor = kernel.make_compound(&[cut, other_cut])?;
        let left = solids(&mut kernel, a)?;
        let right = solids(&mut kernel, b)?;
        let pts = [
            [-1.0, 0.5, 0.5],
            [0.5, 0.5, 0.5],
            [1.5, 0.5, 0.5],
            [2.5, 0.5, 0.5],
            [1.0, 0.5, 0.5],
        ];
        for (op, shape) in [
            (Operation::Difference, cut),
            (Operation::Intersection, common),
            (Operation::Union, union),
            (Operation::Xor, xor),
        ] {
            let result = solids(&mut kernel, shape)?;
            let evidence = assess(&mut kernel, [&left, &right, &result], op, &pts, 1e-7)?;
            assert_eq!(evidence.counts.matched, 4);
            assert_eq!(evidence.counts.boundary_ambiguous, 1);
            assert_eq!(evidence.exit_code(), ExitCode::SUCCESS);
        }
        let wrong = assess(
            &mut kernel,
            [&left, &right, &left],
            Operation::Difference,
            &pts,
            1e-7,
        )?;
        assert_eq!(wrong.counts.mismatched, 1);
        assert_eq!(wrong.exit_code(), ExitCode::FAILURE);
        Ok(())
    }
}

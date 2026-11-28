//! Compare per-face mesh area and signed volume contributions before and
//! after STEP exchange. Match using complete native boundary vertex
//! multisets in world coordinates, never face order or surface type.
//! Ambiguous or unmatched signatures are reported, not guessed.

use anyhow::Result;
use clap::Parser;
use monstertruck_io::step::load::Table;
use monstertruck_mesh::{InnerProduct, Point3, PolygonMesh};
use monstertruck_meshing::{analyzers::CalcVolume, tessellation::MeshableShape};
use monstertruck_modeling::{Curve, Surface};
use monstertruck_topology::compress::CompressedShell;
use serde_json::{Value, json};
use std::{
    collections::BTreeMap,
    fs::File,
    io::{BufReader, BufWriter, Write},
    path::PathBuf,
};

#[derive(Parser)]
struct Args {
    native: PathBuf,
    step: PathBuf,
    output: PathBuf,
    #[arg(long, default_value_t = 0.00125)]
    tolerance: f64,
    /// Quantization used only for diagnostic face association.
    #[arg(long, default_value_t = 1e-8)]
    resolution: f64,
}

type MeshedShell =
    CompressedShell<Point3, monstertruck_mesh::PolylineCurve<Point3>, Option<PolygonMesh>>;
type Signature = Vec<[u64; 3]>;

fn signature(shell: &MeshedShell, index: usize, resolution: f64) -> Option<Signature> {
    let mut points = shell.faces[index]
        .boundaries
        .iter()
        .flatten()
        .flat_map(|edge| {
            let (a, b) = shell.edges[edge.index].vertices;
            [shell.vertices[a], shell.vertices[b]]
        })
        .map(|p| {
            let values = [p.x, p.y, p.z].map(|value| (value / resolution).round());
            values
                .iter()
                .all(|value| value.is_finite())
                .then(|| values.map(|value| if value == 0.0 { 0 } else { value.to_bits() }))
        })
        .collect::<Option<Vec<_>>>()?;
    points.sort_unstable();
    Some(points)
}

fn mesh_area(mesh: &PolygonMesh) -> f64 {
    mesh.face_iter()
        .map(|polygon| {
            let origin = mesh.positions()[polygon[0].pos];
            polygon[1..]
                .windows(2)
                .map(|pair| {
                    let a = mesh.positions()[pair[0].pos] - origin;
                    let b = mesh.positions()[pair[1].pos] - origin;
                    a.cross(&b).norm() * 0.5
                })
                .sum::<f64>()
        })
        .sum()
}

#[derive(Clone)]
struct FaceRecord {
    face: usize,
    shell: Option<u64>,
    area: f64,
    volume: f64,
    polygons: usize,
}

impl FaceRecord {
    fn from_mesh(mesh: &MeshedShell, face: usize, shell: Option<u64>) -> Option<Self> {
        let polygon = mesh.faces[face].surface.as_ref()?;
        let orientation = if mesh.faces[face].orientation {
            1.0
        } else {
            -1.0
        };
        Some(Self {
            face,
            shell,
            area: mesh_area(polygon),
            volume: polygon.volume() * orientation,
            polygons: polygon.faces().len(),
        })
    }
    fn json(&self) -> Value {
        json!({"face":self.face,"shell":self.shell,"area":self.area,"volume":self.volume,"polygons":self.polygons})
    }
}

fn records(
    mesh: &MeshedShell,
    shell: Option<u64>,
    resolution: f64,
) -> BTreeMap<Signature, Vec<FaceRecord>> {
    let mut result = BTreeMap::<Signature, Vec<FaceRecord>>::new();
    for index in 0..mesh.faces.len() {
        if let Some((signature, record)) =
            signature(mesh, index, resolution).zip(FaceRecord::from_mesh(mesh, index, shell))
        {
            result.entry(signature).or_default().push(record);
        }
    }
    result
}

fn main() -> Result<()> {
    let args = Args::parse();
    anyhow::ensure!(
        args.resolution.is_finite()
            && args.resolution > 0.0
            && args.tolerance.is_finite()
            && args.tolerance >= monstertruck_modeling::KERNEL_EPSILON,
        "invalid diagnostic resolution"
    );
    let native: CompressedShell<Point3, Curve, Surface> =
        serde_json::from_reader(BufReader::new(File::open(&args.native)?))?;
    let native_records = records(&native.triangulation(args.tolerance), None, args.resolution);
    let table = Table::from_step(&std::fs::read_to_string(&args.step)?)?;
    let mut imported = BTreeMap::<Signature, Vec<FaceRecord>>::new();
    for (&id, shell) in &table.shell {
        let shell = table.to_compressed_trimmed_shell(shell)?;
        for (signature, records) in records(
            &shell.triangulation(args.tolerance),
            Some(id),
            args.resolution,
        ) {
            imported.entry(signature).or_default().extend(records);
        }
    }
    let mut matched = Vec::new();
    let mut unmatched_native = Vec::new();
    let mut unmatched_step = Vec::new();
    for (signature, original) in &native_records {
        match (original.as_slice(),imported.get(signature).map(Vec::as_slice)){
            ([a],Some([b]))=>matched.push(json!({"native":a.json(),"step":b.json(),"area_difference":b.area-a.area,"volume_difference":b.volume-a.volume})),
            _=>unmatched_native.extend(original.iter().map(FaceRecord::json)),
        }
    }
    for (signature, restored) in &imported {
        if !(restored.len() == 1
            && native_records
                .get(signature)
                .is_some_and(|records| records.len() == 1))
        {
            unmatched_step.extend(restored.iter().map(FaceRecord::json));
        }
    }
    let mut out = BufWriter::new(File::create_new(args.output)?);
    serde_json::to_writer(
        &mut out,
        &json!({"matched":matched,"unmatched_native":unmatched_native,"unmatched_step":unmatched_step}),
    )?;
    out.flush()?;
    Ok(())
}

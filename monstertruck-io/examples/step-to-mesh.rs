//! Parse STEP data, extract shape, and meshing.

use anyhow::{Context, Result, anyhow, bail};
use clap::Parser;
use monstertruck_assembly::assy::NodeEntity;
use monstertruck_io::step::load::{Table, convert::ProductEntity};
use monstertruck_meshing::prelude::*;
use monstertruck_topology::compress::{CompressedEdge, CompressedFace, CompressedShell};
use std::path::Path;

type MeshedCShell = CompressedShell<Point3, PolylineCurve<Point3>, Option<PolygonMesh>>;

#[derive(Parser, Debug)]
struct Args {
    /// name of input step file
    input_step_file: String,
    /// name of output mesh file
    #[arg(default_value = "output.obj")]
    output_mesh_file: String,
    /// output only edge file
    /// Ignored when outputting obj files
    #[arg(long("only-edge"))]
    only_edge: bool,
    /// output the condition of meshes
    /// Ignored when outputting vtk files
    #[arg(long("condition-check"))]
    condition_check: bool,
}

fn main() -> Result<()> {
    let Args {
        input_step_file,
        output_mesh_file,
        only_edge,
        condition_check,
    } = Args::parse();

    println!("reading file...");
    let step_file = std::fs::read_to_string(&input_step_file)
        .with_context(|| format!("failed to read {input_step_file}"))?;
    let table = Table::from_step(&step_file).context("failed to parse STEP")?;
    println!("meshing...");
    let polyshells = step_to_mesh(&table)?;

    let path: &Path = output_mesh_file.as_ref();
    match path.extension().and_then(|extension| extension.to_str()) {
        Some("obj") => output_obj(&polyshells, path, condition_check)?,
        Some("vtu") if only_edge => output_vtk_edge(polyshells, path)?,
        Some("vtu") => output_vtk(polyshells, path)?,
        extension => bail!("unsupported output extension: {extension:?}"),
    }
    Ok(())
}

fn step_to_mesh(table: &Table) -> Result<Vec<MeshedCShell>> {
    let assy = table.step_assy().context("failed to build STEP assembly")?;

    let node_map = |ProductEntity {
                        shape: indices,
                        attrs: name,
                    }: &ProductEntity| {
        let shape = indices
            .iter()
            .filter_map(|idx| {
                let shells = if let Some(step_solid) = table.manifold_solid_brep.get(idx) {
                    table
                        .to_compressed_solid(step_solid)
                        .map_err(|err| eprintln!("failed to convert solid: {err}"))
                        .ok()?
                        .boundaries
                } else {
                    let step_shells = table.shell_based_surface_model.get(idx)?;
                    table
                        .to_compressed_shells(step_shells)
                        .map_err(|err| eprintln!("failed to convert shells: {err}"))
                        .ok()?
                };
                let meshed_shells = shells
                    .into_iter()
                    .map(|shell| {
                        let bdd = shell_bounding_box(&shell);
                        let tol = f64::max(bdd.diameter() * 0.001, KERNEL_EPSILON);
                        shell.robust_triangulation(tol)
                    })
                    .collect::<Vec<_>>();
                Some(meshed_shells)
            })
            .collect::<Vec<_>>();
        NodeEntity {
            shape,
            attrs: name.clone(),
        }
    };

    // Keep the original fallible STEP transforms while the expensive geometry
    // conversion is parallelized. Each path is converted transactionally below.
    let meshed_assy = assy.par_map(node_map, Clone::clone);
    let mut result = Vec::new();

    for top in meshed_assy.top_nodes() {
        for path in meshed_assy.paths_iter(top.index()) {
            let matrix = path.edges().iter().try_fold(
                Matrix4::identity(),
                |matrix, edge| -> Result<Matrix4> {
                    let transform = Matrix4::try_from(&edge.entity().matrix)
                        .map_err(|error| anyhow!("invalid STEP assembly transform: {error}"))?;
                    Ok(matrix * transform)
                },
            )?;
            let shapes = path
                .terminal_node()
                .map_or_else(Vec::new, |node| node.shape().clone());
            for mut shell in shapes.into_iter().flatten() {
                for vertex in &mut shell.vertices {
                    *vertex = matrix.transform_point(*vertex);
                }
                for edge in &mut shell.edges {
                    edge.curve.transform_by(matrix);
                }
                for face in &mut shell.faces {
                    if let Some(surface) = face.surface.as_mut() {
                        surface.transform_by(matrix);
                    }
                }
                result.push(shell);
            }
        }
    }
    Ok(result)
}

/// Computes bounding box from vertices, edge curves, and bounded surface corners.
fn shell_bounding_box<C, S>(shell: &CompressedShell<Point3, C, S>) -> BoundingBox<Point3>
where
    C: ParametricCurve3D + BoundedCurve,
    S: ParametricSurface3D, {
    let mut bdd: BoundingBox<Point3> = shell.vertices.iter().collect();
    for edge in &shell.edges {
        let (t0, t1) = edge.curve.range_tuple();
        // Sample 5 points per edge to capture curvature extent.
        for i in 0..=4 {
            let t = t0 + (t1 - t0) * f64::from(i) / 4.0;
            bdd.push(edge.curve.evaluate(t));
        }
    }
    for face in &shell.faces {
        let (u_range, v_range) = face.surface.try_range_tuple();
        if let (Some((u0, u1)), Some((v0, v1))) = (u_range, v_range) {
            bdd.push(face.surface.evaluate(u0, v0));
            bdd.push(face.surface.evaluate(u1, v0));
            bdd.push(face.surface.evaluate(u0, v1));
            bdd.push(face.surface.evaluate(u1, v1));
            bdd.push(
                face.surface
                    .evaluate(f64::midpoint(u0, u1), f64::midpoint(v0, v1)),
            );
        }
    }
    bdd
}

fn output_obj(polyshells: &[MeshedCShell], path: &Path, condition_check: bool) -> Result<()> {
    let mut polymesh = PolygonMesh::default();
    for shell in polyshells {
        let mut poly = shell.to_polygon();
        poly.put_together_same_attrs(KERNEL_EPSILON * 50.0)
            .remove_degenerate_faces()
            .remove_unused_attrs();
        if condition_check {
            println!("{:?}", poly.shell_condition());
        }
        polymesh.merge(poly);
    }
    let obj_file = std::fs::File::create(path)
        .with_context(|| format!("failed to create {}", path.display()))?;
    obj::write(&polymesh, obj_file)
        .with_context(|| format!("failed to write {}", path.display()))?;
    Ok(())
}

fn output_vtk(polyshells: Vec<MeshedCShell>, path: &Path) -> Result<()> {
    use vtkio::model::{ByteOrder, DataSet, Vtk};

    let mut pieces = Vec::new();
    for CompressedShell {
        vertices,
        edges,
        faces,
        ..
    } in polyshells
    {
        let faces = faces
            .into_iter()
            .filter_map(|face| {
                Some(CompressedFace {
                    boundaries: face.boundaries,
                    orientation: face.orientation,
                    surface: face.surface?,
                })
            })
            .collect();
        let data_set = CompressedShell {
            vertices,
            edges,
            faces,
            vertex_stable_ids: None,
            edge_stable_ids: None,
            face_stable_ids: None,
        }
        .to_data_set();
        let DataSet::UnstructuredGrid {
            pieces: shell_pieces,
            ..
        } = data_set
        else {
            bail!("mesh conversion unexpectedly produced a non-unstructured VTK data set");
        };
        pieces.extend(shell_pieces);
    }

    let vtk = Vtk {
        version: (1, 0).into(),
        title: String::new(),
        byte_order: ByteOrder::LittleEndian,
        file_path: None,
        data: DataSet::UnstructuredGrid { meta: None, pieces },
    };

    let vtk_file = std::fs::File::create(path)
        .with_context(|| format!("failed to create {}", path.display()))?;
    vtk.write_xml(vtk_file)
        .with_context(|| format!("failed to write {}", path.display()))?;
    Ok(())
}

fn output_vtk_edge(polyshells: Vec<MeshedCShell>, path: &Path) -> Result<()> {
    use vtkio::model::{
        Attributes, ByteOrder, CellType, Cells, DataSet, IOBuffer, Piece, UnstructuredGridPiece,
        VertexNumbers, Vtk,
    };

    let mut pieces = Vec::new();
    for CompressedEdge { curve, .. } in polyshells.into_iter().flat_map(|shell| shell.edges) {
        if curve.0.is_empty() {
            continue;
        }
        let points = curve
            .0
            .iter()
            .copied()
            .flat_map(Into::<[f64; 3]>::into)
            .collect::<Vec<_>>();
        let len = u64::try_from(curve.0.len()).context("edge has too many vertices for VTK")?;
        let connectivity = (0..len).chain([0, len - 1]).collect();
        pieces.push(Piece::Inline(Box::new(UnstructuredGridPiece {
            points: IOBuffer::F64(points),
            cells: Cells {
                cell_verts: VertexNumbers::XML {
                    connectivity,
                    offsets: vec![len, len + 1, len + 2],
                },
                types: vec![CellType::PolyLine, CellType::Vertex, CellType::Vertex],
            },
            data: Attributes::default(),
        })));
    }

    let vtk = Vtk {
        version: (1, 0).into(),
        title: String::new(),
        byte_order: ByteOrder::LittleEndian,
        file_path: None,
        data: DataSet::UnstructuredGrid { meta: None, pieces },
    };

    let vtk_file = std::fs::File::create(path)
        .with_context(|| format!("failed to create {}", path.display()))?;
    vtk.write_xml(vtk_file)
        .with_context(|| format!("failed to write {}", path.display()))?;
    Ok(())
}

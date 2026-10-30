#![allow(
    clippy::similar_names,
    reason = "parameter grid and indexed mesh names encode paired geometric axes and topology"
)]

use std::iter;

use super::{
    Cdt, CdtEdge, EuclideanSpace, Faces, FixedVertexHandle, HashMap, InnerFaceIterator, InnerSpace,
    Instant, Invertible, MIN_ALLOWED_VALUE, MetricSpace, Origin, ParametricSurface3D, Point2,
    Point3, PolygonMesh, PreMeshableSurface, SPoint2, StandardAttributes, StandardVertex,
    StructuringFilter, TOLERANCE, TessellationPrimitiveMode, TessellationPrimitiveOptions,
    Triangulation, Vector2, Vector3, VertexIterator, Zero, array, boundary, log_mesh_trace,
};
use boundary::{PolyBoundary, UvKey, boundary_segment_parameter, uv_key};

const MAX_SURFACE_SPAN_INTERVALS: usize = 32;
const SURFACE_SPAN_TOLERANCE_FACTOR: f64 = 8.0;

/// Ensures the mesh winding direction globally matches the stored vertex normals.
///
/// Sums the area-weighted dot product of each face's geometric normal against its
/// vertex normals. If the total is negative the entire mesh is inverted, keeping
/// all faces consistently wound.
///
/// # Why the weight is the raw cross-product sum and not `FaceNormal`
///
/// Ledger class C15, measured 2026-08-01. `FaceNormal::new` NORMALIZES, so a
/// face of zero area answers `0/0 = NaN`. A sphere's parameter grid ends in a
/// row of such faces at each pole, and one NaN term makes the whole `vote`
/// NaN -- at which point `vote < 0.0` is FALSE and the mesh is never corrected.
///
/// `occt-sphere.step` measured **-523.253857** on a `+523.5988` ball for
/// exactly this reason: 158 of its 24,964 mesh faces (one pole strip) voted
/// NaN, the raw grid winding is `-Su x Sv`, and the correction that exists to
/// catch that never ran. The torus next to it in the same fixture set has no
/// degenerate row, voted finite, and was right all along.
///
/// The raw sum is also what the doc line above always CLAIMED -- it is
/// area-weighted, so a degenerate face contributes exactly zero instead of
/// poisoning the sum, and a large face outvotes a sliver.
fn ensure_winding_matches_normals(mesh: &mut PolygonMesh) {
    let positions = mesh.positions();
    let normals = mesh.normals();
    let mut vote = 0.0f64;
    for face in mesh.faces().face_iter() {
        let center = face
            .iter()
            .fold(Vector3::zero(), |sum, v| sum + positions[v.pos].to_vec())
            / crate::count_as_f64(face.len());
        // Twice the area-weighted normal of the (possibly non-planar) polygon.
        // Deliberately NOT normalized: see the note above.
        let geo = face
            .windows(2)
            .chain(iter::once([face[face.len() - 1], face[0]].as_ref()))
            .fold(Vector3::zero(), |sum, v| {
                let vec0 = positions[v[0].pos].to_vec() - center;
                let vec1 = positions[v[1].pos].to_vec() - center;
                sum + vec0.cross(vec1)
            });
        let vtx: Vector3 = face.iter().fold(Vector3::zero(), |sum, v| {
            sum + v.nor.map_or_else(Vector3::zero, |i| normals[i])
        });
        let term = geo.dot(vtx);
        if term.is_finite() {
            vote += term;
        }
    }
    if vote < 0.0 {
        mesh.invert();
    }
}

pub(super) fn spade_round(x: f64) -> f64 {
    if f64::abs(x) < MIN_ALLOWED_VALUE {
        0.0
    } else {
        x
    }
}

fn add_constraint_through_existing_vertices(
    triangulation: &mut Cdt,
    split_vertices: &[(FixedVertexHandle, Point2)],
    front: (FixedVertexHandle, Point2),
    back: (FixedVertexHandle, Point2),
) {
    if front.0 == back.0 {
        return;
    }
    let mut chain = split_vertices
        .iter()
        .filter_map(|(handle, point)| {
            if *handle == front.0 || *handle == back.0 {
                None
            } else {
                boundary_segment_parameter(*point, front.1, back.1)
                    .map(|parameter| (parameter, *handle))
            }
        })
        .collect::<Vec<_>>();
    chain.sort_by(|lhs, rhs| lhs.0.total_cmp(&rhs.0));
    let handles = iter::once(front.0)
        .chain(chain.into_iter().map(|(_, handle)| handle))
        .chain(iter::once(back.0))
        .collect::<Vec<_>>();
    handles.windows(2).for_each(|window| {
        if window[0] != window[1] && triangulation.can_add_constraint(window[0], window[1]) {
            triangulation.add_constraint(window[0], window[1]);
        }
    });
}

fn surface_is_nearly_planar(surface: &impl PreMeshableSurface, boundary: &PolyBoundary) -> bool {
    let u0 = boundary.uv_min.x;
    let u1 = boundary.uv_max.x;
    let v0 = boundary.uv_min.y;
    let v1 = boundary.uv_max.y;
    if [u0, u1, v0, v1].into_iter().all(f64::is_finite) {
        let center = (f64::midpoint(u0, u1), f64::midpoint(v0, v1));
        let normal = surface.normal(center.0, center.1);
        [(u0, v0), (u1, v0), (u0, v1), (u1, v1), center]
            .into_iter()
            .all(|(u, v)| normal.dot(surface.normal(u, v)) > 1.0 - 1.0e-6)
    } else {
        false
    }
}

/// Tessellates a bounded surface without any trimming.
///
/// Generates a structured grid from parameter division, then triangulates
/// each quad cell into two triangles. Skips CDT and inclusion tests entirely.
pub(super) fn untrimmed_tessellation<S>(
    surface: &S,
    range: ((f64, f64), (f64, f64)),
    tolerance: f64,
    primitive_mode: TessellationPrimitiveMode,
) -> PolygonMesh
where
    S: PreMeshableSurface,
{
    let (udiv, vdiv) = surface.parameter_division(range, tolerance);
    let nu = udiv.len();
    let nv = vdiv.len();
    let mut positions = Vec::with_capacity(nu * nv);
    let mut uv_coords = Vec::with_capacity(nu * nv);
    let mut normals = Vec::with_capacity(nu * nv);
    for u in &udiv {
        for v in &vdiv {
            positions.push(surface.evaluate(*u, *v));
            uv_coords.push(Vector2::new(*u, *v));
            normals.push(surface.normal(*u, *v));
        }
    }
    let idx = |i: usize, j: usize| -> usize { i * nv + j };
    let sv = |k: usize| -> StandardVertex { [k, k, k].into() };
    let (tri_faces, quad_faces) = if primitive_mode == TessellationPrimitiveMode::Triangles {
        let tri_faces: Vec<[StandardVertex; 3]> = (0..nu - 1)
            .flat_map(|i| {
                (0..nv - 1).flat_map(move |j| {
                    let a = idx(i, j);
                    let b = idx(i, j + 1);
                    let c = idx(i + 1, j + 1);
                    let d = idx(i + 1, j);
                    [[sv(a), sv(b), sv(c)], [sv(a), sv(c), sv(d)]]
                })
            })
            .collect();
        (tri_faces, Vec::new())
    } else {
        let quad_faces: Vec<[StandardVertex; 4]> = (0..nu - 1)
            .flat_map(|i| {
                (0..nv - 1).map(move |j| {
                    let a = idx(i, j);
                    let b = idx(i, j + 1);
                    let c = idx(i + 1, j + 1);
                    let d = idx(i + 1, j);
                    [sv(a), sv(b), sv(c), sv(d)]
                })
            })
            .collect();
        (Vec::new(), quad_faces)
    };
    let mut mesh = PolygonMesh::new_unchecked(
        StandardAttributes {
            positions,
            uv_coords,
            normals,
        },
        Faces::from_tri_and_quad_faces(tri_faces, quad_faces),
    );
    ensure_winding_matches_normals(&mut mesh);
    mesh
}

/// Tessellates one surface trimmed by polyline.
pub(super) fn trimming_tessellation<S>(
    surface: &S,
    polyboundary: &PolyBoundary,
    tolerance: f64,
    primitive_config: TessellationPrimitiveOptions,
    face_idx: usize,
) -> PolygonMesh
where
    S: PreMeshableSurface,
{
    let mut mesh = if primitive_config.mode == TessellationPrimitiveMode::IsoQuads {
        iso_quad_trimmed_tessellation(surface, polyboundary, tolerance, face_idx).unwrap_or_else(
            || cdt_trimming_tessellation(surface, polyboundary, tolerance, face_idx),
        )
    } else {
        let mut mesh = cdt_trimming_tessellation(surface, polyboundary, tolerance, face_idx);
        apply_primitive_mode(&mut mesh, primitive_config);
        mesh
    };
    ensure_winding_matches_normals(&mut mesh);
    mesh
}

fn cdt_trimming_tessellation<S>(
    surface: &S,
    polyboundary: &PolyBoundary,
    tolerance: f64,
    face_idx: usize,
) -> PolygonMesh
where
    S: PreMeshableSurface,
{
    let mut triangulation = Cdt::new();
    let mut boundary_map = HashMap::<FixedVertexHandle, Point3>::default();
    let insert_boundary_start = Instant::now();
    let (boundary_vertices, added_constraints, skipped_constraints) =
        polyboundary.insert_to(&mut triangulation, &mut boundary_map);
    log_mesh_trace(
        face_idx,
        "cdt_boundary",
        format!(
            "loops={} boundary_vertices={boundary_vertices} constraints={added_constraints} skipped_constraints={skipped_constraints}",
            polyboundary.loops.len(),
        ),
        insert_boundary_start,
    );
    let insert_surface_start = Instant::now();
    insert_surface(
        &mut triangulation,
        surface,
        polyboundary,
        tolerance,
        face_idx,
    );
    log_mesh_trace(
        face_idx,
        "cdt_surface",
        format!("vertices={}", triangulation.num_vertices()),
        insert_surface_start,
    );
    triangulation_into_polymesh(
        triangulation.vertices(),
        triangulation.inner_faces(),
        surface,
        polyboundary,
        &boundary_map,
    )
}

fn iso_quad_full_cells(udiv: &[f64], vdiv: &[f64], polyboundary: &PolyBoundary) -> Vec<bool> {
    let (nu, nv) = (udiv.len(), vdiv.len());
    let idx = |i: usize, j: usize| i * nv + j;
    let inside = udiv
        .iter()
        .flat_map(|u| {
            vdiv.iter()
                .map(move |v| polyboundary.include(Point2::new(*u, *v)))
        })
        .collect::<Vec<_>>();
    (0..nu - 1)
        .flat_map(|i| {
            let inside = &inside;
            (0..nv - 1).map(move |j| {
                let corner_indices = [idx(i, j), idx(i, j + 1), idx(i + 1, j + 1), idx(i + 1, j)];
                if !corner_indices.into_iter().all(|index| inside[index]) {
                    return false;
                }
                let (u0, u1) = (udiv[i], udiv[i + 1]);
                let (v0, v1) = (vdiv[j], vdiv[j + 1]);
                let center = Point2::new(f64::midpoint(u0, u1), f64::midpoint(v0, v1));
                let edge_midpoints = [
                    Point2::new(f64::midpoint(u0, u1), v0),
                    Point2::new(u1, f64::midpoint(v0, v1)),
                    Point2::new(f64::midpoint(u0, u1), v1),
                    Point2::new(u0, f64::midpoint(v0, v1)),
                ];
                polyboundary.include(center)
                    && edge_midpoints
                        .into_iter()
                        .all(|midpoint| polyboundary.include(midpoint))
            })
        })
        .collect()
}

fn iso_quad_interior_mesh<S>(
    surface: &S,
    udiv: &[f64],
    vdiv: &[f64],
    polyboundary: &PolyBoundary,
) -> Option<(PolygonMesh, Vec<bool>)>
where
    S: PreMeshableSurface,
{
    let nv = vdiv.len();
    let idx = |i: usize, j: usize| i * nv + j;
    let vertex = |index: usize| StandardVertex::from([index, index, index]);
    let full_cells = iso_quad_full_cells(udiv, vdiv, polyboundary);
    let nv_cells = nv - 1;
    let interior_quads = full_cells
        .iter()
        .enumerate()
        .filter(|(_, is_full)| **is_full)
        .map(|(cell_index, _)| {
            let i = cell_index / nv_cells;
            let j = cell_index % nv_cells;
            [
                vertex(idx(i, j)),
                vertex(idx(i, j + 1)),
                vertex(idx(i + 1, j + 1)),
                vertex(idx(i + 1, j)),
            ]
        })
        .collect::<Vec<_>>();
    if interior_quads.is_empty() {
        return None;
    }
    let positions = udiv
        .iter()
        .flat_map(|u| vdiv.iter().map(move |v| surface.evaluate(*u, *v)))
        .collect();
    let uv_coords = udiv
        .iter()
        .flat_map(|u| vdiv.iter().map(move |v| Vector2::new(*u, *v)))
        .collect();
    let normals = udiv
        .iter()
        .flat_map(|u| vdiv.iter().map(move |v| surface.normal(*u, *v)))
        .collect();
    let mesh = PolygonMesh::new_unchecked(
        StandardAttributes {
            positions,
            uv_coords,
            normals,
        },
        Faces::from_tri_and_quad_faces(Vec::new(), interior_quads),
    );
    Some((mesh, full_cells))
}

fn remove_full_cells_from_boundary_mesh(
    boundary_mesh: &mut PolygonMesh,
    udiv: &[f64],
    vdiv: &[f64],
    full_cells: &[bool],
) {
    let boundary_triangles = boundary_mesh
        .tri_faces()
        .iter()
        .copied()
        .filter(|triangle| {
            !triangle_is_in_full_cell(*triangle, boundary_mesh.uv_coords(), udiv, vdiv, full_cells)
        })
        .collect();
    let boundary_quads = boundary_mesh
        .quad_faces()
        .iter()
        .copied()
        .filter(|quad| {
            !quad_is_in_full_cell(*quad, boundary_mesh.uv_coords(), udiv, vdiv, full_cells)
        })
        .collect();
    let editor = boundary_mesh.editor();
    *editor.faces = Faces::from_tri_and_quad_faces(boundary_triangles, boundary_quads);
}

fn iso_quad_trimmed_tessellation<S>(
    surface: &S,
    polyboundary: &PolyBoundary,
    tolerance: f64,
    face_idx: usize,
) -> Option<PolygonMesh>
where
    S: PreMeshableSurface,
{
    let range = (
        (polyboundary.uv_min.x, polyboundary.uv_max.x),
        (polyboundary.uv_min.y, polyboundary.uv_max.y),
    );
    let division_start = Instant::now();
    let (udiv, vdiv) = densify_surface_divisions(
        &surface,
        range,
        surface.parameter_division(range, tolerance),
        tolerance,
    );
    log_mesh_trace(
        face_idx,
        "iso_division",
        format!("udiv={} vdiv={}", udiv.len(), vdiv.len()),
        division_start,
    );
    if udiv.len() < 2 || vdiv.len() < 2 {
        return None;
    }
    let (mut interior_mesh, full_cells) =
        iso_quad_interior_mesh(surface, &udiv, &vdiv, polyboundary)?;
    let mut boundary_mesh = cdt_trimming_tessellation(surface, polyboundary, tolerance, face_idx);
    remove_full_cells_from_boundary_mesh(&mut boundary_mesh, &udiv, &vdiv, &full_cells);
    interior_mesh.merge(boundary_mesh);
    Some(interior_mesh)
}

fn triangle_is_in_full_cell(
    triangle: [StandardVertex; 3],
    uv_coords: &[Vector2],
    udiv: &[f64],
    vdiv: &[f64],
    full_cells: &[bool],
) -> bool {
    triangle_center_uv(triangle, uv_coords)
        .is_some_and(|uv| uv_is_in_full_cell(uv, udiv, vdiv, full_cells))
}

fn quad_is_in_full_cell(
    quad: [StandardVertex; 4],
    uv_coords: &[Vector2],
    udiv: &[f64],
    vdiv: &[f64],
    full_cells: &[bool],
) -> bool {
    quad_center_uv(quad, uv_coords).is_some_and(|uv| uv_is_in_full_cell(uv, udiv, vdiv, full_cells))
}

fn triangle_center_uv(triangle: [StandardVertex; 3], uv_coords: &[Vector2]) -> Option<Vector2> {
    let uv0 = uv_coords[*triangle[0].uv.as_ref()?];
    let uv1 = uv_coords[*triangle[1].uv.as_ref()?];
    let uv2 = uv_coords[*triangle[2].uv.as_ref()?];
    Some((uv0 + uv1 + uv2) / 3.0)
}

fn quad_center_uv(quad: [StandardVertex; 4], uv_coords: &[Vector2]) -> Option<Vector2> {
    let uv0 = uv_coords[*quad[0].uv.as_ref()?];
    let uv1 = uv_coords[*quad[1].uv.as_ref()?];
    let uv2 = uv_coords[*quad[2].uv.as_ref()?];
    let uv3 = uv_coords[*quad[3].uv.as_ref()?];
    Some((uv0 + uv1 + uv2 + uv3) / 4.0)
}

fn uv_is_in_full_cell(uv: Vector2, udiv: &[f64], vdiv: &[f64], full_cells: &[bool]) -> bool {
    if vdiv.len() < 2 {
        false
    } else {
        let nv_cells = vdiv.len() - 1;
        if let (Some(i), Some(j)) = (
            parameter_cell_index(udiv, uv.x),
            parameter_cell_index(vdiv, uv.y),
        ) {
            full_cells.get(i * nv_cells + j).copied().unwrap_or(false)
        } else {
            false
        }
    }
}

fn parameter_cell_index(parameters: &[f64], value: f64) -> Option<usize> {
    if parameters.len() < 2 || value < parameters[0] || value > parameters[parameters.len() - 1] {
        None
    } else {
        let upper = parameters.partition_point(|sample| *sample <= value);
        if upper == 0 {
            None
        } else if upper >= parameters.len() {
            Some(parameters.len() - 2)
        } else {
            Some(upper - 1)
        }
    }
}

fn apply_primitive_mode(mesh: &mut PolygonMesh, primitive_config: TessellationPrimitiveOptions) {
    if primitive_config.mode != TessellationPrimitiveMode::Triangles {
        mesh.quadrangulate(
            primitive_config.plane_tolerance,
            primitive_config.score_tolerance,
        );
        if primitive_config.mode == TessellationPrimitiveMode::AllQuads {
            force_all_triangles_to_quads(mesh, primitive_config);
        }
    }
}

fn force_all_triangles_to_quads(
    mesh: &mut PolygonMesh,
    primitive_config: TessellationPrimitiveOptions,
) {
    if mesh.tri_faces().is_empty() {
        return;
    }
    let normal_blend_angle = primitive_config.normal_blend_angle;
    let triangles = mesh.tri_faces().clone();
    let mut quadrangles = mesh.quad_faces().clone();
    let mut midpoint_cache = HashMap::<(VertexKey, VertexKey), StandardVertex>::default();
    for triangle in triangles {
        let [vertex0, vertex1, vertex2] = triangle;
        let point0 = mesh.positions()[vertex0.pos];
        let point1 = mesh.positions()[vertex1.pos];
        let point2 = mesh.positions()[vertex2.pos];
        let point01 = point0.midpoint(point1);
        let point12 = point1.midpoint(point2);
        let point20 = point2.midpoint(point0);
        let center_point =
            Point3::from_vec((point0.to_vec() + point1.to_vec() + point2.to_vec()) / 3.0);
        let split_quads = [
            [point0, point01, center_point, point20],
            [point1, point12, center_point, point01],
            [point2, point20, center_point, point12],
        ];
        let quality_ok = split_quads
            .into_iter()
            .all(|quad| quad_passes_quality_gates(quad, primitive_config));
        if quality_ok {
            let midpoint01 = midpoint_for_edge(
                mesh,
                &mut midpoint_cache,
                vertex0,
                vertex1,
                normal_blend_angle,
            );
            let midpoint12 = midpoint_for_edge(
                mesh,
                &mut midpoint_cache,
                vertex1,
                vertex2,
                normal_blend_angle,
            );
            let midpoint20 = midpoint_for_edge(
                mesh,
                &mut midpoint_cache,
                vertex2,
                vertex0,
                normal_blend_angle,
            );
            let center = create_centroid_vertex(mesh, triangle, normal_blend_angle);
            quadrangles.extend([
                [vertex0, midpoint01, center, midpoint20],
                [vertex1, midpoint12, center, midpoint01],
                [vertex2, midpoint20, center, midpoint12],
            ]);
        } else {
            quadrangles.push([vertex0, vertex1, vertex2, vertex2]);
        }
    }
    {
        let editor = mesh.editor();
        *editor.faces = Faces::from_tri_and_quad_faces(Vec::new(), quadrangles);
    }
    ensure_winding_matches_normals(mesh);
}

type VertexKey = (usize, Option<usize>, Option<usize>);

fn midpoint_for_edge(
    mesh: &mut PolygonMesh,
    midpoint_cache: &mut HashMap<(VertexKey, VertexKey), StandardVertex>,
    vertex0: StandardVertex,
    vertex1: StandardVertex,
    normal_blend_angle: f64,
) -> StandardVertex {
    let key = edge_key(vertex0, vertex1);
    midpoint_cache.get(&key).copied().unwrap_or_else(|| {
        let midpoint = create_midpoint_vertex(mesh, vertex0, vertex1, normal_blend_angle);
        midpoint_cache.insert(key, midpoint);
        midpoint
    })
}

fn edge_key(vertex0: StandardVertex, vertex1: StandardVertex) -> (VertexKey, VertexKey) {
    let key0 = (vertex0.pos, vertex0.uv, vertex0.nor);
    let key1 = (vertex1.pos, vertex1.uv, vertex1.nor);
    if key0 <= key1 {
        (key0, key1)
    } else {
        (key1, key0)
    }
}

fn quad_passes_quality_gates(
    quad: [Point3; 4],
    primitive_config: TessellationPrimitiveOptions,
) -> bool {
    let area = quad_area(quad);
    let convex = is_convex_quad(quad);
    let corner_angles_ok = quad.iter().enumerate().all(|(index, _)| {
        corner_angle(quad, index)
            .is_some_and(|angle| angle <= primitive_config.maximum_corner_angle)
    });
    area >= primitive_config.minimum_area && convex && corner_angles_ok
}

fn quad_area(quad: [Point3; 4]) -> f64 {
    triangle_area(quad[0], quad[1], quad[2]) + triangle_area(quad[0], quad[2], quad[3])
}

fn triangle_area(point0: Point3, point1: Point3, point2: Point3) -> f64 {
    let vector0 = point1.to_vec() - point0.to_vec();
    let vector1 = point2.to_vec() - point0.to_vec();
    vector0.cross(vector1).magnitude() * 0.5
}

fn is_convex_quad(quad: [Point3; 4]) -> bool {
    let diagonal0 = quad[1].to_vec() - quad[0].to_vec();
    let diagonal1 = quad[2].to_vec() - quad[0].to_vec();
    let normal = diagonal0.cross(diagonal1);
    if normal.so_small() {
        false
    } else {
        let (has_positive, has_negative, has_nonzero) = quad
            .iter()
            .enumerate()
            .map(|(index, _)| corner_turn(quad, index, normal))
            .fold(
                (false, false, false),
                |(positive, negative, nonzero), turn| {
                    (
                        positive || turn > TOLERANCE,
                        negative || turn < -TOLERANCE,
                        nonzero || f64::abs(turn) > TOLERANCE,
                    )
                },
            );
        has_nonzero && !(has_positive && has_negative)
    }
}

fn corner_turn(quad: [Point3; 4], index: usize, normal: Vector3) -> f64 {
    let previous = quad[(index + 3) % 4];
    let current = quad[index];
    let next = quad[(index + 1) % 4];
    let edge0 = current.to_vec() - previous.to_vec();
    let edge1 = next.to_vec() - current.to_vec();
    edge0.cross(edge1).dot(normal)
}

fn corner_angle(quad: [Point3; 4], index: usize) -> Option<f64> {
    let previous = quad[(index + 3) % 4];
    let current = quad[index];
    let next = quad[(index + 1) % 4];
    let vector0 = previous.to_vec() - current.to_vec();
    let vector1 = next.to_vec() - current.to_vec();
    let length0 = vector0.magnitude();
    let length1 = vector1.magnitude();
    if length0.so_small() || length1.so_small() {
        None
    } else {
        let cosine = (vector0.dot(vector1) / (length0 * length1)).clamp(-1.0, 1.0);
        Some(cosine.acos())
    }
}

fn create_midpoint_vertex(
    mesh: &mut PolygonMesh,
    vertex0: StandardVertex,
    vertex1: StandardVertex,
    normal_blend_angle: f64,
) -> StandardVertex {
    let point0 = mesh.positions()[vertex0.pos];
    let point1 = mesh.positions()[vertex1.pos];
    let position = point0.midpoint(point1);
    let position_index = mesh.positions().len();
    mesh.push_position(position);

    let uv_coord = match (vertex0.uv, vertex1.uv) {
        (Some(index0), Some(index1)) => {
            Some((mesh.uv_coords()[index0] + mesh.uv_coords()[index1]) / 2.0)
        }
        _ => None,
    };
    let uv_index = uv_coord.map(|uv_coord| {
        let index = mesh.uv_coords().len();
        mesh.push_uv_coord(uv_coord);
        index
    });

    let normal = match (vertex0.nor, vertex1.nor) {
        (Some(index0), Some(index1)) => Some(blend_normals(
            mesh.normals()[index0],
            mesh.normals()[index1],
            normal_blend_angle,
        )),
        _ => None,
    };
    let normal_index = normal.map(|normal| {
        let index = mesh.normals().len();
        mesh.extend_normals([normal]);
        index
    });

    StandardVertex {
        pos: position_index,
        uv: uv_index,
        nor: normal_index,
    }
}

fn create_centroid_vertex(
    mesh: &mut PolygonMesh,
    triangle: [StandardVertex; 3],
    normal_blend_angle: f64,
) -> StandardVertex {
    let [vertex0, vertex1, vertex2] = triangle;
    let point0 = mesh.positions()[vertex0.pos];
    let point1 = mesh.positions()[vertex1.pos];
    let point2 = mesh.positions()[vertex2.pos];
    let position = Point3::from_vec((point0.to_vec() + point1.to_vec() + point2.to_vec()) / 3.0);
    let position_index = mesh.positions().len();
    mesh.push_position(position);

    let uv_coord = match (vertex0.uv, vertex1.uv, vertex2.uv) {
        (Some(index0), Some(index1), Some(index2)) => Some(
            (mesh.uv_coords()[index0] + mesh.uv_coords()[index1] + mesh.uv_coords()[index2]) / 3.0,
        ),
        _ => None,
    };
    let uv_index = uv_coord.map(|uv_coord| {
        let index = mesh.uv_coords().len();
        mesh.push_uv_coord(uv_coord);
        index
    });

    let normal = match (vertex0.nor, vertex1.nor, vertex2.nor) {
        (Some(index0), Some(index1), Some(index2)) => {
            let normal01 = blend_normals(
                mesh.normals()[index0],
                mesh.normals()[index1],
                normal_blend_angle,
            );
            Some(blend_normals(
                normal01,
                mesh.normals()[index2],
                normal_blend_angle,
            ))
        }
        _ => None,
    };
    let normal_index = normal.map(|normal| {
        let index = mesh.normals().len();
        mesh.extend_normals([normal]);
        index
    });

    StandardVertex {
        pos: position_index,
        uv: uv_index,
        nor: normal_index,
    }
}

fn blend_normals(normal0: Vector3, normal1: Vector3, normal_blend_angle: f64) -> Vector3 {
    let cosine = normal0.dot(normal1);
    let min_cosine = f64::cos(normal_blend_angle);
    let blended = if cosine < min_cosine {
        normal0
    } else {
        normal0 + normal1
    };
    let magnitude = blended.magnitude();
    if magnitude.so_small() {
        normal0
    } else {
        blended / magnitude
    }
}

fn surface_axis_point<S: ParametricSurface3D>(
    surface: &S,
    axis: usize,
    axis_value: f64,
    other_value: f64,
) -> Point3 {
    if axis == 0 {
        surface.evaluate(axis_value, other_value)
    } else {
        surface.evaluate(other_value, axis_value)
    }
}

fn division_axis_length<S: ParametricSurface3D>(
    surface: &S,
    divisions: &[f64],
    axis: usize,
    other_value: f64,
) -> f64 {
    divisions
        .windows(2)
        .map(|window| {
            surface_axis_point(surface, axis, window[0], other_value).distance(surface_axis_point(
                surface,
                axis,
                window[1],
                other_value,
            ))
        })
        .sum()
}

fn desired_span_intervals(length: f64, tolerance: f64) -> usize {
    let max_segment_length = f64::max(tolerance * SURFACE_SPAN_TOLERANCE_FACTOR, TOLERANCE);
    if length.is_finite() && max_segment_length.is_finite() && max_segment_length > 0.0 {
        let required = f64::ceil(length / max_segment_length).max(1.0);
        let mut intervals = 1usize;
        while intervals < MAX_SURFACE_SPAN_INTERVALS && crate::count_as_f64(intervals) < required {
            intervals += 1;
        }
        intervals
    } else {
        1
    }
}

fn densify_divisions(
    mut divisions: Vec<f64>,
    range: (f64, f64),
    desired_intervals: usize,
) -> Vec<f64> {
    if divisions.len() < 2 {
        divisions = vec![range.0, range.1];
    }
    if let (Some(start), Some(end)) = (divisions.first().copied(), divisions.last().copied())
        && desired_intervals > divisions.len().saturating_sub(1)
        && start.is_finite()
        && end.is_finite()
        && !(end - start).so_small()
    {
        divisions = (0..=desired_intervals)
            .map(|index| {
                start
                    + (end - start) * crate::count_as_f64(index)
                        / crate::count_as_f64(desired_intervals)
            })
            .collect();
    }
    divisions
}

fn densify_surface_divisions<S: ParametricSurface3D>(
    surface: &S,
    range: ((f64, f64), (f64, f64)),
    divisions: (Vec<f64>, Vec<f64>),
    tolerance: f64,
) -> (Vec<f64>, Vec<f64>) {
    let (udiv, vdiv) = divisions;
    let vmid = f64::midpoint(range.1.0, range.1.1);
    let umid = f64::midpoint(range.0.0, range.0.1);
    let u_length = division_axis_length(surface, &udiv, 0, vmid);
    let v_length = division_axis_length(surface, &vdiv, 1, umid);
    (
        densify_divisions(udiv, range.0, desired_span_intervals(u_length, tolerance)),
        densify_divisions(vdiv, range.1, desired_span_intervals(v_length, tolerance)),
    )
}

fn insert_surface_grid(
    triangulation: &mut Cdt,
    polyline: &PolyBoundary,
    udiv: Vec<f64>,
    vdiv: &[f64],
) -> usize {
    let split_vertices = triangulation
        .vertices()
        .map(|vertex| {
            let point = *vertex.as_ref();
            (vertex.fix(), Point2::new(point.x, point.y))
        })
        .collect::<Vec<_>>();
    let insert_res: Vec<Vec<Option<_>>> = udiv
        .into_iter()
        .map(|u| {
            vdiv.iter()
                .map(|v| {
                    let uv = Point2::new(u, *v);
                    if polyline.include(uv) {
                        triangulation
                            .insert(SPoint2::new(u, *v))
                            .ok()
                            .map(|handle| (handle, uv))
                    } else {
                        None
                    }
                })
                .collect()
        })
        .collect();
    insert_res.windows(2).for_each(|rows| {
        rows[0].windows(2).zip(&rows[1]).for_each(|(pair, next)| {
            if let Some(front) = pair[0] {
                if let Some(back) = pair[1] {
                    add_constraint_through_existing_vertices(
                        triangulation,
                        &split_vertices,
                        front,
                        back,
                    );
                }
                if let Some(next) = next {
                    add_constraint_through_existing_vertices(
                        triangulation,
                        &split_vertices,
                        front,
                        *next,
                    );
                }
            }
        });
        let last = rows[0].len() - 1;
        if let (Some(front), Some(back)) = (rows[0][last], rows[1][last]) {
            add_constraint_through_existing_vertices(triangulation, &split_vertices, front, back);
        }
    });
    insert_res
        .iter()
        .flatten()
        .filter(|vertex| vertex.is_some())
        .count()
}

/// Inserts parameter divisions into triangulation.
fn insert_surface(
    triangulation: &mut Cdt,
    surface: impl PreMeshableSurface,
    polyline: &PolyBoundary,
    tolerance: f64,
    face_idx: usize,
) {
    let range = (
        (polyline.uv_min.x, polyline.uv_max.x),
        (polyline.uv_min.y, polyline.uv_max.y),
    );
    if surface_is_nearly_planar(&surface, polyline) {
        log_mesh_trace(
            face_idx,
            "parameter_division",
            format!(
                "range=({:.6},{:.6})x({:.6},{:.6}) planar_skip=true",
                range.0.0, range.0.1, range.1.0, range.1.1,
            ),
            Instant::now(),
        );
        log_mesh_trace(
            face_idx,
            "insert_surface",
            "inserted_points=0",
            Instant::now(),
        );
        return;
    }
    let division_start = Instant::now();
    let (udiv, vdiv) = densify_surface_divisions(
        &surface,
        range,
        surface.parameter_division(range, tolerance),
        tolerance,
    );
    log_mesh_trace(
        face_idx,
        "parameter_division",
        format!(
            "range=({:.6},{:.6})x({:.6},{:.6}) udiv={} vdiv={}",
            range.0.0,
            range.0.1,
            range.1.0,
            range.1.1,
            udiv.len(),
            vdiv.len(),
        ),
        division_start,
    );
    let insert_start = Instant::now();
    let inserted_points = insert_surface_grid(triangulation, polyline, udiv, &vdiv);
    log_mesh_trace(
        face_idx,
        "insert_surface",
        format!("inserted_points={inserted_points}"),
        insert_start,
    );
}

/// Converts triangulation into `PolygonMesh`.
fn triangulation_into_polymesh<'a>(
    vertices: VertexIterator<'a, SPoint2, (), CdtEdge<()>, ()>,
    triangles: InnerFaceIterator<'a, SPoint2, (), CdtEdge<()>, ()>,
    surface: &impl ParametricSurface3D,
    polyline: &PolyBoundary,
    boundary_map: &HashMap<FixedVertexHandle, Point3>,
) -> PolygonMesh {
    let mut positions = Vec::<Point3>::new();
    let mut uv_coords = Vec::<Vector2>::new();
    let mut normals = Vec::<Vector3>::new();
    let mut surface_point_cache = HashMap::<UvKey, Point3>::default();
    let mut normal_cache = HashMap::<UvKey, Vector3>::default();
    let vmap: HashMap<_, _> = vertices
        .enumerate()
        .map(|(i, v)| {
            let p = *v.as_ref();
            let uv = Point2::new(p.x, p.y);
            let key = uv_key(uv);
            let idx = v.fix();
            let point = boundary_map.get(&idx).copied().unwrap_or_else(|| {
                *surface_point_cache
                    .entry(key)
                    .or_insert_with(|| surface.evaluate(p.x, p.y))
            });
            let normal = *normal_cache
                .entry(key)
                .or_insert_with(|| surface.normal(p.x, p.y));
            positions.push(point);
            uv_coords.push(Vector2::new(p.x, p.y));
            normals.push(normal);
            (idx, i)
        })
        .collect();
    let tri_faces: Vec<[StandardVertex; 3]> = triangles
        .map(|tri| tri.vertices())
        .filter(|tri| {
            const fn sp2cg(p: SPoint2) -> Point2 { Point2::new(p.x, p.y) }
            let tri = array![i => sp2cg(*tri[i].as_ref()); 3];
            let (a, b) = (tri[1] - tri[0], tri[2] - tri[0]);
            let c = tri[0] + (a + b) / 3.0;
            let area = a.x * b.y - a.y * b.x;
            polyline.include(c) && !area.so_small2()
        })
        .map(|tri| {
            let idcs = array![i => vmap[&tri[i].fix()]; 3];
            array![i => [idcs[i], idcs[i], idcs[i]].into(); 3]
        })
        .collect();
    PolygonMesh::new_unchecked(
        StandardAttributes {
            positions,
            uv_coords,
            normals,
        },
        Faces::from_tri_and_quad_faces(tri_faces, Vec::new()),
    )
}

use super::*;
use crate::attributes::{AttributeValue, SolidAttributes};
use crate::compress::{
    CompressedEdge, CompressedEdgeUse, CompressedTrimmedFace, CompressedTrimmedShell,
};
use monstertruck_geometry::prelude::{Line, ParameterCurve, Plane, Point2, Point3};

// A pentagonal prism whose fifth lateral face is a planar chamfer.  Inverting the chamfer removes
// that face and restores the sharp x=2 / y=2 carrier intersection.  This is deliberately a full
// closed solid so the transactional executor is tested against the runtime topology validator,
// not merely against array bookkeeping.
type Trim = ParameterCurve<Line<Point2>, Plane>;
type TestSolid = CompressedTrimmedSolid<Point3, Line<Point3>, Plane, Trim>;

#[derive(Clone, Copy)]
struct FixtureIds {
    vertices: [StableId; 10],
    edges: [StableId; 15],
    faces: [StableId; 7],
}

fn sid(raw: u64) -> StableId { StableId::new(raw) }

fn edge_use(index: usize, orientation: bool) -> CompressedEdgeUse<Trim> {
    CompressedEdgeUse {
        index,
        orientation,
        trim_curve: None,
    }
}

fn face(surface: Plane, uses: &[(usize, bool)]) -> CompressedTrimmedFace<Plane, Trim> {
    CompressedTrimmedFace {
        boundaries: vec![
            uses.iter()
                .map(|&(index, orientation)| edge_use(index, orientation))
                .collect(),
        ],
        orientation: true,
        surface,
    }
}

fn chamfered_prism() -> (TestSolid, FixtureIds) {
    let points = [
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(2.0, 0.0, 0.0),
        Point3::new(2.0, 1.5, 0.0),
        Point3::new(1.5, 2.0, 0.0),
        Point3::new(0.0, 2.0, 0.0),
        Point3::new(0.0, 0.0, 1.0),
        Point3::new(2.0, 0.0, 1.0),
        Point3::new(2.0, 1.5, 1.0),
        Point3::new(1.5, 2.0, 1.0),
        Point3::new(0.0, 2.0, 1.0),
    ];
    let edge_vertices = [
        (0, 1),
        (1, 2),
        (2, 3),
        (3, 4),
        (4, 0),
        (5, 6),
        (6, 7),
        (7, 8),
        (8, 9),
        (9, 5),
        (0, 5),
        (1, 6),
        (2, 7),
        (3, 8),
        (4, 9),
    ];
    let edges = edge_vertices
        .map(|(front, back)| CompressedEdge {
            vertices: (front, back),
            curve: Line(points[front], points[back]),
        })
        .to_vec();

    // Each shared edge is used once in each direction.  The treatment face is face 4.
    let faces = vec![
        face(
            Plane::new(points[0], points[1], points[4]),
            &[(0, true), (1, true), (2, true), (3, true), (4, true)],
        ),
        face(
            Plane::new(points[5], points[9], points[6]),
            &[(9, false), (8, false), (7, false), (6, false), (5, false)],
        ),
        face(
            Plane::new(points[0], points[5], points[1]),
            &[(0, false), (10, true), (5, true), (11, false)],
        ),
        face(
            Plane::new(points[1], points[6], points[2]),
            &[(1, false), (11, true), (6, true), (12, false)],
        ),
        face(
            Plane::new(points[2], points[7], points[3]),
            &[(2, false), (12, true), (7, true), (13, false)],
        ),
        face(
            Plane::new(points[3], points[8], points[4]),
            &[(3, false), (13, true), (8, true), (14, false)],
        ),
        face(
            Plane::new(points[4], points[9], points[0]),
            &[(4, false), (14, true), (9, true), (10, false)],
        ),
    ];

    let ids = FixtureIds {
        vertices: [101, 102, 103, 104, 105, 106, 107, 108, 109, 110].map(sid),
        edges: [
            201, 202, 203, 204, 205, 206, 207, 208, 209, 210, 211, 212, 213, 214, 215,
        ]
        .map(sid),
        faces: [301, 302, 303, 304, 305, 306, 307].map(sid),
    };
    let attributes = fixture_attributes(ids);

    (
        CompressedTrimmedSolid {
            boundaries: vec![CompressedTrimmedShell {
                vertices: points.to_vec(),
                edges,
                faces,
                vertex_stable_ids: Some(ids.vertices.to_vec()),
                edge_stable_ids: Some(ids.edges.to_vec()),
                face_stable_ids: Some(ids.faces.to_vec()),
            }],
            id_allocator: None,
            attributes: Some(attributes),
        },
        ids,
    )
}

fn fixture_attributes(ids: FixtureIds) -> SolidAttributes {
    let mut attributes = SolidAttributes::new();
    attributes.faces.set(
        "kind",
        ids.faces[4],
        AttributeValue::String("treatment".to_owned()),
    );
    attributes
        .edges
        .set("kind", ids.edges[12], AttributeValue::Bool(true));
    attributes
        .vertices
        .set("kind", ids.vertices[2], AttributeValue::Bool(true));
    attributes.faces.set(
        "selection",
        ids.faces[3],
        AttributeValue::IdSet(vec![
            ids.faces[4],
            ids.edges[12],
            ids.vertices[2],
            ids.edges[0],
        ]),
    );
    attributes
}

fn inverse_chamfer_patch(ids: FixtureIds) -> CompressedTrimmedTopologyPatch<Point3, Line<Point3>> {
    let sharp_bottom = Point3::new(2.0, 2.0, 0.0);
    let sharp_top = Point3::new(2.0, 2.0, 1.0);
    let p1 = Point3::new(2.0, 0.0, 0.0);
    let p4 = Point3::new(0.0, 2.0, 0.0);
    let p6 = Point3::new(2.0, 0.0, 1.0);
    let p9 = Point3::new(0.0, 2.0, 1.0);

    CompressedTrimmedTopologyPatch {
        remove_faces: vec![ids.faces[4]],
        new_vertices: vec![sharp_bottom, sharp_top],
        new_edges: vec![NewCompressedEdge {
            vertices: (PatchVertexRef::New(0), PatchVertexRef::New(1)),
            curve: Line(sharp_bottom, sharp_top),
        }],
        rewrite_edges: vec![
            RetainedCompressedEdge {
                edge: ids.edges[1],
                vertices: (
                    PatchVertexRef::Existing(ids.vertices[1]),
                    PatchVertexRef::New(0),
                ),
                curve: Line(p1, sharp_bottom),
            },
            RetainedCompressedEdge {
                edge: ids.edges[3],
                vertices: (
                    PatchVertexRef::New(0),
                    PatchVertexRef::Existing(ids.vertices[4]),
                ),
                curve: Line(sharp_bottom, p4),
            },
            RetainedCompressedEdge {
                edge: ids.edges[6],
                vertices: (
                    PatchVertexRef::Existing(ids.vertices[6]),
                    PatchVertexRef::New(1),
                ),
                curve: Line(p6, sharp_top),
            },
            RetainedCompressedEdge {
                edge: ids.edges[8],
                vertices: (
                    PatchVertexRef::New(1),
                    PatchVertexRef::Existing(ids.vertices[9]),
                ),
                curve: Line(sharp_top, p9),
            },
        ],
        face_splices: vec![
            // Bottom/top terminal edges disappear because the neighboring retained edges are
            // extended to the restored sharp vertices.
            CompressedFaceSplice {
                face: ids.faces[0],
                remove_edges: vec![ids.edges[2]],
                replacement: vec![],
            },
            CompressedFaceSplice {
                face: ids.faces[1],
                remove_edges: vec![ids.edges[7]],
                replacement: vec![],
            },
            // Both carriers share the one restored sharp edge with opposite local orientation.
            CompressedFaceSplice {
                face: ids.faces[3],
                remove_edges: vec![ids.edges[12]],
                replacement: vec![PatchEdgeUse {
                    edge: PatchEdgeRef::New(0),
                    orientation: false,
                }],
            },
            CompressedFaceSplice {
                face: ids.faces[5],
                remove_edges: vec![ids.edges[13]],
                replacement: vec![PatchEdgeUse {
                    edge: PatchEdgeRef::New(0),
                    orientation: true,
                }],
            },
        ],
        face_boundary_replacements: vec![],
        require_unused_edges: vec![ids.edges[2], ids.edges[7], ids.edges[12], ids.edges[13]],
    }
}

#[test]
fn inverse_chamfer_patch_preserves_identity_and_validates_closed_topology() -> anyhow::Result<()> {
    let (mut solid, ids) = chamfered_prism();
    // Verify the fixture itself is valid before exercising surgery.
    let _: crate::trimmed::TrimmedSolid<Point3, Line<Point3>, Plane, Trim> =
        solid.clone().try_into()?;

    let result = solid.apply_topology_patch(&inverse_chamfer_patch(ids))?;
    assert_eq!(result.new_vertex_ids.len(), 2);
    assert_eq!(result.new_edge_ids.len(), 1);
    assert!(
        result
            .new_vertex_ids
            .iter()
            .chain(&result.new_edge_ids)
            .all(|id| id.raw() > ids.faces[6].raw())
    );

    // Surviving topology keeps its imported identity.
    for id in [
        ids.vertices[0],
        ids.vertices[1],
        ids.vertices[4],
        ids.vertices[5],
        ids.vertices[6],
        ids.vertices[9],
        ids.edges[0],
        ids.edges[1],
        ids.edges[3],
        ids.edges[6],
        ids.edges[8],
        ids.faces[0],
        ids.faces[1],
        ids.faces[3],
        ids.faces[5],
    ] {
        assert!(matches!(
            solid.location_by_stable_id(id)?,
            StableIdLocation::Unique { .. }
        ));
    }

    // Treatment topology and now-unused tangent vertices are compacted away.
    for id in [
        ids.faces[4],
        ids.edges[2],
        ids.edges[7],
        ids.edges[12],
        ids.edges[13],
        ids.vertices[2],
        ids.vertices[3],
        ids.vertices[7],
        ids.vertices[8],
    ] {
        assert_eq!(solid.location_by_stable_id(id)?, StableIdLocation::Missing);
    }

    // Every use of geometry changed by this patch has an exact regenerated pcurve.
    let shell = &solid.boundaries[0];
    let changed_ids = [
        ids.edges[1],
        ids.edges[3],
        ids.edges[6],
        ids.edges[8],
        result.new_edge_ids[0],
    ];
    let edge_ids = shell
        .edge_stable_ids
        .as_deref()
        .ok_or_else(|| anyhow::anyhow!("edge StableIds disappeared"))?;
    for face in &shell.faces {
        for edge_use in face.boundaries.iter().flatten() {
            if changed_ids.contains(&edge_ids[edge_use.index]) {
                assert!(edge_use.trim_curve.is_some());
            }
        }
    }

    // Deleted identity cannot linger either as an attribute key or inside an IdSet.
    let attributes = solid
        .attributes
        .as_ref()
        .ok_or_else(|| anyhow::anyhow!("attributes disappeared"))?;
    assert!(attributes.faces.get("kind", ids.faces[4]).is_none());
    assert!(attributes.edges.get("kind", ids.edges[12]).is_none());
    assert!(attributes.vertices.get("kind", ids.vertices[2]).is_none());
    let AttributeValue::IdSet(selection) = attributes
        .faces
        .get("selection", ids.faces[3])
        .ok_or_else(|| anyhow::anyhow!("surviving selection attribute disappeared"))?
    else {
        anyhow::bail!("selection attribute changed type");
    };
    assert_eq!(selection, &vec![ids.edges[0]]);
    Ok(())
}

#[test]
fn complete_face_boundary_replacement_reuses_edges_and_regenerates_trims() -> anyhow::Result<()> {
    let (mut solid, ids) = chamfered_prism();
    let patch = CompressedTrimmedTopologyPatch {
        remove_faces: vec![],
        new_vertices: vec![],
        new_edges: vec![],
        rewrite_edges: vec![],
        face_splices: vec![],
        face_boundary_replacements: vec![CompressedFaceBoundaryReplacement {
            face: ids.faces[0],
            boundaries: vec![vec![
                PatchEdgeUse {
                    edge: PatchEdgeRef::Existing(ids.edges[2]),
                    orientation: true,
                },
                PatchEdgeUse {
                    edge: PatchEdgeRef::Existing(ids.edges[3]),
                    orientation: true,
                },
                PatchEdgeUse {
                    edge: PatchEdgeRef::Existing(ids.edges[4]),
                    orientation: true,
                },
                PatchEdgeUse {
                    edge: PatchEdgeRef::Existing(ids.edges[0]),
                    orientation: true,
                },
                PatchEdgeUse {
                    edge: PatchEdgeRef::Existing(ids.edges[1]),
                    orientation: true,
                },
            ]],
        }],
        require_unused_edges: vec![],
    };

    let result = solid.apply_topology_patch(&patch)?;
    assert_eq!(result.new_vertex_ids, Vec::<StableId>::new());
    assert_eq!(result.new_edge_ids, Vec::<StableId>::new());
    let StableIdLocation::Unique {
        boundary, index, ..
    } = solid.location_by_stable_id(ids.faces[0])?
    else {
        anyhow::bail!("retained face StableId became non-unique");
    };
    let shell = &solid.boundaries[boundary];
    let edge_ids = shell
        .edge_stable_ids
        .as_deref()
        .ok_or_else(|| anyhow::anyhow!("edge StableIds disappeared"))?;
    let wire = shell.faces[index]
        .boundaries
        .first()
        .ok_or_else(|| anyhow::anyhow!("replacement face lost its boundary"))?;
    assert_eq!(
        wire.iter()
            .map(|edge_use| edge_ids[edge_use.index])
            .collect::<Vec<_>>(),
        vec![
            ids.edges[2],
            ids.edges[3],
            ids.edges[4],
            ids.edges[0],
            ids.edges[1]
        ]
    );
    assert!(wire.iter().all(|edge_use| edge_use.trim_curve.is_some()));
    let _: crate::trimmed::TrimmedSolid<Point3, Line<Point3>, Plane, Trim> = solid.try_into()?;
    Ok(())
}

#[test]
fn incomplete_inverse_patch_rolls_back_transactionally() -> anyhow::Result<()> {
    let (mut solid, ids) = chamfered_prism();
    let before = solid.clone();
    let mut patch = inverse_chamfer_patch(ids);
    // Leave the bottom terminal treatment edge in place while still declaring that the normalized
    // inverse operation must eliminate it.  The executor must reject this partial heal.
    patch
        .face_splices
        .retain(|splice| splice.face != ids.faces[0]);

    let error = solid
        .apply_topology_patch(&patch)
        .err()
        .ok_or_else(|| anyhow::anyhow!("incomplete patch unexpectedly succeeded"))?;
    assert!(matches!(
        error,
        CompressedTopologyPatchError::RequiredEdgeStillUsed { edge }
            if edge == ids.edges[2]
    ));
    assert_eq!(solid, before);
    Ok(())
}

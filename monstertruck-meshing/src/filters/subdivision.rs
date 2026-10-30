use super::{EuclideanSpace, Faces, Point3, PolygonMesh, StandardVertex, Vector3, Zero};
use rustc_hash::FxHashMap as HashMap;
use std::f64::consts::PI;

/// subdivision surface algorithms
pub trait Subdivision {
    /// Extended Loop method
    ///
    /// # Remarks
    /// Confirm:
    /// - All faces are triangles.
    /// - `self.shell_condition()` is `Oriented` or `Closed` before use.
    ///
    /// This method does NOT check these conditions.
    fn loop_subdivision(&mut self) -> &mut Self;
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
struct Edge(usize, usize);
#[derive(Clone, Copy, Debug)]
struct EdgeInfo {
    idx: usize,
    first_wing: usize,
    second_wing: Option<usize>,
}

impl Edge {
    #[inline]
    const fn new(v0: usize, v1: usize) -> Self {
        if v0 <= v1 { Self(v0, v1) } else { Self(v1, v0) }
    }
}

impl EdgeInfo {
    #[inline]
    const fn new(idx: usize, first_wing: usize) -> Self {
        Self {
            idx,
            first_wing,
            second_wing: None,
        }
    }
    #[inline]
    const fn add_second_wing(&mut self, second_wing: usize) {
        self.second_wing = Some(second_wing);
    }
}

impl Subdivision for PolygonMesh {
    fn loop_subdivision(&mut self) -> &mut Self {
        let mut edge_map = HashMap::default();
        let mut vertex_adjacency = vec![Vec::new(); self.positions().len()];
        for face in self.tri_faces() {
            if face
                .iter()
                .any(|vertex| vertex.pos >= vertex_adjacency.len())
            {
                return self;
            }
            let face = [face[0].pos, face[1].pos, face[2].pos];
            add_vertex_edge(
                &mut edge_map,
                &mut vertex_adjacency,
                face[1],
                face[2],
                face[0],
            );
            add_vertex_edge(
                &mut edge_map,
                &mut vertex_adjacency,
                face[2],
                face[0],
                face[1],
            );
            add_vertex_edge(
                &mut edge_map,
                &mut vertex_adjacency,
                face[0],
                face[1],
                face[2],
            );
        }

        let Some(mut positions) = vertex_adjacency
            .iter()
            .enumerate()
            .map(|(vertex, adjacency)| {
                calc_new_position(vertex, adjacency, &edge_map, self.positions())
            })
            .collect::<Option<Vec<_>>>()
        else {
            return self;
        };
        let Some(new_edge_positions) =
            edge_positions(&edge_map, &vertex_adjacency, self.positions())
        else {
            return self;
        };
        positions.extend(new_edge_positions);

        let len = self.positions().len();
        let mut tri_faces = Vec::with_capacity(self.tri_faces().len() * 4);
        for v in self.tri_faces() {
            let Some(opposite_edge0) = edge_map.get(&Edge::new(v[1].pos, v[2].pos)) else {
                return self;
            };
            let Some(opposite_edge1) = edge_map.get(&Edge::new(v[2].pos, v[0].pos)) else {
                return self;
            };
            let Some(opposite_edge2) = edge_map.get(&Edge::new(v[0].pos, v[1].pos)) else {
                return self;
            };
            let e: [StandardVertex; 3] = [
                (opposite_edge0.idx + len).into(),
                (opposite_edge1.idx + len).into(),
                (opposite_edge2.idx + len).into(),
            ];
            tri_faces.extend([
                [v[0], e[2], e[1]],
                [e[0], v[2], e[1]],
                [e[0], e[2], v[1]],
                [e[0], e[1], e[2]],
            ]);
        }
        {
            let editor = self.editor();
            editor.attributes.positions = positions;
            *editor.faces = Faces::from_tri_and_quad_faces(tri_faces, Vec::new());
        }
        self
    }
}

fn add_vertex_edge(
    edge_map: &mut HashMap<Edge, EdgeInfo>,
    vertex_adjacency: &mut [Vec<usize>],
    v0: usize,
    v1: usize,
    v2: usize,
) {
    let idx = edge_map.len();
    edge_map
        .entry(Edge::new(v0, v1))
        .and_modify(|tuple| tuple.add_second_wing(v2))
        .or_insert_with(|| {
            vertex_adjacency[v0].push(v1);
            vertex_adjacency[v1].push(v0);
            EdgeInfo::new(idx, v2)
        });
}

#[derive(Clone, Copy, Debug)]
enum VertexBoundaryCondition {
    #[expect(
        dead_code,
        reason = "corner classification remains part of the subdivision boundary model for alternate filters"
    )]
    Corner(usize, usize),
    Boundary(usize, usize),
    Inner,
}

impl VertexBoundaryCondition {
    fn new(v: usize, adjacency: &[usize], edge_map: &HashMap<Edge, EdgeInfo>) -> Option<Self> {
        let mut binfo = (None, None);
        for &w in adjacency {
            let edge = edge_map.get(&Edge::new(v, w))?;
            binfo = match (edge.second_wing, binfo) {
                (None, (None, _)) => (Some(w), None),
                (None, (Some(x), None)) => (Some(x), Some(w)),
                _ => binfo,
            };
        }
        Some(if let (Some(x), Some(y)) = binfo {
            if adjacency.len() == 2 {
                Self::Corner(x, y)
            } else {
                Self::Boundary(x, y)
            }
        } else {
            Self::Inner
        })
    }
}

fn calc_new_position(
    v: usize,
    adjacency: &[usize],
    edge_map: &HashMap<Edge, EdgeInfo>,
    positions: &[Point3],
) -> Option<Point3> {
    use VertexBoundaryCondition::{Boundary, Corner, Inner};
    let point = *positions.get(v)?;
    if adjacency.is_empty() {
        return Some(point);
    }
    Some(
        match VertexBoundaryCondition::new(v, adjacency, edge_map)? {
            Corner(_, _) => point,
            Boundary(w0, w1) => {
                point + ((*positions.get(w0)? - point) + (*positions.get(w1)? - point)) / 8.0
            }
            Inner => {
                let alpha =
                    3.0 / 8.0 + f64::cos(2.0 * PI / crate::count_as_f64(adjacency.len())) / 4.0;
                let alpha = alpha.mul_add(-alpha, 5.0 / 8.0) / crate::count_as_f64(adjacency.len());
                let mut delta = Vector3::zero();
                for &w in adjacency {
                    delta += *positions.get(w)? - point;
                }
                point + delta * alpha
            }
        },
    )
}

fn edge_positions(
    edge_map: &HashMap<Edge, EdgeInfo>,
    vertex_adjacency: &[Vec<usize>],
    positions: &[Point3],
) -> Option<Vec<Point3>> {
    use VertexBoundaryCondition::{Boundary, Corner, Inner};
    let mut res = vec![Point3::origin(); edge_map.len()];
    for (edge, info) in edge_map {
        let Edge(v0, v1) = *edge;
        let EdgeInfo {
            idx,
            first_wing,
            second_wing,
        } = *info;
        let p0 = *positions.get(v0)?;
        let p1 = *positions.get(v1)?;
        res[idx] = match second_wing {
            Some(second_wing) => {
                let gamma = match VertexBoundaryCondition::new(v0, &vertex_adjacency[v0], edge_map)?
                {
                    Boundary(_, _) => {
                        5.0 / 8.0
                            - f64::cos(PI / crate::count_as_f64(vertex_adjacency[v0].len())) / 4.0
                    }
                    Inner => {
                        match VertexBoundaryCondition::new(v1, &vertex_adjacency[v1], edge_map)? {
                            Boundary(_, _) => {
                                3.0 / 8.0
                                    + f64::cos(PI / crate::count_as_f64(vertex_adjacency[v1].len()))
                                        / 4.0
                            }
                            Inner => 0.5,
                            Corner(_, _) => return None,
                        }
                    }
                    Corner(_, _) => return None,
                };
                p0 + (p1 - p0) * gamma
                    + (*positions.get(first_wing)? - p0) / 8.0
                    + (*positions.get(second_wing)? - p1) / 8.0
            }
            None => p0.midpoint(p1),
        };
    }
    Some(res)
}

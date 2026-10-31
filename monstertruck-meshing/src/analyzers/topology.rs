use super::{Faces, PolygonMesh, Vertex};
use itertools::Itertools;
use monstertruck_topology::shell::ShellCondition;
use rustc_hash::{FxHashMap as HashMap, FxHashSet as HashSet};

/// Extracts boundaries, and check shell condition: closed or orientation.
pub trait Topology {
    /// Returns a vector of all boundaries as line strip.
    fn extract_boundaries(&self) -> Vec<Vec<usize>>;
    /// Determines the shell conditions: non-regular, regular, oriented, or closed.\
    /// The complexity increases in proportion to the number of edges.
    ///
    /// Examples for each condition can be found on the page of
    /// [`ShellCondition`](https://docs.rs/monstertruck-topology/0.2.0/monstertruck_topology/shell/enum.ShellCondition.html).
    fn shell_condition(&self) -> ShellCondition;
}

#[derive(Clone, Debug)]
struct Boundaries {
    checked: HashSet<[usize; 2]>,
    boundary: HashMap<[usize; 2], bool>,
    condition: ShellCondition,
}

impl Boundaries {
    #[inline]
    fn new() -> Self {
        Self {
            checked: HashSet::default(),
            boundary: HashMap::default(),
            condition: ShellCondition::Oriented,
        }
    }
    #[inline]
    fn insert(&mut self, edge: [Vertex; 2]) {
        let ori = edge[0].pos < edge[1].pos;
        let edge = if ori {
            [edge[0].pos, edge[1].pos]
        } else {
            [edge[1].pos, edge[0].pos]
        };
        self.condition = self.condition
            & match (self.checked.insert(edge), self.boundary.insert(edge, ori)) {
                (true, None) => ShellCondition::Oriented,
                (false, None) | (true, Some(_)) => ShellCondition::Irregular,
                (false, Some(ori0)) => {
                    self.boundary.remove(&edge);
                    if ori == ori0 {
                        ShellCondition::Regular
                    } else {
                        ShellCondition::Oriented
                    }
                }
            };
    }

    #[inline]
    fn condition(&self) -> ShellCondition {
        if self.condition == ShellCondition::Oriented && self.boundary.is_empty() {
            ShellCondition::Closed
        } else {
            self.condition
        }
    }
}

impl FromIterator<[Vertex; 2]> for Boundaries {
    fn from_iter<I: IntoIterator<Item = [Vertex; 2]>>(iter: I) -> Self {
        let mut boundaries = Self::new();
        for edge in iter {
            boundaries.insert(edge);
        }
        boundaries
    }
}

fn face_edge_iter<T: Copy>(face: &[T]) -> impl Iterator<Item = [T; 2]> + '_ {
    face.iter().circular_tuple_windows().map(|(v, w)| [*v, *w])
}

impl Topology for Faces {
    fn extract_boundaries(&self) -> Vec<Vec<usize>> {
        let boundaries: Boundaries = self.face_iter().flat_map(face_edge_iter).collect();
        let mut vemap: HashMap<usize, usize> = boundaries
            .boundary
            .into_iter()
            .map(|(edge, ori)| {
                let [front, back] = edge;
                if ori { (front, back) } else { (back, front) }
            })
            .collect();

        let mut res = Vec::new();
        while !vemap.is_empty() {
            let mut wire = Vec::new();
            let Some((&front0, &front1)) = vemap.iter().next() else {
                break;
            };
            let front = (front0, front1);
            vemap.remove(&front.0);
            wire.push(front.0);
            let mut cursor = front.1;
            while cursor != front.0 {
                wire.push(cursor);
                cursor = vemap.remove(&cursor).unwrap_or(front.0);
            }
            res.push(wire);
        }
        res
    }
    fn shell_condition(&self) -> ShellCondition {
        self.face_iter()
            .flat_map(face_edge_iter)
            .collect::<Boundaries>()
            .condition()
    }
}

impl Topology for PolygonMesh {
    fn extract_boundaries(&self) -> Vec<Vec<usize>> { self.faces().extract_boundaries() }
    fn shell_condition(&self) -> ShellCondition { self.faces().shell_condition() }
}

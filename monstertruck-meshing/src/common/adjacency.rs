use super::{Faces, Vertex};
use rustc_hash::{FxHashMap as HashMap, FxHashSet as HashSet};

pub trait Adjacency {
    #[allow(
        dead_code,
        reason = "the trait keeps vertex adjacency alongside the currently used face adjacency operation"
    )]
    fn vertex_adjacency(&self, num_of_vertices: usize) -> Vec<Vec<usize>>;
    /// create the adjacency list of the faces
    fn face_adjacency(&self, use_normal: bool) -> Vec<Vec<usize>>;
}

impl Adjacency for Faces {
    fn vertex_adjacency(&self, num_of_vertices: usize) -> Vec<Vec<usize>> {
        let mut already = HashSet::default();
        let mut res = vec![Vec::new(); num_of_vertices];
        for face in self.face_iter() {
            let Some((&first_vertex, rest)) = face.split_first() else {
                continue;
            };
            let Some(&last_vertex) = face.last() else {
                continue;
            };
            rest.iter()
                .scan(first_vertex, |previous, &current| {
                    let edge = (previous.pos, current.pos);
                    *previous = current;
                    Some(edge)
                })
                .chain(std::iter::once((last_vertex.pos, first_vertex.pos)))
                .for_each(|edge| {
                    if edge.0 >= num_of_vertices || edge.1 >= num_of_vertices {
                        return;
                    }
                    let edge = if edge.0 < edge.1 {
                        edge
                    } else {
                        (edge.1, edge.0)
                    };
                    if already.insert(edge) {
                        res[edge.0].push(edge.1);
                        res[edge.1].push(edge.0);
                    }
                });
        }
        res
    }
    fn face_adjacency(&self, use_normal: bool) -> Vec<Vec<usize>> {
        let len = self.len();
        let mut face_adjacency = vec![Vec::<usize>::new(); len];
        let mut edge_face_map = HashMap::default();
        for (i, face) in self.face_iter().enumerate() {
            let Some((&first, rest)) = face.split_first() else {
                continue;
            };
            let Some(&last) = face.last() else {
                continue;
            };
            rest.iter()
                .scan(first, |previous, &current| {
                    let edge = [*previous, current];
                    *previous = current;
                    Some(edge)
                })
                .chain(std::iter::once([last, first]))
                .for_each(|edge| {
                    signup_adjacency(
                        i,
                        edge[0],
                        edge[1],
                        &mut face_adjacency,
                        &mut edge_face_map,
                        use_normal,
                    );
                });
        }
        face_adjacency
    }
}

fn signup_adjacency(
    i: usize,
    v0: Vertex,
    v1: Vertex,
    face_adjacency: &mut [Vec<usize>],
    edge_face_map: &mut HashMap<[(usize, Option<usize>); 2], usize>,
    use_normal: bool,
) {
    let edge = match (v0.pos < v1.pos, use_normal) {
        (true, true) => [(v0.pos, v0.nor), (v1.pos, v1.nor)],
        (false, true) => [(v1.pos, v1.nor), (v0.pos, v0.nor)],
        (true, false) => [(v0.pos, None), (v1.pos, None)],
        (false, false) => [(v1.pos, None), (v0.pos, None)],
    };
    match edge_face_map.get(&edge) {
        Some(j) => {
            face_adjacency[i].push(*j);
            face_adjacency[*j].push(i);
        }
        None => {
            edge_face_map.insert(edge, i);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn malformed_deserialized_faces_do_not_panic_adjacency_building() {
        let faces: Faces = serde_json::from_value(serde_json::json!({
            "tri_faces": [],
            "quad_faces": [],
            "other_faces": [
                [],
                [
                    {"pos": 0, "uv": null, "nor": null},
                    {"pos": 9, "uv": null, "nor": null},
                    {"pos": 1, "uv": null, "nor": null}
                ]
            ]
        }))
        .expect("fixture Faces JSON is valid");

        assert_eq!(faces.vertex_adjacency(2), vec![vec![1], vec![0]]);
        assert_eq!(faces.face_adjacency(false), vec![Vec::<usize>::new(); 2]);
    }
}

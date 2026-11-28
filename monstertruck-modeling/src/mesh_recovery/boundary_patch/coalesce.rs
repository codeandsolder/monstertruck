//! Local convex-cell coalescence for otherwise uncertifiable source facets.
//!
//! Internal diagonals are not geometric constraints unless prescribed.
//! A merged cell must cover exactly the projected union of its source
//! polygons, retain every source vertex, and pay for nonplanarity.

use super::{BoundaryPatchApproximation, BoundaryPatchError, BoundaryPatchOptions};
use crate::{
    Edge, Point3, Shell, Vector3,
    mesh_recovery::{
        BoundaryFacetApproximation, BoundaryFacetError, BoundaryFacetInput, EdgeIncidence,
        approximate_boundary_facet, edge_key, source_edge_incidence,
        trace_oriented_boundary_cycles,
    },
};
use monstertruck_topology::shell::ShellCondition;
use robust::{Coord, orient2d};
use rustc_hash::FxHashMap;
use std::collections::{BTreeMap, BTreeSet};

pub(super) struct AssemblyInput<'a> {
    pub source: &'a [Point3],
    pub polygons: &'a [&'a [usize]],
    pub edges: &'a BTreeMap<[usize; 2], Edge>,
    pub incidence: &'a EdgeIncidence,
    pub protected: &'a BTreeSet<[usize; 2]>,
    pub options: BoundaryPatchOptions<'a>,
    pub vertex_displacement: f64,
    pub expected_closed: bool,
}

enum CellGeometry {
    Original(usize),
    Union {
        members: Vec<usize>,
        boundary: Vec<usize>,
        source_penalty: f64,
    },
}

impl CellGeometry {
    fn members(&self) -> &[usize] {
        match self {
            Self::Original(id) => std::slice::from_ref(id),
            Self::Union { members, .. } => members,
        }
    }

    fn boundary<'a>(&'a self, polygons: &'a [&[usize]]) -> &'a [usize] {
        match self {
            Self::Original(id) => polygons[*id],
            Self::Union { boundary, .. } => boundary,
        }
    }

    const fn penalty(&self) -> f64 {
        match self {
            Self::Original(_) => 0.0,
            Self::Union { source_penalty, .. } => *source_penalty,
        }
    }
}

struct Cell {
    geometry: CellGeometry,
    approximation: Option<BoundaryFacetApproximation>,
}

struct Repair {
    cells: Vec<usize>,
    geometry: CellGeometry,
    approximation: BoundaryFacetApproximation,
}

// Drop the dominant normal coordinate. Orientation predicates act on
// original coordinate pairs, not rounded arbitrary-plane projections.
struct Projection {
    axis: usize,
    reverse: bool,
    normal: Vector3,
    origin: Point3,
}

impl Projection {
    fn new(source: &[Point3], boundary: &[usize]) -> Option<Self> {
        let origin = source[*boundary.first()?];
        let normal = boundary
            .iter()
            .zip(boundary.iter().cycle().skip(1))
            .map(|(&a, &b)| (source[a] - origin).cross(&(source[b] - origin)))
            .fold(Vector3::new(0.0, 0.0, 0.0), |sum, term| sum + term);
        let length = normal.norm();
        if !length.is_finite() || length <= 0.0 {
            return None;
        }
        let normal = normal / length;
        let components = [normal.x, normal.y, normal.z];
        let axis = (0..3).max_by(|&a, &b| components[a].abs().total_cmp(&components[b].abs()))?;
        Some(Self {
            axis,
            reverse: components[axis] < 0.0,
            normal,
            origin,
        })
    }

    fn project(&self, p: Point3) -> Coord<f64> {
        let coords = [p.x, p.y, p.z];
        let mut x = coords[(self.axis + 1) % 3];
        let mut y = coords[(self.axis + 2) % 3];
        if self.reverse {
            std::mem::swap(&mut x, &mut y);
        }
        Coord { x, y }
    }

    // All vertices lie on the left of every oriented edge and at least
    // one lies strictly to the left. Repeated projected positions would
    // admit a degenerate or multiply traversed boundary, so reject them.
    fn convex(&self, source: &[Point3], polygon: &[usize]) -> bool {
        let points = polygon
            .iter()
            .map(|&id| self.project(source[id]))
            .collect::<Vec<_>>();
        if points.len() < 3
            || points.iter().enumerate().any(|(i, a)| {
                points[..i]
                    .iter()
                    .any(|b| (a.x - b.x).abs() == 0.0 && (a.y - b.y).abs() == 0.0)
            })
        {
            return false;
        }
        points
            .iter()
            .zip(points.iter().cycle().skip(1))
            .all(|(&a, &b)| {
                let mut positive = false;
                let valid = points.iter().all(|&p| {
                    let turn = orient2d(a, b, p);
                    positive |= turn > 0.0;
                    turn.is_finite() && turn >= 0.0
                });
                valid && positive
            })
    }

    // The projected polygon fans cover identical sets. Relative to the
    // same affine reference plane, the dropped-coordinate heights of
    // both fans lie in the convex range of the original vertex heights.
    // Charge that whole range, including arithmetic conditioning.
    fn height_width(&self, source: &[Point3], vertices: &[usize]) -> f64 {
        let component = [self.normal.x, self.normal.y, self.normal.z][self.axis].abs();
        let (low, high, scale) = vertices.iter().map(|&id| source[id] - self.origin).fold(
            (0.0_f64, 0.0_f64, 0.0_f64),
            |(lo, hi, scale), p| {
                let h = self.normal.dot(&p);
                (lo.min(h), hi.max(h), scale.max(p.norm()))
            },
        );
        (64.0 * f64::EPSILON).mul_add(scale, high - low) / component
    }
}

impl AssemblyInput<'_> {
    fn same_owner(&self, a: usize, b: usize) -> bool {
        self.options
            .owner_by_polygon
            .is_none_or(|owners| owners[a] == owners[b])
    }

    fn geometry(&self, members: Vec<usize>) -> Option<CellGeometry> {
        let incidence = source_edge_incidence(members.iter().map(|&id| self.polygons[id]));
        let mut next = FxHashMap::default();
        let mut incoming = BTreeSet::new();
        for (key, uses) in incidence {
            match uses.as_slice() {
                [(_, forward)] => {
                    let [a, b] = if *forward { key } else { [key[1], key[0]] };
                    if next.insert(a, b).is_some() || !incoming.insert(b) {
                        return None;
                    }
                }
                [(_, a), (_, b)] if a != b && !self.protected.contains(&key) => {}
                _ => return None,
            }
        }
        let [boundary]: [Vec<usize>; 1] = trace_oriented_boundary_cycles(next)?.try_into().ok()?;
        let vertices = members
            .iter()
            .flat_map(|&id| self.polygons[id].iter().copied())
            .collect::<BTreeSet<_>>();
        // Every original vertex remains a native boundary vertex, even
        // after a diagonal disappears. This preserves the vertex budget
        // and point incidence without a new interior-vertex proof.
        if vertices.len() != boundary.len() {
            return None;
        }
        let projection = Projection::new(self.source, &boundary)?;
        if !projection.convex(self.source, &boundary)
            || !members
                .iter()
                .all(|&id| projection.convex(self.source, self.polygons[id]))
        {
            return None;
        }
        // The positively oriented convex source polygons have cancelled
        // internal edge uses and the same simple convex outer boundary.
        // Their projected winding numbers sum to one inside and zero
        // outside, hence no gaps or overlapping interiors are possible.
        let source_penalty = projection.height_width(self.source, &boundary);
        if !source_penalty.is_finite() || source_penalty >= self.options.budget.chord {
            return None;
        }
        Some(CellGeometry::Union {
            members,
            boundary,
            source_penalty,
        })
    }

    fn boundary_input(&self, geometry: &CellGeometry) -> (Vec<Point3>, Vec<Edge>) {
        let indices = geometry.boundary(self.polygons);
        let points = indices.iter().map(|&id| self.source[id]).collect();
        let edges = indices
            .iter()
            .zip(indices.iter().cycle().skip(1))
            .map(|(&a, &b)| {
                let (key, forward) = edge_key(a, b);
                let edge = &self.edges[&key];
                if forward {
                    edge.clone()
                } else {
                    edge.inverse()
                }
            })
            .collect();
        (points, edges)
    }

    fn certify(
        &self,
        geometry: &CellGeometry,
    ) -> Result<BoundaryFacetApproximation, BoundaryFacetError> {
        let (points, edges) = self.boundary_input(geometry);
        let mut result = approximate_boundary_facet(
            &points,
            &edges,
            self.options.budget.chord - geometry.penalty(),
        )?;
        result.deviation_upper_bound += geometry.penalty();
        if !result.deviation_upper_bound.is_finite()
            || result.deviation_upper_bound > self.options.budget.chord
        {
            return Err(BoundaryFacetError::Uncertified);
        }
        Ok(result)
    }

    fn neighbors(
        &self,
        group: &[usize],
        cells: &[Option<Cell>],
        membership: &[usize],
    ) -> BTreeSet<usize> {
        let mut neighbors = BTreeSet::new();
        for id in group {
            let Some(cell) = &cells[*id] else {
                continue;
            };
            for &member in cell.geometry.members() {
                let polygon = self.polygons[member];
                for (&a, &b) in polygon.iter().zip(polygon.iter().cycle().skip(1)) {
                    let key = edge_key(a, b).0;
                    if self.protected.contains(&key) {
                        continue;
                    }
                    for &(other, _) in &self.incidence[&key] {
                        let owner = membership[other];
                        if self.same_owner(member, other) && group.binary_search(&owner).is_err() {
                            neighbors.insert(owner);
                        }
                    }
                }
            }
        }
        neighbors
    }

    // Best-first connected subsets, ordered by source polygon count.
    // Nonconvex intermediate groups may grow to a valid convex union.
    // Bound both queued hypotheses and attempted certificates.
    fn repair(
        &self,
        seed: usize,
        cells: &[Option<Cell>],
        membership: &[usize],
        remaining: &mut usize,
    ) -> Option<Repair> {
        let seed_count = cells[seed].as_ref()?.geometry.members().len();
        let limit = remaining.saturating_add(1);
        let mut visited = BTreeSet::from([vec![seed]]);
        let mut frontier = BTreeSet::from([(seed_count, vec![seed])]);
        while let Some((count, group)) = frontier.pop_first() {
            if group.len() > 1 {
                *remaining = remaining.checked_sub(1)?;
                let mut members = Vec::with_capacity(count);
                for &id in &group {
                    members.extend_from_slice(cells[id].as_ref()?.geometry.members());
                }
                members.sort_unstable();
                if let Some(geometry) = self.geometry(members)
                    && let Ok(approximation) = self.certify(&geometry)
                {
                    return Some(Repair {
                        cells: group,
                        geometry,
                        approximation,
                    });
                }
            }
            if *remaining == 0 {
                break;
            }
            for other in self.neighbors(&group, cells, membership) {
                if visited.len() >= limit {
                    break;
                }
                let mut next = group.clone();
                let index = next.binary_search(&other).err()?;
                next.insert(index, other);
                if visited.insert(next.clone()) {
                    let count = count + cells[other].as_ref()?.geometry.members().len();
                    frontier.insert((count, next));
                }
            }
        }
        None
    }

    fn failure(&self, geometry: &CellGeometry, error: BoundaryFacetError) -> BoundaryPatchError {
        let (points, edges) = self.boundary_input(geometry);
        BoundaryPatchError::Facet {
            polygon: geometry.members()[0],
            error,
            source_vertex_ids: geometry.boundary(self.polygons).to_vec(),
            input: Box::new(BoundaryFacetInput::capture(
                points,
                &edges,
                self.options.budget.chord - geometry.penalty(),
            )),
        }
    }

    pub(super) fn assemble(self) -> Result<BoundaryPatchApproximation, BoundaryPatchError> {
        use BoundaryPatchError::InvalidShell;
        let mut cells = (0..self.polygons.len())
            .map(|id| {
                Some(Cell {
                    geometry: CellGeometry::Original(id),
                    approximation: None,
                })
            })
            .collect::<Vec<_>>();
        let mut membership = (0..cells.len()).collect::<Vec<_>>();
        let mut remaining = self.options.max_coalescence_trials;
        for seed in 0..cells.len() {
            let Some(cell) = &cells[seed] else {
                continue;
            };
            if cell.approximation.is_some() {
                continue;
            }
            match self.certify(&cell.geometry) {
                Ok(approximation) => {
                    cells[seed].as_mut().ok_or(InvalidShell)?.approximation = Some(approximation);
                }
                Err(error) => {
                    let repair = self
                        .repair(seed, &cells, &membership, &mut remaining)
                        .ok_or_else(|| self.failure(&cell.geometry, error))?;
                    for id in repair.cells {
                        cells[id] = None;
                    }
                    for &member in repair.geometry.members() {
                        membership[member] = seed;
                    }
                    cells[seed] = Some(Cell {
                        geometry: repair.geometry,
                        approximation: Some(repair.approximation),
                    });
                }
            }
        }
        let mut faces = Vec::new();
        let mut source_face_ranges = vec![0..0; self.polygons.len()];
        let mut deviation_upper_bound = 0.0_f64;
        let mut seam_deviation_upper_bound = 0.0_f64;
        let mut coalesced_cells = 0;
        for cell in cells.into_iter().flatten() {
            let approximation = cell.approximation.ok_or(InvalidShell)?;
            deviation_upper_bound = deviation_upper_bound.max(approximation.deviation_upper_bound);
            seam_deviation_upper_bound =
                seam_deviation_upper_bound.max(approximation.seam_deviation_upper_bound);
            let first = faces.len();
            faces.extend(approximation.faces);
            coalesced_cells += usize::from(cell.geometry.members().len() > 1);
            for &member in cell.geometry.members() {
                source_face_ranges[member] = first..faces.len();
            }
        }
        let shell = faces.into_iter().collect::<Shell>();
        let expected = if self.expected_closed {
            ShellCondition::Closed
        } else {
            ShellCondition::Oriented
        };
        if shell.shell_condition() != expected
            || !shell.is_connected()
            || !shell.singular_vertices().is_empty()
            || !shell.is_geometric_consistent()
        {
            return Err(InvalidShell);
        }
        Ok(BoundaryPatchApproximation {
            shell,
            source_face_ranges,
            coalesced_cells,
            vertex_displacement: self.vertex_displacement,
            deviation_upper_bound,
            seam_deviation_upper_bound,
        })
    }
}

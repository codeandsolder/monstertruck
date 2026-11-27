//! Exact local reproduction of a boundary-constrained facet certificate.

use super::{BoundaryFacetApproximation, BoundaryFacetError, approximate_boundary_facet};
use crate::{Curve, Edge, Point3, builder};
use serde::{Deserialize, Serialize};

/// Complete geometric input for a single facet approximation.
///
/// Native vertices and oriented curves are listed in source boundary order.
/// Replaying creates one shared native vertex per corner, and a single edge
/// between each consecutive pair. No mesh discovery or carrier fitting is
/// needed.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BoundaryFacetInput {
    /// Original polygon points in boundary order.
    pub source: Vec<Point3>,
    /// Reconciled native corner positions, in the same order.
    pub native: Vec<Point3>,
    /// Exact oriented boundary curves, one per polygon side.
    pub curves: Vec<Curve>,
    /// Continuous facet deviation budget.
    pub budget: f64,
}

impl BoundaryFacetInput {
    pub(in crate::mesh_recovery) fn capture(
        source: Vec<Point3>,
        edges: &[Edge],
        budget: f64,
    ) -> Self {
        Self {
            source,
            native: edges.iter().map(|edge| edge.front().point()).collect(),
            curves: edges.iter().map(Edge::oriented_curve).collect(),
            budget,
        }
    }

    /// Reconstruct the original oriented geometric input and apply the same
    /// continuous certificate used by patch recovery.
    ///
    /// # Errors
    /// Rejects invalid incidence and the same geometric or numerical
    /// failures as [`approximate_boundary_facet`].
    pub fn approximate(&self) -> Result<BoundaryFacetApproximation, BoundaryFacetError> {
        use BoundaryFacetError::InvalidBoundary;
        if self.source.len() < 3
            || self.source.len() != self.native.len()
            || self.source.len() != self.curves.len()
            || !self
                .native
                .iter()
                .all(|p| [p.x, p.y, p.z].into_iter().all(f64::is_finite))
        {
            return Err(InvalidBoundary);
        }
        let vertices = self
            .native
            .iter()
            .copied()
            .map(builder::vertex)
            .collect::<Vec<_>>();
        let edges = vertices
            .iter()
            .zip(vertices.iter().cycle().skip(1))
            .zip(self.curves.iter().cloned())
            .map(|((a, b), curve)| Edge::new(a, b, curve))
            .collect::<Result<Vec<_>, _>>()
            .map_err(|_| InvalidBoundary)?;
        approximate_boundary_facet(&self.source, &edges, self.budget)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        BsplineCurve, EuclideanSpace, KnotVector, NurbsCurve, Vector4,
        mesh_recovery::{
            BoundaryPatchError, BoundaryPatchOptions, PrescribedSourceEdge, SurfaceRegionBudget,
            approximate_boundary_patch,
        },
    };

    #[test]
    fn serialized_input_replays_curved_fan_success_and_failure()
    -> Result<(), Box<dyn std::error::Error>> {
        let source = vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(2.0, 0.0, 0.0),
            Point3::new(0.0, 2.0, 0.0),
        ];
        let vertices = source
            .iter()
            .copied()
            .map(builder::vertex)
            .collect::<Vec<_>>();
        let curve = Curve::NurbsCurve(NurbsCurve::new(BsplineCurve::new(
            KnotVector::bezier_knot(2),
            vec![
                source[0].to_vec().extend(1.0),
                Vector4::new(0.2, 0.5, 0.0, 1.0),
                source[1].to_vec().extend(1.0),
            ],
        )?));
        let edges = [
            Edge::new(&vertices[0], &vertices[1], curve)?,
            builder::line(&vertices[1], &vertices[2]),
            builder::line(&vertices[2], &vertices[0]),
        ];
        for (budget, expected) in [(0.3, true), (0.1, false)] {
            let direct = approximate_boundary_facet(&source, &edges, budget);
            assert_eq!(direct.is_ok(), expected);
            let input = BoundaryFacetInput::capture(source.clone(), &edges, budget);
            let decoded: BoundaryFacetInput = serde_json::from_slice(&serde_json::to_vec(&input)?)?;
            let replay = decoded.approximate();
            match (direct, replay) {
                (Ok(a), Ok(b)) => {
                    assert_eq!(a.faces.len(), b.faces.len());
                    assert_eq!(
                        a.deviation_upper_bound.to_bits(),
                        b.deviation_upper_bound.to_bits()
                    );
                    assert_eq!(
                        a.seam_deviation_upper_bound.to_bits(),
                        b.seam_deviation_upper_bound.to_bits()
                    );
                }
                (Err(a), Err(b)) => assert_eq!(a, b),
                _ => return Err("replay changed certificate outcome".into()),
            }
        }
        let prescribed = [PrescribedSourceEdge {
            vertices: [0, 1],
            edge: edges[0].clone(),
        }];
        let polygons: [&[usize]; 1] = [&[0, 1, 2]];
        let failure = approximate_boundary_patch(
            &source,
            &polygons,
            &prescribed,
            BoundaryPatchOptions::new(SurfaceRegionBudget {
                vertex: 0.0,
                chord: 0.1,
            }),
        )
        .err()
        .ok_or("expected a failed patch")?;
        let BoundaryPatchError::Facet {
            polygon,
            error,
            source_vertex_ids,
            input,
        } = failure
        else {
            return Err("expected a replayable facet failure".into());
        };
        assert_eq!(polygon, 0);
        assert_eq!(source_vertex_ids, [0, 1, 2]);
        assert!(matches!(input.approximate(), Err(replayed) if replayed == error));
        Ok(())
    }
}

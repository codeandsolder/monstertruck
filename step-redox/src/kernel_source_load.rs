use crate::SourceEntity;
use anyhow::Result;
use serde::Serialize;

/// Parser-neutral source-load diagnostics exposed by a CAD-kernel adapter.
///
/// This boundary deliberately owns its vocabulary in Step Redox. Backend-native
/// report types must be translated before they cross it, so Complete IR does not
/// acquire a dependency on any particular kernel's STEP loader API.
pub trait SourceStepLoadInspector {
    /// Inspect how the kernel backend would ingest the original STEP source.
    ///
    /// # Errors
    /// Returns an error when the previously parsed source cannot be resolved or when the
    /// backend's diagnostic load path itself fails unexpectedly.
    /// One parsed exchange and its kernel-owned typed entity tables.
    type Document: AsRef<ruststep::ast::Exchange>;

    /// Build the typed resolution graph once per source document.
    ///
    /// # Errors
    /// Returns an error if the backend cannot initialize the source document.
    fn resolve_source_exchange(&self, exchange: ruststep::ast::Exchange) -> Result<Self::Document>;

    /// Inspect the already resolved source document.
    ///
    /// # Errors
    /// Returns an error if the backend cannot import or describe a source solid.
    fn inspect_source_document(
        &self,
        document: &Self::Document,
    ) -> Result<KernelSourceLoadEvidence>;
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize)]
pub struct KernelSourceLoadEvidence {
    pub sections: Vec<KernelSourceSectionLoadEvidence>,
    pub solids: Vec<KernelSourceSolidLoadEvidence>,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize)]
pub struct KernelSourceSectionLoadEvidence {
    pub data_section: usize,
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub entity_losses: Vec<KernelEntityLoadLoss>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum KernelEntityLoadLossKind {
    DeserializeRefused,
    UnexpectedParameterShape,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize)]
pub struct KernelEntityLoadLoss {
    pub kind: KernelEntityLoadLossKind,
    pub entity_type: String,
    pub tally: KernelLoadTally,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize)]
pub struct KernelSourceSolidLoadEvidence {
    pub source: SourceEntity,
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub boundaries: Vec<KernelBoundaryLoadEvidence>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub topology_provenance: Option<KernelTopologyProvenance>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub failure: Option<KernelSourceSolidLoadFailure>,
}

/// Opaque kernel topology identity scoped to one imported solid.
///
/// This is deliberately a distinct namespace from STEP entity IDs. Values are
/// only comparable within the [`KernelSourceSolidLoadEvidence`] that owns them.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize)]
#[serde(transparent)]
pub struct KernelTopologyId(u64);

impl KernelTopologyId {
    #[cfg(any(test, feature = "cad-kernel-monstertruck"))]
    pub(crate) const fn from_raw(raw: u64) -> Self { Self(raw) }

    #[must_use]
    pub const fn raw(self) -> u64 { self.0 }
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize)]
pub struct KernelTopologyEntityProvenance {
    pub source: SourceEntity,
    pub topology_id: KernelTopologyId,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize)]
pub struct KernelFaceTopologyProvenance {
    #[serde(skip_serializing_if = "Option::is_none")]
    pub source: Option<SourceEntity>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub shell_member: Option<SourceEntity>,
    pub topology_id: KernelTopologyId,
}

/// Source-to-kernel identity journal for topology retained during STEP import.
#[derive(Debug, Clone, PartialEq, Eq, Serialize)]
pub struct KernelTopologyProvenance {
    pub vertices: Vec<KernelTopologyEntityProvenance>,
    pub edges: Vec<KernelTopologyEntityProvenance>,
    pub faces: Vec<KernelFaceTopologyProvenance>,
}

/// Distinguish absence from conflicting retained topology identities.
///
/// These states are not interchangeable: an ambiguous shell-member mapping must
/// never trigger fallback to the underlying face's source identity.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum UniqueTopologyLookup {
    Missing,
    Unique(KernelTopologyId),
    Ambiguous,
}

impl UniqueTopologyLookup {
    const fn into_unique(self) -> Option<KernelTopologyId> {
        match self {
            Self::Unique(id) => Some(id),
            Self::Missing | Self::Ambiguous => None,
        }
    }
}

impl KernelTopologyProvenance {
    fn unique_id(ids: impl Iterator<Item = KernelTopologyId>) -> UniqueTopologyLookup {
        let mut first = None;
        for id in ids {
            match first {
                None => first = Some(id),
                Some(existing) if existing == id => {}
                Some(_) => return UniqueTopologyLookup::Ambiguous,
            }
        }
        first.map_or(UniqueTopologyLookup::Missing, UniqueTopologyLookup::Unique)
    }

    pub(crate) fn vertex_for_source_id(
        &self,
        data_section: usize,
        source_id: u64,
    ) -> Option<KernelTopologyId> {
        Self::unique_id(
            self.vertices
                .iter()
                .filter(|entry| {
                    entry.source.data_section == data_section && entry.source.id == source_id
                })
                .map(|entry| entry.topology_id),
        )
        .into_unique()
    }

    pub(crate) fn edge_for_source_id(
        &self,
        data_section: usize,
        source_id: u64,
    ) -> Option<KernelTopologyId> {
        Self::unique_id(
            self.edges
                .iter()
                .filter(|entry| {
                    entry.source.data_section == data_section && entry.source.id == source_id
                })
                .map(|entry| entry.topology_id),
        )
        .into_unique()
    }

    pub(crate) fn face_for_source_id(
        &self,
        data_section: usize,
        source_id: u64,
    ) -> Option<KernelTopologyId> {
        let shell_member = Self::unique_id(
            self.faces
                .iter()
                .filter(|entry| {
                    entry.shell_member.is_some_and(|source| {
                        source.data_section == data_section && source.id == source_id
                    })
                })
                .map(|entry| entry.topology_id),
        );
        match shell_member {
            UniqueTopologyLookup::Unique(id) => Some(id),
            UniqueTopologyLookup::Ambiguous => None,
            UniqueTopologyLookup::Missing => Self::unique_id(
                self.faces
                    .iter()
                    .filter(|entry| {
                        entry.source.is_some_and(|source| {
                            source.data_section == data_section && source.id == source_id
                        })
                    })
                    .map(|entry| entry.topology_id),
            )
            .into_unique(),
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize)]
#[serde(tag = "kind", rename_all = "snake_case")]
pub enum KernelSourceSolidLoadFailure {
    TableEntryUnavailable,
    ConversionFailed { detail: String },
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum KernelBoundaryRole {
    Outer,
    Void,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize)]
pub struct KernelBoundaryLoadEvidence {
    pub boundary_index: usize,
    pub role: KernelBoundaryRole,
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub category_counts: Vec<KernelLoadCategoryCount>,
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub losses: Vec<KernelLoadLoss>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum KernelLoadCategory {
    Face,
    Wire,
    EdgeUse,
    Edge,
    Vertex,
    Representation,
    Solid,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize)]
pub struct KernelLoadCategoryCount {
    pub category: KernelLoadCategory,
    pub listed: usize,
    pub kept: usize,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum KernelLoadLossReason {
    FaceUnresolved,
    SurfaceUnresolved,
    SurfaceRefused,
    BoundUnresolved,
    DegenerateVertexLoop,
    BoundNotALoopWeImplement,
    EdgeUseUnresolved,
    EdgeUnresolved,
    EdgeGeometryMissing,
    EdgeCurveRefused,
    EdgeEndpointUnresolved,
    VertexUnresolved,
    RelatedRepresentationUnresolved,
    RepresentationItemNotASolid,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize)]
pub struct KernelLoadLoss {
    pub category: KernelLoadCategory,
    pub reason: KernelLoadLossReason,
    pub tally: KernelLoadTally,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize)]
pub struct KernelLoadTally {
    pub count: usize,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub first_source: Option<SourceEntity>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub first_detail: Option<String>,
}

#[cfg(feature = "cad-kernel-monstertruck")]
mod monstertruck_adapter {
    use super::{
        KernelBoundaryLoadEvidence, KernelBoundaryRole, KernelEntityLoadLoss,
        KernelEntityLoadLossKind, KernelFaceTopologyProvenance, KernelLoadCategory,
        KernelLoadCategoryCount, KernelLoadLoss, KernelLoadLossReason, KernelLoadTally,
        KernelSourceLoadEvidence, KernelSourceSectionLoadEvidence, KernelSourceSolidLoadEvidence,
        KernelSourceSolidLoadFailure, KernelTopologyEntityProvenance, KernelTopologyId,
        KernelTopologyProvenance, SourceEntity, SourceStepLoadInspector,
    };
    use crate::brep::{
        self, CircleSupport, ConeSupport, CurveSupport, CylinderSupport, FaceLoop, LineSupport,
        OrientedEdgeUse, PlaneSupport, SphereSupport, SplineExtrusionSupport, SurfaceSupport,
        TorusSupport,
    };
    use crate::cad_kernel::monstertruck::MonstertruckKernel;
    use crate::kernel_topology::{
        KernelEditedSourceEvidence, KernelFaceBoundaryReplacement, KernelFaceSplice,
        KernelPatchCurve, KernelPatchEdgeRef, KernelPatchVertexRef, KernelTopologyPatch,
        SourceTopologyEditor,
    };
    use crate::reduced_brep::{ReducedFace, ReducedSolid};
    use crate::step_document::StepSection;
    use crate::step_graph;
    use anyhow::Result;
    use monstertruck_geometry::analytic_carrier::{AnalyticCarrier, AnalyticSupport};
    use monstertruck_geometry::analytic_line_segment::{
        SegmentLocus, certify_bspline_segment_locus, certify_nurbs_segment_locus,
    };
    #[cfg(test)]
    use monstertruck_io::step::load::ParsedStepSource;
    use monstertruck_io::step::load::convert::StepCompressedTrimmedSolid;
    use monstertruck_io::step::load::step_geometry::{
        BsplineSurface, Conic3D, Curve3D, DEFAULT_NUMERICAL_TOLERANCES, ElementarySurface,
        ExactParameterBoundary2D, KERNEL_EPSILON, Line, Matrix4, Point3, Processor, Surface,
        SweepSurface, Transform, TrimmedCurve, TryIntoHomogeneousBsplineSurface, UnitCircle,
        Vector3, Vector4,
    };
    use monstertruck_io::step::load::{
        EntityLoadReport, LossCategory, LossReason, LossTally, ResolvedStepDocument,
        ShellLoadReport, StepSolidImport, StepTopologyProvenance,
    };
    use monstertruck_io::step::save::{self, CompleteStepDisplay};
    use monstertruck_modeling::InnerSpace;
    use monstertruck_topology::StableId;
    use monstertruck_topology::compress::StableIdEntity;
    use monstertruck_topology::compressed_surgery::{
        CompressedFaceBoundaryReplacement, CompressedFaceSplice, CompressedTrimmedTopologyPatch,
        NewCompressedEdge, PatchEdgeRef, PatchEdgeUse, PatchVertexRef, RetainedCompressedEdge,
    };
    use monstertruck_topology::repair::unify_same_domain_compressed_trimmed_faces;

    impl SourceStepLoadInspector for MonstertruckKernel {
        type Document = ResolvedStepDocument;

        fn resolve_source_exchange(
            &self,
            exchange: ruststep::ast::Exchange,
        ) -> Result<Self::Document> {
            Ok(ResolvedStepDocument::from_exchange(exchange))
        }

        fn inspect_source_document(
            &self,
            document: &Self::Document,
        ) -> Result<KernelSourceLoadEvidence> {
            let exchange = document.exchange();
            let mut evidence = KernelSourceLoadEvidence {
                sections: Vec::with_capacity(exchange.data.len()),
                solids: Vec::new(),
            };

            for (data_section, section) in exchange.data.iter().enumerate() {
                let index = step_graph::build_index(&section.entities);
                let source = StepSection::new(data_section, &section.entities, &index);
                let table = document.table(data_section).ok_or_else(|| {
                    anyhow::anyhow!("typed source document omitted DATA section {data_section}")
                })?;

                evidence.sections.push(KernelSourceSectionLoadEvidence {
                    data_section,
                    entity_losses: entity_losses(&table.entity_report, source),
                });

                for entity in source
                    .iter()
                    .filter(|entity| entity.record_named("MANIFOLD_SOLID_BREP").is_some())
                {
                    let source_entity = entity.source();
                    if !table.manifold_solid_brep.contains_key(&source_entity.id) {
                        evidence.solids.push(KernelSourceSolidLoadEvidence {
                            source: source_entity,
                            boundaries: Vec::new(),
                            topology_provenance: None,
                            failure: Some(KernelSourceSolidLoadFailure::TableEntryUnavailable),
                        });
                        continue;
                    }

                    match document.import_solid(data_section, source_entity.id) {
                        Ok(import) => evidence.solids.push(KernelSourceSolidLoadEvidence {
                            source: source_entity,
                            boundaries: import
                                .shell_reports
                                .iter()
                                .enumerate()
                                .map(|(boundary_index, report)| KernelBoundaryLoadEvidence {
                                    boundary_index,
                                    role: if boundary_index == 0 {
                                        KernelBoundaryRole::Outer
                                    } else {
                                        KernelBoundaryRole::Void
                                    },
                                    category_counts: category_counts(report),
                                    losses: shell_losses(report, source),
                                })
                                .collect(),
                            topology_provenance: Some(map_topology_provenance(
                                &import.provenance,
                                source,
                            )?),
                            failure: None,
                        }),
                        Err(error) => evidence.solids.push(KernelSourceSolidLoadEvidence {
                            source: source_entity,
                            boundaries: Vec::new(),
                            topology_provenance: None,
                            failure: Some(KernelSourceSolidLoadFailure::ConversionFailed {
                                detail: error.to_string(),
                            }),
                        }),
                    }
                }
            }

            evidence.solids.sort_by_key(|solid| {
                (
                    solid.source.data_section,
                    solid.source.ordinal,
                    solid.source.id,
                )
            });
            Ok(evidence)
        }
    }

    fn map_topology_provenance(
        provenance: &StepTopologyProvenance,
        source: StepSection<'_>,
    ) -> Result<KernelTopologyProvenance> {
        let map_entity = |entry: &monstertruck_io::step::load::StepEntityProvenance| -> Result<KernelTopologyEntityProvenance> {
            let source_id = entry.source_id;
            let stable_id = entry.stable_id;
            let source = source
                .entity(source_id)
                .ok_or_else(|| {
                    anyhow::anyhow!("STEP provenance references missing entity #{source_id}")
                })?
                .source();
            Ok(KernelTopologyEntityProvenance {
                source,
                topology_id: KernelTopologyId::from_raw(stable_id.raw()),
            })
        };

        let vertices = provenance
            .vertices
            .iter()
            .map(&map_entity)
            .collect::<Result<Vec<_>>>()?;
        let edges = provenance
            .edges
            .iter()
            .map(&map_entity)
            .collect::<Result<Vec<_>>>()?;
        let faces = provenance
            .faces
            .iter()
            .map(|entry| {
                let map_optional = |id: Option<u64>| -> Result<Option<SourceEntity>> {
                    id.map(|source_id| {
                        source
                            .entity(source_id)
                            .ok_or_else(|| {
                                anyhow::anyhow!(
                                    "STEP face provenance references missing entity #{source_id}"
                                )
                            })
                            .map(crate::step_document::StepEntity::source)
                    })
                    .transpose()
                };
                Ok(KernelFaceTopologyProvenance {
                    source: map_optional(entry.source.face_surface_id)?,
                    shell_member: map_optional(entry.source.shell_member_id)?,
                    topology_id: KernelTopologyId::from_raw(entry.stable_id.raw()),
                })
            })
            .collect::<Result<Vec<_>>>()?;

        Ok(KernelTopologyProvenance {
            vertices,
            edges,
            faces,
        })
    }

    /// A deterministic radial frame for axisymmetric support geometry.
    fn native_radial_direction(axis: Vector3) -> [f64; 3] {
        let candidate = if axis.x.abs() < 0.8 {
            Vector3::unit_x()
        } else {
            Vector3::unit_y()
        };
        let radial = candidate
            - axis
                * candidate
                    .x
                    .mul_add(axis.x, candidate.y.mul_add(axis.y, candidate.z * axis.z));
        let length = radial
            .x
            .mul_add(radial.x, radial.y.mul_add(radial.y, radial.z * radial.z))
            .sqrt();
        [radial.x / length, radial.y / length, radial.z / length]
    }

    /// Classify the full geometric carrier independently of trim and STEP
    /// entity spelling. Line-generated sweep surfaces use the same analytic
    /// identities as elementary cylinders/cones.
    fn native_surface_certificate(surface: &Surface) -> Option<AnalyticCarrier> {
        match surface {
            Surface::Elementary(ElementarySurface::Plane(plane)) => AnalyticCarrier::plane(plane),
            Surface::Elementary(ElementarySurface::Sphere(sphere)) => AnalyticCarrier::sphere_pose(
                sphere.entity().0.center(),
                sphere.entity().0.radius(),
                *sphere.transform(),
                sphere.orientation(),
            ),
            Surface::Elementary(ElementarySurface::CylindricalSurface(cylinder)) => {
                AnalyticCarrier::linear_revolution(
                    *cylinder.entity().entity_curve(),
                    cylinder.entity().origin(),
                    cylinder.entity().axis(),
                    *cylinder.transform(),
                    cylinder.orientation(),
                )
            }
            Surface::Elementary(ElementarySurface::ConicalSurface(cone)) => {
                AnalyticCarrier::linear_revolution(
                    *cone.entity().entity_curve(),
                    cone.entity().origin(),
                    cone.entity().axis(),
                    *cone.transform(),
                    cone.orientation(),
                )
            }
            Surface::Elementary(ElementarySurface::ToroidalSurface(torus)) => {
                AnalyticCarrier::torus(torus)
            }
            Surface::Sweep(SweepSurface::RevolutionSurface(revolution)) => {
                match revolution.entity().entity_curve() {
                    Curve3D::Line(line) => AnalyticCarrier::linear_revolution(
                        *line,
                        revolution.entity().origin(),
                        revolution.entity().axis(),
                        *revolution.transform(),
                        revolution.orientation(),
                    ),
                    _ => None,
                }
            }
            _ => None,
        }
    }

    /// A transformed *planar* unit circle needs equal, orthogonal X/Y
    /// columns. The perpendicular Z column is irrelevant to its 3D locus:
    /// actual STEP placements commonly scale X/Y by the radius but leave Z
    /// unit length, so a 3D-similarity requirement would discard them.
    fn native_circle_support(
        circle: &Processor<TrimmedCurve<UnitCircle<Point3>>, Matrix4>,
    ) -> Option<CircleSupport> {
        let transform = *circle.transform();
        if ![
            transform.x.x,
            transform.x.y,
            transform.x.z,
            transform.x.w,
            transform.y.x,
            transform.y.y,
            transform.y.z,
            transform.y.w,
            transform.z.x,
            transform.z.y,
            transform.z.z,
            transform.z.w,
            transform.w.x,
            transform.w.y,
            transform.w.z,
            transform.w.w,
        ]
        .into_iter()
        .all(f64::is_finite)
            || transform.x.w.abs() > 1.0e-12
            || transform.y.w.abs() > 1.0e-12
            || (transform.w.w - 1.0).abs() > 1.0e-12
        {
            return None;
        }
        let x = transform.x.truncate();
        let y = transform.y.truncate();
        let radius_mm = x.magnitude();
        let tolerance = 1.0e-10 * radius_mm.max(1.0);
        if radius_mm <= KERNEL_EPSILON
            || (y.magnitude() - radius_mm).abs() > tolerance
            || x.dot(y).abs() > tolerance * radius_mm
        {
            return None;
        }
        let normal = x.cross(y).normalize();
        let center = transform.transform_point(Point3::new(0.0, 0.0, 0.0));
        let radial = x / radius_mm;
        Some(CircleSupport {
            center_mm: [center.x, center.y, center.z],
            normal: [normal.x, normal.y, normal.z],
            x_direction: [radial.x, radial.y, radial.z],
            radius_mm,
        })
    }

    /// Homogeneous spline control weights are defined only up to a shared
    /// positive scale. Normalize edges and faces identically before comparing
    /// their exact source profiles.
    fn native_spline_weights(raw: &[f64]) -> Option<Vec<f64>> {
        let scale = *raw.first()?;
        if !scale.is_finite() || scale <= 0.0 {
            return None;
        }
        let normalized = raw.iter().map(|weight| weight / scale).collect::<Vec<_>>();
        if normalized
            .iter()
            .any(|weight| !weight.is_finite() || *weight <= 0.0)
        {
            return None;
        }
        Some(normalized)
    }

    /// Extract an actual Euclidean control point, not a homogeneous pole.
    fn native_net_point(point: Vector4) -> Option<Point3> {
        if ![point.x, point.y, point.z, point.w]
            .into_iter()
            .all(f64::is_finite)
            || point.w <= 0.0
        {
            return None;
        }
        let euclidean = point.truncate() / point.w;
        Some(Point3::new(euclidean.x, euclidean.y, euclidean.z))
    }

    /// A positive-weight B-spline lies exactly in a plane when every one of
    /// its Euclidean control points lies in that plane. Preserve the U x V
    /// normal rather than canonicalizing away the native chart orientation.
    fn native_planar_net(net: &BsplineSurface<Vector4>) -> Option<PlaneSupport> {
        let rows = net.control_points();
        let first = rows.first()?;
        let origin = native_net_point(*first.first()?)?;
        let u = rows.iter().skip(1).find_map(|row| {
            let delta = native_net_point(*row.first()?)? - origin;
            (delta.magnitude() > KERNEL_EPSILON).then_some(delta)
        })?;
        let v = first.iter().skip(1).find_map(|pole| {
            let delta = native_net_point(*pole)? - origin;
            (delta.magnitude() > KERNEL_EPSILON).then_some(delta)
        })?;
        let cross = u.cross(v);
        let length = cross.magnitude();
        if !length.is_finite() || length <= KERNEL_EPSILON {
            return None;
        }
        let normal = cross / length;
        let mut max_residual_mm = 0.0_f64;
        for pole in rows.iter().flatten() {
            let point = native_net_point(*pole)?;
            let residual = (point - origin).dot(normal).abs();
            if !residual.is_finite() || residual > 1.0e-7 {
                return None;
            }
            max_residual_mm = max_residual_mm.max(residual);
        }
        Some(PlaneSupport {
            origin_mm: [origin.x, origin.y, origin.z],
            normal: [normal.x, normal.y, normal.z],
            max_residual_mm,
        })
    }

    /// Recognize exact ruled spline surfaces independently of which parameter
    /// runs along the profile. The native B-rep lets us transpose this *view*
    /// without modifying STEP parameter curves or their trims.
    fn native_spline_extrusion(net: &BsplineSurface<Vector4>) -> Option<SplineExtrusionSupport> {
        let control = net.control_points();
        let (profile_degree, profile_knots, sweep_knots, rows) = if net.vdegree() == 1
            && net.udegree() > 0
            && control.len() >= 2
            && control.iter().all(|row| row.len() == 2)
        {
            (
                net.udegree(),
                net.knot_vector_u(),
                net.knot_vector_v(),
                control
                    .iter()
                    .map(|row| [row[0], row[1]])
                    .collect::<Vec<_>>(),
            )
        } else if net.udegree() == 1
            && net.vdegree() > 0
            && control.len() == 2
            && control[0].len() >= 2
            && control[0].len() == control[1].len()
        {
            (
                net.vdegree(),
                net.knot_vector_v(),
                net.knot_vector_u(),
                (0..control[0].len())
                    .map(|index| [control[0][index], control[1][index]])
                    .collect::<Vec<_>>(),
            )
        } else {
            return None;
        };
        if sweep_knots.len() != 4
            || !sweep_knots.iter().all(|knot| knot.is_finite())
            || (sweep_knots[0] - sweep_knots[1]).abs() > 1.0e-12
            || (sweep_knots[2] - sweep_knots[3]).abs() > 1.0e-12
            || sweep_knots[1] >= sweep_knots[2]
            || profile_knots.iter().any(|knot| !knot.is_finite())
        {
            return None;
        }
        let mut profile_points = Vec::with_capacity(rows.len());
        let mut profile_weights = Vec::with_capacity(rows.len());
        let mut sweep = None::<Vector3>;
        let mut max_residual_mm = 0.0_f64;
        for [first, last] in rows {
            if ![first.w, last.w].into_iter().all(f64::is_finite)
                || first.w <= 0.0
                || last.w <= 0.0
                || (first.w - last.w).abs() > 1.0e-12 * first.w.abs().max(last.w.abs())
            {
                return None;
            }
            let p0 = native_net_point(first)?;
            let p1 = native_net_point(last)?;
            let delta = p1 - p0;
            let length = delta.magnitude();
            if !length.is_finite() || length <= 1.0e-7 {
                return None;
            }
            if let Some(expected) = sweep {
                let residual = (delta - expected).magnitude();
                if residual > 1.0e-7 {
                    return None;
                }
                max_residual_mm = max_residual_mm.max(residual);
            } else {
                sweep = Some(delta);
            }
            profile_points.push([p0.x, p0.y, p0.z]);
            profile_weights.push(first.w);
        }
        let normalized_weights = native_spline_weights(&profile_weights)?;
        let weights = normalized_weights
            .iter()
            .any(|weight| (weight - 1.0).abs() > 1.0e-12)
            .then_some(normalized_weights);
        let delta = sweep?;
        Some(SplineExtrusionSupport {
            profile: brep::BSplineSupport {
                degree: profile_degree,
                control_points_mm: profile_points,
                knots: profile_knots.iter().copied().collect(),
                weights,
            },
            extrusion_mm: [delta.x, delta.y, delta.z],
            max_residual_mm,
        })
    }

    /// Classify a true native sweep without requiring STEP's surface name:
    /// line-generated sweeps are planar and circular sweeps are cylindrical.
    fn native_exact_sweep_surface(surface: &Surface) -> Option<SurfaceSupport> {
        let Surface::Sweep(SweepSurface::ExtrusionSurface(sweep)) = surface else {
            return None;
        };
        if let Curve3D::Line(Line(start, end)) = sweep.entity_curve() {
            let normal = (*end - *start).cross(sweep.extruding_vector());
            let magnitude = normal.magnitude();
            if magnitude.is_finite() && magnitude > KERNEL_EPSILON {
                let normal = normal / magnitude;
                return Some(SurfaceSupport::Plane(PlaneSupport {
                    origin_mm: [start.x, start.y, start.z],
                    normal: [normal.x, normal.y, normal.z],
                    max_residual_mm: 0.0,
                }));
            }
        }
        if let Curve3D::Conic(Conic3D::Ellipse(circle)) = sweep.entity_curve()
            && let Some(circle) = native_circle_support(circle)
        {
            let vector = sweep.extruding_vector();
            let norm = vector.x.hypot(vector.y).hypot(vector.z);
            let normal = Vector3::new(circle.normal[0], circle.normal[1], circle.normal[2]);
            if norm.is_finite()
                && norm > KERNEL_EPSILON
                && (normal.dot(vector / norm).abs() - 1.0).abs() <= 1.0e-10
            {
                let axis = vector / norm;
                return Some(SurfaceSupport::Cylinder(CylinderSupport {
                    axis_origin_mm: circle.center_mm,
                    axis: [axis.x, axis.y, axis.z],
                    x_direction: circle.x_direction,
                    radius_mm: circle.radius_mm,
                }));
            }
        }
        None
    }

    /// Convert full checked analytic carriers from the edited native B-rep.
    /// Unknown geometry stays unknown, with no fitted or sampled replacement.
    fn native_surface_support(surface: &Surface, id: u64) -> SurfaceSupport {
        if let Some(support) = native_exact_sweep_surface(surface) {
            return support;
        }
        let certificate = native_surface_certificate(surface);
        match certificate.map(AnalyticCarrier::support) {
            Some(AnalyticSupport::Plane { origin, .. }) => {
                let Surface::Elementary(ElementarySurface::Plane(plane)) = surface else {
                    return SurfaceSupport::Other { entity_id: id };
                };
                // A carrier certificate canonicalizes an unoriented plane.
                // Preserve the native chart normal for face-sense recovery.
                let normal = plane.normal();
                SurfaceSupport::Plane(PlaneSupport {
                    origin_mm: [origin.x, origin.y, origin.z],
                    normal: [normal.x, normal.y, normal.z],
                    max_residual_mm: 0.0,
                })
            }
            Some(AnalyticSupport::Cylinder {
                origin,
                axis,
                radius,
            }) => SurfaceSupport::Cylinder(CylinderSupport {
                axis_origin_mm: [origin.x, origin.y, origin.z],
                axis: [axis.x, axis.y, axis.z],
                x_direction: native_radial_direction(axis),
                radius_mm: radius,
            }),
            Some(AnalyticSupport::Cone {
                apex,
                axis,
                slope,
                positive_nappe,
            }) => {
                let directed = axis * if positive_nappe { 1.0 } else { -1.0 };
                let reference = apex + directed;
                SurfaceSupport::Cone(ConeSupport {
                    reference_origin_mm: [reference.x, reference.y, reference.z],
                    axis: [directed.x, directed.y, directed.z],
                    x_direction: native_radial_direction(directed),
                    reference_radius_mm: slope,
                    semi_angle_rad: slope.atan(),
                })
            }
            Some(AnalyticSupport::Sphere { center, radius }) => {
                let axis = Vector3::unit_z();
                SurfaceSupport::Sphere(SphereSupport {
                    center_mm: [center.x, center.y, center.z],
                    axis: [axis.x, axis.y, axis.z],
                    x_direction: native_radial_direction(axis),
                    radius_mm: radius,
                })
            }
            Some(AnalyticSupport::Torus {
                center,
                axis,
                major_radius,
                minor_radius,
            }) => SurfaceSupport::Torus(TorusSupport {
                center_mm: [center.x, center.y, center.z],
                axis: [axis.x, axis.y, axis.z],
                x_direction: native_radial_direction(axis),
                major_radius_mm: major_radius,
                minor_radius_mm: minor_radius,
            }),
            None => {
                let Some(net) = surface.try_into_homogeneous_bspline_surface() else {
                    return SurfaceSupport::Other { entity_id: id };
                };
                if let Some(plane) = native_planar_net(&net) {
                    return SurfaceSupport::Plane(plane);
                }
                native_spline_extrusion(&net).map_or(
                    SurfaceSupport::Other { entity_id: id },
                    SurfaceSupport::SplineExtrusion,
                )
            }
        }
    }

    /// Extract analytic curve support directly from the kernel's exact 3D
    /// carrier. Polynomial and rational splines are promoted to a line only
    /// when a whole-curve control-hull certificate covers the complete chord;
    /// unknown geometry is never replaced by sampled polyline segments.
    fn native_curve_support(curve: &Curve3D, start: Point3, end: Point3, id: u64) -> CurveSupport {
        let line = || {
            let vector = end - start;
            CurveSupport::Line(LineSupport {
                origin_mm: [start.x, start.y, start.z],
                direction: [vector.x, vector.y, vector.z],
            })
        };
        let is_line = match curve {
            Curve3D::Line(_) => true,
            Curve3D::BsplineCurve(spline) => {
                certify_bspline_segment_locus(spline, start, end, 1.0e-7)
                    == SegmentLocus::WithinTolerance
            }
            Curve3D::NurbsCurve(spline) => {
                certify_nurbs_segment_locus(spline, start, end, 1.0e-7)
                    == SegmentLocus::WithinTolerance
            }
            Curve3D::SurfaceCurve(surface) => {
                return native_curve_support(surface.leader(), start, end, id);
            }
            Curve3D::IntersectionCurve(intersection) => {
                return native_curve_support(intersection.leader(), start, end, id);
            }
            _ => false,
        };
        if is_line {
            return line();
        }
        match curve {
            Curve3D::Conic(Conic3D::Ellipse(circle)) => native_circle_support(circle)
                .map_or(CurveSupport::Other { entity_id: id }, CurveSupport::Circle),
            Curve3D::BsplineCurve(spline) => CurveSupport::BSpline(brep::BSplineSupport {
                degree: spline.degree(),
                control_points_mm: spline
                    .control_points()
                    .iter()
                    .map(|p| [p.x, p.y, p.z])
                    .collect(),
                knots: spline.knot_vector().iter().copied().collect(),
                weights: None,
            }),
            Curve3D::NurbsCurve(spline) => {
                let points = spline.control_points();
                if points.iter().any(|p| p.w <= 0.0 || !p.w.is_finite()) {
                    return CurveSupport::Other { entity_id: id };
                }
                CurveSupport::BSpline(brep::BSplineSupport {
                    degree: spline.degree(),
                    control_points_mm: points
                        .iter()
                        .map(|p| [p.x / p.w, p.y / p.w, p.z / p.w])
                        .collect(),
                    knots: spline.knot_vector().iter().copied().collect(),
                    weights: {
                        let Some(weights) =
                            native_spline_weights(&points.iter().map(|p| p.w).collect::<Vec<_>>())
                        else {
                            return CurveSupport::Other { entity_id: id };
                        };
                        weights
                            .iter()
                            .any(|weight| (weight - 1.0).abs() > 1.0e-12)
                            .then_some(weights)
                    },
                })
            }
            _ => CurveSupport::Other { entity_id: id },
        }
    }

    /// Exact-trim compressed solid and original STEP provenance after a
    /// transactional edit. Stable IDs and allocator state remain alive until a
    /// downstream consumer explicitly requests a STEP snapshot.
    pub struct MonstertruckEditedSource {
        source: SourceEntity,
        imported: StepSolidImport<StepCompressedTrimmedSolid>,
        source_topology: KernelTopologyProvenance,
    }

    impl SourceTopologyEditor for MonstertruckKernel {
        type EditedSource = MonstertruckEditedSource;

        fn edit_source_topology(
            &self,
            document: &Self::Document,
            source: SourceEntity,
            patches: &[KernelTopologyPatch],
        ) -> Result<Self::EditedSource> {
            let section = document
                .exchange()
                .data
                .get(source.data_section)
                .ok_or_else(|| {
                    anyhow::anyhow!(
                        "source solid references missing DATA section {}",
                        source.data_section
                    )
                })?;
            let index = step_graph::build_index(&section.entities);
            let section_source = StepSection::new(source.data_section, &section.entities, &index);
            let actual_source = section_source
                .entity(source.id)
                .ok_or_else(|| anyhow::anyhow!("source solid #{} is missing", source.id))?
                .source();
            if actual_source != source {
                anyhow::bail!(
                    "source solid identity changed: requested {source:?}, found {actual_source:?}"
                );
            }

            let table = document.table(source.data_section).ok_or_else(|| {
                anyhow::anyhow!(
                    "typed source document omitted DATA section {}",
                    source.data_section
                )
            })?;
            if !table.manifold_solid_brep.contains_key(&source.id) {
                anyhow::bail!("source solid #{} is unavailable to Monstertruck", source.id);
            }
            let imported = document
                .import_solid(source.data_section, source.id)?
                .clone();
            let source_topology = map_topology_provenance(&imported.provenance, section_source)?;
            let mut edited = MonstertruckEditedSource {
                source,
                imported,
                source_topology,
            };
            self.apply_edited_source_patches(&mut edited, patches)?;
            Ok(edited)
        }

        fn apply_edited_source_patches(
            &self,
            edited: &mut Self::EditedSource,
            patches: &[KernelTopologyPatch],
        ) -> Result<()> {
            let mut candidate = edited.imported.solid.clone();
            for patch in patches {
                candidate.apply_topology_patch(&translate_patch(patch))?;
            }
            edited.imported.solid = candidate;
            Ok(())
        }

        fn normalize_edited_source(&self, edited: &mut Self::EditedSource) -> Result<usize> {
            let mut candidate = edited.imported.solid.clone();
            let mut accepted = 0;
            for shell in &mut candidate.boundaries {
                let result = unify_same_domain_compressed_trimmed_faces(
                    shell,
                    |left, right| {
                        monstertruck_io::step::load::step_geometry::same_surface_domain(
                            left,
                            right,
                            KERNEL_EPSILON,
                            DEFAULT_NUMERICAL_TOLERANCES.same_domain_similarity,
                        )
                    },
                    ExactParameterBoundary2D::exact_parameter_boundary_2d,
                )?;
                accepted += result.accepted;
            }
            // All shells must have succeeded before modifying the native shape.
            edited.imported.solid = candidate;
            Ok(accepted)
        }

        fn edited_source_reduced_brep(&self, edited: &Self::EditedSource) -> Result<ReducedSolid> {
            let mut faces = Vec::new();
            // Dense synthetic identities are deliberately independent of the
            // original source journal. They identify a *current* exact shape.
            let mut next_id = 1_u64;
            let mut allocate = || {
                let id = next_id;
                next_id += 1;
                id
            };
            for shell in &edited.imported.solid.boundaries {
                let vertex_ids = shell
                    .vertices
                    .iter()
                    .map(|_| allocate())
                    .collect::<Vec<_>>();
                let edge_ids = shell.edges.iter().map(|_| allocate()).collect::<Vec<_>>();
                for face in &shell.faces {
                    let face_id = allocate();
                    let surface = native_surface_support(&face.surface, face_id);
                    let loops = face
                        .boundaries
                        .iter()
                        .enumerate()
                        .map(|(wire_index, wire)| {
                            let bound_id = allocate();
                            let loop_id = allocate();
                            let edges = wire
                                .iter()
                                .map(|use_| {
                                    let edge = &shell.edges[use_.index];
                                    let (front, back) = if use_.orientation {
                                        edge.vertices
                                    } else {
                                        (edge.vertices.1, edge.vertices.0)
                                    };
                                    let start = shell.vertices[front];
                                    let end = shell.vertices[back];
                                    let support = native_curve_support(
                                        &edge.curve,
                                        start,
                                        end,
                                        edge_ids[use_.index],
                                    );
                                    OrientedEdgeUse {
                                        oriented_edge_id: allocate(),
                                        edge_id: edge_ids[use_.index],
                                        curve_id: edge_ids[use_.index],
                                        curve_same_sense: true,
                                        parameter_forward: use_.orientation,
                                        start_vertex: vertex_ids[front],
                                        end_vertex: vertex_ids[back],
                                        start_mm: [start.x, start.y, start.z],
                                        end_mm: [end.x, end.y, end.z],
                                        support,
                                    }
                                })
                                .collect();
                            FaceLoop {
                                bound_id,
                                loop_id,
                                outer: wire_index == 0,
                                orientation: true,
                                edges,
                            }
                        })
                        .collect();
                    faces.push(ReducedFace {
                        id: face_id,
                        surface,
                        surface_forward: face.orientation,
                        loops,
                    });
                }
            }
            Ok(ReducedSolid::from_native_faces(edited.source.id, faces))
        }

        fn inspect_edited_source(
            &self,
            edited: &Self::EditedSource,
        ) -> Result<KernelEditedSourceEvidence> {
            let solid = &edited.imported.solid;
            // Build one native StableId index across all boundary shells.
            // Duplicated IDs are not unique topology identities and must not
            // be reported as unambiguously retained source entities.
            let mut native_ids = std::collections::HashMap::<u64, (StableIdEntity, usize)>::new();
            for boundary in &solid.boundaries {
                for (ids, entity) in [
                    (
                        boundary.vertex_stable_ids.as_deref(),
                        StableIdEntity::Vertex,
                    ),
                    (boundary.edge_stable_ids.as_deref(), StableIdEntity::Edge),
                    (boundary.face_stable_ids.as_deref(), StableIdEntity::Face),
                ] {
                    for id in ids.into_iter().flatten().filter(|id| id.is_assigned()) {
                        native_ids
                            .entry(id.raw())
                            .and_modify(|entry| entry.1 += 1)
                            .or_insert((entity, 1));
                    }
                }
            }
            let retained = |id: KernelTopologyId, class| {
                native_ids
                    .get(&id.raw())
                    .is_some_and(|&(entity, count)| entity == class && count == 1)
            };
            let source_topology = &edited.source_topology;
            Ok(KernelEditedSourceEvidence {
                source: edited.source,
                boundary_shells: solid.boundaries.len(),
                faces: solid.boundaries.iter().map(|shell| shell.faces.len()).sum(),
                retained_source_topology: KernelTopologyProvenance {
                    vertices: source_topology
                        .vertices
                        .iter()
                        .filter(|entry| retained(entry.topology_id, StableIdEntity::Vertex))
                        .cloned()
                        .collect(),
                    edges: source_topology
                        .edges
                        .iter()
                        .filter(|entry| retained(entry.topology_id, StableIdEntity::Edge))
                        .cloned()
                        .collect(),
                    faces: source_topology
                        .faces
                        .iter()
                        .filter(|entry| retained(entry.topology_id, StableIdEntity::Face))
                        .cloned()
                        .collect(),
                },
            })
        }

        fn edited_source_to_step(&self, edited: &Self::EditedSource) -> Result<String> {
            Ok(CompleteStepDisplay::new(
                save::StepModel::from(&edited.imported.solid),
                save::StepHeaderDescriptor {
                    organization_system: "step-redox Monstertruck source-topology surgery".into(),
                    ..Default::default()
                },
            )
            .to_string())
        }
    }

    fn translate_patch(
        patch: &KernelTopologyPatch,
    ) -> CompressedTrimmedTopologyPatch<Point3, Curve3D> {
        CompressedTrimmedTopologyPatch {
            remove_faces: patch
                .remove_faces
                .iter()
                .map(|id| StableId::new(id.raw()))
                .collect(),
            new_vertices: patch.new_vertices.iter().copied().map(point3).collect(),
            new_edges: patch
                .new_edges
                .iter()
                .map(|edge| NewCompressedEdge {
                    vertices: (
                        translate_vertex_ref(edge.vertices.0),
                        translate_vertex_ref(edge.vertices.1),
                    ),
                    curve: translate_curve(&edge.curve),
                })
                .collect(),
            rewrite_edges: patch
                .rewrite_edges
                .iter()
                .map(|edge| RetainedCompressedEdge {
                    edge: StableId::new(edge.edge.raw()),
                    vertices: (
                        translate_vertex_ref(edge.vertices.0),
                        translate_vertex_ref(edge.vertices.1),
                    ),
                    curve: translate_curve(&edge.curve),
                })
                .collect(),
            face_splices: patch
                .face_splices
                .iter()
                .map(translate_face_splice)
                .collect(),
            face_boundary_replacements: patch
                .face_boundary_replacements
                .iter()
                .map(translate_face_boundary_replacement)
                .collect(),
            require_unused_edges: patch
                .require_unused_edges
                .iter()
                .map(|id| StableId::new(id.raw()))
                .collect(),
        }
    }

    const fn translate_vertex_ref(reference: KernelPatchVertexRef) -> PatchVertexRef {
        match reference {
            KernelPatchVertexRef::Existing(id) => PatchVertexRef::Existing(StableId::new(id.raw())),
            KernelPatchVertexRef::New(index) => PatchVertexRef::New(index),
        }
    }

    const fn translate_edge_ref(reference: KernelPatchEdgeRef) -> PatchEdgeRef {
        match reference {
            KernelPatchEdgeRef::Existing(id) => PatchEdgeRef::Existing(StableId::new(id.raw())),
            KernelPatchEdgeRef::New(index) => PatchEdgeRef::New(index),
        }
    }

    fn translate_face_boundary_replacement(
        replacement: &KernelFaceBoundaryReplacement,
    ) -> CompressedFaceBoundaryReplacement {
        CompressedFaceBoundaryReplacement {
            face: StableId::new(replacement.face.raw()),
            boundaries: replacement
                .boundaries
                .iter()
                .map(|wire| {
                    wire.iter()
                        .map(|edge_use| PatchEdgeUse {
                            edge: translate_edge_ref(edge_use.edge),
                            orientation: edge_use.orientation,
                        })
                        .collect()
                })
                .collect(),
        }
    }

    fn translate_face_splice(splice: &KernelFaceSplice) -> CompressedFaceSplice {
        CompressedFaceSplice {
            face: StableId::new(splice.face.raw()),
            remove_edges: splice
                .remove_edges
                .iter()
                .map(|id| StableId::new(id.raw()))
                .collect(),
            replacement: splice
                .replacement
                .iter()
                .map(|edge_use| PatchEdgeUse {
                    edge: translate_edge_ref(edge_use.edge),
                    orientation: edge_use.orientation,
                })
                .collect(),
        }
    }

    fn translate_curve(curve: &KernelPatchCurve) -> Curve3D {
        match *curve {
            KernelPatchCurve::Line { start_mm, end_mm } => {
                Curve3D::Line(Line(point3(start_mm), point3(end_mm)))
            }
            KernelPatchCurve::CircleArc {
                center_mm,
                normal,
                x_direction,
                radius_mm,
                start_angle_rad,
                end_angle_rad,
            } => {
                let y_direction = [
                    normal[1].mul_add(x_direction[2], -normal[2] * x_direction[1]),
                    normal[2].mul_add(x_direction[0], -normal[0] * x_direction[2]),
                    normal[0].mul_add(x_direction[1], -normal[1] * x_direction[0]),
                ];
                let transform = Matrix4::from_cols(
                    Vector4::new(
                        x_direction[0] * radius_mm,
                        x_direction[1] * radius_mm,
                        x_direction[2] * radius_mm,
                        0.0,
                    ),
                    Vector4::new(
                        y_direction[0] * radius_mm,
                        y_direction[1] * radius_mm,
                        y_direction[2] * radius_mm,
                        0.0,
                    ),
                    Vector4::new(normal[0], normal[1], normal[2], 0.0),
                    Vector4::new(center_mm[0], center_mm[1], center_mm[2], 1.0),
                );
                Curve3D::Conic(Conic3D::Ellipse(Processor::with_transform(
                    TrimmedCurve::new(
                        UnitCircle::<Point3>::new(),
                        (start_angle_rad, end_angle_rad),
                    ),
                    transform,
                )))
            }
        }
    }

    const fn point3(point: [f64; 3]) -> Point3 { Point3::new(point[0], point[1], point[2]) }

    fn entity_losses(
        report: &EntityLoadReport,
        source: StepSection<'_>,
    ) -> Vec<KernelEntityLoadLoss> {
        report
            .refused()
            .map(|(entity_type, tally)| KernelEntityLoadLoss {
                kind: KernelEntityLoadLossKind::DeserializeRefused,
                entity_type: entity_type.to_owned(),
                tally: map_tally(tally, source),
            })
            .chain(
                report
                    .shape_unexpected()
                    .map(|(entity_type, tally)| KernelEntityLoadLoss {
                        kind: KernelEntityLoadLossKind::UnexpectedParameterShape,
                        entity_type: entity_type.to_owned(),
                        tally: map_tally(tally, source),
                    }),
            )
            .collect()
    }

    fn category_counts(report: &ShellLoadReport) -> Vec<KernelLoadCategoryCount> {
        LossCategory::ALL
            .into_iter()
            .filter_map(|category| {
                let count = report.count(category);
                (count.listed != 0 || count.kept != 0).then_some(KernelLoadCategoryCount {
                    category: map_category(category),
                    listed: count.listed,
                    kept: count.kept,
                })
            })
            .collect()
    }

    fn shell_losses(report: &ShellLoadReport, source: StepSection<'_>) -> Vec<KernelLoadLoss> {
        report
            .losses()
            .map(|(category, reason, tally)| KernelLoadLoss {
                category: map_category(category),
                reason: map_reason(reason),
                tally: map_tally(tally, source),
            })
            .collect()
    }

    fn map_tally(tally: &LossTally, source: StepSection<'_>) -> KernelLoadTally {
        KernelLoadTally {
            count: tally.count,
            first_source: tally
                .first_id
                .and_then(|id| source.entity(id))
                .map(crate::step_document::StepEntity::source),
            first_detail: tally.first_detail.clone(),
        }
    }

    const fn map_category(category: LossCategory) -> KernelLoadCategory {
        match category {
            LossCategory::Face => KernelLoadCategory::Face,
            LossCategory::Wire => KernelLoadCategory::Wire,
            LossCategory::EdgeUse => KernelLoadCategory::EdgeUse,
            LossCategory::Edge => KernelLoadCategory::Edge,
            LossCategory::Vertex => KernelLoadCategory::Vertex,
            LossCategory::Representation => KernelLoadCategory::Representation,
            LossCategory::Solid => KernelLoadCategory::Solid,
        }
    }

    const fn map_reason(reason: LossReason) -> KernelLoadLossReason {
        match reason {
            LossReason::FaceUnresolved => KernelLoadLossReason::FaceUnresolved,
            LossReason::SurfaceUnresolved => KernelLoadLossReason::SurfaceUnresolved,
            LossReason::SurfaceRefused => KernelLoadLossReason::SurfaceRefused,
            LossReason::BoundUnresolved => KernelLoadLossReason::BoundUnresolved,
            LossReason::DegenerateVertexLoop => KernelLoadLossReason::DegenerateVertexLoop,
            LossReason::BoundNotALoopWeImplement => KernelLoadLossReason::BoundNotALoopWeImplement,
            LossReason::EdgeUseUnresolved => KernelLoadLossReason::EdgeUseUnresolved,
            LossReason::EdgeUnresolved => KernelLoadLossReason::EdgeUnresolved,
            LossReason::EdgeGeometryMissing => KernelLoadLossReason::EdgeGeometryMissing,
            LossReason::EdgeCurveRefused => KernelLoadLossReason::EdgeCurveRefused,
            LossReason::EdgeEndpointUnresolved => KernelLoadLossReason::EdgeEndpointUnresolved,
            LossReason::VertexUnresolved => KernelLoadLossReason::VertexUnresolved,
            LossReason::RelatedRepresentationUnresolved => {
                KernelLoadLossReason::RelatedRepresentationUnresolved
            }
            LossReason::RepresentationItemNotASolid => {
                KernelLoadLossReason::RepresentationItemNotASolid
            }
        }
    }

    #[cfg(test)]
    mod tests {
        use super::*;

        fn inspect_bytes(input: &[u8]) -> Result<KernelSourceLoadEvidence> {
            let kernel = MonstertruckKernel::project_standard();
            let document =
                kernel.resolve_source_exchange(ParsedStepSource::parse(input)?.exchange)?;
            kernel.inspect_source_document(&document)
        }

        #[test]
        fn maps_table_swallow_to_parser_neutral_source() -> Result<()> {
            let input = b"ISO-10303-21;\nHEADER;\nFILE_DESCRIPTION(('x'),'1');\nFILE_NAME('a','b',(''),(''),'x','y','');\nFILE_SCHEMA(('AUTOMOTIVE_DESIGN'));\nENDSEC;\nDATA;\n#11=ORIENTED_EDGE('bad');\nENDSEC;\nEND-ISO-10303-21;\n";
            let evidence = inspect_bytes(input)?;
            assert_eq!(evidence.sections.len(), 1);
            assert_eq!(evidence.sections[0].entity_losses.len(), 1);
            let loss = &evidence.sections[0].entity_losses[0];
            assert_eq!(
                loss.kind,
                KernelEntityLoadLossKind::UnexpectedParameterShape
            );
            assert_eq!(loss.entity_type, "ORIENTED_EDGE");
            assert_eq!(loss.tally.count, 1);
            assert_eq!(loss.tally.first_source.map(|source| source.id), Some(11));
            assert_eq!(
                loss.tally.first_source.map(|source| source.data_section),
                Some(0)
            );
            Ok(())
        }

        #[test]
        fn distinguishes_unretained_source_solid_from_conversion_failure() -> Result<()> {
            let unretained = b"ISO-10303-21;\nHEADER;\nFILE_DESCRIPTION(('x'),'1');\nFILE_NAME('a','b',(''),(''),'x','y','');\nFILE_SCHEMA(('AUTOMOTIVE_DESIGN'));\nENDSEC;\nDATA;\n#11=MANIFOLD_SOLID_BREP('bad');\nENDSEC;\nEND-ISO-10303-21;\n";
            let evidence = inspect_bytes(unretained)?;
            assert_eq!(evidence.solids.len(), 1);
            assert_eq!(
                evidence.solids[0].failure,
                Some(KernelSourceSolidLoadFailure::TableEntryUnavailable)
            );
            assert!(evidence.solids[0].topology_provenance.is_none());

            let conversion_failure = b"ISO-10303-21;\nHEADER;\nFILE_DESCRIPTION(('x'),'1');\nFILE_NAME('a','b',(''),(''),'x','y','');\nFILE_SCHEMA(('AUTOMOTIVE_DESIGN'));\nENDSEC;\nDATA;\n#11=MANIFOLD_SOLID_BREP('bad',#12);\nENDSEC;\nEND-ISO-10303-21;\n";
            let evidence = inspect_bytes(conversion_failure)?;
            assert_eq!(evidence.solids.len(), 1);
            let Some(KernelSourceSolidLoadFailure::ConversionFailed { detail }) =
                evidence.solids[0].failure.as_ref()
            else {
                anyhow::bail!("source solid conversion failure was not typed separately");
            };
            assert_ne!(detail, "");
            assert!(evidence.solids[0].topology_provenance.is_none());
            Ok(())
        }

        #[test]
        fn exact_native_sweep_and_spline_nets_supply_constructive_side_support() -> Result<()> {
            use monstertruck_io::step::load::step_geometry::{
                KnotVector, NurbsSurface, StepExtrusionSurface,
            };
            let cylinder_transform = Matrix4::from_cols(
                Vector4::new(3.0, 0.0, 0.0, 0.0),
                Vector4::new(0.0, 3.0, 0.0, 0.0),
                Vector4::new(0.0, 0.0, 1.0, 0.0),
                Vector4::new(1.0, -2.0, 0.0, 1.0),
            );
            let circle = Curve3D::Conic(Conic3D::Ellipse(Processor::with_transform(
                TrimmedCurve::new(UnitCircle::<Point3>::new(), (0.0, std::f64::consts::TAU)),
                cylinder_transform,
            )));
            let circle_wall = Surface::Sweep(SweepSurface::ExtrusionSurface(
                StepExtrusionSurface::by_extrusion(circle, Vector3::new(0.0, 0.0, 4.0)),
            ));
            let SurfaceSupport::Cylinder(cylinder) = native_surface_support(&circle_wall, 1) else {
                anyhow::bail!("native radius-three circular sweep did not become a cylinder");
            };
            assert!((cylinder.radius_mm - 3.0).abs() < 1.0e-12);
            assert_eq!(cylinder.axis_origin_mm, [1.0, -2.0, 0.0]);

            let plane_wall = Surface::Sweep(SweepSurface::ExtrusionSurface(
                StepExtrusionSurface::by_extrusion(
                    Curve3D::Line(Line(Point3::new(0.0, 0.0, 0.0), Point3::new(3.0, 0.0, 0.0))),
                    Vector3::new(0.0, 0.0, 4.0),
                ),
            ));
            assert!(matches!(
                native_surface_support(&plane_wall, 2),
                SurfaceSupport::Plane(_)
            ));

            let knots = (KnotVector::bezier_knot(2), KnotVector::bezier_knot(1));
            let points = [
                Point3::new(0.0, 0.0, 0.0),
                Point3::new(1.0, 1.0, 0.0),
                Point3::new(2.0, 0.0, 0.0),
            ];
            let polynomial_net = BsplineSurface::new(
                knots.clone(),
                points
                    .iter()
                    .map(|point| vec![*point, *point + Vector3::new(0.0, 0.0, 4.0)])
                    .collect(),
            )?;
            let SurfaceSupport::SplineExtrusion(polynomial) =
                native_surface_support(&Surface::Bspline(polynomial_net), 3)
            else {
                anyhow::bail!("polynomial native net did not supply spline extrusion support");
            };
            assert_eq!(polynomial.profile.degree, 2);
            assert_eq!(polynomial.profile.weights, None);
            assert_eq!(polynomial.extrusion_mm, [0.0, 0.0, 4.0]);

            // Homogeneous source weights may be globally rescaled without
            // changing the represented curve; recover a normalized profile.
            let weights = [4.0, 3.0, 4.0];
            let rational_net = BsplineSurface::new(
                knots,
                points
                    .iter()
                    .zip(weights)
                    .map(|(point, weight)| {
                        [
                            Vector4::new(
                                point.x * weight,
                                point.y * weight,
                                point.z * weight,
                                weight,
                            ),
                            Vector4::new(point.x * weight, point.y * weight, 4.0 * weight, weight),
                        ]
                        .to_vec()
                    })
                    .collect(),
            )?;
            let SurfaceSupport::SplineExtrusion(rational) =
                native_surface_support(&Surface::Nurbs(NurbsSurface::new(rational_net)), 4)
            else {
                anyhow::bail!("rational native net did not supply spline extrusion support");
            };
            assert_eq!(rational.profile.weights, Some(vec![1.0, 0.75, 1.0]));
            assert_eq!(rational.extrusion_mm, [0.0, 0.0, 4.0]);
            Ok(())
        }

        #[test]
        fn topology_editor_roundtrips_exact_selected_source_solid() -> Result<()> {
            let input = include_bytes!("../validation/fixtures/production_radial_slot.step");
            let kernel = MonstertruckKernel::project_standard();
            let document =
                kernel.resolve_source_exchange(ParsedStepSource::parse(input)?.exchange)?;
            let evidence = kernel.inspect_source_document(&document)?;
            let source = evidence
                .solids
                .first()
                .ok_or_else(|| anyhow::anyhow!("fixture has no source solid"))?
                .source;
            let edited = kernel.edit_source_topology(&document, source, &[])?;
            let snapshot = kernel.edited_source_to_step(&edited)?;
            let exchange = ParsedStepSource::parse(snapshot.as_bytes())?.exchange;
            assert_eq!(exchange.data.len(), 1);
            let document = ResolvedStepDocument::from_exchange(exchange);
            let section = document
                .table(0)
                .ok_or_else(|| anyhow::anyhow!("kernel-edited STEP omitted DATA section 0"))?;
            assert!(!section.manifold_solid_brep.is_empty());
            Ok(())
        }

        #[test]
        fn reports_source_solid_boundaries_positionally() -> Result<()> {
            let input = include_bytes!("../validation/fixtures/production_radial_slot.step");
            let evidence = inspect_bytes(input)?;
            assert_ne!(evidence.solids.len(), 0);
            for solid in &evidence.solids {
                assert!(solid.failure.is_none());
                assert_ne!(solid.boundaries.len(), 0);
                assert_eq!(solid.boundaries[0].boundary_index, 0);
                assert_eq!(solid.boundaries[0].role, KernelBoundaryRole::Outer);
                let provenance = solid
                    .topology_provenance
                    .as_ref()
                    .ok_or_else(|| anyhow::anyhow!("successful source solid omitted provenance"))?;
                assert_ne!(provenance.vertices.len(), 0);
                assert_ne!(provenance.edges.len(), 0);
                assert_ne!(provenance.faces.len(), 0);
                assert!(provenance.vertices.iter().all(|entry| {
                    entry.source.data_section == solid.source.data_section
                        && entry.topology_id.raw() != 0
                }));
                assert!(provenance.edges.iter().all(|entry| {
                    entry.source.data_section == solid.source.data_section
                        && entry.topology_id.raw() != 0
                }));
                assert!(provenance.faces.iter().all(|entry| {
                    entry.topology_id.raw() != 0
                        && entry
                            .source
                            .is_none_or(|source| source.data_section == solid.source.data_section)
                        && entry
                            .shell_member
                            .is_none_or(|source| source.data_section == solid.source.data_section)
                }));
            }
            Ok(())
        }
    }
}

#[cfg(test)]
mod provenance_lookup_tests {
    use super::*;

    const fn source(data_section: usize, id: u64) -> SourceEntity {
        SourceEntity {
            data_section,
            ordinal: 0,
            id,
        }
    }

    const fn face(
        source_id: Option<SourceEntity>,
        shell_member: Option<SourceEntity>,
        kernel_id: u64,
    ) -> KernelFaceTopologyProvenance {
        KernelFaceTopologyProvenance {
            source: source_id,
            shell_member,
            topology_id: KernelTopologyId::from_raw(kernel_id),
        }
    }

    fn provenance(faces: Vec<KernelFaceTopologyProvenance>) -> KernelTopologyProvenance {
        KernelTopologyProvenance {
            vertices: Vec::new(),
            edges: Vec::new(),
            faces,
        }
    }

    #[test]
    fn unique_shell_membership_has_priority_over_underlying_face() {
        let journal = provenance(vec![
            face(Some(source(0, 41)), Some(source(0, 9)), 100),
            face(Some(source(0, 9)), None, 200),
        ]);
        assert_eq!(
            journal.face_for_source_id(0, 9),
            Some(KernelTopologyId::from_raw(100)),
        );
    }

    #[test]
    fn underlying_face_lookup_is_used_only_when_shell_membership_is_absent() {
        let journal = provenance(vec![
            face(Some(source(0, 9)), None, 200),
            face(Some(source(1, 9)), Some(source(1, 9)), 300),
        ]);
        assert_eq!(
            journal.face_for_source_id(0, 9),
            Some(KernelTopologyId::from_raw(200)),
        );
        assert_eq!(
            journal.face_for_source_id(1, 9),
            Some(KernelTopologyId::from_raw(300)),
        );
        assert_eq!(journal.face_for_source_id(2, 9), None);
    }

    #[test]
    fn conflicting_shell_membership_never_falls_back_to_a_unique_underlying_face() {
        // In the old Option-based lookup, the two shell members became None,
        // and the distinct underlying face silently became an editable target.
        let journal = provenance(vec![
            face(None, Some(source(0, 9)), 100),
            face(None, Some(source(0, 9)), 101),
            face(Some(source(0, 9)), None, 200),
        ]);
        assert_eq!(journal.face_for_source_id(0, 9), None);
    }

    #[test]
    fn repeated_same_kernel_identity_is_unique_but_distinct_fallbacks_are_not() {
        let repeated = provenance(vec![
            face(None, Some(source(0, 9)), 100),
            face(None, Some(source(0, 9)), 100),
        ]);
        assert_eq!(
            repeated.face_for_source_id(0, 9),
            Some(KernelTopologyId::from_raw(100)),
        );
        let ambiguous_fallback = provenance(vec![
            face(Some(source(0, 9)), None, 100),
            face(Some(source(0, 9)), None, 101),
        ]);
        assert_eq!(ambiguous_fallback.face_for_source_id(0, 9), None);
    }

    #[test]
    fn vertices_and_edges_remain_fail_closed_on_conflicting_kernel_ids() {
        let mut journal = provenance(Vec::new());
        journal.vertices = vec![
            KernelTopologyEntityProvenance {
                source: source(0, 9),
                topology_id: KernelTopologyId::from_raw(100),
            },
            KernelTopologyEntityProvenance {
                source: source(0, 9),
                topology_id: KernelTopologyId::from_raw(101),
            },
        ];
        journal.edges = journal.vertices.clone();
        assert_eq!(journal.vertex_for_source_id(0, 9), None);
        assert_eq!(journal.edge_for_source_id(0, 9), None);
        assert_eq!(journal.vertex_for_source_id(1, 9), None);
    }
}

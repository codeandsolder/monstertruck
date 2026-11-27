//! Contact-constrained toroidal fillet proposals. Near tangency, fitting
//! independent surface residuals does not determine the intersection locus.
//! Recover the meridian center from a plane and coaxial revolution, and
//! its radius from the independently witnessed transverse tube.
use super::{
    JointFitCandidate, JointFitRejection, Neighborhood, SurfaceCarrier, SurfaceRegionPartition,
    satisfies,
};
use crate::{InnerSpace, KERNEL_EPSILON, Point3, Vector3};

#[derive(Clone, Copy)]
struct Torus {
    center: Point3,
    axis: Vector3,
    major: f64,
    minor: f64,
}

impl Torus {
    fn from_carrier(carrier: SurfaceCarrier) -> Option<Self> {
        let SurfaceCarrier::CircularRevolution {
            axial_center,
            axis,
            radial_center,
            profile_radius,
        } = carrier
        else {
            return None;
        };
        (radial_center > KERNEL_EPSILON && profile_radius > KERNEL_EPSILON).then_some(Self {
            center: axial_center,
            axis,
            major: radial_center,
            minor: profile_radius,
        })
    }

    fn carrier(self) -> SurfaceCarrier {
        SurfaceCarrier::CircularRevolution {
            axial_center: self.center,
            axis: self.axis,
            radial_center: self.major,
            profile_radius: self.minor,
        }
    }

    // Rigid placement plus meridian parameter variation bounds the
    // normalized signed distance change on the entire source convex hull.
    fn change_bound(self, candidate: Self, extent: f64) -> f64 {
        let angle = self
            .axis
            .cross(candidate.axis)
            .magnitude()
            .atan2(self.axis.dot(candidate.axis));
        (self.center - candidate.center).magnitude()
            + angle * extent
            + (self.major - candidate.major).abs()
            + (self.minor - candidate.minor).abs()
    }
}

#[derive(Clone, Copy)]
struct Revolution {
    origin: Point3,
    axis: Vector3,
    radius: f64,
    slope: f64,
}
impl Revolution {
    fn from_carrier(carrier: SurfaceCarrier) -> Option<Self> {
        let SurfaceCarrier::LinearRevolution {
            origin,
            axis,
            radius_at_origin,
            slope,
        } = carrier
        else {
            return None;
        };
        Some(Self {
            origin,
            axis,
            radius: radius_at_origin,
            slope,
        })
    }
    fn carrier(self) -> SurfaceCarrier {
        SurfaceCarrier::LinearRevolution {
            origin: self.origin,
            axis: self.axis,
            radius_at_origin: self.radius,
            slope: self.slope,
        }
    }
    fn radius_at(self, p: Point3) -> f64 {
        self.slope
            .mul_add((p - self.origin).dot(self.axis), self.radius)
    }
}

#[derive(Clone, Copy)]
struct Circle {
    center: Point3,
    axis: Vector3,
    radius: f64,
}
impl Circle {
    fn project(self, point: Point3) -> Option<Point3> {
        let d = point - self.center;
        let radial = d - self.axis * d.dot(self.axis);
        let length = radial.magnitude();
        (length > 0.0 && self.radius > 0.0).then(|| self.center + radial * (self.radius / length))
    }
}

fn anchored_torus(
    original: Torus,
    plane: (Point3, Vector3),
    coaxial: Revolution,
    tube: Revolution,
    angular: f64,
) -> Option<Torus> {
    let (origin, mut axis) = plane;
    if original.axis.dot(axis) < 0.0 {
        axis = -axis;
    }
    if original.axis.cross(axis).magnitude() > angular
        || coaxial.axis.cross(axis).magnitude() > angular
        || tube.axis.dot(axis).abs() > angular
    {
        return None;
    }
    let side = (original.center - origin).dot(axis).signum();
    let alignment = coaxial.axis.dot(axis);
    let at_plane =
        coaxial.origin - coaxial.axis * ((coaxial.origin - origin).dot(axis) / alignment);
    let shift = coaxial.axis * (side / alignment);
    // The transverse section radius is affine in the meridian center,
    // and that center is affine in the unknown minor radius. Solve the
    // resulting scalar equation instead of iterating a cone/circle fit.
    let denominator = (-tube.slope).mul_add(shift.dot(tube.axis), 1.0);
    if denominator.abs() <= f64::EPSILON {
        return None;
    }
    let minor = tube.radius_at(at_plane) / denominator;
    let center = at_plane + shift * minor;
    let radius = coaxial.radius_at(center);
    // Signed side of the affine meridian: circle center lies one normal
    // radius away from the adjacent straight meridian.
    let major = (original.major - radius)
        .signum()
        .mul_add(minor * 1.0_f64.hypot(coaxial.slope), radius);
    let candidate = Torus {
        center,
        axis,
        major,
        minor,
    };
    (minor > KERNEL_EPSILON && major > KERNEL_EPSILON).then_some(candidate)
}

fn contact_circle(
    torus: Torus,
    other: SurfaceCarrier,
    angular: f64,
    tolerance: f64,
) -> Option<Circle> {
    if let SurfaceCarrier::Plane { origin, normal } = other {
        if normal.cross(torus.axis).magnitude() > angular {
            return None;
        }
        let z = (origin - torus.center).dot(torus.axis);
        if (z.abs() - torus.minor).abs() > tolerance {
            return None;
        }
        return Some(Circle {
            center: torus.center + torus.axis * (z.signum() * torus.minor),
            axis: torus.axis,
            radius: torus.major,
        });
    }
    let line = Revolution::from_carrier(other)?;
    let relative = torus.center - line.origin;
    let nearest = line.origin + line.axis * relative.dot(line.axis);
    if line.axis.cross(torus.axis).magnitude() <= angular {
        let radius = line.radius_at(torus.center);
        let slope = line.slope * line.axis.dot(torus.axis);
        let distance = torus.major - radius;
        let scale = 1.0_f64.hypot(slope);
        if (distance.abs() / scale - torus.minor).abs() > tolerance
            || (nearest - torus.center).magnitude() > tolerance
        {
            return None;
        }
        let z = slope * distance / scale.powi(2);
        return Some(Circle {
            center: torus.center + torus.axis * z,
            axis: torus.axis,
            radius: torus.major - distance / scale.powi(2),
        });
    }
    if line.axis.dot(torus.axis).abs() > angular {
        return None;
    }
    let radial = nearest - torus.center;
    let radial = radial - torus.axis * radial.dot(torus.axis);
    if (radial.magnitude() - torus.major).abs() > tolerance
        || (line.radius_at(nearest) - torus.minor).abs() > tolerance
        || (nearest - torus.center).dot(torus.axis).abs() > tolerance
    {
        return None;
    }
    if radial.magnitude() <= f64::EPSILON {
        return None;
    }
    let direction = radial.normalize();
    let circle = Circle {
        center: torus.center + direction * torus.major,
        axis: torus.axis.cross(direction),
        radius: torus.minor,
    };
    // The source-facing radial branch and the complete native trim are
    // checked separately; the circle is a constructive candidate.
    Some(circle)
}

fn carriers_at(
    partition: &SurfaceRegionPartition,
    region: usize,
    neighbors: &[usize],
) -> Vec<SurfaceCarrier> {
    neighbors
        .iter()
        .copied()
        .filter(|&id| id != region)
        .filter_map(|id| partition.carriers[id])
        .collect()
}

fn proposals(original: Torus, neighbors: &[SurfaceCarrier], angular: f64) -> Vec<Torus> {
    let revolutions = neighbors
        .iter()
        .filter_map(|&c| Revolution::from_carrier(c))
        .collect::<Vec<_>>();
    let mut candidates = Vec::new();
    for &plane in neighbors {
        let SurfaceCarrier::Plane { origin, normal } = plane else {
            continue;
        };
        for &coaxial in &revolutions {
            for &tube in &revolutions {
                let Some(torus) =
                    anchored_torus(original, (origin, normal), coaxial, tube, angular)
                else {
                    continue;
                };
                let tolerance = KERNEL_EPSILON * 0.01;
                if [plane, coaxial.carrier(), tube.carrier()]
                    .into_iter()
                    .all(|other| contact_circle(torus, other, angular, tolerance).is_some())
                {
                    candidates.push(torus);
                }
            }
        }
    }
    candidates
}

pub(super) fn eligible(
    partition: &SurfaceRegionPartition,
    region: usize,
    neighbors: &[usize],
    angular: f64,
) -> bool {
    partition.carriers[region]
        .and_then(Torus::from_carrier)
        .is_some_and(|original| {
            !proposals(
                original,
                &carriers_at(partition, region, neighbors),
                angular,
            )
            .is_empty()
        })
}

fn nearest_junction(
    source: Point3,
    prior: Point3,
    local: &[SurfaceCarrier],
    circles: &[Circle],
    allowance: f64,
    angular: f64,
) -> Option<Point3> {
    let tolerance = KERNEL_EPSILON * 0.01;
    let admissible = |point: Point3| {
        (point - source).magnitude() <= allowance
            && satisfies(point, local, tolerance)
            && circles.iter().all(|circle| {
                circle
                    .project(point)
                    .is_some_and(|p| (p - point).magnitude() <= tolerance)
            })
    };
    circles
        .iter()
        .filter_map(|c| c.project(source))
        .filter(|&p| admissible(p))
        .min_by(|&a, &b| {
            (a - source)
                .magnitude2()
                .total_cmp(&(b - source).magnitude2())
        })
        .or_else(|| {
            super::fixed_sample_position(source, prior, None, local, allowance, angular)
                .filter(|&p| admissible(p))
        })
}

fn candidate_positions(
    neighborhood: &Neighborhood<'_>,
    region: usize,
    torus: Torus,
    prior: &[Point3],
) -> Result<(Vec<Point3>, Vec<f64>, f64), JointFitRejection> {
    let carrier = |id| {
        if id == region {
            Some(torus.carrier())
        } else {
            neighborhood.partition.carriers[id]
        }
    };
    let mut contacts = std::collections::BTreeMap::new();
    let angular = neighborhood.options.angular_resolution;
    let tolerance = KERNEL_EPSILON * 0.01;
    for &id in &neighborhood.context.region_neighbors[region] {
        if let Some(circle) = carrier(id).and_then(|c| contact_circle(torus, c, angular, tolerance))
        {
            contacts.insert(id, circle);
        }
    }
    let mut positions = prior.to_vec();
    let mut allowances = vec![0.0; prior.len()];
    let mut maximum = 0.0_f64;
    for sample in &neighborhood.samples {
        let source = neighborhood.evidence.positions[sample.id];
        let rejection = JointFitRejection::Sample {
            vertex: sample.id,
            allowance: sample.allowance,
        };
        let local = sample
            .regions
            .iter()
            .map(|&id| carrier(id))
            .collect::<Option<Vec<_>>>()
            .ok_or(rejection)?;
        allowances[sample.id] = sample.allowance;
        if !sample.boundary {
            if !satisfies(source, &local, sample.allowance) {
                return Err(rejection);
            }
            continue;
        }
        let circles = sample
            .regions
            .iter()
            .filter_map(|id| contacts.get(id).copied())
            .collect::<Vec<_>>();
        let point = nearest_junction(
            source,
            prior[sample.id],
            &local,
            &circles,
            sample.allowance,
            angular,
        )
        .ok_or(rejection)?;
        maximum = maximum.max((point - source).magnitude());
        positions[sample.id] = point;
    }
    if !super::preserves_polygon_orientation(
        neighborhood.evidence,
        &positions,
        &neighborhood.supports,
    ) {
        return Err(JointFitRejection::Orientation);
    }
    Ok((positions, allowances, maximum))
}

pub(super) fn fit(
    neighborhood: &Neighborhood<'_>,
    region: usize,
    prior: &[Point3],
) -> Result<JointFitCandidate, JointFitRejection> {
    let original = Torus::from_carrier(
        neighborhood.partition.carriers[region].ok_or(JointFitRejection::Neighborhood)?,
    )
    .ok_or(JointFitRejection::Unchanged)?;
    let neighbors = carriers_at(
        neighborhood.partition,
        region,
        &neighborhood.context.region_neighbors[region],
    );
    let extent = neighborhood
        .samples
        .iter()
        .map(|s| (neighborhood.evidence.positions[s.id] - original.center).magnitude())
        .fold(0.0_f64, f64::max);
    let budget = neighborhood
        .samples
        .iter()
        .map(|s| s.allowance)
        .fold(0.0_f64, f64::max)
        - neighborhood.context.carrier_change[region];
    let mut candidates = proposals(
        original,
        &neighbors,
        neighborhood.options.angular_resolution,
    )
    .into_iter()
    .filter_map(|torus| {
        let bound = original.change_bound(torus, extent);
        (bound.is_finite() && bound > 0.0 && bound <= budget).then_some((bound, torus))
    })
    .collect::<Vec<_>>();
    candidates.sort_by(|a, b| a.0.total_cmp(&b.0));
    let mut last = JointFitRejection::Unchanged;
    for (bound, torus) in candidates {
        match candidate_positions(neighborhood, region, torus, prior) {
            Ok((positions, allowances, maximum_displacement)) => {
                let mut partition = neighborhood.partition.clone();
                partition.carriers[region] = Some(torus.carrier());
                return Ok(JointFitCandidate {
                    partition,
                    positions,
                    allowances,
                    maximum_displacement,
                    refitted_regions: 1,
                    max_additional_carrier_residual: bound,
                    carrier_changes: vec![(region, bound)],
                });
            }
            Err(error) => last = error,
        }
    }
    Err(last)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn fillet_meridian_is_constrained_by_independent_neighbors()
    -> Result<(), Box<dyn std::error::Error>> {
        let original = Torus {
            center: Point3::new(0.0, 0.0, -0.000_000_08),
            axis: Vector3::unit_z(),
            major: 0.600_000_1,
            minor: 0.050_000_12,
        };
        let coaxial = Revolution {
            origin: Point3::new(0.0, 0.0, 1.0),
            axis: Vector3::unit_z(),
            radius: 0.55,
            slope: 0.0,
        };
        let tube = Revolution {
            origin: Point3::new(0.6, 0.0, 0.0),
            axis: Vector3::unit_y(),
            radius: 0.05,
            slope: 0.0,
        };
        let plane = SurfaceCarrier::Plane {
            origin: Point3::new(0.0, 0.0, 0.05),
            normal: Vector3::unit_z(),
        };
        let torus = anchored_torus(
            original,
            (Point3::new(0.0, 0.0, 0.05), Vector3::unit_z()),
            coaxial,
            tube,
            1e-8,
        )
        .ok_or("anchored torus")?;
        assert!(original.change_bound(torus, 1.0) < 1e-6);
        for (carrier, point) in [
            (plane, Point3::new(0.600_06, 0.0, 0.05)),
            (
                SurfaceCarrier::LinearRevolution {
                    origin: tube.origin,
                    axis: tube.axis,
                    radius_at_origin: tube.radius,
                    slope: 0.0,
                },
                Point3::new(0.6, 0.0, 0.05),
            ),
            (
                SurfaceCarrier::LinearRevolution {
                    origin: coaxial.origin,
                    axis: coaxial.axis,
                    radius_at_origin: coaxial.radius,
                    slope: 0.0,
                },
                Point3::new(0.55, 0.0, 0.0),
            ),
        ] {
            let projected = contact_circle(torus, carrier, 1e-8, 1e-8)
                .and_then(|c| c.project(point))
                .ok_or("contact circle")?;
            assert!(satisfies(projected, &[torus.carrier(), carrier], 1e-10));
        }
        Ok(())
    }
    // Closed annular-sector boundary is not necessary to test the joint
    // transaction; every torus boundary edge has a second source carrier.
    // The open outer boundaries of the four strips remain unchanged.
    fn append_strip(
        positions: &mut Vec<Point3>,
        faces: &mut Vec<Vec<usize>>,
        labels: &mut Vec<usize>,
        boundary: &[usize],
        label: usize,
        offset: impl Fn(Point3) -> Point3,
    ) {
        let start = positions.len();
        positions.extend(
            boundary
                .iter()
                .map(|&id| offset(positions[id]))
                .collect::<Vec<_>>(),
        );
        for (i, pair) in boundary.windows(2).enumerate() {
            let [a, b] = [pair[0], pair[1]];
            let [c, d] = [start + i, start + i + 1];
            // Reverse the retained torus boundary orientation.
            faces.extend([vec![b, a, c], vec![b, c, d]]);
            labels.extend([label, label]);
        }
    }

    fn fillet_mesh() -> Result<(super::super::MeshEvidence, Vec<usize>), Box<dyn std::error::Error>>
    {
        use crate::mesh_recovery::{EvidenceOptions, analyze_mesh_evidence};
        use monstertruck_mesh::{Faces, PolygonMesh, StandardAttributes};
        use std::f64::consts::FRAC_PI_2;
        let count = 8_usize;
        let stride = count + 1;
        let mut positions = Vec::new();
        let mut faces = Vec::new();
        let mut labels = Vec::new();
        for i in 0..=count {
            let u = FRAC_PI_2 * crate::usize_to_f64(i) / crate::usize_to_f64(count);
            let (sin_u, cos_u) = u.sin_cos();
            for j in 0..=count {
                let v = FRAC_PI_2 + FRAC_PI_2 * crate::usize_to_f64(j) / crate::usize_to_f64(count);
                let (sin_v, cos_v) = v.sin_cos();
                // Small radial source noise on the tangent plane rim,
                // but not at the two triply constrained end junctions.
                let noise = if j == 0 && i != 0 && i != count {
                    0.000_03 * (2.0 * u).sin()
                } else {
                    0.0
                };
                let rho = 0.05_f64.mul_add(cos_v, 0.6) + noise;
                positions.push(Point3::new(rho * cos_u, rho * sin_u, 0.05 * sin_v));
            }
        }
        for i in 0..count {
            for j in 0..count {
                let a = i * stride + j;
                let [b, c, d] = [a + stride, a + stride + 1, a + 1];
                faces.extend([vec![a, b, c], vec![a, c, d]]);
                labels.extend([0, 0]);
            }
        }
        // Each loop is directed along the source torus's boundary.
        let bottom = (0..=count).map(|i| i * stride).collect::<Vec<_>>();
        let top = (0..=count)
            .rev()
            .map(|i| i * stride + count)
            .collect::<Vec<_>>();
        let left = (0..=count).rev().collect::<Vec<_>>();
        let right = (0..=count).map(|j| count * stride + j).collect::<Vec<_>>();
        append_strip(&mut positions, &mut faces, &mut labels, &bottom, 1, |p| {
            Point3::new(p.x * 1.1, p.y * 1.1, p.z)
        });
        append_strip(&mut positions, &mut faces, &mut labels, &top, 2, |p| {
            p - Vector3::unit_z() * 0.1
        });
        append_strip(&mut positions, &mut faces, &mut labels, &left, 3, |p| {
            p - Vector3::unit_y() * 0.1
        });
        append_strip(&mut positions, &mut faces, &mut labels, &right, 4, |p| {
            p - Vector3::unit_x() * 0.1
        });
        let mesh = PolygonMesh::new(
            StandardAttributes {
                positions,
                ..Default::default()
            },
            Faces::from_iter(faces),
        )?;
        let evidence = analyze_mesh_evidence(&mesh, EvidenceOptions::default())?;
        let support_labels = evidence
            .planar_supports
            .iter()
            .map(|support| {
                let first = *support.source_face_indices.first().ok_or("empty support")?;
                let label = labels[evidence.provenance.original_polygon_ids[first]];
                if !support
                    .source_face_indices
                    .iter()
                    .all(|&id| labels[evidence.provenance.original_polygon_ids[id]] == label)
                {
                    return Err("support crosses fixture carrier labels");
                }
                Ok(label)
            })
            .collect::<Result<Vec<_>, _>>()?;
        Ok((evidence, support_labels))
    }

    #[test]
    fn constrained_fillet_transaction_repairs_shared_tangent_rim()
    -> Result<(), Box<dyn std::error::Error>> {
        use super::super::{JointFitContext, JointFitOptions};
        use crate::mesh_recovery::{NativeTrimOptions, RegionBoundaryIndex, SharedTrimOptions};
        let (evidence, labels) = fillet_mesh()?;
        let index = RegionBoundaryIndex::new(&evidence)?;
        let cylinder = |origin, axis, radius| SurfaceCarrier::LinearRevolution {
            origin,
            axis,
            radius_at_origin: radius,
            slope: 0.0,
        };
        let carriers = [
            Torus {
                center: Point3::new(0.0, 0.0, -0.000_000_08),
                axis: Vector3::unit_z(),
                major: 0.600_000_1,
                minor: 0.050_000_12,
            }
            .carrier(),
            SurfaceCarrier::Plane {
                origin: Point3::new(0.0, 0.0, 0.05),
                normal: Vector3::unit_z(),
            },
            cylinder(Point3::new(0.0, 0.0, 0.0), Vector3::unit_z(), 0.55),
            cylinder(Point3::new(0.6, 0.0, 0.0), Vector3::unit_y(), 0.05),
            cylinder(Point3::new(0.0, 0.6, 0.0), Vector3::unit_x(), 0.05),
        ];
        let partition = index.repartition_connected(&labels, |id| carriers.get(id).copied())?;
        let torus = partition.region_by_polygon[0];
        let policy = JointFitOptions {
            maximum_adjustment: 0.001,
            declared_uncertainty: 0.0001,
            angular_resolution: 1e-8,
        };
        let context = JointFitContext::new(&evidence, &partition, policy).ok_or("context")?;
        assert!(context.anchored_fillet_possible(&partition, torus));
        let neighborhood =
            Neighborhood::new(&evidence, &partition, &[torus], &context).ok_or("neighborhood")?;
        let candidate = fit(&neighborhood, torus, &evidence.positions)?;
        assert_eq!(candidate.refitted_regions, 1);
        assert!(candidate.maximum_displacement <= policy.declared_uncertainty);
        assert!(candidate.max_additional_carrier_residual < 1e-6);
        for sample in &neighborhood.samples {
            let source = evidence.positions[sample.id];
            let p = candidate.positions[sample.id];
            assert!((p - source).magnitude() <= sample.allowance);
            if sample.boundary {
                let local = sample
                    .regions
                    .iter()
                    .filter_map(|&id| candidate.partition.carriers[id])
                    .collect::<Vec<_>>();
                assert!(satisfies(p, &local, 1e-8));
            }
        }
        let trims = index.assemble_shared_trims(
            &candidate.partition,
            &candidate.positions,
            SharedTrimOptions {
                absolute_tolerance: 1e-7,
                relative_tolerance: 1e-4,
                max_interval_edges: 128,
            },
        )?;
        let native = index.construct_native_trims(
            &candidate.partition,
            &trims,
            &candidate.positions,
            NativeTrimOptions {
                max_vertex_displacement: 1e-7,
                intersection_tolerance: 1e-8,
                angular_resolution: 1e-8,
            },
        )?;
        let region = native
            .regions
            .iter()
            .find(|r| r.region == torus)
            .ok_or("torus region")?;
        assert!(region.faces.is_ok(), "{:?}", region.faces.as_ref().err());
        assert!(!region.requires_chart_cut);
        for (trim, edges) in trims.trims.iter().zip(native.trims) {
            if trim.regions.1.is_some()
                && (trim.regions.0 == torus || trim.regions.1 == Some(torus))
            {
                assert!(edges.is_ok(), "{:?}", edges.err());
            }
        }
        Ok(())
    }
}

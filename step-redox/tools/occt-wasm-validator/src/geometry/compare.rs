use super::{GeometryReport, Solid};
use anyhow::{Result, ensure};
use clap::ValueEnum;
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Copy, Serialize, Deserialize, ValueEnum, Default)]
#[serde(rename_all = "snake_case")]
pub enum Alignment {
    #[default]
    Identity,
    BboxCenter,
    CenterOfMass,
}

#[derive(Debug, Default, Deserialize, Serialize)]
#[serde(default, deny_unknown_fields)]
pub struct Thresholds {
    topology_exact: bool,
    volume_rel: Option<f64>,
    area_rel: Option<f64>,
    edge_length_rel: Option<f64>,
    solid_volume_sum_rel: Option<f64>,
    solid_area_sum_rel: Option<f64>,
    bbox_extent_abs_mm: Option<f64>,
    center_abs_mm: Option<f64>,
    com_abs_mm: Option<f64>,
    inertia_rel: Option<f64>,
    solid_topology_exact: bool,
    solid_center_abs_mm: Option<f64>,
    solid_bbox_extent_abs_mm: Option<f64>,
    solid_volume_rel: Option<f64>,
    solid_area_rel: Option<f64>,
    // The independent raster stage consumes these, not the metric comparator.
    pub render_mean_abs: Option<f64>,
    pub render_changed_fraction: Option<f64>,
    pub silhouette_xor: Option<f64>,
    pub strict_mesh_hash: bool,
}

#[derive(Debug, Serialize)]
#[serde(untagged)]
enum Number {
    Count(i64),
    Measure(f64),
}

#[derive(Debug, Serialize)]
pub struct Check {
    name: String,
    value: Number,
    limit: Number,
    #[serde(skip_serializing_if = "Option::is_none")]
    unit: Option<&'static str>,
    #[serde(skip_serializing_if = "Option::is_none")]
    reference: Option<u32>,
    #[serde(skip_serializing_if = "Option::is_none")]
    candidate: Option<u32>,
    pass: bool,
}

impl Check {
    fn measure(name: &str, value: f64, limit: f64, unit: &'static str) -> Result<Self> {
        ensure!(
            value.is_finite() && value >= 0.0,
            "{name}: invalid measured error"
        );
        ensure!(
            limit.is_finite() && limit >= 0.0,
            "{name}: invalid threshold"
        );
        Ok(Self {
            name: name.into(),
            value: Number::Measure(value),
            limit: Number::Measure(limit),
            unit: Some(unit),
            reference: None,
            candidate: None,
            pass: value <= limit,
        })
    }
    fn count(name: &str, value: i64, reference: Option<u32>, candidate: Option<u32>) -> Self {
        Self {
            name: name.into(),
            value: Number::Count(value),
            limit: Number::Count(0),
            unit: None,
            reference,
            candidate,
            pass: value == 0,
        }
    }
    fn exact(name: &str, reference: u32, candidate: u32) -> Self {
        Self::count(
            name,
            i64::from(candidate) - i64::from(reference),
            Some(reference),
            Some(candidate),
        )
    }
}

#[derive(Debug, Serialize)]
pub struct Comparison {
    alignment: Alignment,
    pub alignment_translation_mm: [f64; 3],
    checks: Vec<Check>,
    pub metrics_pass: bool,
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub(crate) deferred_checks: Vec<&'static str>,
}

fn difference<const N: usize>(a: [f64; N], b: [f64; N]) -> [f64; N] {
    std::array::from_fn(|i| a[i] - b[i])
}
fn translated_delta(a: [f64; 3], b: [f64; 3], shift: [f64; 3]) -> [f64; 3] {
    difference(a, std::array::from_fn(|i| b[i] + shift[i]))
}
fn norm(values: impl IntoIterator<Item = f64>) -> f64 {
    values.into_iter().map(|x| x * x).sum::<f64>().sqrt()
}
fn max_abs(values: impl IntoIterator<Item = f64>) -> f64 {
    values.into_iter().map(f64::abs).fold(0.0, f64::max)
}
fn rel_error(a: f64, b: f64) -> f64 { (a - b).abs() / a.abs().max(b.abs()).max(1e-30) }

#[derive(Debug)]
struct Match {
    index: usize,
    topology_same: bool,
    center: f64,
    extent: f64,
    volume: f64,
    area: f64,
}
impl Match {
    fn new(reference: &Solid, candidate: &Solid, index: usize, shift: [f64; 3]) -> Self {
        Self {
            index,
            topology_same: reference.topology == candidate.topology,
            center: norm(translated_delta(
                reference.bbox_center,
                candidate.bbox_center,
                shift,
            )),
            extent: max_abs(difference(reference.bbox_extent, candidate.bbox_extent)),
            volume: rel_error(reference.volume, candidate.volume),
            area: rel_error(reference.area, candidate.area),
        }
    }
    fn score(&self) -> f64 {
        (if self.topology_same { 0.0 } else { 1e9 })
            + self.center * 1e6
            + self.extent * 1e5
            + self.volume * 1e4
            + self.area * 1e4
    }
}

fn append_measures(
    checks: &mut Vec<Check>,
    measures: &[(&str, f64, Option<f64>, &'static str)],
) -> Result<()> {
    for &(name, value, limit, unit) in measures {
        if let Some(limit) = limit {
            checks.push(Check::measure(name, value, limit, unit)?);
        }
    }
    Ok(())
}

fn solid_checks(
    reference: &[Solid],
    candidate: &[Solid],
    shift: [f64; 3],
    t: &Thresholds,
    checks: &mut Vec<Check>,
) -> Result<()> {
    if !t.solid_topology_exact
        && [
            t.solid_center_abs_mm,
            t.solid_bbox_extent_abs_mm,
            t.solid_volume_rel,
            t.solid_area_rel,
        ]
        .iter()
        .all(Option::is_none)
    {
        return Ok(());
    }
    if reference.len() != candidate.len() {
        checks.push(Check::exact(
            "solid_inventory_count",
            u32::try_from(reference.len())?,
            u32::try_from(candidate.len())?,
        ));
        return Ok(());
    }
    // Stable sort preserves source order for indistinguishable references,
    // while the candidate index breaks equal-score ties exactly as the oracle.
    let mut order: Vec<_> = reference.iter().collect();
    order.sort_by(|a, b| {
        b.topology
            .faces
            .cmp(&a.topology.faces)
            .then_with(|| b.volume.total_cmp(&a.volume))
            .then_with(|| b.area.total_cmp(&a.area))
    });
    let mut matched = vec![false; candidate.len()];
    let (mut center, mut extent, mut volume, mut area) = (0.0_f64, 0.0_f64, 0.0_f64, 0.0_f64);
    let mut mismatches = 0_i64;
    for reference in order {
        let best = candidate
            .iter()
            .enumerate()
            .filter(|(i, _)| !matched[*i])
            .map(|(i, candidate)| Match::new(reference, candidate, i, shift))
            .min_by(|a, b| {
                a.score()
                    .total_cmp(&b.score())
                    .then_with(|| a.index.cmp(&b.index))
            })
            .ok_or_else(|| anyhow::anyhow!("solid inventory matcher exhausted candidates"))?;
        ensure!(best.score().is_finite(), "nonfinite solid match score");
        matched[best.index] = true;
        mismatches += i64::from(!best.topology_same);
        center = center.max(best.center);
        extent = extent.max(best.extent);
        volume = volume.max(best.volume);
        area = area.max(best.area);
    }
    if t.solid_topology_exact {
        checks.push(Check::count(
            "solid_topology_mismatches",
            mismatches,
            None,
            None,
        ));
    }
    append_measures(
        checks,
        &[
            ("solid_center_abs_mm", center, t.solid_center_abs_mm, "mm"),
            (
                "solid_bbox_extent_abs_mm",
                extent,
                t.solid_bbox_extent_abs_mm,
                "mm",
            ),
            ("solid_volume_rel", volume, t.solid_volume_rel, ""),
            ("solid_area_rel", area, t.solid_area_rel, ""),
        ],
    )
}

pub fn compare(
    reference: &GeometryReport,
    candidate: &GeometryReport,
    alignment: Alignment,
    t: &Thresholds,
) -> Result<Comparison> {
    let ref_bounds = &reference.bounds;
    let cand_bounds = &candidate.bounds;
    let shift = match alignment {
        Alignment::Identity => [0.0; 3],
        Alignment::BboxCenter => difference(ref_bounds.bbox_center, cand_bounds.bbox_center),
        Alignment::CenterOfMass => difference(reference.center_of_mass, candidate.center_of_mass),
    };
    ensure!(shift.iter().all(|x| x.is_finite()), "nonfinite alignment");
    let mut checks = vec![Check::measure(
        "brep_valid",
        u8::from(!candidate.brep_valid).into(),
        0.0,
        "",
    )?];
    if t.topology_exact {
        for (name, reference, candidate) in [
            ("solids_exact", reference.solids, candidate.solids),
            ("shells_exact", reference.shells, candidate.shells),
            (
                "faces_exact",
                reference.topology.faces,
                candidate.topology.faces,
            ),
            (
                "edges_exact",
                reference.topology.edges,
                candidate.topology.edges,
            ),
            (
                "vertices_exact",
                reference.topology.vertices,
                candidate.topology.vertices,
            ),
        ] {
            checks.push(Check::exact(name, reference, candidate));
        }
    }
    let inertia_ref = reference.inertia.as_flattened();
    let inertia_cand = candidate.inertia.as_flattened();
    let inertia_delta = norm(inertia_ref.iter().zip(inertia_cand).map(|(a, b)| a - b))
        / norm(inertia_ref.iter().copied())
            .max(norm(inertia_cand.iter().copied()))
            .max(1e-30);
    append_measures(
        &mut checks,
        &[
            (
                "volume_rel",
                rel_error(reference.volume, candidate.volume),
                t.volume_rel,
                "",
            ),
            (
                "area_rel",
                rel_error(reference.area, candidate.area),
                t.area_rel,
                "",
            ),
            (
                "edge_length_rel",
                rel_error(reference.edge_length, candidate.edge_length),
                t.edge_length_rel,
                "",
            ),
            (
                "solid_volume_sum_rel",
                rel_error(reference.solid_volume_sum, candidate.solid_volume_sum),
                t.solid_volume_sum_rel,
                "",
            ),
            (
                "solid_area_sum_rel",
                rel_error(reference.solid_area_sum, candidate.solid_area_sum),
                t.solid_area_sum_rel,
                "",
            ),
            (
                "bbox_extent_abs_mm",
                max_abs(difference(ref_bounds.bbox_extent, cand_bounds.bbox_extent)),
                t.bbox_extent_abs_mm,
                "mm",
            ),
            (
                "center_abs_mm",
                norm(translated_delta(
                    ref_bounds.bbox_center,
                    cand_bounds.bbox_center,
                    shift,
                )),
                t.center_abs_mm,
                "mm",
            ),
            (
                "com_abs_mm",
                norm(translated_delta(
                    reference.center_of_mass,
                    candidate.center_of_mass,
                    shift,
                )),
                t.com_abs_mm,
                "mm",
            ),
            ("inertia_rel", inertia_delta, t.inertia_rel, ""),
        ],
    )?;
    solid_checks(
        &reference.solid_inventory,
        &candidate.solid_inventory,
        shift,
        t,
        &mut checks,
    )?;
    let deferred_checks = [
        ("render_mean_abs", t.render_mean_abs.is_some()),
        (
            "render_changed_fraction",
            t.render_changed_fraction.is_some(),
        ),
        ("silhouette_xor", t.silhouette_xor.is_some()),
        ("strict_mesh_hash", t.strict_mesh_hash),
    ]
    .into_iter()
    .filter_map(|(name, enabled)| enabled.then_some(name))
    .collect();
    Ok(Comparison {
        alignment,
        alignment_translation_mm: shift,
        deferred_checks,
        metrics_pass: checks.iter().all(|check| check.pass),
        checks,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::{Bounds, Topology};

    fn topology() -> Topology {
        Topology {
            faces: 6,
            edges: 24,
            vertices: 48,
        }
    }

    fn report(offset: f64) -> GeometryReport {
        let solid = |x| Solid {
            volume: 1.0,
            area: 6.0,
            bbox_center: [x + offset, 0.5, 0.5],
            bbox_extent: [1.0; 3],
            topology: topology(),
        };
        GeometryReport {
            path: "example.step".into(),
            bytes: 42,
            roots: 1,
            transferred: 1,
            brep_valid: true,
            volume: 2.0,
            area: 12.0,
            edge_length: 48.0,
            center_of_mass: [1.5 + offset, 0.5, 0.5],
            inertia: [[1.0, 0.0, 0.0], [0.0, 3.0, 0.0], [0.0, 0.0, 3.0]],
            bounds: Bounds {
                bbox: [offset, 0.0, 0.0, offset + 3.0, 1.0, 1.0],
                bbox_center: [1.5 + offset, 0.5, 0.5],
                bbox_extent: [3.0, 1.0, 1.0],
            },
            topology: Topology {
                faces: 12,
                edges: 48,
                vertices: 96,
            },
            solids: 2,
            shells: 2,
            solid_volume_sum: 2.0,
            solid_area_sum: 12.0,
            solid_inventory: vec![solid(0.5), solid(2.5)],
        }
    }

    fn all_thresholds() -> Thresholds {
        Thresholds {
            topology_exact: true,
            volume_rel: Some(0.0),
            area_rel: Some(0.0),
            edge_length_rel: Some(0.0),
            solid_volume_sum_rel: Some(0.0),
            solid_area_sum_rel: Some(0.0),
            bbox_extent_abs_mm: Some(0.0),
            center_abs_mm: Some(0.0),
            com_abs_mm: Some(0.0),
            inertia_rel: Some(0.0),
            solid_topology_exact: true,
            solid_center_abs_mm: Some(0.0),
            solid_bbox_extent_abs_mm: Some(0.0),
            solid_volume_rel: Some(0.0),
            solid_area_rel: Some(0.0),
            render_mean_abs: Some(0.001),
            strict_mesh_hash: true,
            ..Thresholds::default()
        }
    }

    #[test]
    fn alignment_and_reordered_solids() -> Result<()> {
        let reference = report(0.0);
        let mut candidate = report(17.0);
        candidate.solid_inventory.reverse();
        let thresholds = all_thresholds();
        for alignment in [Alignment::BboxCenter, Alignment::CenterOfMass] {
            let comparison = compare(&reference, &candidate, alignment, &thresholds)?;
            assert!(comparison.metrics_pass);
            assert_eq!(comparison.checks.len(), 20);
            assert_eq!(
                comparison.deferred_checks,
                ["render_mean_abs", "strict_mesh_hash"]
            );
        }
        assert!(!compare(&reference, &candidate, Alignment::Identity, &thresholds)?.metrics_pass);
        Ok(())
    }

    #[test]
    fn invalid_geometry_and_thresholds_fail() -> Result<()> {
        let reference = report(0.0);
        let mut candidate = report(0.0);
        let mut thresholds = all_thresholds();
        candidate.solid_inventory.pop();
        let comparison = compare(&reference, &candidate, Alignment::Identity, &thresholds)?;
        assert!(!comparison.metrics_pass);
        assert!(
            comparison
                .checks
                .iter()
                .any(|c| c.name == "solid_inventory_count" && !c.pass)
        );
        candidate = report(0.0);
        candidate.brep_valid = false;
        assert!(!compare(&reference, &candidate, Alignment::Identity, &thresholds)?.metrics_pass);
        thresholds.volume_rel = Some(-1.0);
        assert!(compare(&reference, &candidate, Alignment::Identity, &thresholds).is_err());
        assert!(serde_json::from_str::<Thresholds>(r#"{"volme_rel":0.01}"#).is_err());
        Ok(())
    }

    #[test]
    fn geometry_report_round_trips_native_schema() -> Result<()> {
        let original = report(0.0);
        let bytes = serde_json::to_vec(&original)?;
        let restored = serde_json::from_slice(&bytes)?;
        assert!(
            compare(&original, &restored, Alignment::Identity, &all_thresholds())?.metrics_pass
        );
        Ok(())
    }
}

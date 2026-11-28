use crate::{
    complete_ir, dominant_extrusions, formed_sheet, periodic_chains, solid_extrusions,
    solid_revolutions,
};
use anyhow::Result;
use std::collections::HashMap;
use std::ops::Range;
use std::sync::OnceLock;
use step_p21::ast::{EntityInstance, Exchange};

type EntityIndex = HashMap<u64, usize>;
type SourceTolerances = HashMap<u64, f64>;

/// Read-only STEP evidence, initialized only when a detector needs it.
/// Caches are section-local: STEP IDs are not unique across DATA sections.
#[derive(Default)]
struct SectionCache {
    index: OnceLock<EntityIndex>,
    source_tolerances: OnceLock<SourceTolerances>,
}

impl SectionCache {
    fn index(&self, entities: &[EntityInstance]) -> &EntityIndex {
        self.index
            .get_or_init(|| crate::step_graph::build_index(entities))
    }

    fn source_tolerances(&self, entities: &[EntityInstance]) -> &SourceTolerances {
        self.source_tolerances.get_or_init(|| {
            crate::source_accuracy::source_tolerance_by_representation_item(
                entities,
                self.index(entities),
            )
        })
    }
}

/// A single owned signature tree with the exact output span of each source
/// DATA section. Counts of source solids are insufficient: the detector skips
/// solids whose face/topology evidence cannot be read.
struct SurfaceSignatureCache {
    signatures: Vec<solid_revolutions::SolidSurfaceSignature>,
    section_ranges: Vec<Range<usize>>,
}

/// Owned STEP analysis context with reusable, lazily initialized per-section
/// entity indexes and representation-context uncertainty evidence.
///
/// All detectors remain read-only. Public standalone detectors still accept
/// arbitrary entity slices without relying on a session.
pub struct AnalysisSession {
    exchange: Exchange,
    sections: Vec<SectionCache>,
    // Surface signatures are deep owned trees; retain them once and lend
    // references rather than cloning a large tree for every diagnostic call.
    surface_signatures: OnceLock<SurfaceSignatureCache>,
}

impl AnalysisSession {
    fn from_exchange(exchange: Exchange) -> Self {
        let sections = std::iter::repeat_with(SectionCache::default)
            .take(exchange.data.len())
            .collect();
        Self {
            exchange,
            sections,
            surface_signatures: OnceLock::new(),
        }
    }

    /// Decode and parse one supported STEP exchange.
    ///
    /// # Errors
    /// Returns an error when the input cannot be decoded or parsed, or when it
    /// contains optional exchange sections that step-redox does not support.
    pub fn from_bytes(input: &[u8]) -> Result<Self> {
        let exchange = monstertruck_io::step::load::ParsedStepSource::parse(input)?.exchange;
        Ok(Self::from_exchange(exchange))
    }

    fn indexed_sections(
        &self,
    ) -> impl Iterator<Item = (&[EntityInstance], &EntityIndex, &SectionCache)> {
        self.exchange
            .data
            .iter()
            .zip(&self.sections)
            .map(|(section, cache)| {
                let entities = section.entities.as_slice();
                (entities, cache.index(entities), cache)
            })
    }

    #[must_use = "The computed formed sheet evidence collection is returned and must be retained or inspected."]
    pub fn formed_sheet_evidence(&self) -> Vec<formed_sheet::FormedSheetEvidence> {
        self.indexed_sections()
            .flat_map(|(entities, index, _)| {
                formed_sheet::detect_formed_sheet_evidence_with_index(entities, index)
            })
            .collect()
    }

    /// Borrow the session's immutable surface/loop signatures. The first call
    /// builds them from each DATA section in source order; subsequent calls
    /// return the same evidence without geometry traversal or deep cloning.
    ///
    /// Call `.to_vec()` explicitly only when an owned snapshot is necessary.
    #[must_use = "The solid surface signatures accessor exposes its value only through the returned reference."]
    pub fn solid_surface_signatures(&self) -> &[solid_revolutions::SolidSurfaceSignature] {
        &self
            .surface_signatures
            .get_or_init(|| {
                let mut signatures = Vec::new();
                let mut section_ranges = Vec::with_capacity(self.sections.len());
                for (entities, index, _) in self.indexed_sections() {
                    let start = signatures.len();
                    signatures.extend(
                        solid_revolutions::detect_solid_surface_signatures_with_index(
                            entities, index,
                        ),
                    );
                    section_ranges.push(start..signatures.len());
                }
                SurfaceSignatureCache {
                    signatures,
                    section_ranges,
                }
            })
            .signatures
    }

    #[must_use = "The computed dominant extrusion candidates collection is returned and must be retained or inspected."]
    pub fn dominant_extrusion_candidates(
        &self,
    ) -> Vec<dominant_extrusions::DominantExtrusionCandidate> {
        dominant_extrusions::detect_dominant_extrusion_candidates(self.solid_surface_signatures())
    }

    #[must_use = "The computed solid revolutions collection is returned and must be retained or inspected."]
    pub fn solid_revolutions(&self) -> Vec<solid_revolutions::RecoveredSolidRevolution> {
        self.indexed_sections()
            .flat_map(|(entities, index, cache)| {
                solid_revolutions::detect_solid_revolutions_with_tolerances(
                    entities,
                    index,
                    cache.source_tolerances(entities),
                    &crate::brep::manifold_solid_ids(entities),
                )
            })
            .collect()
    }

    #[must_use = "The computed radial slot revolutions collection is returned and must be retained or inspected."]
    pub fn radial_slot_revolutions(&self) -> Vec<solid_revolutions::RecoveredRadialSlotRevolution> {
        self.indexed_sections()
            .flat_map(|(entities, index, cache)| {
                solid_revolutions::detect_radial_slot_revolutions_with_tolerances(
                    entities,
                    index,
                    cache.source_tolerances(entities),
                    &crate::brep::manifold_solid_ids(entities),
                )
            })
            .collect()
    }

    #[must_use = "The computed solid extrusions collection is returned and must be retained or inspected."]
    pub fn solid_extrusions(&self) -> Vec<solid_extrusions::RecoveredSolidExtrusion> {
        self.indexed_sections()
            .flat_map(|(entities, index, _)| {
                solid_extrusions::detect_solid_extrusions_with_index(
                    entities,
                    index,
                    &crate::brep::manifold_solid_ids(entities),
                )
            })
            .collect()
    }

    #[must_use = "The computed periodic chains collection is returned and must be retained or inspected."]
    pub fn periodic_chains(&self) -> Vec<periodic_chains::PeriodicChainPattern> {
        self.indexed_sections()
            .flat_map(|(entities, index, _)| {
                periodic_chains::detect_periodic_chains_with_index(entities, index)
            })
            .collect()
    }

    /// Build complete semantic CAD IR directly from the already parsed exchange.
    ///
    /// Unlike explicit Compact STEP cleanup, this analysis does not rewrite,
    /// serialize, or reparse the input solely to measure rewrite idempotence.
    ///
    /// # Errors
    /// Returns an error if detector output cannot be converted into exact CAD IR.
    pub fn complete_ir(&self) -> Result<complete_ir::CompleteIr> {
        // An already-requested signature census can also supply Complete IR's
        // compact per-source summaries. Do not initialize the full signature
        // tree for a cold Complete IR request: that would retain heavyweight
        // diagnostics that are only needed transiently in that call.
        let signatures = self.surface_signatures.get();
        let context = self
            .indexed_sections()
            .enumerate()
            .map(
                |(data_section, (entities, index, cache))| complete_ir::SectionEvidence {
                    index,
                    source_tolerances: cache.source_tolerances(entities),
                    surface_signatures: signatures.map(|cached| {
                        &cached.signatures[cached.section_ranges[data_section].clone()]
                    }),
                },
            )
            .collect::<Vec<_>>();
        complete_ir::recover_complete_ir_exchange_with_context(&self.exchange, &context)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use serde::Serialize;
    use std::hint::black_box;
    use std::time::Instant;

    fn json_equal<T: Serialize, U: Serialize>(cached: &T, standalone: &U) -> Result<()> {
        assert_eq!(
            serde_json::to_value(cached)?,
            serde_json::to_value(standalone)?
        );
        Ok(())
    }

    #[test]
    fn session_reuses_lazy_per_section_evidence_without_changing_detectors() -> Result<()> {
        let bytes = include_bytes!("../validation/fixtures/production_radial_slot.step");
        let session = AnalysisSession::from_bytes(bytes)?;
        let [section] = session.exchange.data.as_slice() else {
            anyhow::bail!("expected one DATA section");
        };
        let entities = &section.entities;
        let [cache] = session.sections.as_slice() else {
            anyhow::bail!("expected one cache");
        };
        assert!(cache.index.get().is_none());
        assert!(cache.source_tolerances.get().is_none());

        json_equal(
            &session.solid_surface_signatures(),
            &solid_revolutions::detect_solid_surface_signatures(entities),
        )?;
        let index = cache
            .index
            .get()
            .ok_or_else(|| anyhow::anyhow!("index missing"))?;
        assert!(cache.source_tolerances.get().is_none());
        json_equal(
            &session.formed_sheet_evidence(),
            &formed_sheet::detect_formed_sheet_evidence(entities),
        )?;
        json_equal(
            &session.solid_extrusions(),
            &solid_extrusions::detect_solid_extrusions(entities),
        )?;
        json_equal(
            &session.periodic_chains(),
            &periodic_chains::detect_periodic_chains(entities),
        )?;
        assert!(std::ptr::eq(
            index,
            cache
                .index
                .get()
                .ok_or_else(|| anyhow::anyhow!("index invalidated"))?,
        ));
        assert!(cache.source_tolerances.get().is_none());

        json_equal(
            &session.solid_revolutions(),
            &solid_revolutions::detect_solid_revolutions(entities),
        )?;
        let tolerances = cache
            .source_tolerances
            .get()
            .ok_or_else(|| anyhow::anyhow!("tolerances missing"))?;
        assert!(!tolerances.is_empty());
        assert_eq!(
            tolerances,
            &crate::source_accuracy::source_tolerance_by_representation_item(entities, index),
            "session uncertainty evidence must match standalone source parsing"
        );
        json_equal(
            &session.radial_slot_revolutions(),
            &solid_revolutions::detect_radial_slot_revolutions(entities),
        )?;
        assert!(std::ptr::eq(
            tolerances,
            cache
                .source_tolerances
                .get()
                .ok_or_else(|| anyhow::anyhow!("tolerance map invalidated"))?,
        ));
        assert!(std::ptr::eq(
            index,
            cache
                .index
                .get()
                .ok_or_else(|| anyhow::anyhow!("index invalidated"))?,
        ));
        Ok(())
    }

    #[test]
    fn duplicate_entity_ids_in_distinct_data_sections_keep_independent_indexes() -> Result<()> {
        let bytes = include_bytes!("../validation/fixtures/native_conical_frustum.step");
        let mut exchange = monstertruck_io::step::load::ParsedStepSource::parse(bytes)?.exchange;
        let [section] = exchange.data.as_slice() else {
            anyhow::bail!("expected one original DATA section");
        };
        let duplicate = section.clone();
        exchange.data.push(duplicate);
        let session = AnalysisSession::from_exchange(exchange);
        assert!(
            session
                .sections
                .iter()
                .all(|cache| cache.index.get().is_none())
        );
        let signatures = session.solid_surface_signatures();
        assert_eq!(signatures.len(), 2);
        assert_eq!(signatures[0], signatures[1]);
        let [first, second] = session.sections.as_slice() else {
            anyhow::bail!("expected exactly two independent caches");
        };
        let first_index = first
            .index
            .get()
            .ok_or_else(|| anyhow::anyhow!("first index missing"))?;
        let second_index = second
            .index
            .get()
            .ok_or_else(|| anyhow::anyhow!("second index missing"))?;
        assert_eq!(first_index, second_index);
        assert!(!std::ptr::eq(first_index, second_index));
        let revolutions = session.solid_revolutions();
        assert_eq!(revolutions.len(), 2);
        assert_eq!(revolutions[0], revolutions[1]);
        let first_tolerances = first
            .source_tolerances
            .get()
            .ok_or_else(|| anyhow::anyhow!("first uncertainty missing"))?;
        let second_tolerances = second
            .source_tolerances
            .get()
            .ok_or_else(|| anyhow::anyhow!("second uncertainty missing"))?;
        assert_eq!(first_tolerances, second_tolerances);
        assert!(!std::ptr::eq(first_tolerances, second_tolerances));
        // Cross-section source identity must survive full CAD report generation.
        let report = session.complete_ir()?;
        assert_eq!(report.body_count, 2);
        assert_eq!(session.sections.len(), 2);
        Ok(())
    }

    #[test]
    fn surface_signatures_are_borrowed_and_shared_with_dominant_extrusions() -> Result<()> {
        let bytes = include_bytes!("../validation/fixtures/native_split_ring_hollow_cylinder.step");
        let session = AnalysisSession::from_bytes(bytes)?;
        assert!(session.surface_signatures.get().is_none());
        // Other detectors should not eagerly retain the surface signature tree.
        let _ = session.formed_sheet_evidence();
        assert!(session.surface_signatures.get().is_none());

        let original =
            solid_revolutions::detect_solid_surface_signatures(&session.exchange.data[0].entities);
        let first = session.solid_surface_signatures();
        assert_eq!(first, original.as_slice());
        assert!(!first.is_empty());
        assert!(std::ptr::eq(first, session.solid_surface_signatures()));

        let expected = dominant_extrusions::detect_dominant_extrusion_candidates(&original);
        assert_eq!(session.dominant_extrusion_candidates(), expected);
        assert!(std::ptr::eq(first, session.solid_surface_signatures()));

        // The reverse call order primes the same cache exactly once.
        let reverse = AnalysisSession::from_bytes(bytes)?;
        assert!(reverse.surface_signatures.get().is_none());
        assert_eq!(reverse.dominant_extrusion_candidates(), expected);
        assert_eq!(reverse.solid_surface_signatures(), original.as_slice());
        Ok(())
    }

    #[test]
    fn surface_signature_cache_is_shared_across_readers_and_data_sections() -> Result<()> {
        let bytes = include_bytes!("../validation/fixtures/native_conical_frustum.step");
        let mut exchange = monstertruck_io::step::load::ParsedStepSource::parse(bytes)?.exchange;
        let [section] = exchange.data.as_slice() else {
            anyhow::bail!("expected exactly one source DATA section");
        };
        exchange.data.push(section.clone());
        let session = AnalysisSession::from_exchange(exchange);
        assert!(session.surface_signatures.get().is_none());

        let signatures = std::thread::scope(|scope| {
            let handles = (0..4)
                .map(|_| scope.spawn(|| session.solid_surface_signatures()))
                .collect::<Vec<_>>();
            let mut results = Vec::with_capacity(handles.len());
            for handle in handles {
                results.push(
                    handle
                        .join()
                        .map_err(|_| anyhow::anyhow!("signature reader panicked"))?,
                );
            }
            Ok::<_, anyhow::Error>(results)
        })?;
        assert!(
            signatures
                .windows(2)
                .all(|pair| std::ptr::eq(pair[0], pair[1]))
        );
        let signatures = session.solid_surface_signatures();
        assert_eq!(signatures.len(), 2);
        assert_eq!(signatures[0], signatures[1]);
        assert!(std::ptr::eq(signatures, session.solid_surface_signatures()));
        assert_eq!(session.sections.len(), 2);
        // Equal IDs must not merge distinct DATA sections or reorder results.
        let standalone = session
            .exchange
            .data
            .iter()
            .flat_map(|section| {
                solid_revolutions::detect_solid_surface_signatures(&section.entities)
            })
            .collect::<Vec<_>>();
        assert_eq!(signatures, standalone.as_slice());
        Ok(())
    }

    #[test]
    fn complete_ir_borrows_preexisting_signatures_but_does_not_eagerly_cache_them() -> Result<()> {
        for bytes in [
            include_bytes!("../validation/fixtures/production_radial_slot.step").as_slice(),
            include_bytes!("../validation/fixtures/native_conical_frustum.step").as_slice(),
        ] {
            let cold = AnalysisSession::from_bytes(bytes)?;
            assert!(cold.surface_signatures.get().is_none());
            let cold_report = cold.complete_ir()?;
            assert!(cold.surface_signatures.get().is_none());
            let warm = AnalysisSession::from_bytes(bytes)?;
            let signatures = warm.solid_surface_signatures();
            let cached = warm
                .surface_signatures
                .get()
                .ok_or_else(|| anyhow::anyhow!("missing expected surface cache"))?;
            assert_eq!(cached.section_ranges.len(), 1);
            assert_eq!(cached.section_ranges[0], 0..signatures.len());
            let warm_report = warm.complete_ir()?;
            json_equal(&warm_report, &cold_report)?;
            assert!(std::ptr::eq(signatures, warm.solid_surface_signatures()));
            let warm_report_again = warm.complete_ir()?;
            json_equal(&warm_report_again, &cold_report)?;
        }
        Ok(())
    }

    #[test]
    fn complete_ir_reuses_exact_per_data_ranges_including_skipped_sections() -> Result<()> {
        let bytes = include_bytes!("../validation/fixtures/native_conical_frustum.step");
        let mut exchange = monstertruck_io::step::load::ParsedStepSource::parse(bytes)?.exchange;
        let [valid] = exchange.data.as_slice() else {
            anyhow::bail!("expected one original DATA section");
        };
        let mut empty = valid.clone();
        empty.entities.clear();
        exchange.data = vec![empty.clone(), valid.clone(), empty, valid.clone()];
        let cold = AnalysisSession::from_exchange(exchange.clone());
        let report_without_cache = cold.complete_ir()?;
        assert!(cold.surface_signatures.get().is_none());

        let warm = AnalysisSession::from_exchange(exchange);
        let signatures = warm.solid_surface_signatures();
        assert_eq!(signatures.len(), 2);
        assert_eq!(signatures[0], signatures[1]);
        let cache = warm
            .surface_signatures
            .get()
            .ok_or_else(|| anyhow::anyhow!("surface cache disappeared"))?;
        assert_eq!(cache.section_ranges, [0..0, 0..1, 1..1, 1..2]);
        let report = warm.complete_ir()?;
        json_equal(&report, &report_without_cache)?;
        assert_eq!(report.surface_signature_summaries.len(), 2);
        assert_eq!(
            report
                .surface_signature_summaries
                .iter()
                .map(|summary| summary.source.data_section)
                .collect::<Vec<_>>(),
            [1, 3]
        );
        assert!(std::ptr::eq(signatures, warm.solid_surface_signatures()));
        Ok(())
    }

    #[test]
    fn signature_section_ranges_follow_successful_detection_not_source_solid_counts() -> Result<()>
    {
        let bytes = include_bytes!("../validation/fixtures/native_conical_frustum.step");
        let mut exchange = monstertruck_io::step::load::ParsedStepSource::parse(bytes)?.exchange;
        let [valid] = exchange.data.as_slice() else {
            anyhow::bail!("expected one original DATA section");
        };
        let mut incomplete = valid.clone();
        // Keep the source solid record, but remove its face evidence. A detector
        // must skip this incomplete solid rather than reserve a fake cache slot.
        incomplete.entities.retain(|entity| {
            crate::step_graph::simple_record(entity)
                .is_none_or(|record| record.name != "ADVANCED_FACE")
        });
        assert_eq!(
            solid_revolutions::detect_solid_surface_signatures(&incomplete.entities).len(),
            0
        );
        exchange.data = vec![incomplete, valid.clone(), valid.clone()];
        let session = AnalysisSession::from_exchange(exchange);
        let signatures = session.solid_surface_signatures();
        let cache = session
            .surface_signatures
            .get()
            .ok_or_else(|| anyhow::anyhow!("surface cache disappeared"))?;
        assert_eq!(cache.section_ranges, [0..0, 0..1, 1..2]);
        assert_eq!(signatures.len(), 2);
        assert_eq!(signatures[0], signatures[1]);
        Ok(())
    }

    /// Manual end-to-end diagnostic: compare a cold Complete IR session with
    /// a session whose surface-signature tree was previously requested.
    /// The benchmark excludes STEP parsing and does not assert elapsed time.
    #[test]
    #[ignore = "manual Complete IR signature reuse benchmark"]
    fn complete_ir_signature_reuse_benchmark() -> Result<()> {
        let bytes = include_bytes!("../validation/fixtures/production_radial_slot.step");
        let cold = AnalysisSession::from_bytes(bytes)?;
        let warm = AnalysisSession::from_bytes(bytes)?;
        let start = Instant::now();
        let precomputed = warm.solid_surface_signatures();
        let cache_initialization = start.elapsed();
        let serialized_bytes = serde_json::to_vec(precomputed)?.len();
        let runs = 64;
        let mut cold_samples = Vec::with_capacity(runs);
        let mut cached_samples = Vec::with_capacity(runs);
        // Alternate execution order so unrelated CPU load doesn't skew
        // two long consecutive blocks differently.
        for iteration in 0..runs {
            if iteration % 2 == 0 {
                let start = Instant::now();
                black_box(cold.complete_ir()?);
                cold_samples.push(start.elapsed());
                let start = Instant::now();
                black_box(warm.complete_ir()?);
                cached_samples.push(start.elapsed());
            } else {
                let start = Instant::now();
                black_box(warm.complete_ir()?);
                cached_samples.push(start.elapsed());
                let start = Instant::now();
                black_box(cold.complete_ir()?);
                cold_samples.push(start.elapsed());
            }
        }
        assert_eq!(
            serde_json::to_value(cold.complete_ir()?)?,
            serde_json::to_value(warm.complete_ir()?)?
        );
        cold_samples.sort_unstable();
        cached_samples.sort_unstable();
        let cold_median = cold_samples[runs / 2];
        let cached_median = cached_samples[runs / 2];
        let cold_total = cold_samples.iter().sum::<std::time::Duration>();
        let cached_total = cached_samples.iter().sum::<std::time::Duration>();
        println!(
            "Complete IR paired signature evidence: runs={runs}, payload_bytes={serialized_bytes}, signatures={}, cache_initialize={cache_initialization:?}, cold_total={cold_total:?}, cached_total={cached_total:?}, cold_median={cold_median:?}, cached_median={cached_median:?}, median_ratio={:.3}x",
            precomputed.len(),
            cold_median.as_secs_f64() / cached_median.as_secs_f64()
        );
        Ok(())
    }

    /// Profile the surface-signature census, a deep owned clone, and borrowed
    /// dominant-extrusion analysis on an actual STEP file. Intentionally not a
    /// timing assertion; compare this before and after session caching.
    #[test]
    #[ignore = "manual profile of signature retention and borrowing"]
    fn surface_signature_cost_benchmark() -> Result<()> {
        let session = AnalysisSession::from_bytes(include_bytes!(
            "../validation/fixtures/native_split_ring_hollow_cylinder.step"
        ))?;
        let (entities, index, _) = session
            .indexed_sections()
            .next()
            .ok_or_else(|| anyhow::anyhow!("no DATA section"))?;
        let original =
            solid_revolutions::detect_solid_surface_signatures_with_index(entities, index);
        let faces = original
            .iter()
            .map(|solid| solid.faces.len())
            .sum::<usize>();
        let serialized_bytes = serde_json::to_vec(&original)?.len();
        let runs = 36;
        let start = Instant::now();
        for _ in 0..runs {
            black_box(
                solid_revolutions::detect_solid_surface_signatures_with_index(entities, index),
            );
        }
        let recompute = start.elapsed();
        let start = Instant::now();
        for _ in 0..runs {
            black_box(original.clone());
        }
        let deep_clone = start.elapsed();
        let start = Instant::now();
        for _ in 0..runs {
            black_box(dominant_extrusions::detect_dominant_extrusion_candidates(
                black_box(&original),
            ));
        }
        let dominant_only = start.elapsed();
        let start = Instant::now();
        let from_cache = session.solid_surface_signatures();
        let cache_build = start.elapsed();
        assert_eq!(from_cache, original.as_slice());
        let start = Instant::now();
        for _ in 0..runs {
            black_box(session.solid_surface_signatures());
            black_box(session.dominant_extrusion_candidates());
        }
        let borrowed_reuse = start.elapsed();
        println!(
            "signature census: solids={}, faces={}, serialized_bytes={}, runs={}, recompute={:?}, deep_clone={:?}, dominant_only={:?}, cache_build={:?}, borrowed_signature_and_dominant={:?}",
            original.len(),
            faces,
            serialized_bytes,
            runs,
            recompute,
            deep_clone,
            dominant_only,
            cache_build,
            borrowed_reuse
        );
        Ok(())
    }

    /// Run manually to compare the exact indexing/context work eliminated by
    /// multi-detector sessions. The wall time is diagnostic, not an assertion.
    #[test]
    #[ignore = "manual microbenchmark"]
    fn index_and_context_reuse_benchmark() -> Result<()> {
        let session = AnalysisSession::from_bytes(include_bytes!(
            "../validation/fixtures/native_split_ring_hollow_cylinder.step"
        ))?;
        let [section] = session.exchange.data.as_slice() else {
            anyhow::bail!("expected one DATA section");
        };
        let entities = section.entities.as_slice();
        let [cache] = session.sections.as_slice() else {
            anyhow::bail!("expected one cache");
        };
        let iterations = 300;
        let start = Instant::now();
        for _ in 0..iterations {
            let index = black_box(crate::step_graph::build_index(entities));
            black_box(
                crate::source_accuracy::source_tolerance_by_representation_item(entities, &index),
            );
        }
        let rebuilt = start.elapsed();
        let start = Instant::now();
        for _ in 0..iterations {
            black_box(cache.index(entities));
            black_box(cache.source_tolerances(entities));
        }
        let reused = start.elapsed();
        println!(
            "cached index+source tolerance: {iterations} lookups, rebuilt={rebuilt:?}, cached={reused:?}, ratio={:.1}x",
            rebuilt.as_secs_f64() / reused.as_secs_f64()
        );

        // Real mixed-detector workload: this includes geometry recognition,
        // so the difference is not inflated by measuring only index creation.
        let session = AnalysisSession::from_bytes(include_bytes!(
            "../validation/fixtures/production_radial_slot.step"
        ))?;
        let entities = session.exchange.data[0].entities.as_slice();
        let iterations = 96;
        let start = Instant::now();
        for _ in 0..iterations {
            black_box(formed_sheet::detect_formed_sheet_evidence(entities));
            black_box(solid_revolutions::detect_solid_surface_signatures(entities));
            black_box(solid_revolutions::detect_solid_revolutions(entities));
            black_box(solid_revolutions::detect_radial_slot_revolutions(entities));
            black_box(solid_extrusions::detect_solid_extrusions(entities));
            black_box(periodic_chains::detect_periodic_chains(entities));
        }
        let standalone = start.elapsed();
        let start = Instant::now();
        for _ in 0..iterations {
            black_box(session.formed_sheet_evidence());
            black_box(session.solid_surface_signatures());
            black_box(session.solid_revolutions());
            black_box(session.radial_slot_revolutions());
            black_box(session.solid_extrusions());
            black_box(session.periodic_chains());
        }
        let cached_detectors = start.elapsed();
        println!(
            "six full detectors x {iterations}: standalone={standalone:?}, cached={cached_detectors:?}, ratio={:.2}x",
            standalone.as_secs_f64() / cached_detectors.as_secs_f64()
        );
        Ok(())
    }
}

# Deduplicated body corpus (v3)

The whole-model manifest remains the authority for array and assembly recovery.
Use a body corpus to iterate on feature recovery without repeatedly parsing
complete assemblies. The Rust commands operate directly on the same scanner,
canonicalization, and extraction implementation as the standalone body tool.

Build `step-body-tool` from the workspace:

```sh
cargo build -p step-redox --bin step-body-tool --release
```

With a UTF-8 TSV containing columns `uuid` and `path` (additional columns
are allowed, relative paths are resolved against the manifest's directory):

```sh
step-body-tool corpus build --manifest sources.tsv --out corpus-v3 \
  --cache scan-cache --jobs 4
step-body-tool corpus validate --corpus corpus-v3 --jobs 4
step-body-tool corpus census --manifest sources.tsv --out supports.jsonl \
  --face-count 6 --contains-support plane --jobs 4
step-body-tool corpus find-gap --manifest sources.tsv --start 0 --jobs 4
```

Use `--no-extract` for metadata-only builds and `--require-canonical` to
reject fallback identities or skipped solids. Census supports repeated
`--contains-support` filters (all required). Gap search returns the first
failing source in manifest order, with a zero-based row index.

## Identity and cache

Version 3 deliberately replaces the Python manifest format. It does not
modify or read older corpus outputs. Canonical bodies with known length
scales group by canonical geometry fingerprint and exact unit-scale bits.
Unknown-scale identities are source-local; unsupported bodies are identified
by source UUID and DATA-section/solid ID and are never geometrically merged.
Representatives prefer extractable bodies, then smaller closures and stable
source/entity ordering. Unit information comes from the production scanner,
not whole-file regular expressions.

The reusable scan cache is separate from the corpus. Its namespace contains
the exact scanner executable SHA-256 and its entries are keyed by source
content SHA-256. Concurrent equal-content reads share a parse within one
process. Cache entries are atomically replaced and checked for structural
invariants on load. Changing input bytes or the scanner invalidates old
results. The cache is disposable; original source STEP bytes are authoritative.

## Publication and validation

A build publishes a new directory (never overwrites an existing corpus).
The staging directory resides on the same filesystem so publication is one
rename. Reports record all observed source and extraction errors; a completed
but failing build publishes `FAILED`, whereas successful builds publish
`COMPLETE`.

- `summary.json`: schema, executable/manifest digests, counters and result.
- `model-scans.jsonl`: ordered, content-identified scan results or failures.
- `occurrences.jsonl`: individual source/body provenance records.
- `bodies.jsonl`: representatives and occurrence counts.
- `extraction-reports.json`: materialized object byte digests and failures.
- `objects/<prefix>/<fingerprint>.step`: one object per extractable identity
  (absent for metadata-only builds).

Every isolated STEP body is published as a new standalone product with its
own standard product-definition graph referencing its actual isolated shape
representation. This matters for input assemblies with orphan instance
representations: a geometrically complete STEP entity closure can parse
successfully while OCCT still transfers **zero** product roots. The validator
requires a root-linked product definition for every materialized object and
rejects older unrooted output. Regenerate the affected new corpus rather than
modifying it in place; the source-content scan cache remains reusable through
its scanner-executable revision namespace.

The validator rebuilds the identity/representative inventory from saved scans,
checks its correspondence with occurrence and body manifests, checks markers
and counts, and rescans materialized objects after verifying their exact byte
digests. Source-local fallback extraction is checked for topology and
provenance rather than claimed as canonical geometry equivalence. These checks
are corpus integrity checks, not independent geometric recovery certificates.

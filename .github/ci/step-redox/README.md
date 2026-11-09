# STEP CAD corpus and OpenCascade CI tooling

This directory is CI-only Python tooling, not a geometry implementation or
runtime dependency. It uses uv with the local pyproject.toml and uv.lock.

- harness.py and semantic_probe.py validate the committed STEP fixtures.
- occt_* utilities provide independent OpenCascade geometry comparisons and
  corpus diagnostics. They do not participate in constructive recovery.
- Core STEP parsing, feature detection, geometry and topology belong to Rust.

From the repository root, call uv run --frozen --project .github/ci/step-redox
python .github/ci/step-redox/harness.py, passing fixture manifest, executable
paths and report directories. Run full corpus validation on release/manual CI.
Do not relax repository-wide strict Python policy to accommodate this tooling.

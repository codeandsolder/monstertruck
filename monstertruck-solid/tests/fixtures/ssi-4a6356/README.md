# 4a6356 constructive STEP — upper-round Boolean operands

These are native `monstertruck_topology::Solid` JSON-compressed B-reps,
serialized from the existing source-derived six-stage CAD replay immediately
after its fourth (void) stage, and from its analytic upper-round cutter.
Both were evaluated under the default numerical policy (not an experimental
mesh sampling override); the kernel modeling budget was 0.00001 mm while the
classic Boolean numerical length was independently 0.000001 mm.

They capture a real regression: two exact cylinder-cylinder intersection
coedges had matching topological endpoints and interior geometry but
parameterized SSI leaders whose endpoint inverse queries were inaccurate.
The generated Boolean output had 47 faces but two unpaired two-edge loops.

This fixture tests operation correctness without recomputing earlier CAD
stages or roundtripping through STEP strings. The complete source-derived
CAD replay remains a separate end-to-end test in `step-redox`.

# Source-derived production radial-slot exact B-rep operands

`host.json` and `cutter.json` are native `monstertruck_modeling::Solid`
serde JSON snapshots, obtained directly from the two child nodes of the
production source-recovery model in `step-redox` before its Boolean difference.
The authoritative input STEP fixture is
`step-redox/validation/fixtures/production_radial_slot.step`.

They capture a general two-stage numerical robustness regression:

- A cylinder/plane SSI leader slightly overruns its original finite trim. The
  corrected exact endpoint is valid on both carrier surfaces, but the sampled
  polygon boundary is ~0.00027 mm away (normal 0.001 mm sampling radius).
  Discovering a sampled polyline edge at the exact numerical tolerance fails.
- The slot cutter extends 0.0001 mm beyond the original bore. Reattaching an
  already certified source-trim endpoint to the cutter cap would move that
  vertex off the original exact source surface trim, breaking downstream face
  subdivision. Both original-trim endpoint provenance and numerical movement
  must remain certified independently of mesh candidate discovery.
- Two semicylindrical source faces cross the cutter's upper cap; their SSI
  arcs must partition that cap into distinct material and nonmaterial regions.

No experimental tolerance overrides are used. The fast tests use these saved
native B-reps rather than rerunning source STEP parsing or creating a second
copy of CAD geometry logic. The existing production STEP end-to-end test still
exercises the original parser, recovery model, Boolean and STEP writer.

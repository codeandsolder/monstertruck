#[cfg(feature = "cad-kernel-monstertruck")]
use anyhow::{Context, Result, bail};
#[cfg(feature = "cad-kernel-monstertruck")]
#[path = "support/extrusion_core_args.rs"]
mod extrusion_core_args;

#[cfg(feature = "cad-kernel-monstertruck")]
fn main() -> Result<()> {
    use step_redox::cad_kernel::{CadKernel, monstertruck::MonstertruckKernel};

    let usage =
        "usage: step-extrusion-core-healed-rebuild INPUT.step SOLID_ID CAP_A CAP_B OUTPUT.step";
    let args = extrusion_core_args::ExtrusionCoreArgs::parse(usage)?;
    let output = &args.output;
    let bytes = args.read_input()?;
    let candidate = args.recover(
        &bytes,
        "treatment-healed extrusion core",
        step_redox::recover_treatment_healed_extrusion_core_from_cap_faces_bytes,
    )?;
    let fragment = step_redox::cad_recovery::recover_solid_extrusion_fragment(&candidate.core)?;

    let kernel = MonstertruckKernel::project_standard();
    let evaluated = kernel.evaluate(&fragment.model, fragment.root)?;
    let summary = kernel.summarize(&evaluated);
    if !summary.geometrically_consistent || summary.shells != 1 {
        bail!("Monstertruck rebuilt an invalid treatment-healed extrusion host: {summary:?}");
    }
    std::fs::write(output, kernel.to_step(&evaluated)?)
        .with_context(|| format!("write {}", output.display()))?;
    Ok(())
}

#[cfg(not(feature = "cad-kernel-monstertruck"))]
fn main() -> anyhow::Result<()> {
    anyhow::bail!("step-extrusion-core-healed-rebuild requires --features cad-kernel-monstertruck")
}

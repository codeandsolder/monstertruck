use anyhow::{Context, Result, ensure};
use std::{path::PathBuf, time::Instant};
use step_redox::{cad_kernel::{CadKernel,monstertruck::MonstertruckKernel}, complete_ir::recover_complete_ir_with_kernel_bytes, compare_analytic_solid_boundaries_bytes, detect_solid_surface_signatures_bytes};
fn main()->Result<()> {
    let mut args=std::env::args_os().skip(1).map(PathBuf::from);
    let input=args.next().context("input")?;
    let output=args.next().context("output directory")?;
    let name=input.file_name().context("input basename")?.to_string_lossy().chars().take(8).collect::<String>();
    let now=Instant::now();
    let bytes=std::fs::read(&input)?;
    let kernel=MonstertruckKernel::project_standard();
    let report=recover_complete_ir_with_kernel_bytes(&bytes,&kernel)?;
    std::fs::write(output.join(format!("{name}.json")),serde_json::to_vec(&report)?)?;
    println!("{name} constructive {} strict {} recovery {:?}",report.constructively_recovered_solid_count,report.source_boundary_verified_solids.len(),now.elapsed());
    if report.constructively_recovered_solid_count!=1 {return Ok(());}
    let fragments=report.selected_fragments().collect::<Vec<_>>();
    let [selected]=fragments.as_slice() else {anyhow::bail!("expected one fragment")};
    let evaluated=kernel.evaluate(&selected.fragment.model,selected.fragment.root)?;
    let summary=kernel.summarize(&evaluated);
    ensure!(summary.is_single_solid(),"{summary:?}");
    let rebuilt=kernel.to_step(&evaluated)?;
    std::fs::write(output.join(format!("{name}.rebuilt.step")),rebuilt.as_bytes())?;
    let signatures=detect_solid_surface_signatures_bytes(rebuilt.as_bytes())?;
    let [sig]=signatures.as_slice() else {anyhow::bail!("expected single exported solid")};
    let comparison=compare_analytic_solid_boundaries_bytes(&bytes,report.body_scan.solids[0].source.id,rebuilt.as_bytes(),sig.solid_id)?;
    println!("{name} faces {} nodes {} independent equal {} supported {}/{} reasons {:?}/{:?}",summary.faces,selected.fragment.model.nodes.len(),comparison.equivalent,comparison.first.supported,comparison.second.supported,comparison.first.reason,comparison.second.reason);
    std::fs::write(output.join(format!("{name}.compare.json")),serde_json::to_vec(&comparison)?)?;
    Ok(())
}

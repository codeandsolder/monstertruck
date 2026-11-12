//! Opt-in CUDA kernel compilation. Ordinary portable builds never invoke nvcc.
use std::{env, path::PathBuf, process::Command};

fn main() {
    println!("cargo:rerun-if-changed=cuda/seed_kernel.cu");
    if env::var_os("CARGO_FEATURE_CUDA_PRESEARCH").is_none() {
        return;
    }
    let nvcc = env::var_os("NVCC")
        .map(PathBuf::from)
        .unwrap_or_else(|| PathBuf::from("/opt/cuda/bin/nvcc"));
    let out = PathBuf::from(env::var_os("OUT_DIR").expect("Cargo supplies OUT_DIR"));
    let object = out.join("seed_kernel.o");
    let archive = out.join("libmonstertruck_seed_cuda.a");
    let status = Command::new(nvcc)
        .args([
            "-std=c++17",
            "-O3",
            "--fmad=false",
            "-arch=sm_75",
            "-Xcompiler=-fPIC",
            "-c",
        ])
        .arg("cuda/seed_kernel.cu")
        .arg("-o")
        .arg(&object)
        .status()
        .expect("NVCC required for opt-in cuda-presearch feature");
    assert!(status.success(), "NVCC failed compiling seed_kernel.cu");
    let status = Command::new("ar")
        .arg("crs")
        .arg(&archive)
        .arg(&object)
        .status()
        .expect("ar required for CUDA static archive");
    assert!(status.success(), "ar failed creating CUDA static archive");
    println!("cargo:rustc-link-search=native={}", out.display());
    println!("cargo:rustc-link-lib=static=monstertruck_seed_cuda");
    println!("cargo:rustc-link-search=native=/opt/cuda/lib64");
    println!("cargo:rustc-link-lib=dylib=cudart");
    println!("cargo:rustc-link-lib=dylib=stdc++");
}

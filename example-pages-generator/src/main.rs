#![cfg_attr(not(debug_assertions), deny(warnings))]
#![deny(clippy::all, rust_2018_idioms)]

use std::fmt::Write as _;
use std::io::{self, BufRead, BufReader};
use std::process::{Child, Command};
use std::thread::JoinHandle;

const EXAMPLES: &[&str] = &[
    "wgsl-sandbox",
    "bsp-animation",
    "collision-sphere",
    "material-samples",
    "rotate-objects",
    "simple-obj-viewer",
    "simple-shape-viewer",
    "textured-cube",
];

type OutputThread = JoinHandle<io::Result<()>>;

fn output_threads(child: &mut Child) -> io::Result<[OutputThread; 2]> {
    let stdout = child
        .stdout
        .take()
        .ok_or_else(|| io::Error::other("child stdout is not piped"))?;
    let stderr = child
        .stderr
        .take()
        .ok_or_else(|| io::Error::other("child stderr is not piped"))?;

    let stdout = std::thread::spawn(move || {
        for line in BufReader::new(stdout).lines() {
            println!("{}", line?);
        }
        Ok(())
    });
    let stderr = std::thread::spawn(move || {
        for line in BufReader::new(stderr).lines() {
            eprintln!("{}", line?);
        }
        Ok(())
    });
    Ok([stdout, stderr])
}

fn wait_for_child(
    child: &mut Child,
    threads: [OutputThread; 2],
    description: &str,
) -> io::Result<()> {
    let status = child.wait()?;
    for thread in threads {
        thread
            .join()
            .map_err(|_| io::Error::other(format!("{description} output thread panicked")))??;
    }
    if status.success() {
        Ok(())
    } else {
        Err(io::Error::other(format!(
            "{description} exited with status {status}"
        )))
    }
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut child = Command::new("cargo")
        .args([
            "build",
            "--target",
            "wasm32-unknown-unknown",
            "--examples",
            "--release",
            "-p",
            "monstertruck-gpu",
            "-p",
            "monstertruck-render",
        ])
        .stdout(std::process::Stdio::piped())
        .stderr(std::process::Stdio::piped())
        .spawn()?;
    let threads = output_threads(&mut child)?;
    wait_for_child(&mut child, threads, "cargo build")?;

    let mut sum = String::new();
    for dir in EXAMPLES {
        let output_dir = format!("dist/{dir}");
        std::fs::create_dir_all(&output_dir)?;
        let wasm_path = format!("target/wasm32-unknown-unknown/release/examples/{dir}.wasm");
        let mut child = Command::new("wasm-bindgen")
            .args(["--target", "web", "--out-dir", &output_dir, &wasm_path])
            .stdout(std::process::Stdio::piped())
            .stderr(std::process::Stdio::piped())
            .spawn()?;
        let threads = output_threads(&mut child)?;
        wait_for_child(&mut child, threads, &format!("wasm-bindgen for {dir}"))?;

        std::fs::write(
            format!("{output_dir}/index.html"),
            include_str!("example-index.html").replace("{example}", dir),
        )?;
        write!(&mut sum, "<li><a href=\"{dir}/index.html\">{dir}</a></li>")?;
    }
    std::fs::write(
        "dist/index.html",
        include_str!("index.html").replace("<!-- index -->", &sum),
    )?;
    Ok(())
}

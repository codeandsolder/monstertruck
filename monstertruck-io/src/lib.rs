//! Read and write CAD exchange formats as monstertruck B-rep.
//!
//! One crate, one feature per format. A caller reaches every format monstertruck
//! supports through a single dependency instead of one crate per format, and
//! compiles only the formats it asks for.
//!
//! # Formats
//!
//! | format | feature | default | state |
//! | --- | --- | --- | --- |
//! | ISO 10303-21 (STEP) | `step` | yes | read and write, ours |
//! | IGES 4.0 to 5.3 | `iges` | no | read only, via cadmpeg 0.6 |
//!
//! ```
//! use monstertruck_io::step::load::{*, step_geometry::*};
//!
//! let bytes = include_bytes!(concat!(
//!     env!("CARGO_MANIFEST_DIR"),
//!     "/../resources/step/occt-cube.step",
//! ));
//! let table = Table::from_step_bytes(bytes).unwrap();
//! let step_shell = table.shell.values().next().unwrap();
//! let compressed = table.to_compressed_shell(step_shell).unwrap();
//! ```
//!
//! # STEP
//!
//! [`step`] is monstertruck's native STEP reader and writer. Import through
//! [`step::load`] and export through [`step::save`].
//!
//! It deliberately does **not** go through cadmpeg. It is measured against a
//! corpus, routes analytic surfaces onto closed forms so booleans stay exact,
//! and reports what a file lost. The cadmpeg 0.6 bridge is a separate,
//! opt-in alternative used for comparative import testing, not our STEP writer.
//! Earlier cadmpeg losses on circular edges (cadmpeg/cadmpeg#79) were fixed
//! upstream (cadmpeg/cadmpeg#83).
//!
//! # IGES: the conversion works, the decoder is the limit
//!
//! Decoding is [`cadmpeg`](https://github.com/cadmpeg/cadmpeg)'s, which
//! monstertruck has no reason to reimplement. Turning the recovered geometry into
//! monstertruck's B-rep is ours, and it is written: see [`cadmpeg`] for what
//! reaches which surface variant, and why analytic carriers stay analytic.
//!
//! Measured against `cadmpeg_ir::examples::unit_cube` -- cadmpeg's own document,
//! not a hand-built one -- a solid body converts to a shell of 8 vertices, 12
//! edges and 6 faces that [`monstertruck_topology::Shell::extract`] accepts and
//! reports CLOSED.
//!
//! ## What still stops a real file
//!
//! **`cadmpeg-codec-iges` 0.6 decodes fixed-format ASCII IGES 4.0 through 5.3.**
//! It returns decoder losses with their original provenance, even when the
//! decoder succeeds. Unsupported body types remain named refusals.
//!
//! Carriers this crate refuses on its own are refused BY NAME
//! ([`Error::UnsupportedSurfaceKind`], [`Error::UnsupportedCurveKind`]): ellipses
//! and the other conics, composite curves, procedural carriers without exact
//! caches, pole loops, and coedge-local edges. Periodic NURBS are supported.
//!
//! ## Nothing is ever silently empty
//!
//! An empty document yields [`Error::NoGeometry`]. The normal IGES reader returns
//! [`cadmpeg::Import`] with convertible bodies, explicitly named refusals, and
//! decoder losses. Callers needing all-or-nothing behavior can request a strict
//! conversion, which fails on the first unsupported body.
//!
//! Adding a format is cheap because the cadmpeg codecs all decode into ONE
//! intermediate representation, `cadmpeg_ir::CadIr`: the per-format work is only
//! which decoder to call, and the part that carries the risk -- mapping recovered
//! geometry onto [`monstertruck_topology`] -- is shared, in [`cadmpeg`].

#![cfg_attr(not(debug_assertions), deny(warnings))]
#![deny(clippy::all, rust_2018_idioms)]

mod error;

pub use error::{Error, Result};

#[cfg(feature = "cadmpeg")]
pub mod cadmpeg;
#[cfg(feature = "iges")]
pub mod iges;
#[cfg(feature = "step")]
pub mod step;

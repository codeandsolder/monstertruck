//! Fillet operations for [`Shell`](monstertruck_topology::Shell) edges.
//!
//! Provides rolling-ball fillet operations: single-edge fillets,
//! fillets with side face updates, and fillets along open or closed wire chains.
//! The [`fillet_edges`] function provides a high-level API that automatically
//! resolves face adjacency from edge IDs.
//!
//! Fillets are POST-CSG and kernel-independent -- nothing here references a
//! boolean backend, which is why the crate sits below both the published
//! `monstertruck-solid` marching kernel and any external SSI boolean backend
//! rather than inside either.

#![cfg_attr(not(debug_assertions), deny(warnings))]
#![deny(clippy::all, rust_2018_idioms)]

pub(crate) use ahash::AHashSet as HashSet;

#[expect(
    clippy::cast_precision_loss,
    reason = "fillet sampling and subdivision counts are bounded by iterable geometry and converted only for normalized f64 parameters"
)]
pub(crate) const fn count_as_f64(value: usize) -> f64 { value as f64 }

#[expect(
    clippy::cast_precision_loss,
    reason = "extended fillet sampling indices range only from -1 through a bounded subdivision count plus one"
)]
pub(crate) const fn signed_index_as_f64(value: isize) -> f64 { value as f64 }

mod edge_select;
mod ops;

mod convert;
mod error;
mod geometry;
mod params;
mod topology;
mod types;

#[cfg(test)]
mod tests;

pub use convert::{FilletIntersectionCurve, FilletableCurve, FilletableSurface};
pub use edge_select::{fillet_edges, fillet_edges_generic};
pub use error::FilletError;
pub use ops::{fillet, fillet_along_wire, fillet_with_side};
pub use params::{FilletOptions, FilletProfile, RadiusSpec};
pub use types::ParameterCurveLinear;

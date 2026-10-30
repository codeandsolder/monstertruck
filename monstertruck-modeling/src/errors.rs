use thiserror::Error;

/// Modeling errors.
#[derive(Debug, PartialEq, Eq, Error)]
pub enum Error {
    /// Wrapper of topological error.
    #[error(transparent)]
    FromTopology(#[from] monstertruck_topology::errors::Error),
    /// Tried to attach a plane to a wire that was not on one plane.
    /// cf. [`builder::try_attach_plane`](../builder/fn.try_attach_plane.html)
    #[error("cannot attach a plane to a wire that is not on one plane.")]
    WireNotInOnePlane,
    /// Tried to create homotopy for two wires with different numbers of edges.
    /// cf. [`builder::try_wire_homotopy`](../builder/fn.try_wire_homotopy.html)
    #[error("The wires must contain the same number of edges to create a homotopy.")]
    NotSameNumberOfEdges,
    /// One or more wires are not closed.
    #[error("all wires must be closed for profile construction.")]
    OpenWire,
    /// Ambiguous nesting: a loop is not clearly inside or outside another.
    #[error("ambiguous nesting between loops; loops may overlap or touch.")]
    AmbiguousNesting,
    /// No outer loop found among the provided wires.
    #[error("no outer loop found; at least one wire must have positive signed area.")]
    NoOuterLoop,
    /// Projected profile geometry contains non-finite coordinates or area.
    #[error("profile geometry contains non-finite projected values.")]
    NonFiniteProfileGeometry,
    /// A circular-arc tangent vector is zero or near zero, so it cannot define a direction.
    #[error("circular-arc tangent is degenerate (zero or near-zero).")]
    DegenerateCircularArcTangent,
    /// A circular-arc tangent is parallel (or anti-parallel) to the chord between the endpoints,
    /// so the plane of the arc is under-determined.
    #[error("circular-arc tangent is parallel to the chord between the arc endpoints.")]
    CircularArcTangentParallelToChord,
    /// Three points supplied for a circular arc are coincident, collinear, or non-finite.
    #[error("circular-arc points do not define a finite circle.")]
    DegenerateCircularArcPoints,
    /// A wire revolution requested a non-unit/non-finite axis, invalid angle, or invalid subdivision count.
    #[error("wire revolution parameters are invalid.")]
    InvalidRevolveParameters,
    /// The generic sweep produced topology that cannot be reduced to the axis-aware wire form.
    #[error("wire revolution produced unexpected intermediate topology.")]
    UnexpectedRevolveTopology,
    /// A composite sweep was requested without any path segments.
    #[error("composite sweep requires at least one segment.")]
    EmptyCompositeSweep,
    /// A composite sweep segment has zero/invalid motion or an invalid rotation subdivision.
    #[error("composite sweep contains a degenerate or invalid segment.")]
    InvalidCompositeSweepSegment,
    /// An internally generated sweep segment did not have the expected single-shell, capped topology.
    #[error("composite sweep produced unexpected intermediate topology.")]
    UnexpectedCompositeSweepTopology,
    /// A closed composite sweep does not return its final section to the starting section.
    #[error("closed composite sweep does not return to its starting section.")]
    CompositeSweepNotClosed,
}

#[test]
fn print_messages() {
    use std::io::Write;
    writeln!(
        &mut std::io::stderr(),
        "****** test of the expressions of error messages ******\n"
    )
    .unwrap();
    writeln!(
        &mut std::io::stderr(),
        "{}\n",
        Error::FromTopology(monstertruck_topology::errors::Error::SameVertex)
    )
    .unwrap();
    writeln!(&mut std::io::stderr(), "{}\n", Error::WireNotInOnePlane).unwrap();
    writeln!(&mut std::io::stderr(), "{}\n", Error::OpenWire).unwrap();
    writeln!(&mut std::io::stderr(), "{}\n", Error::AmbiguousNesting).unwrap();
    writeln!(&mut std::io::stderr(), "{}\n", Error::NoOuterLoop).unwrap();
    writeln!(
        &mut std::io::stderr(),
        "*******************************************************"
    )
    .unwrap();
}

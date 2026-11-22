"""Independent signed boundary integration with shared reference locations."""

import math
from dataclasses import dataclass

from OCP.Bnd import Bnd_Box
from OCP.BRep import BRep_Tool
from OCP.BRepBndLib import BRepBndLib
from OCP.BRepCheck import BRepCheck_Analyzer
from OCP.BRepGProp import BRepGProp_Domain, BRepGProp_Face, BRepGProp_VinertGK
from OCP.gp import gp_Pnt
from OCP.TopAbs import (
    TopAbs_FACE,
    TopAbs_FORWARD,
    TopAbs_REVERSED,
    TopAbs_SHELL,
    TopAbs_SOLID,
)
from OCP.TopExp import TopExp_Explorer
from OCP.TopoDS import TopoDS, TopoDS_Face, TopoDS_Shape


@dataclass(frozen=True)
class Integral:
    volume: float
    summed_face_error_estimate: float


@dataclass(frozen=True)
class ReferenceSample:
    reference: tuple[float, float, float]
    a: Integral
    b: Integral
    delta: float
    coarse_delta: float
    a_quadrature_drift: float
    b_quadrature_drift: float


@dataclass(frozen=True)
class PairedVolume:
    method: str
    quadrature_eps: tuple[float, float]
    absolute_tolerance: float
    max_absolute_delta: float
    max_quadrature_drift: float
    within_absolute_tolerance: bool
    reference_spread: dict[str, float]
    samples: list[ReferenceSample]


def closed_faces(shape: TopoDS_Shape) -> list[TopoDS_Face]:
    solids = TopExp_Explorer(shape, TopAbs_SOLID)
    if not solids.More() or not BRepCheck_Analyzer(shape).IsValid():
        raise ValueError("paired integration requires a valid solid")
    solid = solids.Current()
    solids.Next()
    if solids.More():
        raise ValueError("paired integration requires exactly one solid per operand")
    faces = []
    shells = TopExp_Explorer(solid, TopAbs_SHELL)
    while shells.More():
        shell = shells.Current()
        if not BRep_Tool.IsClosed_s(shell):
            raise ValueError("paired integration requires closed shells")
        explorer = TopExp_Explorer(shell, TopAbs_FACE)
        while explorer.More():
            face = TopoDS.Face(explorer.Current())
            if face.Orientation() in (TopAbs_FORWARD, TopAbs_REVERSED):
                if BRep_Tool.Surface_s(face) is None:
                    raise ValueError("face has no exact surface")
                faces.append(face)
            explorer.Next()
        shells.Next()
    if not faces:
        raise ValueError("solid has no oriented surface faces")
    return faces


def shared_references(
    a: TopoDS_Shape, b: TopoDS_Shape
) -> list[tuple[float, float, float]]:
    # One symmetric geometry-derived frame, independent of vertex multiplicity.
    # All four comparisons are required: never choose the best volume residual.
    box = Bnd_Box()
    for shape in (a, b):
        BRepBndLib.AddOptimal_s(shape, box, False, False)
    if box.IsVoid() or box.IsOpen():
        raise ValueError("paired integration requires finite bounds")
    xmin, ymin, zmin = box.CornerMin().Coord()
    xmax, ymax, zmax = box.CornerMax().Coord()
    cx, cy, cz = (
        (lo / 2 + hi / 2) for lo, hi in ((xmin, xmax), (ymin, ymax), (zmin, zmax))
    )
    return [(cx, cy, cz), (xmax, cy, cz), (cx, ymax, cz), (cx, cy, zmax)]


def integrate(
    faces: list[TopoDS_Face], reference: tuple[float, float, float], eps: float
) -> Integral:
    # Same adaptive Gauss-Kronrod, spline-span splitting, and signed face
    # orientations as OCCT's whole-solid method, but without choosing a new
    # location for each operand or accumulating in traversal-dependent order.
    integrator = BRepGProp_VinertGK()
    integrator.SetLocation(gp_Pnt(*reference))
    surface = BRepGProp_Face(True)
    domain = BRepGProp_Domain()
    masses = []
    errors = []
    for face in faces:
        surface.Load(face)
        if face.NbChildren():
            domain.Init(face)
            error = integrator.Perform(surface, domain, eps, False, False)
        else:
            error = integrator.Perform(surface, eps, False, False)
        mass = integrator.Mass()
        if error < 0 or not math.isfinite(error) or not math.isfinite(mass):
            raise ValueError("OCCT face integration failed")
        masses.append(mass)
        errors.append(error)
    return Integral(math.fsum(masses), math.fsum(errors))


def paired_volume(
    a: TopoDS_Shape, b: TopoDS_Shape, absolute_tolerance: float
) -> PairedVolume:
    if not math.isfinite(absolute_tolerance) or absolute_tolerance <= 0:
        raise ValueError("volume tolerance must be finite and positive")
    lhs, rhs = closed_faces(a), closed_faces(b)
    references = shared_references(a, b)
    samples = []
    for reference in references:
        coarse_a, coarse_b = (
            integrate(faces, reference, 1e-12) for faces in (lhs, rhs)
        )
        fine_a, fine_b = (integrate(faces, reference, 1e-13) for faces in (lhs, rhs))
        samples.append(
            ReferenceSample(
                reference=reference,
                a=fine_a,
                b=fine_b,
                delta=fine_a.volume - fine_b.volume,
                coarse_delta=coarse_a.volume - coarse_b.volume,
                a_quadrature_drift=abs(fine_a.volume - coarse_a.volume),
                b_quadrature_drift=abs(fine_b.volume - coarse_b.volume),
            )
        )
    max_delta = max(
        abs(value)
        for sample in samples
        for value in (sample.delta, sample.coarse_delta)
    )
    max_drift = max(
        value
        for sample in samples
        for value in (sample.a_quadrature_drift, sample.b_quadrature_drift)
    )
    # A paired result can agree even when each individual integral drifts with
    # origin (e.g. inherited trim closure residual). Expose that separately.
    volumes_a = [sample.a.volume for sample in samples]
    volumes_b = [sample.b.volume for sample in samples]
    return PairedVolume(
        method="shared-reference signed Gauss-Kronrod with spline spans",
        quadrature_eps=(1e-12, 1e-13),
        absolute_tolerance=absolute_tolerance,
        max_absolute_delta=max_delta,
        max_quadrature_drift=max_drift,
        within_absolute_tolerance=max_delta <= absolute_tolerance
        and max_drift <= absolute_tolerance,
        reference_spread={
            "a": max(volumes_a) - min(volumes_a),
            "b": max(volumes_b) - min(volumes_b),
        },
        samples=samples,
    )

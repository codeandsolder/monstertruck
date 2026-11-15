"""Independent OpenCascade Boolean material-occupancy contract.

Use source STEP operands and a candidate STEP result, NOT the candidate
kernel's internal point classifier. Source operands may contain multiple
disjoint solids. Boundary-ambiguous witnesses are reported, never coerced.
"""

import argparse
import json
import random
import sys
from itertools import starmap
from pathlib import Path

from OCP.BRepCheck import BRepCheck_Analyzer
from OCP.BRepClass3d import BRepClass3d_SolidClassifier
from OCP.gp import gp_Pnt
from OCP.IFSelect import IFSelect_RetDone
from OCP.STEPControl import STEPControl_Reader
from OCP.TopAbs import TopAbs_IN, TopAbs_ON, TopAbs_SOLID
from OCP.TopExp import TopExp_Explorer
from OCP.TopoDS import TopoDS


def read_classifiers(path: Path) -> list[BRepClass3d_SolidClassifier]:
    reader = STEPControl_Reader()
    if reader.ReadFile(str(path)) != IFSelect_RetDone or not reader.TransferRoots():
        raise ValueError(f"cannot import STEP: {path}")
    explorer = TopExp_Explorer(reader.OneShape(), TopAbs_SOLID)
    result = []
    while explorer.More():
        solid = TopoDS.Solid(explorer.Current())
        if not BRepCheck_Analyzer(solid).IsValid():
            raise ValueError(f"input has invalid BRep solid: {path}")
        result.append(BRepClass3d_SolidClassifier(solid))
        explorer.Next()
    return result


def classify(
    solids: list[BRepClass3d_SolidClassifier],
    point: tuple[float, float, float],
    tol: float,
) -> bool | None:
    p = gp_Pnt(*point)
    ambiguous = False
    for classifier in solids:
        classifier.Perform(p, tol)
        state = classifier.State()
        if state == TopAbs_IN:
            return True
        if state == TopAbs_ON:
            ambiguous = True
    return None if ambiguous else False


def material_occupancy(operation: str, left: bool, right: bool) -> bool:
    match operation:
        case "difference":
            return left and not right
        case "intersection":
            return left and right
        case "union":
            return left or right
        case "xor":
            return left != right
        case _:
            raise ValueError(f"unsupported operation: {operation}")


def read_points(path: Path) -> list[tuple[float, float, float]]:
    raw = json.loads(path.read_text(encoding="utf-8"))
    if isinstance(raw, dict):
        raw = raw["points"]
    if not isinstance(raw, list):
        raise TypeError(f"expected a list of 3D points in {path}")
    result = []
    for point in raw:
        if len(point) != 3:
            raise ValueError(f"expected a 3D point, got {point!r}")
        result.append(tuple(float(value) for value in point))
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--operation",
        required=True,
        choices=["difference", "intersection", "union", "xor"],
    )
    parser.add_argument("--left", required=True, type=Path)
    parser.add_argument("--right", required=True, type=Path)
    parser.add_argument("--result", required=True, type=Path)
    parser.add_argument(
        "--point",
        nargs=3,
        type=float,
        action="append",
        default=[],
        metavar=("X", "Y", "Z"),
    )
    parser.add_argument(
        "--points-json", type=Path, help='JSON list of [x,y,z], or {"points": [...]}'
    )
    parser.add_argument(
        "--box",
        nargs=6,
        type=float,
        metavar=("XMIN", "XMAX", "YMIN", "YMAX", "ZMIN", "ZMAX"),
    )
    parser.add_argument("--samples", type=int, default=0)
    parser.add_argument("--seed", type=int, default=20261115)
    parser.add_argument("--tolerance", type=float, default=1e-7)
    parser.add_argument("--report", type=Path)
    args = parser.parse_args()
    if not (0.0 < args.tolerance < float("inf")):
        parser.error("--tolerance must be finite and positive")
    if args.samples < 0 or (args.samples > 0) != (args.box is not None):
        parser.error("--box and positive --samples must appear together")
    points = [tuple(point) for point in args.point]
    if args.points_json:
        points.extend(read_points(args.points_json))
    if args.box:
        bounds = list(zip(args.box[::2], args.box[1::2], strict=True))
        if any(not (low < high) for low, high in bounds):
            parser.error("--box needs ordered finite coordinate bounds")
        rng = random.Random(args.seed)
        points.extend(tuple(starmap(rng.uniform, bounds)) for _ in range(args.samples))
    if not points:
        parser.error("supply --point, --points-json, or --box with --samples")
    if any(
        any(not (-float("inf") < v < float("inf")) for v in point) for point in points
    ):
        parser.error("points must have finite coordinates")

    left = read_classifiers(args.left)
    right = read_classifiers(args.right)
    result = read_classifiers(args.result)

    counts = {"total": 0, "matched": 0, "mismatched": 0, "boundary_ambiguous": 0}
    violations = []
    for index, point in enumerate(points):
        actual = classify(result, point, args.tolerance)
        a = classify(left, point, args.tolerance)
        b = classify(right, point, args.tolerance)
        counts["total"] += 1
        if None in (a, b, actual):
            counts["boundary_ambiguous"] += 1
            continue
        expected = material_occupancy(args.operation, a, b)
        if expected == actual:
            counts["matched"] += 1
        else:
            counts["mismatched"] += 1
            violations.append(
                {
                    "index": index,
                    "point": point,
                    "left": a,
                    "right": b,
                    "expected": expected,
                    "actual": actual,
                }
            )
    report = {
        "operation": args.operation,
        "left": str(args.left),
        "right": str(args.right),
        "result": str(args.result),
        "tolerance": args.tolerance,
        "counts": counts,
        "violations": violations,
    }
    if args.report:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))
    if counts["mismatched"]:
        return 1
    if not counts["matched"]:
        print("No unambiguous occupancy witnesses were available", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())

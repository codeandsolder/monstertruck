import unittest

from occt_paired_volume import closed_faces, paired_volume
from OCP.BRepAlgoAPI import BRepAlgoAPI_Cut, BRepAlgoAPI_Fuse
from OCP.BRepPrimAPI import BRepPrimAPI_MakeBox
from OCP.gp import gp_Pnt
from OCP.TopAbs import TopAbs_FACE
from OCP.TopExp import TopExp_Explorer
from OCP.TopoDS import TopoDS_Shape

TOLERANCE = 1e-10


def box_at(
    x: float, y: float, z: float, dx: float, dy: float, dz: float
) -> TopoDS_Shape:
    return BRepPrimAPI_MakeBox(gp_Pnt(x, y, z), dx, dy, dz).Shape()


class PairedVolumeTests(unittest.TestCase):
    def test_subdivision_and_operand_order(self) -> None:
        whole = box_at(0.0, 0.0, 0.0, 2.0, 3.0, 4.0)
        split = BRepAlgoAPI_Fuse(
            box_at(0.0, 0.0, 0.0, 1.0, 3.0, 4.0), box_at(1.0, 0.0, 0.0, 1.0, 3.0, 4.0)
        ).Shape()
        self.assertGreater(len(closed_faces(split)), len(closed_faces(whole)))
        forward = paired_volume(whole, split, TOLERANCE)
        reverse = paired_volume(split, whole, TOLERANCE)
        self.assertTrue(forward.within_absolute_tolerance)
        self.assertTrue(reverse.within_absolute_tolerance)
        for lhs, rhs in zip(forward.samples, reverse.samples, strict=True):
            self.assertEqual(lhs.reference, rhs.reference)
            self.assertAlmostEqual(lhs.a.volume, 24.0, delta=TOLERANCE)
            self.assertAlmostEqual(lhs.b.volume, 24.0, delta=TOLERANCE)
            self.assertAlmostEqual(lhs.delta, -rhs.delta, delta=TOLERANCE)

    def test_volume_change_is_not_hidden_by_reference(self) -> None:
        whole = box_at(0.0, 0.0, 0.0, 2.0, 3.0, 4.0)
        for height_change in (1e-6, 1e-10):
            with self.subTest(height_change=height_change):
                changed = box_at(0.0, 0.0, 0.0, 2.0, 3.0, 4.0 + height_change)
                result = paired_volume(whole, changed, TOLERANCE)
                self.assertFalse(result.within_absolute_tolerance)
                for sample in result.samples:
                    self.assertAlmostEqual(
                        sample.delta, -6 * height_change, delta=TOLERANCE
                    )

    def test_oriented_inner_shell_subtracts_material(self) -> None:
        outer = box_at(0.0, 0.0, 0.0, 3.0, 3.0, 3.0)
        cavity = box_at(1.0, 1.0, 1.0, 1.0, 1.0, 1.0)
        hollow = BRepAlgoAPI_Cut(outer, cavity).Shape()
        result = paired_volume(hollow, hollow, TOLERANCE)
        self.assertTrue(result.within_absolute_tolerance)
        for sample in result.samples:
            self.assertAlmostEqual(sample.a.volume, 26.0, delta=TOLERANCE)

    def test_open_face_rejected(self) -> None:
        solid = box_at(0.0, 0.0, 0.0, 2.0, 3.0, 4.0)
        face = TopExp_Explorer(solid, TopAbs_FACE).Current()
        with self.assertRaises(ValueError):
            paired_volume(solid, face, TOLERANCE)


if __name__ == "__main__":
    unittest.main()

//! Exercise 2D/3D transform derives on generic tuple structs and enum variants.
//! Both derive macros share one implementation but must select distinct matrix types.
#![cfg(feature = "derive")]

use monstertruck_core::math::{Matrix3, Matrix4, Point2, Point3, Vector2, Vector3};
use monstertruck_traits::{Transformed, TransformedM3, TransformedM4};

#[derive(Clone, TransformedM3)]
struct Wrapped2<T>(T)
where T: Transformed<Matrix3>;

#[derive(Clone, TransformedM4)]
struct Wrapped3<T>(T)
where T: Transformed<Matrix4>;

#[derive(Clone, TransformedM3)]
enum Choice2<T>
where T: Transformed<Matrix3> {
    Tuple(T),
    Named { inner: T },
}

#[derive(Clone, TransformedM4)]
enum Choice3<T>
where T: Transformed<Matrix4> {
    Tuple(T),
    Named { inner: T },
}

#[test]
fn generic_tuple_wrappers_preserve_in_place_and_copied_transform() {
    let move_2d = Matrix3::new_translation(&Vector2::new(3.0, -2.0));
    let mut original_2d = Wrapped2(Point2::new(1.0, 2.0));
    let translated = original_2d.transformed(move_2d);
    assert_eq!(translated.0, Point2::new(4.0, 0.0));
    assert_eq!(original_2d.0, Point2::new(1.0, 2.0));
    original_2d.transform_by(move_2d);
    assert_eq!(original_2d.0, translated.0);

    let move_3d = Matrix4::new_translation(&Vector3::new(2.0, -5.0, 3.0));
    let mut original_3d = Wrapped3(Point3::new(1.0, 2.0, 4.0));
    let translated = original_3d.transformed(move_3d);
    assert_eq!(translated.0, Point3::new(3.0, -3.0, 7.0));
    assert_eq!(original_3d.0, Point3::new(1.0, 2.0, 4.0));
    original_3d.transform_by(move_3d);
    assert_eq!(original_3d.0, translated.0);
}

#[test]
fn generic_enum_variants_dispatch_both_transformation_methods() {
    let move_2d = Matrix3::new_translation(&Vector2::new(3.0, -2.0));
    for named in [false, true] {
        let mut original = if named {
            Choice2::Named {
                inner: Point2::new(1.0, 2.0),
            }
        } else {
            Choice2::Tuple(Point2::new(1.0, 2.0))
        };
        let copy = original.transformed(move_2d);
        match copy {
            Choice2::Tuple(point) | Choice2::Named { inner: point } => {
                assert_eq!(point, Point2::new(4.0, 0.0));
            }
        }
        original.transform_by(move_2d);
        match original {
            Choice2::Tuple(point) | Choice2::Named { inner: point } => {
                assert_eq!(point, Point2::new(4.0, 0.0));
            }
        }
    }

    let move_3d = Matrix4::new_translation(&Vector3::new(2.0, -5.0, 3.0));
    for named in [false, true] {
        let mut original = if named {
            Choice3::Named {
                inner: Point3::new(1.0, 2.0, 4.0),
            }
        } else {
            Choice3::Tuple(Point3::new(1.0, 2.0, 4.0))
        };
        let copy = original.transformed(move_3d);
        match copy {
            Choice3::Tuple(point) | Choice3::Named { inner: point } => {
                assert_eq!(point, Point3::new(3.0, -3.0, 7.0));
            }
        }
        original.transform_by(move_3d);
        match original {
            Choice3::Tuple(point) | Choice3::Named { inner: point } => {
                assert_eq!(point, Point3::new(3.0, -3.0, 7.0));
            }
        }
    }
}

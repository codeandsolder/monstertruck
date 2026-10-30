/// This macro replaces the prelude module and loads the main public structures
/// in `monstertruck-topology`.
///
/// It redefines topological elements so they are related to geometric elements.
/// See [`imported`] for the exact structures and aliases loaded.
///
/// # Examples
/// ```
/// struct Point;
/// struct Curve;
/// struct Surface;
/// monstertruck_topology::prelude!(Point, Curve, Surface);
///
/// use std::any::type_name;
/// assert_eq!(
///     type_name::<Face>(),
///     type_name::<monstertruck_topology::Face<Point, Curve, Surface>>(),
/// );
/// ```
///
/// By adding the `pub` option to the end, one can also add `pub` to the generated code.
/// ```
/// mod reexport_topology {
///     monstertruck_topology::prelude!((), (), (), pub);
/// }
///
/// use std::any::type_name;
/// assert_eq!(
///     type_name::<reexport_topology::Shell>(),
///     type_name::<monstertruck_topology::Shell<(), (), ()>>(),
/// );
/// ```
/// ```compile_fail
/// mod first {
///     pub mod second {
///         monstertruck_topology::prelude!((), (), (), pub(super));
///     }
///     fn inside_public_range() {
///         let _ = second::Vertex::new(());
///     }
/// }
/// // compile fail
/// fn outside_public_range() {
///     let _ = first::second::Vertex::new(());
/// }
/// ```
///
/// [`imported`]: crate::imported
#[macro_export]
macro_rules! prelude {
    ($point: ty, $curve: ty, $surface: ty) => {
        $crate::prelude!(@impl []
            [expect(unused_imports, reason = "prelude consumers intentionally choose only a subset of convenience imports")]
            [expect(dead_code, reason = "prelude consumers intentionally choose only a subset of convenience aliases")]
            $point, $curve, $surface
        );
    };
    ($point: ty, $curve: ty, $surface: ty, $pub: tt $($super: tt)?) => {
        $crate::prelude!(@impl [$pub $($super)?] [] [] $point, $curve, $surface);
    };
    (@impl [$($vis: tt)*] [$($use_attr: meta)?] [$($alias_attr: meta)?]
        $point: ty, $curve: ty, $surface: ty
    ) => {
        $(#[$use_attr])?
        $($vis)* use $crate::{
            compress::CompressedEdgeIndex,
            shell::ShellCondition,
            VertexDisplayFormat,
            EdgeDisplayFormat,
            WireDisplayFormat,
            FaceDisplayFormat,
            ShellDisplayFormat,
            SolidDisplayFormat,
            wire,
            shell,
        };

        /// Vertex, the minimum topological unit.
        $(#[$alias_attr])?
        $($vis)* type Vertex = $crate::Vertex<$point>;
        /// Edge, which consists two vertices.
        $(#[$alias_attr])?
        $($vis)* type Edge = $crate::Edge<$point, $curve>;
        /// Wire, a path or cycle which consists some edges.
        $(#[$alias_attr])?
        $($vis)* type Wire = $crate::Wire<$point, $curve>;
        /// Face, attached to a simple and closed wire.
        $(#[$alias_attr])?
        $($vis)* type Face = $crate::Face<$point, $curve, $surface>;
        /// Shell, a connected compounded faces.
        $(#[$alias_attr])?
        $($vis)* type Shell = $crate::Shell<$point, $curve, $surface>;
        /// Solid, attached to a closed shells.
        $(#[$alias_attr])?
        $($vis)* type Solid = $crate::Solid<$point, $curve, $surface>;

        /// The id of vertex. `Copy` trait is implemented.
        $(#[$alias_attr])?
        $($vis)* type VertexId = $crate::VertexId<$point>;
        /// The id that does not depend on the direction of the edge.
        $(#[$alias_attr])?
        $($vis)* type EdgeId = $crate::EdgeId<$curve>;
        /// The id that does not depend on the direction of the face.
        $(#[$alias_attr])?
        $($vis)* type FaceId = $crate::FaceId<$surface>;

        /// Serialized compressed edge
        $(#[$alias_attr])?
        $($vis)* type CompressedEdge = $crate::compress::CompressedEdge<$curve>;
        /// Serialized compressed face
        $(#[$alias_attr])?
        $($vis)* type CompressedFace = $crate::compress::CompressedFace<$surface>;
        /// Serialized compressed shell
        $(#[$alias_attr])?
        $($vis)* type CompressedShell = $crate::compress::CompressedShell<$point, $curve, $surface>;
        /// Serialized compressed solid
        $(#[$alias_attr])?
        $($vis)* type CompressedSolid = $crate::compress::CompressedSolid<$point, $curve, $surface>;
    };
}
#[doc(hidden)]
pub mod empty_geometries {
    #![expect(
        missing_debug_implementations,
        reason = "empty marker geometries intentionally carry no state worth formatting"
    )]
    pub struct Point;
    pub struct Curve;
    pub struct Surface;
}
#[doc(hidden)]
pub use empty_geometries::*;
prelude!(Point, Curve, Surface, pub);

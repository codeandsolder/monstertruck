# Native geometry fixtures

These shape fixtures contain the same geometry and topology as the `resources4truck` submodule at `9bf9de1426c5fc2f2b8d63f501dca0d3c53ebd91`. Only the representation of points, vectors and matrices changed: coordinate objects become arrays in coordinate order, and matrix column objects become flat column-major arrays. Numeric values and all non-math fields are unchanged. The `just create-shape-json` recipe regenerates them using the native geometry types.

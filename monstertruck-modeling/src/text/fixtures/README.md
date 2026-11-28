# Synthetic OpenType fixtures

All four tiny fonts here contain original, deliberately simple geometric outlines
(.notdef, space, A and O), generated with fontTools 4.61.1 for this project.
No external typeface or font artwork is incorporated. They may be used under
the repository's Apache-2.0 license.

The design uses 1000 font units per em, 600 units of A/O advance and 250
units of space advance. A contains a quadratic outline in the TrueType
font, a cubic outline in the CFF OpenType font; O has two oppositely
oriented rectangular contours. Collection face 1 and the variable font's
wght=900 master expand the outline width from 500 to 750 units and the
advance from 600 to 900. The default weight is 400. The tests exercise
the parser, scaling, variation, collection selection and uniform placement
without relying on host fonts or rasterization.

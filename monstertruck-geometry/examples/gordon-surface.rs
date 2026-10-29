//! Demonstrates `BsplineSurface::gordon` -- creating a surface from two
//! families of curves using the boolean sum formula.

use monstertruck_geometry::prelude::{
    BsplineCurve, BsplineSurface, KnotVector, ParametricSurface, Point3,
};

fn main() -> std::result::Result<(), Box<dyn std::error::Error>> {
    // A 3x3 grid of intersection points on a curved surface.
    let grid: Vec<Vec<Point3>> = (0..3)
        .map(|i| {
            let y = f64::from(i);
            (0..3)
                .map(|j| {
                    let x = f64::from(j);
                    let z = (x * 0.5).sin() * (y * 0.5).cos();
                    Point3::new(x, y, z)
                })
                .collect()
        })
        .collect();

    // U-curves: connect grid points along the x-direction (j varies).
    let u_curves: Vec<BsplineCurve<Point3>> = grid
        .iter()
        .map(|row| BsplineCurve::new_unchecked(KnotVector::bezier_knot(2), row.clone()))
        .collect();

    // V-curves: connect grid points along the y-direction (i varies).
    let v_curves: Vec<BsplineCurve<Point3>> = (0..3)
        .map(|j| {
            let pts: Vec<Point3> = (0..3).map(|i| grid[i][j]).collect();
            BsplineCurve::new_unchecked(KnotVector::bezier_knot(2), pts)
        })
        .collect();

    let surface = BsplineSurface::gordon(u_curves, v_curves, &grid)?;

    println!("Gordon surface from 3x3 curve network:");
    println!("  u-degree = {}", surface.udegree());
    println!("  v-degree = {}", surface.vdegree());

    // Verify the surface interpolates the grid points.
    println!("\nGrid point verification (should match input):");
    for (v, row) in [0.0, 0.5, 1.0].into_iter().zip(&grid) {
        for (u, &expected) in [0.0, 0.5, 1.0].into_iter().zip(row) {
            let pt = surface.subs(u, v);
            let err = (pt.z - expected.z)
                .mul_add(
                    pt.z - expected.z,
                    (pt.x - expected.x).mul_add(pt.x - expected.x, (pt.y - expected.y).powi(2)),
                )
                .sqrt();
            println!(
                "  ({u:.1}, {v:.1}): ({:.3}, {:.3}, {:.3})  expected ({:.3}, {:.3}, {:.3})  err={err:.2e}",
                pt.x, pt.y, pt.z, expected.x, expected.y, expected.z,
            );
        }
    }
    Ok(())
}

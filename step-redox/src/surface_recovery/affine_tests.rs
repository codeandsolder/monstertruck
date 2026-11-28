use super::*;
use crate::brep::{SurfaceSupport, surface_support};
use std::fmt::Write as _;

fn fixture(
    degree: usize,
    coordinates: &[f64],
    multiplicities: &str,
    knots: &str,
    rational: bool,
) -> anyhow::Result<Vec<EntityInstance>> {
    let mut data = String::new();
    let mut rows = Vec::new();
    for (row, (x, z)) in [(0.0, 0.0), (1.0, 0.5), (2.0, 0.0)].into_iter().enumerate() {
        let mut ids = Vec::new();
        for (column, &v) in coordinates.iter().enumerate() {
            let id = row * coordinates.len() + column + 1;
            ids.push(format!("#{id}"));
            writeln!(
                data,
                "#{id}=CARTESIAN_POINT('',({x:?},{:?},{z:?}));",
                2.0 * v
            )?;
        }
        rows.push(format!("({})", ids.join(",")));
    }
    let grid = format!("({})", rows.join(","));
    if rational {
        let weights = ["1.0", "0.8", "1.0"]
            .map(|w| format!("({})", vec![w; coordinates.len()].join(",")))
            .join(",");
        writeln!(
            data,
            "#100=(BOUNDED_SURFACE() B_SPLINE_SURFACE(2,{degree},{grid},.UNSPECIFIED.,.F.,.F.,.F.) B_SPLINE_SURFACE_WITH_KNOTS((3,3),{multiplicities},(0.,1.),{knots},.UNSPECIFIED.) GEOMETRIC_REPRESENTATION_ITEM() RATIONAL_B_SPLINE_SURFACE(({weights})) REPRESENTATION_ITEM('') SURFACE());"
        )?;
    } else {
        writeln!(
            data,
            "#100=B_SPLINE_SURFACE_WITH_KNOTS('',2,{degree},{grid},.UNSPECIFIED.,.F.,.F.,.F.,(3,3),{multiplicities},(0.,1.),{knots},.UNSPECIFIED.);"
        )?;
    }
    data.push_str("#101=ADVANCED_FACE('',(),#100,.T.);");
    let input = format!(
        "ISO-10303-21;\nHEADER;\nFILE_DESCRIPTION(('x'),'1');\nFILE_NAME('a','b',(''),(''),'x','y','');\nFILE_SCHEMA(('AUTOMOTIVE_DESIGN'));\nENDSEC;\nDATA;\n{data}\nENDSEC;\nEND-ISO-10303-21;\n"
    );
    let exchange = step_p21::parser::parse(&input)?;
    Ok(exchange
        .data
        .into_iter()
        .next()
        .ok_or_else(|| anyhow::anyhow!("missing section"))?
        .entities)
}

#[test]
fn recovers_degree_elevated_and_nonuniform_multispan_affine_extrusions() -> anyhow::Result<()> {
    for rational in [false, true] {
        for (degree, coordinates, multiplicities, knots) in [
            (3, &[0.0, 1.0 / 3.0, 2.0 / 3.0, 1.0][..], "(4,4)", "(0.,1.)"),
            (
                2,
                &[0.0, 0.125, 0.5, 0.875, 1.0][..],
                "(3,1,1,3)",
                "(0.,0.25,0.75,1.)",
            ),
        ] {
            let mut entities = fixture(degree, coordinates, multiplicities, knots, rational)?;
            let index = build_index(&entities);
            let evidence = analyze_v_extrusion_surface(100, &entities, &index, 1.0e-7)
                .ok_or_else(|| anyhow::anyhow!("affine extrusion not recognized"))?;
            assert!(distance(evidence.extrusion_mm, [0.0, 2.0, 0.0]) < 1.0e-12);
            assert!(evidence.max_residual_mm < 1.0e-12);
            let SurfaceSupport::SplineExtrusion(support) = surface_support(100, &entities, &index)
            else {
                anyhow::bail!("native B-rep classifier did not consume the affine evidence");
            };
            assert!(distance(support.extrusion_mm, evidence.extrusion_mm) < 1.0e-12);
            assert_eq!(
                recover_v_extrusion_surfaces(&mut entities).surfaces_recovered,
                1
            );
            let index = build_index(&entities);
            assert_eq!(
                simple_record(&entities[index[&100]]).map(|record| record.name.as_str()),
                Some("SURFACE_OF_LINEAR_EXTRUSION")
            );
        }
    }
    Ok(())
}

#[test]
fn refuses_nonaffine_interior_even_when_boundary_profiles_agree() -> anyhow::Result<()> {
    let mut entities = fixture(
        2,
        &[0.0, 0.2, 0.5, 0.875, 1.0],
        "(3,1,1,3)",
        "(0.,0.25,0.75,1.)",
        false,
    )?;
    assert_eq!(
        recover_v_extrusion_surfaces(&mut entities).surfaces_recovered,
        0
    );
    Ok(())
}

#[test]
fn refuses_discontinuous_and_unsafe_rational_control_nets() -> anyhow::Result<()> {
    let mut discontinuous = fixture(1, &[0.0, 0.5, 0.5, 1.0], "(2,2,2)", "(0.,0.5,1.)", false)?;
    assert_eq!(
        recover_v_extrusion_surfaces(&mut discontinuous).surfaces_recovered,
        0
    );
    let mut entities = fixture(2, &[0.0, 0.5, 1.0], "(3,3)", "(0.,1.)", true)?;
    let index = build_index(&entities);
    let EntityInstance::Complex { subsuper, .. } = &mut entities[index[&100]] else {
        anyhow::bail!("expected rational surface");
    };
    let rational = subsuper
        .0
        .iter_mut()
        .find(|r| r.name == "RATIONAL_B_SPLINE_SURFACE")
        .ok_or_else(|| anyhow::anyhow!("missing rational weights"))?;
    rational.parameter = Parameter::List(vec![Parameter::List(vec![
        Parameter::List(vec![
            Parameter::Real(1.0e-20),
            Parameter::Real(2.0e-20),
            Parameter::Real(1.0e-20),
        ]),
        Parameter::List(vec![Parameter::Real(1.0e-20); 3]),
        Parameter::List(vec![Parameter::Real(1.0e-20); 3]),
    ])]);
    assert_eq!(
        recover_v_extrusion_surfaces(&mut entities).surfaces_recovered,
        0
    );
    Ok(())
}

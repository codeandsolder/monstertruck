use monstertruck_topology::{shell::ShellCondition, *};

// Cases where truck version <= 0.4.0 fails
#[test]
fn singular_vertices_in_certain_boundaries() {
    let v = Vertex::from_points([(); 12]);
    let edge = [
        Edge::new_unchecked(&v[0], &v[1], ()),
        Edge::new_unchecked(&v[1], &v[2], ()),
        Edge::new_unchecked(&v[0], &v[7], ()),
        Edge::new_unchecked(&v[1], &v[8], ()),
        Edge::new_unchecked(&v[2], &v[9], ()),
        Edge::new_unchecked(&v[3], &v[4], ()),
        Edge::new_unchecked(&v[4], &v[3], ()),
        Edge::new_unchecked(&v[5], &v[6], ()),
        Edge::new_unchecked(&v[6], &v[5], ()),
        Edge::new_unchecked(&v[7], &v[8], ()),
        Edge::new_unchecked(&v[8], &v[9], ()),
        Edge::new_unchecked(&v[7], &v[10], ()),
        Edge::new_unchecked(&v[9], &v[11], ()),
        Edge::new_unchecked(&v[10], &v[11], ()),
    ];
    let shell: Shell<(), (), ()> = vec![
        Face::new_unchecked(
            vec![
                wire![
                    edge[9].inverse(),
                    edge[2].inverse(),
                    edge[0].clone(),
                    edge[3].clone(),
                ],
                wire![edge[5].clone(), edge[6].clone()],
            ],
            (),
        ),
        Face::new_unchecked(
            vec![
                wire![
                    edge[3].inverse(),
                    edge[1].clone(),
                    edge[4].clone(),
                    edge[10].inverse(),
                ],
                wire![edge[7].clone(), edge[8].clone()],
            ],
            (),
        ),
        Face::new_unchecked(
            vec![wire![
                edge[10].clone(),
                edge[12].clone(),
                edge[13].inverse(),
                edge[11].inverse(),
                edge[9].clone(),
            ]],
            (),
        ),
    ]
    .into();
    assert_eq!(shell.shell_condition(), ShellCondition::Oriented);
    assert!(shell.singular_vertices().is_empty());
}

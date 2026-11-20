// SPDX-License-Identifier: LGPL-3.0-or-later
use super::*;

#[test]
fn simple_cubic_argon_agrees_with_molecular_shell_sum() {
    let side = 16.0;
    let model = Model::v2_0_1(Conventions::Gxtb);
    let molecular = |positions: &[V]| {
        let mol = Molecule::new(&vec![18; positions.len()], positions, 0., None).unwrap();
        let d4 = RevD4::new(&model, &mol).unwrap();
        // At this separation even the nearest Ar neighbours give zero CN.
        assert!(d4.cn.iter().all(|&cn| cn == 0.0));
        (
            d4.two_body(&vec![0.; positions.len()]).0,
            d4.three_body(&model, &mol, None).0,
        )
    };
    let origin = [0.; 3];
    let cell = [[side, 0., 0.], [0., side, 0.], [0., 0., side]];
    let periodic = Geometry::new(&[18], &[origin], cell, [true; 3]).unwrap();

    // Number of simple-cubic neighbours at squared distances 1..9 in units
    // of side^2. The 50-bohr pair cutoff lies below the tenth shell.
    let multiplicities = [6, 12, 8, 6, 24, 24, 0, 12, 30];
    let mut expected_pair = 0.;
    for (shell, count) in multiplicities.iter().enumerate() {
        let r = side * ((shell + 1) as f64).sqrt();
        expected_pair += f64::from(*count) * molecular(&[origin, [r, 0., 0.]]).0 / 2.;
    }
    // The 25-bohr ATM cutoff allows 12 right face triangles and the eight
    // faces of two inscribed regular tetrahedra per conventional cube.
    let right = molecular(&[origin, cell[0], cell[1]]).1;
    let equilateral = molecular(&[origin, [side, side, 0.], [side, 0., side]]).1;
    let mut pair = [0.];
    periodic.electronic(&[0.], &mut pair, &mut [0.]);
    assert!((pair[0] - expected_pair).abs() < 1e-16);
    assert!((periodic.atm[0] - (12. * right + 8. * equilateral)).abs() < 1e-16);
}

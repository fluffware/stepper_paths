const GAUSS_QUAD_WEIGHTS_2: &[f64] = &[1.0, 1.0];
const GAUSS_QUAD_ABSCISSAS_2: &[f64] = &[-0.5773502692, 0.5773502692];

const GAUSS_QUAD_WEIGHTS_3: &[f64] = &[0.5555555556, 0.88888888889, 0.5555555556];
const GAUSS_QUAD_ABSCISSAS_3: &[f64] = &[-0.7745966692, 0.0, 0.7745966692];

const GAUSS_QUAD_WEIGHTS_4: &[f64] = &[
    0.3478548451374538,
    0.6521451548625461,
    0.6521451548625461,
    0.3478548451374538,
];
const GAUSS_QUAD_ABSCISSAS_4: &[f64] = &[
    -0.8611363115940526,
    -0.3399810435848563,
    0.3399810435848563,
    0.8611363115940526,
];

const GAUSS_QUAD_WEIGHTS_5: &[f64] = &[
    0.2369268850561891,
    0.4786286704993665,
    0.5688888888888889,
    0.4786286704993665,
    0.2369268850561891,
];
const GAUSS_QUAD_ABSCISSAS_5: &[f64] = &[
    -0.906179845938664,
    -0.5384693101056831,
    0.0,
    0.5384693101056831,
    0.906179845938664,
];

const GAUSS_QUAD_WEIGHTS_6: &[f64] = &[
    0.1713244923791704,
    0.3607615730481386,
    0.467913934572691,
    0.467913934572691,
    0.3607615730481386,
    0.1713244923791704,
];
const GAUSS_QUAD_ABSCISSAS_6: &[f64] = &[
    -0.932469514203152,
    -0.6612093864662645,
    -0.2386191860831969,
    0.2386191860831969,
    0.6612093864662645,
    0.932469514203152,
];

const GAUSS_QUAD_WEIGHTS_7: &[f64] = &[
    0.1294849661688697,
    0.2797053914892766,
    0.3818300505051189,
    0.4179591836734694,
    0.3818300505051189,
    0.2797053914892766,
    0.1294849661688697,
];
const GAUSS_QUAD_ABSCISSAS_7: &[f64] = &[
    -0.9491079123427585,
    -0.7415311855993945,
    -0.4058451513773972,
    0.0,
    0.4058451513773972,
    0.7415311855993945,
    0.9491079123427585,
];

const GAUSS_QUAD_WEIGHTS_8: &[f64] = &[
    0.1012285362903763,
    0.2223810344533745,
    0.3137066458778873,
    0.362683783378362,
    0.362683783378362,
    0.3137066458778873,
    0.2223810344533745,
    0.1012285362903763,
];
const GAUSS_QUAD_ABSCISSAS_8: &[f64] = &[
    -0.9602898564975363,
    -0.7966664774136267,
    -0.525532409916329,
    -0.1834346424956498,
    0.1834346424956498,
    0.525532409916329,
    0.7966664774136267,
    0.9602898564975363,
];

const WEIGHTS_N: [&[f64]; 7] = [
    GAUSS_QUAD_WEIGHTS_2,
    GAUSS_QUAD_WEIGHTS_3,
    GAUSS_QUAD_WEIGHTS_4,
    GAUSS_QUAD_WEIGHTS_5,
    GAUSS_QUAD_WEIGHTS_6,
    GAUSS_QUAD_WEIGHTS_7,
    GAUSS_QUAD_WEIGHTS_8,
];
const ABSCISSAS_N: [&[f64]; 7] = [
    GAUSS_QUAD_ABSCISSAS_2,
    GAUSS_QUAD_ABSCISSAS_3,
    GAUSS_QUAD_ABSCISSAS_4,
    GAUSS_QUAD_ABSCISSAS_5,
    GAUSS_QUAD_ABSCISSAS_6,
    GAUSS_QUAD_ABSCISSAS_7,
    GAUSS_QUAD_ABSCISSAS_8,
];

#[allow(clippy::many_single_char_names)]
pub fn gauss_quad_integrate<F>(a: f64, b: f64, n: usize, f: F) -> f64
where
    F: Fn(f64) -> f64,
{
    let weights = WEIGHTS_N[n - 2];
    let abscissas = ABSCISSAS_N[n - 2];
    let mut sum = 0f64;
    let k = (b - a) / 2.0;
    let m = (a + b) / 2.0;
    for i in 0..weights.len() {
        sum += weights[i] * f(k * abscissas[i] + m);
    }
    k * sum
}

#[cfg(test)]
fn assert_close(err: f64, got: f64, expected: f64) {
    println!("{} == {}", got, expected);
    assert!(((expected - got) / expected).abs() <= err);
}

#[test]
fn test_gauss_integrate() {
    const ERR: f64 = 0.0001;
    assert_close(ERR, gauss_quad_integrate(-3.0, 4.0, 2, |x| x), 3.5);
    assert_close(ERR, gauss_quad_integrate(-3.0, 4.0, 3, |x| x), 3.5);
    assert_close(ERR, gauss_quad_integrate(-3.0, 4.0, 4, |x| x), 3.5);
    assert_close(ERR, gauss_quad_integrate(-3.0, 4.0, 5, |x| x), 3.5);
    assert_close(ERR, gauss_quad_integrate(-3.0, 4.0, 6, |x| x), 3.5);
    assert_close(ERR, gauss_quad_integrate(-3.0, 4.0, 7, |x| x), 3.5);
    assert_close(ERR, gauss_quad_integrate(-3.0, 4.0, 8, |x| x), 3.5);

    assert_close(
        ERR,
        gauss_quad_integrate(-3.0, 4.0, 2, |x| x * x * x + x * x + x + 1.0),
        (64.0 + 64.0 / 3.0 + 8.0 + 4.0) - (81.0 / 4.0 - 27.0 / 3.0 + 9.0 / 2.0 - 3.0),
    );

    assert_close(
        0.001,
        gauss_quad_integrate(0.01, 4.0, 8, |x| x.sqrt()),
        5.332666666666666,
    );
}

use super::polynom::Polynom as PolynomT;
use coords::Point;
use super::sqrt_derivative::Terms as SqrtTerms;
use super::inverse_derivative::Terms as InvTerms;

type Polynom = PolynomT<f64>;

const N_COEF: usize = 8;

pub fn length_to_param(c1: Point, c2: Point, p2: Point) -> Polynom {
    let x_diff = Polynom::from(
        [
            3.0 * c1.x,
            6.0 * (c2.x - 2.0 * c1.x),
            3.0 * (p2.x + 3.0 * (c1.x - c2.x)),
        ]
        .as_slice(),
    );
    let y_diff = Polynom::from(
        [
            3.0 * c1.y,
            6.0 * (c2.y - 2.0 * c1.y),
            3.0 * (p2.y + 3.0 * (c1.y - c2.y)),
        ]
        .as_slice(),
    );
    println!("x_diff: {:?} y_diff: {:?}", x_diff, y_diff);
    let base = &(&x_diff * &x_diff) + &(&y_diff * &y_diff);
    println!("bas: {:?}", base);
    let mut sqrt_diffs = Vec::new();
    let t0 = 0.0;
    let mut diff = SqrtTerms::new(base);
    sqrt_diffs.push(diff.value(t0));
    for _ in 1..N_COEF {
	diff = diff.diff();
	sqrt_diffs.push(diff.value(t0));
    }
    println!("sqrt_diffs: {:?}", sqrt_diffs);
    let mut taylor =  Vec::new();
    taylor.push(0.0);
    let mut fac = 1.0;
    for (i, d) in sqrt_diffs.iter().enumerate() {
	taylor.push(d/fac);
	fac *= (i + 2) as f64;
    }
    taylor.reverse();
    println!("taylor: {:?}", taylor);

    let mut inv_diffs = Vec::new();
    let mut diff = InvTerms::first();
    inv_diffs.push(diff.value(&sqrt_diffs));
    for _ in 1..N_COEF {
	diff = diff.diff();
	inv_diffs.push(diff.value(&sqrt_diffs));
    }
    println!("inv_diffs: {:?}", inv_diffs);
    let mut taylor =  Vec::new();
    taylor.push(0.0);
    let mut fac = 1.0;
    for (i, d) in inv_diffs.iter().enumerate() {
	taylor.push(d/fac);
	fac *= (i + 2) as f64;
    }
    taylor.reverse();
    println!("taylor: {:?}", taylor);


    Polynom::from(vec![1.0])
}

#[test]
fn test_length_to_param() {
    length_to_param(
        Point { x: 10.0, y: 0.0 },
        Point { x: 1.1, y: 24.0 },
        Point { x: 23.0, y: 15.0 },
    );
}

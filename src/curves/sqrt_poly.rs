use std::mem;

// Represents c* s^(e/2) + P(t)
#[derive(Debug)]
struct SqrtPoly {
    c: f64,
    e: i32,
    polynom: Vec<f64>, // Polynom, lowest degree coefficient first
}

pub fn polyval(polynom: &[f64], t: f64) -> f64 {
    let mut i = polynom.len();
    let mut sum = 0.0;
    while i >= 1 {
        i -= 1;
        sum = sum * t + polynom[i];
    }
    sum
}

fn polyder(polynom: &[f64]) -> Vec<f64> {
    let len = polynom.len();
    let mut deriv = Vec::with_capacity(len - 1);
    for (i,p) in polynom[1..len].iter().enumerate() {
        deriv.push(p * ((i+1) as f64));
    }
    deriv
}

fn polyint(polynom: &[f64]) -> Vec<f64> {
    let len = polynom.len();
    let mut int = Vec::with_capacity(len + 1);
    int.push(0.0);
    for (i,p) in  polynom.iter().enumerate() {
        int.push(*p / ((i + 1) as f64));
    }
    int
}

pub fn polymul(a: &[f64], b: &[f64]) -> Vec<f64> {
    let alen = a.len();
    let blen = b.len();
    let mut prod = Vec::with_capacity(alen + blen - 1);
    for i in 0..(alen + blen - 1) {
        let mut sum = 0.0;
        let low_b = if i >= alen { i - alen + 1 } else { 0 };
        let high_b = i.min(blen - 1) + 1;
        for j in low_b..high_b {
            sum += a[i - j] * b[j];
        }
        prod.push(sum);
    }
    prod
}

pub fn polyadd<'a>(mut a: &'a [f64], mut b: &'a [f64]) -> Vec<f64> {
    if a.len() > b.len() {
        mem::swap(&mut a, &mut b)
    }
    let mut sum = b.to_vec();
    for i in 0..a.len() {
        sum[i] += a[i]
    }
    sum
}

#[derive(Debug)]
struct DiffSqrtPoly {
    t: f64,         // Value to differentiate at
    base: Vec<f64>, // Start polynom
    base_val: f64,  // base evaluated at t
}

impl DiffSqrtPoly {
    pub fn new(base: &[f64], t: f64) -> DiffSqrtPoly {
        DiffSqrtPoly {
            t,
            base_val: polyval(base, t),
            base: base.to_vec(),
        }
    }

    pub fn eval(&self, sp: &SqrtPoly) -> f64 {
        self.base_val.powf(f64::from(sp.e) / 2.0) * sp.c * polyval(&sp.polynom, self.t)
    }
}

fn diff_single_sqrt_poly(diff: &DiffSqrtPoly, sp: &SqrtPoly) -> (SqrtPoly, SqrtPoly) {
    let diff_sqrt = SqrtPoly {
        c: sp.c * 0.5 * f64::from(sp.e),
        e: sp.e - 2,
        polynom: polymul(&sp.polynom, &polyder(&diff.base)),
    };
    let diff_poly = SqrtPoly {
        c: sp.c,
        e: sp.e,
        polynom: polyder(&sp.polynom),
    };
    (diff_sqrt, diff_poly)
}

fn diff_sqrt_poly(diff: &DiffSqrtPoly, spa: &[SqrtPoly]) -> (Vec<SqrtPoly>, f64) {
    let mut res = Vec::new();
    let mut sum = 0.0;
    for sp in spa {
        let (d1, d2) = diff_single_sqrt_poly(diff, sp);
        sum += diff.eval(&d1);
        res.push(d1);
        if !d2.polynom.is_empty() {
            sum += diff.eval(&d2);
            res.push(d2);
        }
    }
    (res, sum)
}

/// Calculate the next row in Pascal's Triangle
fn binomial_step(row: &mut [u16]) {
    let len = row.len();
    row[0] = 1;
    let mut prev = 1;
    for i in 0..(len - 1) {
        mem::swap(&mut prev, &mut row[i]);
        prev += row[i + 1];
    }
    row[len - 1] = 1;
}
/// Add a constant to the parameter
fn polyoffset(polynom: &[f64], offset: f64) -> Vec<f64> {
    let len = polynom.len();
    let mut res = Vec::with_capacity(len);
    res.resize(len, 0.0f64);
    let mut bin = Vec::new();
    bin.resize(len, 0);
    for i in 0..len {
        binomial_step(&mut bin[..i + 1]);
        let mut f = 1.0;
        for j in 0..(i + 1) {
            res[i - j] += polynom[i] * f64::from(bin[j]) * f;
            f *= offset;
        }
    }
    res
}

fn polyweight(a: &[f64], ta: f64, b: &[f64], tb: f64) -> Vec<f64> {
    polyadd(
        &polymul(&[-tb / (ta - tb), 1.0 / (ta - tb)], a),
        &polymul(&[-ta / (tb - ta), 1.0 / (tb - ta)], b),
    )
}

/// Get values of the n first Taylor coefficients for sqrt(polynom)
/// The polynom takes t - t0 as parameter
pub fn taylor_sequence(polynom: &[f64], t0: f64, n: u32) -> Vec<f64> {
    let diff = DiffSqrtPoly::new(polynom, t0);
    let mut sp = vec![new_sqrt_poly()];
    let mut factorial = 1.0f64;
    let mut coeffs = Vec::with_capacity(n as usize);
    let mut value = polyval(polynom, t0).sqrt();
    coeffs.push(value);
    for i in 1..n {
        (sp, value) = diff_sqrt_poly(&diff, &sp);
        factorial *= f64::from(i);
        coeffs.push(value / factorial);
    }
    coeffs
}

fn new_sqrt_poly() -> SqrtPoly {
    SqrtPoly {
        c: 1.0,
        e: 1,
        polynom: vec![1.0],
    }
}

#[derive(Debug)]
pub struct IntegralSegment {
    start: f64,    // Parameter value at the start of the segment
    int: Vec<f64>, // Polynom for the integral starting at start
}

pub fn build_segments(polynom: &[f64], low: f64, high: f64, n: u32) -> Vec<IntegralSegment> {
    let step = (high - low) / f64::from(n + 1);
    let mut int_seq = Vec::new();
    let mut prev_seq = taylor_sequence(polynom, step * 0.5, 7);
    let mut int = polyint(&polyoffset(&prev_seq, -step * 0.5));

    let mut int_sum = 0.0;
    let end_val = polyval(&int, step / 2.0);
    int[0] = int_sum;
    int_sum += end_val;
    int_seq.push(IntegralSegment { start: low, int });
    for i in 1..n + 1 {
        let t0 = low + step * (f64::from(i) + 0.5);
        let t_prev = t0 - step;
        let seq = taylor_sequence(polynom, t0, 7);
        let mut int = polyint(&polyweight(&prev_seq, 0.0, &polyoffset(&seq, -step), step));
        let end_val = polyval(&int, step);
        int[0] = int_sum;
        int_sum += end_val;
        int_seq.push(IntegralSegment { start: t_prev, int });
        prev_seq = seq;
    }
    let mut int = polyint(&prev_seq);
    let end_val = polyval(&int, step * 0.5);
    int[0] = int_sum;
    int_sum += end_val;
    int_seq.push(IntegralSegment {
        start: high - step / 2.0,
        int,
    });
    let mut int = polyint(&polyoffset(&prev_seq, step * 0.5));
    int[0] = int_sum;
    int_seq.push(IntegralSegment { start: high, int });
    int_seq
}

pub fn eval_int(ints: &[IntegralSegment], t: f64) -> f64 {
    match ints.binary_search_by(|int| int.start.partial_cmp(&t).unwrap()) {
        Ok(i) => ints[i].int[0],
        Err(i) => polyval(&ints[i - 1].int, t - ints[i - 1].start),
    }
}

pub fn find_param(ints: &[IntegralSegment], l: f64) -> f64 {
    match ints.binary_search_by(|seg| seg.int[0].partial_cmp(&l).unwrap()) {
        Ok(i) => ints[i].start,
        Err(i) => {
            assert!(i < ints.len());
            let mut low = 0.0;
            let mut high = ints[i].start - ints[i - 1].start;
            let p = &ints[i - 1].int;
            while high - low > f64::EPSILON {
                let mid = (high + low) * 0.5;
                let v = polyval(p, mid);
                if v > l {
                    high = mid;
                } else {
                    low = mid;
                }
            }
	    low + ints[i - 1].start
        }
    }
}

#[test]
fn test_polymul() {
    assert_eq!(
        polymul(&[1.0, 5.0], &[3.0, 2.0, 6.0]),
        &[3.0, 17.0, 16.0, 30.0]
    );
    assert_eq!(
        polymul(&[3.0, 2.0, 6.0], &[1.0, 5.0]),
        &[3.0, 17.0, 16.0, 30.0]
    );
}

#[test]
fn test_diff_single_sqrt_poly() {
    let p = [36.0, -360.0, 1548.0, -2160.0, 954.0];
    let diff = DiffSqrtPoly::new(&p, 0.5);
    let sp = new_sqrt_poly();
    let (d1, _) = diff_single_sqrt_poly(&diff, &sp);
    assert_eq!(d1.e, -1);
    assert_eq!(d1.c, 0.5);
    assert_eq!(&d1.polynom, &[-360.0, 3096.0, -6480.0, 3816.0]);
}

#[test]
fn test_diff_sqrt_poly() {
    let p = [36.0, -360.0, 1548.0, -2160.0, 954.0];
    let diff = DiffSqrtPoly::new(&p, 0.5);
    let sp = new_sqrt_poly();
    let (res1, value) = diff_sqrt_poly(&diff, &[sp]);
    assert_eq!(res1.len(), 1);
    let d1 = &res1[0];
    assert_eq!(d1.e, -1);
    assert_eq!(d1.c, 0.5);
    assert_eq!(&d1.polynom, &[-360.0, 3096.0, -6480.0, 3816.0]);
    assert_relative_eq!(value, 3.939192985791677);

    let (res2, value) = diff_sqrt_poly(&diff, &res1);
    assert_eq!(res2.len(), 2);
    let d1 = &res2[0];
    assert_eq!(d1.e, -3);
    assert_eq!(d1.c, -0.25);
    assert_eq!(
        &d1.polynom,
        &polymul(
            &[-360.0, 3096.0, -6480.0, 3816.0],
            &[-360.0, 3096.0, -6480.0, 3816.0]
        )
    );
    let d2 = &res2[1];
    assert_eq!(d2.e, -1);
    assert_eq!(d2.c, 0.5);
    assert_eq!(&d2.polynom, &[3096.0, -12960.0, 11448.0]);
    assert_eq!(value, -48.41132345297081);

    let (_res3, value) = diff_sqrt_poly(&diff, &res2);
    assert_relative_eq!(value, -32.19552545438488, epsilon = 1e-6);
}

#[test]
fn test_taylor_sequence() {
    let p = [36.0, -360.0, 1548.0, -2160.0, 954.0];
    let n = 7;
    let t = 0.5;
    let res = taylor_sequence(&p, t, n as u32);

    let mut ps = format!("{}", res[0]);
    for i in 1..n {
        ps += &format!(" + {}*t.^{}", res[i], i);
    }
    println!("{}", ps);
    assert_relative_eq!(res[0], 77706389.0 / 13604465.0, epsilon = 1e-6);
    assert_relative_eq!(res[1], 21837779.0 / 5543719.0, epsilon = 1e-6);
    assert_relative_eq!(res[2], -77115559.0 / 3185848.0, epsilon = 1e-6);
    assert_relative_eq!(res[3], -704987637.0 / 131382413.0, epsilon = 1e-6);
    assert_relative_eq!(res[4], 291003320.0 / 8100977.0, epsilon = 1e-6);
    assert_relative_eq!(res[5], -334588303.0 / 7041954.0, epsilon = 1e-6);
    assert_relative_eq!(res[6], 246779675.0 / 1352379.0, epsilon = 1e-5);
}
#[test]
fn test_binomial_step() {
    let mut c = [0u16; 5];
    binomial_step(&mut c[..1]);
    assert_eq!(c[..1], [1]);
    binomial_step(&mut c[..2]);
    assert_eq!(c[..2], [1, 1]);
    binomial_step(&mut c[..3]);
    assert_eq!(c[..3], [1, 2, 1]);
    binomial_step(&mut c[..4]);
    assert_eq!(c[..4], [1, 3, 3, 1]);
    binomial_step(&mut c[..5]);
    assert_eq!(c[..5], [1, 4, 6, 4, 1]);
}

#[test]
fn test_polyoffset() {
    let p = [2.0, 3.0, 1.0];
    assert_eq!(&polyoffset(&p, 2.0), &[12.0, 7.0, 1.0]);
    let p = [3.0, 4.0, 1.0, 2.0, 6.0];
    assert_eq!(&polyoffset(&p, -2.0), &[79.0, -168.0, 133.0, -46.0, 6.0]);
}
#[test]
fn test_polyweight() {
    let p = [5.0, 0.0, 1.0];
    assert_eq!(&polyweight(&p, 0.2, &p, 0.3), &[5.0, 0.0, 1.0, 0.0]);
}
#[test]
fn test_build_segments() {
    let p = [4.0, 0.0, 0.0, 0.0, 0.0];
    let ints = build_segments(&p, 0.0, 1.0, 9);
    assert_relative_eq!(polyval(&ints[10].int, 1.0 - ints[10].start), 2.0);
    let p = [1.0, 2.0, 1.0, 0.0, 0.0];
    let ints = build_segments(&p, 0.0, 1.0, 9);
    assert_eq!(polyval(&ints[10].int, 1.0 - ints[10].start), 1.5);
    let p = [36.0, -360.0, 1548.0, -2160.0, 954.0];
    let ints = build_segments(&p, 0.0, 1.0, 9);
    assert_relative_eq!(
        polyval(&ints[10].int, 1.0 - ints[10].start),
        4.693604306325399,
        epsilon = 1e-4
    );
}

#[test]
fn test_eval_int() {
    let p = [4.0, 0.0, 0.0, 0.0, 0.0];
    let ints = build_segments(&p, 0.0, 1.0, 9);
    for i in 0..14 {
        let t = f64::from(i) / 13.0;
        let v = eval_int(&ints, t);
        println!("{t} -> {v}");
        assert_relative_eq!(v, t * 2.0, epsilon = 1e-6);
    }

    let p = [4.0, 4.0, -11.0, -6.0, 9.0]; // (-3 t^2 + t + 2)^2
    let ints = build_segments(&p, 0.0, 1.0, 9);
    for i in 0..14 {
        let t = f64::from(i) / 13.0;
        let v = eval_int(&ints, t);
        println!("{t} -> {v}");
        let r = ((-t + 0.5) * t + 2.0) * t;
        assert_relative_eq!(v, r, epsilon = 1e-8);
    }

    for i in 0..10 {
        let t = f64::from(i) / 10.0 + 0.05;
        let v = eval_int(&ints, t);
        println!("{t} -> {v}");
        let r = ((-t + 0.5) * t + 2.0) * t;
        assert_relative_eq!(v, r, epsilon = 1e-8);
    }
}
#[test]
fn test_find_param() {
    let p = [4.0, -12.0, 9.0]; // (-3 t + 2)^2
    let ints = build_segments(&p, 0.0, 1.0, 9);
    let length = eval_int(&ints, 1.0);
    for i in 0..13 {
        let l = f64::from(i) * length / 12.0;
	let t = find_param(&ints, l);
	assert_relative_eq!(l, eval_int(&ints, t), epsilon=1e-9);
    }

    let p = [4.0, -12.0, 17.0, -12.0, 4.0]; // (2 t^2 -3 t + 2)^2
    let ints = build_segments(&p, 0.0, 1.0, 9);
    let length = eval_int(&ints, 1.0);
    for i in 0..13 {
        let l = f64::from(i) * length / 12.0;
	let t = find_param(&ints, l);
	assert_relative_eq!(l, eval_int(&ints, t), epsilon=1e-9);
    }

    
}

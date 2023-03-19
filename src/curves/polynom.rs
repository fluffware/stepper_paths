use std::fmt::Debug;
use std::mem;
use std::ops::RangeInclusive;
use std::ops::{Add, Div, Mul, Neg};

pub struct Polynom<const N: usize, T> {
    // c[i] contains the coefficient for t^(i)
    c: [T; N],
}

fn value<const N: usize, T>(poly: &[T; N], n: usize, t: T) -> T
where
    T: Mul<Output = T> + Add<Output = T> + Copy,
{
    let mut sum = poly[0];
    let mut f = t;
    for p in poly[1..n].iter() {
        sum = sum + f * *p;
        f = f * t;
    }
    sum
}

fn deriv<const N: usize, T>(poly: &[T; N]) -> [T; N]
where
    T: Mul<Output = T> + From<i8> + Copy,
{
    let mut deriv = [T::from(0i8); N];
    for i in 0..N - 1 {
        deriv[i] = T::from((i + 1) as i8) * poly[i + 1];
    }
    deriv[N - 1] = 0i8.into();
    deriv
}

pub trait PolynomCoef:
    Add<Output = Self>
    + Mul<Output = Self>
    + Div<Output = Self>
    + Neg<Output = Self>
    + From<i8>
    + Copy
    + PartialEq
    + PartialOrd
    + Debug
{
    fn close_to_zero(a: Self) -> bool;
}

impl PolynomCoef for f64 {
    fn close_to_zero(a: f64) -> bool {
        (-10.0 * f64::EPSILON..=10.0 * f64::EPSILON).contains(&a)
    }
}

/// Find a single root for polynom
/// There must be at most one root in the given range
fn find_single_root<const N: usize, T>(
    poly: &[T; N],
    n: usize,
    range: &RangeInclusive<T>,
) -> Option<T>
where
    T: PolynomCoef,
{
    let mut low = *range.start();
    let mut high = *range.end();
    let low_value = value(poly, n, low);
    let high_value = value(poly, n, high);
    if high_value * low_value > T::from(0i8) {
        return None;
    }
    if low_value > high_value {
        mem::swap(&mut low, &mut high);
    }

    while !T::close_to_zero(high + -low) {
        let mid = (high + low) / T::from(2i8);
        let mid_value = value(poly, n, mid);
        if mid_value < T::from(0i8) {
            low = mid;
        } else {
            high = mid;
        }
    }
    Some(low)
}

fn find_roots<const N: usize, T>(
    poly: &[T; N],
    mut n: usize,
    range: &RangeInclusive<T>,
) -> ([T; N], usize)
where
    T: PolynomCoef,
{
    // Reduce degree if possible
    while n > 0 && poly[n - 1] == 0i8.into() {
        n -= 1;
    }
    let mut roots = [T::from(0i8); N];
    let n_roots = match n {
        0 => {
            roots[0] = T::from(0i8);
            1
        }
        1 => 0,
        2 => {
            // Linear
            let v1 = poly[0] + poly[1] * *range.start();
            let v2 = poly[0] + poly[1] * *range.end();
            if v1 * v2 < T::from(0i8) {
                roots[0] = -(poly[0] / poly[1]);
                1
            } else {
                0
            }
        }
        _ => {
            let (mut extremes, n_extremes) = find_roots(&deriv(poly), n - 1, range);
            extremes[0..n_extremes].sort_by(|a, b| a.partial_cmp(b).unwrap());
            println!("extremes: {:?}", &extremes[0..n_extremes]);
            let mut low = *range.start();
            let mut n_roots = 0;
            extremes[n_extremes] = *range.end();
            for high in &extremes[0..n_extremes + 1] {
                if *high > *range.end() {
                    break;
                }
                if *high > low {
                    if let Some(root) = find_single_root(poly, n, &(low..=*high)) {
                        println!("Root: {root:?}");
                        if n_roots == 0 || !T::close_to_zero(roots[n_roots - 1] + -root) {
                            roots[n_roots] = root;
                            n_roots += 1;
                        }
                    }
                }
                low = *high;
            }
            n_roots
        }
    };
    (roots, n_roots)
}

impl<const N: usize, T> Polynom<N, T>
where
    T: Mul<Output = T> + From<i8> + Copy,
{
    pub fn value(&self, t: T) -> T
    where
        T: Mul<Output = T> + Add<Output = T>,
    {
        value(&self.c, N, t)
    }

    pub fn deriv(&self) -> [T; N] {
        deriv(&self.c)
    }

    pub fn real_roots(&self, range: &RangeInclusive<T>) -> ([T; N], usize)
    where
        T: PolynomCoef,
    {
        find_roots(&self.c, N, range)
    }
}

#[test]
fn test_find_roots() {
    let p = [0.0, 0.0];
    assert_eq!(find_roots(&p, 2, &(0.0..=1.0)), ([0.0, 0.0], 1));
    let p = [0.1, 0.0];
    assert_eq!(find_roots(&p, 2, &(0.0..=1.0)), ([0.0, 0.0], 0));
    let p = [1.0, -4.0];
    assert_eq!(find_roots(&p, 2, &(0.0..=1.0)), ([0.25, 0.0], 1));

    let p = [-10.0, 2.0, 3.5];
    let (roots, n_roots) = find_roots(&p, 3, &(-10.0..=10.0));
    assert_eq!(n_roots, 2);
    assert_relative_eq!(roots[0], -2.0, epsilon = 1e-8);
    assert_relative_eq!(roots[1], 10.0 / 7.0, epsilon = 1e-8);

    let p = [0.0, 0.0, 0.1];
    let (roots, n_roots) = find_roots(&p, 3, &(-10.0..=10.0));
    assert_eq!(n_roots, 1);
    assert_relative_eq!(roots[0], 0.0);

    let p = [10.0, 2.0, 0.1];
    let (roots, n_roots) = find_roots(&p, 3, &(-11.0..=10.0));
    assert_eq!(n_roots, 1);
    assert_relative_eq!(roots[0], -10.0);

    let p = [0.0, -50.0, 0.1, 2.0];
    let (roots, n_roots) = find_roots(&p, 4, &(-11.0..=10.0));
    assert_eq!(n_roots, 3);
    assert_relative_eq!(roots[0], -(40001.0f64.sqrt() + 1.0) / 40.0, epsilon = 1e-8);
    assert_relative_eq!(roots[1], 0.0, epsilon = 1e-8);
    assert_relative_eq!(roots[2], (40001.0f64.sqrt() - 1.0) / 40.0, epsilon = 1e-8);
}

#[test]
fn test_find_single_root() {
    let p = [0.0, 0.0, 0.1];
    assert_eq!(find_single_root(&p, 3, &(0.0..=1.0)), Some(0.0));
    let p = [-4.0, 0.0, 1.0];
    assert_relative_eq!(
        find_single_root(&p, 3, &(0.0..=2.3)).unwrap(),
        2.0,
        epsilon = 1e-8
    );
    assert_relative_eq!(
        find_single_root(&p, 3, &(-4.0..=0.0)).unwrap(),
        -2.0,
        epsilon = 1e-8
    );
}

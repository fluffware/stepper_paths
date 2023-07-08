use std::fmt::Debug;
use std::mem;
use std::ops::RangeInclusive;
use std::ops::{Add, AddAssign, Div, Mul, Neg, Index, IndexMut};

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


#[macro_export]
macro_rules! poly {
    ($($coef:expr),+) => {
	{
	    let p = vec![$($coef),+];
	    Polynom::from(p)
	}
    }
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

// [i] contains the coefficient for t^(i)
#[derive(Debug, PartialEq, Clone)]
pub struct Polynom<T>(Vec<T>);

impl<T> Polynom<T> {
    pub fn value<V>(&self, t: V) -> V
    where
        T: Copy,
        V: From<T> + Mul<Output = V> + Add<Output = V> + Copy + From<i8>,
    {
	if self.0.is_empty() {
	    return  V::from(0i8);
	}
        let mut sum = V::from(self.0[0]);
        let mut f = t;
        for p in &self.0[1..] {
            sum = sum + f * V::from(*p);
            f = f * t;
        }
        sum
    }

    pub fn derive(&self) -> Polynom<T>
    where
        T: Mul<Output = T> + From<i8> + Copy,
    {
        let n = self.0.len();
        let mut deriv = Vec::with_capacity(n - 1);
        for i in 1..n {
            deriv.push(T::from((i) as i8) * self.0[i]);
        }
        Polynom(deriv)
    }
    pub fn integrate(&self) -> Polynom<T>
    where
        T: Div<Output = T> + From<i8> + Copy,
    {
        let n = self.0.len();
        let mut int = Vec::with_capacity(n + 1);
	int.push(T::from(0i8));
        for i in 0..n {
            int.push(self.0[i] / T::from((i+1) as i8));
        }
        Polynom(int)
    }

   

    /// Find a single root for polynom
    /// There must be at most one root in the given range
    fn find_single_root(&self, range: &RangeInclusive<T>) -> Option<T>
    where
        T: PolynomCoef,
    {
        let mut low = *range.start();
        let mut high = *range.end();
        let low_value = self.value(low);
        let high_value = self.value(high);
        if high_value * low_value > T::from(0i8) {
            return None;
        }
        if low_value > high_value {
            mem::swap(&mut low, &mut high);
        }

        while !T::close_to_zero(high + -low) {
            let mid = (high + low) / T::from(2i8);
            let mid_value = self.value(mid);
            if mid_value < T::from(0i8) {
                low = mid;
            } else {
                high = mid;
            }
        }
        Some(low)
    }

    fn find_roots(&self, range: &RangeInclusive<T>) -> Vec<T>
    where
        T: PolynomCoef,
    {
        let n = self.0.len();
        let mut roots = Vec::new();
        match n {
            0 => {
                roots.push(T::from(0i8));
            }
            1 => {}
            2 => {
                // Linear
                let v1 = self.0[0] + self.0[1] * *range.start();
                let v2 = self.0[0] + self.0[1] * *range.end();
                if v1 * v2 < T::from(0i8) {
                    roots.push(-(self.0[0] / self.0[1]));
                }
            }
            _ => {
                let mut extremes = self.derive().find_roots(range);
                extremes.sort_by(|a, b| a.partial_cmp(b).unwrap());
                println!("extremes: {:?}", &extremes);
                let mut low = *range.start();
                extremes.push(*range.end());
                for high in &extremes {
                    if *high > *range.end() {
                        break;
                    }
                    if *high > low {
                        if let Some(root) = self.find_single_root(&(low..=*high)) {
                            println!("Root: {root:?}");
                            if roots
                                .last()
                                .map(|last| !T::close_to_zero(*last + -root))
                                .unwrap_or(true)
                            {
                                roots.push(root);
                            }
                        }
                    }
                    low = *high;
                }
            }
        };
        roots
    }

    pub fn real_roots(&self, range: &RangeInclusive<T>) -> Vec<T>
    where
        T: PolynomCoef,
    {
        self.find_roots(range)
    }

    pub fn len(&self) -> usize
    {
	self.0.len()
    }
    
    /// Add a constant to the parameter
    pub fn offset(&self, offset: T) -> Polynom<T>
	where T: From<u16> + Mul<T, Output=T> + Copy + AddAssign<T>
    {
	let len = self.0.len();
	let mut res = vec![T::from(0u16);len];
	let mut bin = vec!{0u16;len};
	for i in 0..len {
            binomial_step(&mut bin[..i + 1]);
            let mut f = T::from(1u16);
            for j in 0..(i + 1) {
		res[i - j] += self.0[i] * T::from(bin[j]) * f;
		f = f * offset;
            }
	}
	Polynom(res)
    }
}

impl<T> Add<&Polynom<T>> for &Polynom<T>
where
    T: Add<T, Output = T> + Copy,
{
    type Output = Polynom<T>;
    fn add(self, other: &Polynom<T>) -> Polynom<T> {
        let out_len = self.0.len().max(other.0.len());
        let mut res = Vec::with_capacity(out_len);
        for (a, b) in std::iter::zip(&self.0, &other.0) {
            res.push(*a + *b);
        }
        if self.0.len() > other.0.len() {
            res.extend_from_slice(&self.0[other.0.len()..]);
        } else if self.0.len() < other.0.len() {
            res.extend_from_slice(&other.0[self.0.len()..]);
        }
        Polynom(res)
    }
}

impl<T> Mul<&Polynom<T>> for &Polynom<T>
where
    T: AddAssign<T> + Mul<T, Output = T> + Copy + num::traits::Zero,
{
    type Output = Polynom<T>;
    fn mul(self, other: &Polynom<T>) -> Polynom<T> {
        let out_len = self.0.len() + other.0.len() - 1;
        let mut res = vec![T::zero(); out_len];
        for (i, a) in self.0.iter().enumerate() {
            for (j, b) in other.0.iter().enumerate() {
                res[i + j] += *a * *b;
            }
        }
        Polynom(res)
    }
}

impl<T> Mul<T> for &Polynom<T>
where
    T: AddAssign<T> + Mul<T, Output = T> + Copy + num::traits::Zero,
{
    type Output = Polynom<T>;
    fn mul(self, other: T) -> Polynom<T> {
        let mut res = Vec::with_capacity(self.len());
        for a in self.0.iter() {
            res.push(*a*other);
        }
        Polynom(res)
    }
}

impl<T> From<&[T]> for Polynom<T>
    where T: Copy
{
    fn from(coef: &[T]) -> Polynom<T>
    {
	Polynom(Vec::from(coef))
    }
}

impl<T> From<Vec<T>> for Polynom<T>
    where T: Copy
{
    fn from(coef: Vec<T>) -> Polynom<T>
    {
	Polynom(coef)
    }
}

impl<T> Index<usize> for Polynom<T> {
    type Output = T;
    fn index(&self, i: usize) -> &Self::Output
    {
	self.0.index(i)
    }
}

impl<T> IndexMut<usize> for Polynom<T> {
    fn index_mut(&mut self, i: usize) -> &mut Self::Output
    {
	self.0.index_mut(i)
    }
}

#[test]
fn test_find_roots() {
    let p = Polynom(Vec::new());
    assert_eq!(p.find_roots(&(0.0..=1.0)), &[0.0]);
    let p = poly![0.1];
    assert_eq!(p.find_roots(&(0.0..=1.0)), &[0.0; 0]);
    let p = poly![1.0, -4.0];
    assert_eq!(p.find_roots(&(0.0..=1.0)), &[0.25]);

    let p = poly![-10.0, 2.0, 3.5];
    let roots = p.find_roots(&(-10.0..=10.0));
    assert_eq!(roots.len(), 2);
    assert_relative_eq!(roots[0], -2.0, epsilon = 1e-8);
    assert_relative_eq!(roots[1], 10.0 / 7.0, epsilon = 1e-8);

    let p = poly![0.0, 0.0, 0.1];
    let roots = p.find_roots(&(-10.0..=10.0));
    assert_eq!(roots.len(), 1);
    assert_relative_eq!(roots[0], 0.0);

    let p = poly![10.0, 2.0, 0.1];
    let roots = p.find_roots(&(-11.0..=10.0));
    assert_eq!(roots.len(), 1);
    assert_relative_eq!(roots[0], -10.0);
}

#[test]
fn test_find_single_root() {
    let p = poly![0.0, 0.0, 0.1];
    assert_eq!(p.find_single_root(&(0.0..=1.0)), Some(0.0));
    let p = poly![-4.0, 0.0, 1.0];
    assert_relative_eq!(
        p.find_single_root(&(0.0..=2.3)).unwrap(),
        2.0,
        epsilon = 1e-8
    );
    assert_relative_eq!(
        p.find_single_root(&(-4.0..=0.0)).unwrap(),
        -2.0,
        epsilon = 1e-8
    );
}

#[test]
fn test_value() {
    assert_eq!(poly!(2, 5, 3).value(2), 2 + 10 + 12);
    assert_eq!(poly!(2, 5, 3).value(3.0), 2.0 + 15.0 + 27.0);
    assert_eq!(poly!(2.0, 5.0, 3.0).value(2.0), 24.0);
}

#[test]
fn test_deriv() {
    assert_eq!(poly!(2, 5, 3).derive(), poly!(5, 6));
    assert_eq!(poly!(1.0, 4.0, 2.0, 3.0).derive(), poly!(4.0, 4.0, 9.0));
}

#[test]
fn test_integrate() {
    assert_eq!(poly!(2.0, 5.0, 3.0).integrate(), poly!(0.0, 2.0, 2.5, 1.0));
}
#[test]
fn test_add() {
    assert_eq!(
        &poly!(1, 8, 6) + &poly!(9, 3, 1, 4, 2),
        poly!(10, 11, 7, 4, 2)
    );
    assert_eq!(
        &poly!(1.0, 8.0, 6.7, 2.0) + &poly!(9.0, 4.3, 2.0),
        poly!(10.0, 12.3, 8.7, 2.0)
    );
}

#[test]
fn test_mul() {
    assert_eq!(
        &poly!(1, 8, 6) * &poly!(9, 3, 1, 4, 2),
        poly!(9, 75, 79, 30, 40, 40, 12)
    );
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

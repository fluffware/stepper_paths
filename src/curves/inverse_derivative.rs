use super::tex::ToTeX;
use num::traits::Pow;
use std::fmt::Write as _;
use std::ops::Mul;

#[derive(Clone, Debug, PartialEq)]
struct Term {
    c: f64,
    exp: Vec<i32>, // An exponent for each derivate, starting with the first derivative
}

fn grow_and_add(a: &mut Vec<i32>, i: usize, v: i32) {
    if a.len() <= i {
        a.resize(i + 1, 0);
    }
    a[i] += v;
}

impl Term {
    pub fn value(&self, diff: &[f64]) -> f64 {
        let mut prod = self.c;
        for (e, f) in self.exp.iter().zip(diff) {
            if *e != 0 {
                prod *= f.pow(e);
            }
        }
        prod
    }

    fn diff(&self) -> Terms {
        let mut terms = Vec::new();
        for i in 0..self.exp.len() {
            let exp_i = self.exp[i];
            if exp_i != 0 {
                let mut exp = self.exp.clone();
                grow_and_add(&mut exp, i + 1, 1);
                let c = self.c * f64::from(exp[i]);
                grow_and_add(&mut exp, i, -1);
                terms.push(Term { c, exp });
            }
        }
        Terms(terms)
    }
}

impl Mul<&Term> for &Term {
    type Output = Term;
    fn mul(self, other: &Term) -> Term {
        let c = self.c * other.c;
        let mut exp = Vec::new();
        for (a, b) in std::iter::zip(&self.exp, &other.exp) {
            exp.push(*a + *b);
        }
        if self.exp.len() > other.exp.len() {
            exp.extend_from_slice(&self.exp[other.exp.len()..]);
        } else if self.exp.len() < other.exp.len() {
            exp.extend_from_slice(&other.exp[self.exp.len()..]);
        }
        Term { c, exp }
    }
}

impl<W> ToTeX<W> for Term
where
    W: std::io::Write,
{
    fn to_tex(&self, w: &mut W) -> std::io::Result<()> {
        let mut num = String::new();
        let mut den = String::new();
        for (d, exp) in self.exp.iter().enumerate() {
            if *exp > 0 {
                write!(num, " f^{{({})}}(x)", d + 1).unwrap();
                if *exp != 1 {
                    write!(num, "^{exp}").unwrap();
                }
            }
            if *exp < 0 {
                write!(den, " f^{{({})}}(x)", d + 1).unwrap();
                if *exp != -1 {
                    write!(den, "^{}", -exp).unwrap();
                }
            }
        }
        if self.c == -1.0 {
            write!(w, "-")?;
        } else if self.c != 1.0 {
            write!(w, "{}", self.c)?;
        }
        if den.is_empty() {
            write!(w, "{num}")?;
        } else {
            write!(
                w,
                "\\frac{{{}}}{{{}}}",
                if num.is_empty() { "1".to_string() } else { num },
                den
            )?;
        }
        Ok(())
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct Terms(Vec<Term>);

impl Terms {
    pub fn first() -> Terms {
        Terms(vec![Term {
            c: 1.0,
            exp: vec![-1],
        }])
    }

    pub fn diff(&self) -> Terms {
        let mut terms = Vec::new();
        for t in &self.0 {
            terms.append(&mut t.diff().0);
        }

        // Merge terms with the same exponents
        terms.sort_unstable_by(|a, b| a.exp.cmp(&b.exp));
        terms.dedup_by(|a, b| {
            let equal = a.exp == b.exp;
            if equal {
                b.c += a.c
            };
            equal
        });
        Terms(terms)
    }
    
    pub fn value(&self, diff: &[f64]) -> f64 {
	let mut sum = 0.0;
        for t in &self.0 {
            sum += t.value(diff);
        }
	sum
    }
}
impl Mul<&Term> for &Terms {
    type Output = Terms;
    fn mul(self, other: &Term) -> Terms {
        let mut terms = Vec::new();
        for t in &self.0 {
            terms.push(t * other);
        }
        Terms(terms)
    }
}

impl<W> ToTeX<W> for Terms
where
    W: std::io::Write,
{
    fn to_tex(&self, w: &mut W) -> std::io::Result<()> {
        let mut term_iter = self.0.iter();
        let Some(term) = term_iter.next() else {return Ok(())};
        term.to_tex(w)?;
        while let Some(term) = term_iter.next() {
            write!(w, " + ")?;
            term.to_tex(w)?;
        }
        Ok(())
    }
}

pub fn inverse_diff(terms: &Terms) -> Terms {
    &terms.diff()
        * &Term {
            c: 1.0,
            exp: vec![-1],
        }
}

#[cfg(test)]
use std::io::Write;

#[test]
fn term_diff_test() {
    let term = Term {
        c: 1.0,
        exp: vec![1],
    };
    assert_eq!(
        term.diff(),
        Terms(vec![Term {
            c: 1.0,
            exp: vec![0, 1]
        }])
    );
    let term = Term {
        c: 1.0,
        exp: vec![2],
    };
    let diff1 = term.diff();
    assert_eq!(
        diff1,
        Terms(vec![Term {
            c: 2.0,
            exp: vec![1, 1]
        }])
    );
    let diff2 = diff1.diff();
    assert_eq!(
        diff2,
        Terms(vec![
            Term {
                c: 2.0,
                exp: vec![0, 2]
            },
            Term {
                c: 2.0,
                exp: vec![1, 0, 1]
            }
        ])
    );

    let diff3 = diff2.diff();
    assert_eq!(
        diff3,
        Terms(vec![
            Term {
                c: 6.0,
                exp: vec![0, 1, 1]
            },
            Term {
                c: 2.0,
                exp: vec![1, 0, 0, 1]
            }
        ])
    );

    let diff1 = Terms::first();

    let w = &mut std::io::stdout();
    diff1.to_tex(w).unwrap();
    write!(w, "\n").unwrap();
    let diff2 = inverse_diff(&diff1);
    diff2.to_tex(w).unwrap();
    write!(w, "\n").unwrap();
    let diff3 = inverse_diff(&diff2);
    diff3.to_tex(w).unwrap();
    write!(w, "\n").unwrap();
    let diff4 = inverse_diff(&diff3);
    diff4.to_tex(w).unwrap();
    write!(w, "\n").unwrap();
    let diff5 = inverse_diff(&diff4);
    let diff5 = diff5;
    diff5.to_tex(w).unwrap();
    write!(w, "\n").unwrap();
}

#[test]
fn term_value_test() {}

use super::polynom::Polynom as PolynomT;
type Polynom = PolynomT<f64>;

// Represents c* s^(e+1/2) * P(t)
#[derive(Debug, PartialEq, Clone)]
struct Term {
    e: i32,
    polynom: Polynom,
}


impl Term {
    fn diff(&self, base: &Polynom) -> Vec<Term> {
        let mut terms = Vec::new();
        terms.push(Term {
            e: self.e - 1,
            polynom: &self.polynom * &((&base.derive()) * (f64::from(self.e) + 0.5)),
        });
        if (self.polynom.len() == 1 && self.polynom[0] != 1.0) || self.polynom.len() > 1 {
            terms.push(Term {
                e: self.e,
                polynom: self.polynom.derive(),
            });
        }
        terms
    }

    fn value(&self, base: &Polynom, x: f64) -> f64
    {
	let base_value = base.value(x);
	base_value.powf(self.e as f64 + 0.5)* self.polynom.value(x) 
    }
}

#[derive(Debug, PartialEq, Clone)]
pub struct Terms {
    base: Polynom,
    terms: Vec<Term>,
}

impl Terms {
    pub fn new(base: Polynom) -> Terms
    {
	Terms{base, terms: vec!{Term{e:0, polynom: Polynom::from([1.0].as_slice())}}}
    }
    
    pub fn diff(&self) -> Terms {
        let mut terms = Vec::new();
        for t in &self.terms {
            terms.append(&mut t.diff(&self.base));
        }
        // Merge terms with the same exponents
        terms.sort_unstable_by(|a, b| a.e.cmp(&b.e));
        terms.dedup_by(|a, b| {
            let equal = a.e == b.e;
            if equal {
                b.polynom = &b.polynom + &a.polynom
            };
            equal
        });
        Terms {
            base: self.base.clone(),
            terms,
        }
    }

    pub fn value(&self, x: f64) -> f64
    {
	let mut sum = 0.0;
	for t in &self.terms {
            sum += t.value(&self.base, x);
        }
	sum
    }

    
}
#[cfg(test)]
use ::poly;

#[test]
fn test_diff_term() {
    let t = Term {
        e: 0,
        polynom: poly!(1.0),
    };
    println!("t= {t:?}");
    let base = poly!(1.0, 2.0, 1.0);
    assert_eq!(
        t.diff(&base),
        vec!(Term {
            e: -1,
            polynom: poly!(1.0, 1.0)
        })
    );
}

#[test]
fn test_diff_terms() {
    let t = Terms {
        base: poly!(2.0, 1.0, 1.0, 3.0),
        terms: vec![Term {
            e: 0,
            polynom: poly! {1.0},
        }],
    };
    let d1 = t.diff();
    assert_eq!(
        d1,
        Terms {
            base: poly!(2.0, 1.0, 1.0, 3.0),
            terms: vec! {Term{e:-1, polynom:poly!{0.5,1.0, 4.5}}}
        }
    );
    assert_eq!(d1.value(0.5),1.2020815280171309);
    let d2 = d1.diff();
    assert_eq!(
        d2,
        Terms {
            base: poly!(2.0, 1.0, 1.0, 3.0),
            terms: vec! {Term{e:-2, polynom: poly!(-1.0/4.0,-4.0/4.0,-22.0/4.0, -36.0/4.0,-81.0/4.0)}, Term{e:-1, polynom:poly!{01.0,9.0}} }
        }
    );
    assert_eq!(d2.value(0.5),2.2938543981691604);
}

use std::fmt::Debug;
use std::fmt::Error;
use std::fmt::Formatter;

/// A single range segment
#[derive(PartialEq)]
struct Range<T: Sized + PartialOrd> {
    low: T,
    high: T,
}

impl<T: Debug + PartialOrd> Debug for Range<T> {
    fn fmt(&self, f: &mut Formatter) -> Result<(), Error> {
        write!(f, "[{:?} - {:?}]", self.low, self.high)
    }
}

fn min<T: PartialOrd>(a: T, b: T) -> T {
    if a < b {
        a
    } else {
        b
    }
}

fn max<T: PartialOrd>(a: T, b: T) -> T {
    if a > b {
        a
    } else {
        b
    }
}

/// An ordered sequence of disjoint ranges
#[derive(PartialEq)]
pub struct MultiRange<T: Sized + PartialOrd> {
    ranges: Vec<Range<T>>,
}

impl<T: Sized + PartialOrd + Copy + Debug> MultiRange<T> {
    /// Returns an empty sequence
    pub fn new() -> MultiRange<T> {
        MultiRange { ranges: Vec::new() }
    }

    /// Returns a sequence with one range
    ///
    /// # Arguments
    /// * `low` - Low limit of range, inclusive
    /// * `high` - High limit of range, exclusive

    pub fn range(low: T, high: T) -> MultiRange<T> {
        MultiRange {
            ranges: vec![Range::<T> { low, high }],
        }
    }

    /// Add a slice of ranges
    ///
    /// # Arguments
    /// * `a` - Slice containing (low, high) tuples

    pub fn or_slice(&mut self, a: &[(T, T)]) -> &mut Self {
        for &(low, high) in a {
            self.or(low, high);
        }
        self
    }

    /// Add a single range
    ///
    /// # Arguments
    /// * `low` - Low limit of range, inclusive
    /// * `high` - High limit of range, exclusive

    pub fn or(&mut self, low: T, high: T) -> &mut Self {
        assert!(low <= high);
        let mut first = 0;
        while first < self.ranges.len() {
            if low <= self.ranges[first].high {
                break;
            }
            first += 1;
        }
        let mut last = first;
        while last < self.ranges.len() {
            if self.ranges[last].low > high {
                break;
            }
            last += 1;
        }
        if first == last {
            self.ranges.insert(first, Range { low, high });
        } else {
            let low = min(self.ranges[first].low, low);
            let high = max(self.ranges[last - 1].high, high);
            self.ranges[first] = Range { low, high };
            if first + 1 < last {
                self.ranges.drain(first + 1..last);
            }
        }
        self
    }

    /// Intersect with range.
    /// Returns ranges or parts of ranges between `low` and `high`.
    ///
    /// # Arguments
    /// * `low` - Prune ranges lower than this
    /// * `high` - Prune ranges greater or equal to this

    pub fn and(&self, low: T, high: T) -> Self {
        let mut res = MultiRange::new();
        let mut i = self.ranges.iter();
        while let Some(r) = i.next() {
            if r.high > low {
                let mut rn = r;
                while rn.low < high {
                    res.ranges.push(Range {
                        low: max(low, rn.low),
                        high: min(high, rn.high),
                    });
                    rn = match i.next() {
                        Some(r) => r,
                        None => break,
                    }
                }
                break;
            }
        }
        res
    }
    /// Intersect with multi range.
    /// Returns ranges that's part of both `self` and `other`
    ///
    /// # Arguments
    /// * `other` - `MultiRange` to intersect with.

    pub fn and_range(&self, other: &Self) -> Self {
        let mut r2_iter = other.ranges.iter().peekable();
        let mut res = MultiRange::new();
        let mut r2 = match r2_iter.next() {
            Some(r) => r,
            None => return res,
        };
        'match_ranges:
        // look for a the first range in other that overlaps with r1
        for r1 in &self.ranges {
            while r2.high <= r1.low {
                r2 = match r2_iter.next() {
                    Some(r) => r,
                    None => break 'match_ranges
                }
            }
            if r2.low < r1.high {
                res.ranges.push(Range {low: max(r1.low, r2.low),
                                  high: min(r1.high, r2.high)});
            }
            loop {
                {
                    let rn = match r2_iter.peek() {
                        Some(r) => r,
                        // No more ranges to compare with
                        None => break
                    };
                    if rn.low >= r1.high {
                        break;
                    }
                    r2 = rn;
                    res.ranges.push(Range {low: max(r1.low, rn.low),
                                           high: min(r1.high, rn.high)});
                }
                r2_iter.next();
            }
        }
        res
    }

    /// Returns true if there are no ranges in the sequence.

    pub fn is_empty(&self) -> bool {
        self.ranges.is_empty()
    }

    /// Returns lowest and highest limits of the sequence

    pub fn bounds(&self) -> Option<(T, T)> {
        match (self.ranges.first(), self.ranges.last()) {
            (Some(&Range { low, .. }), Some(&Range { high, .. })) => Some((low, high)),
            _ => None,
        }
    }
}

impl<T: Sized + PartialOrd + Copy + Debug> Default for MultiRange<T> {
    fn default() -> Self {
        Self::new()
    }
}

impl<T: Debug + PartialOrd> Debug for MultiRange<T> {
    fn fmt(&self, f: &mut Formatter) -> Result<(), Error> {
        let mut iter = self.ranges.iter();
        if let Some(r) = iter.next() {
            r.fmt(f)?;
            for r in iter {
                f.write_str(", ")?;
                r.fmt(f)?;
            }
        }
        Ok(())
    }
}

#[test]
fn or_test() {
    let mut r = MultiRange::<i32>::new();
    r.or(10, 34);

    r.or(2, 7);
    assert_eq!(r, *MultiRange::<i32>::new().or_slice(&[(2, 7), (10, 34)]));
    r.or(8, 9);
    assert_eq!(
        r,
        *MultiRange::<i32>::new().or_slice(&[(2, 7), (8, 9), (10, 34)])
    );
    r.or(7, 8);
    assert_eq!(r, *MultiRange::<i32>::new().or_slice(&[(2, 9), (10, 34)]));
    r.or(1, 33);
    assert_eq!(r, *MultiRange::<i32>::new().or_slice(&[(1, 34)]));
    let mut r = MultiRange::<i32>::new();
    r.or(1, 7);
    r.or(8, 23);
    assert_ne!(r, *MultiRange::<i32>::new().or_slice(&[(1, 23)]));
    let mut r = MultiRange::<i32>::new();
    r.or(1, 7);
    r.or(7, 23);
    assert_eq!(r, *MultiRange::<i32>::new().or_slice(&[(1, 23)]));
    r.or(1, 23);
    assert_eq!(r, *MultiRange::<i32>::new().or_slice(&[(1, 23)]));
    r.or(-5, 27);
    assert_eq!(r, *MultiRange::<i32>::new().or_slice(&[(-5, 27)]));
}
#[cfg(test)]
fn check_and<T: PartialOrd + Copy + Debug>(
    a: &MultiRange<T>,
    b: &MultiRange<T>,
    res: &MultiRange<T>,
) {
    assert_eq!(a.and_range(b), *res);
    assert_eq!(b.and_range(a), *res);
    if a.ranges.len() == 1 {
        assert_eq!(b.and(a.ranges[0].low, a.ranges[0].high), *res);
    }
    if b.ranges.len() == 1 {
        assert_eq!(a.and(b.ranges[0].low, b.ranges[0].high), *res);
    }
}

#[cfg(test)]
fn new_from_slice<T: PartialOrd + Copy + Debug>(a: &[(T, T)]) -> MultiRange<T> {
    let mut r = MultiRange::new();
    r.or_slice(a);
    r
}

#[test]
fn and_test() {
    let empty = MultiRange::<i32>::new();
    let mut r = MultiRange::<i32>::new();
    r.or_slice(&[(1, 9), (10, 34)]);
    r = r.and_range(MultiRange::<i32>::new().or_slice(&[(8, 12)]));
    assert_eq!(r, *MultiRange::<i32>::new().or_slice(&[(8, 9), (10, 12)]));

    let mut r = MultiRange::<i32>::new();
    r.or_slice(&[(8, 12)]);
    r = r.and_range(MultiRange::<i32>::new().or_slice(&[(1, 9), (10, 34)]));
    assert_eq!(r, *MultiRange::<i32>::new().or_slice(&[(8, 9), (10, 12)]));

    let mut r = MultiRange::<i32>::new();
    check_and(&r, &new_from_slice(&[]), &empty);
    check_and(&r, &new_from_slice(&[(8, 12)]), &empty);

    r.or_slice(&[(2, 5), (10, 23)]);
    check_and(&r, &new_from_slice(&[(5, 10)]), &empty);
    check_and(
        &r,
        &new_from_slice(&[(5, 11)]),
        &new_from_slice(&[(10, 11)]),
    );
    check_and(&r, &new_from_slice(&[(4, 10)]), &new_from_slice(&[(4, 5)]));
    check_and(&r, &new_from_slice(&[(0, 4)]), &new_from_slice(&[(2, 4)]));
    check_and(
        &r,
        &new_from_slice(&[(20, 40)]),
        &new_from_slice(&[(20, 23)]),
    );
    check_and(&r, &new_from_slice(&[(0, 40)]), &r);
}

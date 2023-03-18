use std::fmt::Debug;
use std::fmt::Error;
use std::fmt::Formatter;

pub trait CurveApprox {
    /// Split a segment of a curve and return the lengths.
    /// Lengths are only returned if the segment is sufficently short to
    /// estimate the length with requred precision.
    /// A suitable split point is chosen by the function, preferably one
    /// that splits the curve in approximately half the length.
    /// Returns (split point, length of lower segment, length of upper segment)
    ///
    /// # Arguments
    ///
    /// * `start` - Parameter value for start of segment
    /// * `end` - Parameter value for end of segment
    fn split_length(&self, start: f64, end: f64) -> (f64, Option<f64>, Option<f64>);
}

struct CurveNodeChildren {
    low: Box<CurveNode>,
    high: Box<CurveNode>,
}

struct CurveNode {
    t_start: f64, // Where this segment starts on the original curve
    t_end: f64,   // Where this segment ends on the original
    length: f64,  // Approximate length of this segment
    children: Option<CurveNodeChildren>,
}

fn indent(f: &mut Formatter, indent: u32) -> Result<(), Error> {
    const FILL: &str = "                    ";
    f.write_str(&FILL[0..FILL.len().min(indent as usize)])
}

fn fmt_node(node: &CurveNode, f: &mut Formatter, ind: u32) -> Result<(), Error> {
    indent(f, ind)?;
    write!(f, "{{{:?} - {:?}, ", node.t_start, node.t_end)?;
    writeln!(f, "{:?}", node.length)?;
    if let Some(children) = &node.children {
        fmt_node(&children.low, f, ind + 2)?;
        fmt_node(&children.high, f, ind + 2)?;
    }
    indent(f, ind)?;
    writeln!(f, "}}")
}

impl Debug for CurveNode {
    fn fmt(&self, f: &mut Formatter) -> Result<(), Error> {
        fmt_node(self, f, 0)
    }
}

impl CurveNode {
    fn new(start: f64, end: f64, len: f64) -> CurveNode {
        CurveNode {
            t_start: start,
            t_end: end,
            length: len,
            children: None,
        }
    }
}

pub struct LinearCurve {
    curve: Box<dyn CurveApprox>,
    node: Box<CurveNode>,
}

fn node_length(curve: &dyn CurveApprox, node: &mut Box<CurveNode>) -> f64 {
    let (low_length, high_length) = match &node.children {
        None => {
            let (split, low_len, high_len) = curve.split_length(node.t_start, node.t_end);
            node.children = Some(CurveNodeChildren {
                low: Box::new(CurveNode::new(node.t_start, split, low_len.unwrap_or(0.0))),
                high: Box::new(CurveNode::new(split, node.t_end, high_len.unwrap_or(0.0))),
            });
            (low_len, high_len)
        }
        Some(CurveNodeChildren { low, high }) => (Some(low.length), Some(high.length)),
    };
    let children = node.children.as_mut().unwrap();
    node.length = low_length.unwrap_or_else(|| node_length(curve, &mut children.low))
        + high_length.unwrap_or_else(|| node_length(curve, &mut children.high));
    node.length
}

fn node_pos(curve: &dyn CurveApprox, node: &mut Box<CurveNode>, rel_pos: f64, error: f64) -> f64 {
    //println!("{:?} - {:?} {}", node.t_start, node.t_end, rel_pos);
    if rel_pos <= error {
        return node.t_start;
    } else if node.length - rel_pos <= error {
        return node.t_end;
    }
    if node.children.is_none() {
        let (split, low_len, high_len) = curve.split_length(node.t_start, node.t_end);
        let mut children = CurveNodeChildren {
            low: Box::new(CurveNode::new(node.t_start, split, low_len.unwrap_or(0.0))),
            high: Box::new(CurveNode::new(split, node.t_end, high_len.unwrap_or(0.0))),
        };
        if low_len.is_none() {
            node_length(curve, &mut children.low);
        };
        if high_len.is_none() {
            node_length(curve, &mut children.high);
        };
        node.children = Some(children);
    };

    let children = node.children.as_mut().unwrap();
    if children.low.length > rel_pos {
        node_pos(curve, &mut children.low, rel_pos, error)
    } else {
        node_pos(
            curve,
            &mut children.high,
            rel_pos - children.low.length,
            error,
        )
    }
}

impl LinearCurve {
    pub fn new(curve: Box<dyn CurveApprox>, start: f64, end: f64) -> LinearCurve {
        let mut top_node = Box::new(CurveNode {
            t_start: start,
            t_end: end,
            length: 0.0,
            children: None,
        });
        node_length(curve.as_ref(), &mut top_node);
        LinearCurve {
            curve,
            node: top_node,
        }
    }

    /// Returns the approximated length of the curve
    pub fn length(&self) -> f64 {
        self.node.length
    }

    /// Returns the corresponding parameter given a position along the curve
    ///
    /// #Arguments
    /// * `pos`: Position along the length of the curve. Must be positive
    pub fn position(&mut self, pos: f64) -> f64 {
        node_pos(self.curve.as_ref(), &mut self.node, pos, 0.00001)
    }
}

#[cfg(test)]
struct TestCurve {
    split_factor: f64,
    max_length: f64,
}

#[cfg(test)]
impl CurveApprox for TestCurve {
    fn split_length(&self, start: f64, end: f64) -> (f64, Option<f64>, Option<f64>) {
        let low = (end - start) * self.split_factor;
        let high = (end - start) * (1.0 - self.split_factor);
        (
            start * (1.0 - self.split_factor) + end * self.split_factor,
            if low <= self.max_length {
                Some(low)
            } else {
                None
            },
            if high <= self.max_length {
                Some(high)
            } else {
                None
            },
        )
    }
}

#[test]
fn test_linear_curve_setup() {
    let curve = Box::new(TestCurve {
        split_factor: 0.7,
        max_length: 0.001,
    });
    let _linear = LinearCurve::new(curve, 0.0, 1.0);
}

#[test]
fn test_linear_curve_length() {
    let curve = Box::new(TestCurve {
        split_factor: 0.7,
        max_length: 0.001,
    });

    let linear = LinearCurve::new(curve, 0.0, 1.0);
    let len = linear.length();
    //println!("{:?}",linear.node);
    assert_relative_eq!(len, 1.0, max_relative = 0.0000001);
}

#[test]
fn test_linear_curve_pos() {
    let curve = Box::new(TestCurve {
        split_factor: 0.7,
        max_length: 0.001,
    });
    let mut linear = LinearCurve::new(curve, 0.0, 1.0);
    //println!("{:?}",linear.node);
    let pos = linear.position(1.0);
    assert_relative_eq!(pos, 1.0, max_relative = 0.00001);

    let pos = linear.position(0.0);
    assert_relative_eq!(pos, 0.0, max_relative = 0.00001);

    for p in 0..1000 {
        let p = p as f64 / 1000.0;
        let pos = linear.position(p);
        assert_abs_diff_eq!(pos, p, epsilon = 0.00001);
    }
}

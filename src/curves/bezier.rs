use super::sqrt_poly::{build_segments, eval_int, find_param, polyadd, polymul, IntegralSegment};
use coords::{Point, Vector};
use curve_approx::CurveInfo;
#[derive(Debug)]
pub struct Bezier {
    segments: Vec<IntegralSegment>,
    // First point is asumed to be (0,0)
    c1: Point,
    c2: Point,
    p2: Point,
}

impl Bezier {
    /// Create a new Bezier curve
    /// End points are (0,0) and p2
    /// c1 is relative to (0,0)
    /// c2 is relative to (0,0)
    pub fn new(c1: Point, c2: Point, p2: Point) -> Bezier {
        let cx = coefficients_from_control(&[c1.x, c2.x, p2.x]);
        let cy = coefficients_from_control(&[c1.y, c2.y, p2.y]);
        let dx = deriv_from_coefficients(&cx);
        let dy = deriv_from_coefficients(&cy);
        let dl2 = polyadd(&polymul(&dx, &dx), &polymul(&dy, &dy));
        let segments = build_segments(&dl2, 0.0, 1.0, 9);

        Bezier {
            segments,
            c1,
            c2,
            p2,
        }
    }

    /// Calculates the point and gradient for a given parameter value
    pub fn bezier_value(&self, t: f64) -> (Point, Vector) {
        let p00 = self.c1 * t;
        let p01 = self.c1 + (self.c2 - self.c1) * t;
        let p02 = self.c2 + (self.p2 - self.c2) * t;
        let p10 = p00 + (p01 - p00) * t;
        let p11 = p01 + (p02 - p01) * t;
        let p20 = p10 + (p11 - p10) * t;
        let mut d;
        if t < 0.5 {
            d = p11 - p20;
            if d.length() < 1e-9 {
                d = self.c2
            }
        } else {
            d = p20 - p10;
            if d.length() < 1e-9 {
                d = self.p2 - self.c1
            }
        }
        (p20, d.unit())
    }
}

impl CurveInfo for Bezier {
    fn length(&self) -> f64 {
        eval_int(&self.segments, 1.0)
    }

    fn value(&self, pos: f64) -> (Point, Vector) {
        let t = if pos <= 0.0 {
            0.0
        } else if pos >= self.length() {
            1.0
        } else {
            find_param(&self.segments, pos)
        };
        self.bezier_value(t)
    }
}

/// Calculate polynom coefficients.
/// Returns coefficients for the polynom c[2]*t^3 + c[1]*t^2 + c[0]*t
///
/// Start point for curve is assumed to be (0,0)
/// # Arguments
/// * `p` - [first control point, second control point, end point]

pub fn coefficients_from_control(p: &[f64; 3]) -> [f64; 3] {
    [
        3.0 * p[0],
        3.0 * (-2.0 * p[0] + p[1]),
        3.0 * (p[0] - p[1]) + p[2],
    ]
}

pub fn deriv_from_coefficients(c: &[f64; 3]) -> [f64; 3] {
    [c[0], 2.0 * c[1], 3.0 * c[2]]
}

use curve_approx::CurveInfo;
use std::cell::RefCell;
use std::cmp::Ordering;
use coords::Point;
use linear_curve::LinearCurve;
use linear_curve::CurveApprox;


pub struct BezierApprox
{
    cx: [f64;3], // coefficient for x coordinate polynom
    cy: [f64;3], // coefficient for y coordinate polynom
    splits: Vec<f64>
}

pub struct Bezier
{
    // First control point is asumed to be (0,0)
    cx: [f64;3], // coefficient for x coordinate polynom
    cy: [f64;3], // coefficient for y coordinate polynom
    dx: [f64;3], // coefficient for x derivative polynom
    dy: [f64;3], // coefficient for y derivative polynom
    linearise: RefCell<LinearCurve>
}

impl Bezier
{
    pub fn new(c1: Point, c2: Point, p2: Point) 
               -> Bezier
    {
        let cx = coefficients_from_control(&[c1.x, c2.x, p2.x]);
        let cy = coefficients_from_control(&[c1.y,  c2.y, p2.y]);
        let dx = deriv_from_coefficients(&cx);
        let dy = deriv_from_coefficients(&cy);
        let splits = find_extremes(&dx,&dy);
         println!("Splits: {:?}", splits);
        let approx = Box::new(BezierApprox{cx,cy,splits});
        Bezier{cx, cy,
               dx, dy,
               linearise: RefCell::new(LinearCurve::new(approx,0.0,1.0))}
        
    }
}

impl CurveApprox for BezierApprox
{
     fn split_length(&self, start: f64, end: f64) 
                     -> (f64, Option<f64>, Option<f64>)
    {
        let split = match self.splits.iter().find(|&&p| p < end && p >= start)
        {
            Some(&p) => p,
            None => (start + end) / 2.0
        };
        println!("Split: {}", split);
        let p_start = Point{x: bezier_poly_value(&self.cx, start),
                            y: bezier_poly_value(&self.cy, start)};
        let p_spilt = Point{x: bezier_poly_value(&self.cx, split),
                            y: bezier_poly_value(&self.cy, split)};
        let p_end = Point{x: bezier_poly_value(&self.cx, end),
                            y: bezier_poly_value(&self.cy, end)};
        (split,Some(1.0),Some(3.0))
    }
}

impl CurveInfo for Bezier
{
    
    
    fn length(&self) -> f64
    {
        // TODO
        self.linearise.borrow().length()
    }
    
    fn value(&self, pos: f64) -> ([f64;2], [f64;2])
    {
        let t = pos;
        ([bezier_poly_value(&self.cx,t), bezier_poly_value(&self.cy,t)],
         [bezier_deriv(&self.cx,t), bezier_deriv(&self.cy,t)])
    }
}

/// Calculates value of bezier polynom
/// Returns c[0]*t^3 + c[1]*t^2 + c[2]*t
pub fn bezier_poly_value(c: &[f64; 3], t: f64) -> f64
{
    ((((c[0] * t) + c[1]) * t) + c[2]) * t
}
        
/// Calculates value of derivative polynom
/// Returns d[0]*t^2 + d[1]*t + d[2]
pub fn bezier_deriv(d: &[f64; 3], t: f64) -> f64
        {
    (d[0] * t + d[1]) * t + d[2]
}
      
/// Calculate polynom coefficients.
/// Returns coefficients for the polynom c[0]*t^3 + c[1]*t^2 + c[2]*t
///
/// Start point for curve is assumed to be (0,0)
/// # Arguments
/// * `p` - [first control point, second control point, end point]
        
pub fn coefficients_from_control(p: &[f64;3]) -> [f64;3]
{
    [
        3.0 * (p[0] - p[1]) + p[2],
        3.0 * (-2.0*p[0] + p[1]),
        3.0 * p[0]
    ]
}

pub fn deriv_from_coefficients(c: &[f64;3]) -> [f64;3]
{
    [
        3.0*c[0],
        2.0*c[1],
        c[2]
    ]
}

#[derive(Clone)]
pub struct BezierSegment
{
    t_start: f64, // Where this segment starts on the original curve
    t_end: f64, // Where this segment endss on the original curve
    px: [f64;4], // x coordinates for control points
    py: [f64;4], // y coordinates for control points
    length: f64 // Approximate length of this segment
}

impl BezierSegment
{
    pub fn as_svg_path(&self) -> String
    {
        format!("M {},{} C {},{} {},{} {},{}",
                self.px[0], self.py[0], self.px[1], self.py[1], 
                self.px[2], self.py[2],  self.px[3], self.py[3])
    }
}

pub trait AsSvgPath
{
    fn as_svg_path(&self) -> String;
}

impl AsSvgPath for Vec<BezierSegment>
{
    fn as_svg_path(&self) -> String
    {
        let mut path = String::new();
        let p :Option<(f64,f64)> = None;
        for s in self.iter() {
            if p != Some((s.px[0], s.py[0])) {
                path += &format!("M {},{} ", s.px[0], s.py[0]);
            }
            path += &format!("C {},{} {},{} {},{} ",
                             s.px[1], s.py[1], s.px[2], s.py[2],
                             s.px[3], s.py[3]);
        }
        return path
    }
}

#[derive(Debug, PartialEq)]
enum Solution2ndDegree
{
    None,
    One(f64),
    Two(f64,f64)
}
/// Find solutions for c[0]*x^2 + c[1]*x + c[2] = 0

fn solve_2nd_degree(c: &[f64;3]) ->Solution2ndDegree
{
    if c[0] == 0.0 {
        if c[1] == 0.0 {return Solution2ndDegree::None}
        return Solution2ndDegree::One(-c[2] / c[1]);
    }
    let p = c[1] / (2.0 * c[0]);
    let q = c[2] / c[0];
    let w = p*p - q;
    if w < 0.0 {return Solution2ndDegree::None}
    if w == 0.0 {return Solution2ndDegree::One(-p)}
    let sqrt_w = w.sqrt();
    Solution2ndDegree::Two(-p + sqrt_w, -p - sqrt_w)
}

fn find_extremes(dx:&[f64;3], dy: &[f64;3]) -> Vec<f64>
{
    let mut split = Vec::<f64>::new();
    // Split at extreme points
    println!("dx= {} t^2 + {} t + {}", dx[0], dx[1], dx[2]);
    println!("dy = {} t^2 + {} t + {}", dy[0], dy[1], dy[2]);
    match solve_2nd_degree(&dx) {
        Solution2ndDegree::None => {},
        Solution2ndDegree::One(x1) => {
            if x1 > 0.0 && x1 < 1.0 {
                split.push(x1)
            }
        },
        Solution2ndDegree::Two(x1,x2) => {
            if x1 > 0.0 && x1 < 1.0 {
                split.push(x1);
            }
            if x2 > 0.0 && x2 < 1.0 {
                split.push(x2)
            }
        }
    }
    match solve_2nd_degree(&dy) {
        Solution2ndDegree::None => {},
        Solution2ndDegree::One(x1) => {
            if x1 > 0.0 && x1 < 1.0 {
                split.push(x1)
            }
        },
        Solution2ndDegree::Two(x1,x2) => {
            if x1 > 0.0 && x1 < 1.0 {
                split.push(x1);
            }
            if x2 > 0.0 && x2 < 1.0 {
                split.push(x2)
            }
        }
    }
    split.sort_by(|a,b| a.partial_cmp(b).unwrap());
    split
}

/* Tests */


#[cfg(test)]
fn check_solution(c:&[f64;3])
{
    let s = solve_2nd_degree(c);
    //println!("{:?}", s);
    match s {
        Solution2ndDegree::None => {
            if c[0] == 0.0 {
                if c[1] != 0.0 {
                    panic!("Linear equation has solution");
                }
            } else {
                let p = c[1] / (2.0*c[0]);
                if p*p - c[2] / c[0] >= 0.0 {
                    panic!("Quadratic equation has solution");
                }
            }
        },
        Solution2ndDegree::One(r) => {
            if bezier_deriv(c, r) >= 1e-10 {
                panic!("Single solution gives non-zero result: {}",r);
            }
        },
        Solution2ndDegree::Two(a,b) => {
            if bezier_deriv(c, a) >= 1e-10 {
                panic!("First solution gives non-zero result");
            }
            if bezier_deriv(c, b) >= 1e-10 {
                panic!("Second solution gives non-zero result");
            }
        }
        
    }
}
        
                
#[test]
fn test_solve_2nd_degree()
{
    check_solution(&[0.0, 3.2, 8.3]);
    
    check_solution(&[23.0, 3.2, -8.3]);
    check_solution(&[23.0, 0.0,-8.3]);
    check_solution(&[23.0, 0.0, 0.0]);
}

#[test]
fn test_bezier()
{
    let c = Bezier::new(Point{x: 1.0, y: 1.0},
                        Point{x: 2.0, y: -1.0}, 
                        Point{x: 3.0, y: 1.0});
    println!("Length: {}", c.length());
}

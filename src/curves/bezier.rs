use curve_approx::CurveInfo;
use coords::{Point,Vector};




#[derive(Debug)]
pub struct Bezier
{
    // First point is asumed to be (0,0)
    cx: [f64;3], // coefficient for x coordinate polynom
    cy: [f64;3], // coefficient for y coordinate polynom
    dx: [f64;3], // coefficient for x derivative polynom
    dy: [f64;3], // coefficient for y derivative polynom
    segments: Vec<BezierSegment>
}

impl Bezier
{
    pub fn new(c1: Point, c2: Point, p2: Point) 
               -> Bezier
    {
        assert_ne!(c1.length(), 0.0);
        assert_ne!(c2.length(), 0.0);
        let cx = coefficients_from_control(&[c1.x, c2.x, p2.x]);
        let cy = coefficients_from_control(&[c1.y,  c2.y, p2.y]);
        let dx = deriv_from_coefficients(&cx);
        let dy = deriv_from_coefficients(&cy);
        let splits = find_extremes(&dx,&dy);
        println!("Splits: {:?}", splits);
        let mut c1r = c1;
        let mut c2r = c2;
        let mut p2r = p2;
        rotate_segment(&mut c1r, &mut c2r, &mut p2r);
        let mut seg = BezierSegment {t_end: 1.0,
                                 p: p2r.x,
                                 c1: c1r,
                                     c2: c2r};
        let mut curved = Vec::new();
        let mut prev_split = 0.0;
        for s in splits {
            let split_pos = (s - prev_split) / (seg.t_end - prev_split);
            println!("Split pos: {}" , split_pos);
            let p1 = split_segment(&mut seg, 
                                   split_pos,
                                   prev_split);
            prev_split = split_pos;
            curved.push(p1);
        }
        curved.push(seg);
        curved.reverse();
        let mut segments = flatten_segments(&mut curved, 1.0001);
        let mut sum = 0f64;
        for s in &mut segments {
            s.p += sum;
            sum = s.p;
        }
        Bezier{cx, cy,
               dx, dy,
               segments
        }
        
    }

    /// Calculates the point and gradient for a given parameter value
    pub fn bezier_value(&self, t: f64) -> (Point, Vector)
    {
        let deriv  = Vector{x: bezier_deriv(&self.dx,t),
                            y: bezier_deriv(&self.dy,t)};
        let scale = 1.0 / deriv.length();
        (Point{x: bezier_poly_value(&self.cx,t), 
               y: bezier_poly_value(&self.cy,t)},
         deriv * scale)
    }
}

#[derive(Debug)]
struct BezierSegment
{
    t_end: f64, // parameter value at end of this segment
    p: f64, /* x coordinate of point, either relative to prevoius
    point or start of whole curve, y is always 0 */
    c1: Point, // first control point
    c2: Point // second point
}
    
/// Splits one coordinate of a bezier curve at `pos`
/// First point is assumed to be (0,0)
///
/// Returns (c11,c12, p_split, c21, c22)
fn split_bezier_coord(c1: f64, c2: f64, p2: f64, pos:f64)
                      -> (f64,f64,f64,f64,f64)
{
    let q0 = c1 * pos;
    let q1 = c1 + (c2 - c1) * pos;
    let q2 = c2 + (p2 - c2) * pos;
    let r0 = q0 + (q1 - q0) * pos;
    let r1 = q1 + (q2 - q1) * pos;
    let p_split = r0 + (r1 - r0) * pos;
    (q0, r0, p_split, r1, q2)
}

fn rotate_point(rs: f64, rc: f64, p: &mut Point)
{
    let x = p.x * rc - p.y * rs;
    let y = p.x * rs + p.y * rc;
    p.x = x;
    p.y = y;
}

fn rotate_segment(c1: &mut Point, c2: &mut Point, p2: & mut Point)
{
    let dist =  (p2.x * p2.x + p2.y * p2.y).sqrt();
    // Rotate segment so that p2.y is 0
    let rs = -p2.y / dist; // Sine of rotation
    let rc = p2.x / dist; // Cosine of rotation
    p2.x = dist;
    p2.y = 0.0;
    rotate_point(rs,rc, c1);
    rotate_point(rs,rc, c2);
}

        
fn split_segment(p2: &mut BezierSegment, pos: f64, t_start: f64) 
                 -> BezierSegment
{
    // split x
    let (c11x, c12x, p1x, c21x, c22x) = 
        split_bezier_coord(p2.c1.x, p2.c2.x, p2.p, pos);
    // split y
    let (c11y, c12y, p1y, c21y, c22y) = 
        split_bezier_coord(p2.c1.y, p2.c2.y, 0.0, pos);
    p2.c1 = Point{x: c21x-p1x, y: c21y-p1y};
    p2.c2 = Point{x: c22x-p1x, y: c22y-p1y};
    let mut p2p = Point{x: p2.p - p1x, y:  -p1y};
    rotate_segment(&mut p2.c1, &mut p2.c2, &mut p2p);
    p2.p = p2p.x;
    
    let mut p1 = BezierSegment{
        t_end: (p2.t_end - t_start) * pos + t_start,
        p: 0.0,
        c1: Point{x: c11x, y: c11y},
        c2: Point{x: c12x, y: c12y}
    };
    let mut p1p = Point{x: p1x, y: p1y};
    rotate_segment(&mut p1.c1, &mut p1.c2, &mut p1p);
    p1.p = p1p.x;
    p1
}

fn segment_is_flat(seg: &BezierSegment, precision: f64) -> bool
{
    let max = seg.c1.length() 
        + (seg.c2 - seg.c1).length() 
        + (Point{x:seg.p, y: 0.0} - seg.c2).length();
    max <= seg.p * precision
}

fn flatten_segments(curved: &mut Vec<BezierSegment>, precision: f64)
                    -> Vec<BezierSegment>
{
    let mut flat = Vec::<BezierSegment>::new();
    let mut t_start = 0f64;
    while let Some(mut s) = curved.pop() {
        if segment_is_flat(&s,precision) {
            t_start = s.t_end;
            flat.push(s);
        } else {
            let s0 = split_segment(&mut s, 0.5, t_start);
            curved.push(s);
            curved.push(s0);
        }
    }
    flat
}
        

impl CurveInfo for Bezier
{
    
    
    fn length(&self) -> f64
    {
       self.segments.last().map_or(0.0, |s| s.p)
    }
    
    fn value(&self, pos: f64) -> (Point, Vector)
    {
        let t = 
            if pos <= 0.0 {
                0.0
            } else if pos >= self.length() {
                1.0
            } else {
                match self.segments.binary_search_by(|s| {
                    s.p.partial_cmp(&pos).unwrap()}) {
                    Ok(i) => self.segments[i].t_end,
                    Err(i) => {
                        let (t_start, p_start) = if i >= 1 {
                            let s = &self.segments[i-1];
                            (s.t_end, s.p)
                        } else {
                            (0.0, 0.0)
                        };
                        let s = &self.segments[i];
                        assert!(pos <= s.p && pos >= p_start);
                        let p = s.p - p_start;
                        let a3 = 3.0 * (s.c1.x - s.c2.x) + p;
                        let a2 = -6.0 * s.c1.x + 3.0 * s.c2.x;
                        let a1 = 3.0 *  s.c1.x;
                        let a0 = -(pos - p_start);

                        /* The curve should have no extremes otherwise 
                        the mapping won't be 1:1. This means the derivative
                        has no zeros. */
                        match solve_2nd_degree(&[3.0 * a3, 2.0 * a2, a1]) {
                            Solution2ndDegree::None => {},
                            Solution2ndDegree::One(x1) 
                                if x1 > 1.0 || x1 < 0.0 => {},
                             Solution2ndDegree::Two(x1,x2) 
                                if (x1 > 1.0 || x1 < 0.0) 
                                && (x2 > 1.0 || x2 < 0.0) => {},
                            
                            _ => panic!("The mapping from position to parameter is not 1:1")
                        }
                        let ts = solve_3rd_degree(&[a3,a2,a1,a0], 0.5);
                        if ts < -0.001 || ts > 1.001 {
                            panic!("ts: {}", ts);
                        }

                        t_start + (s.t_end - t_start) * ts
                    }
                }
            };
        self.bezier_value(t)
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



fn polynom_value(coeffs: &[f64], t: f64) -> f64
{
    let mut sum = coeffs[0];
    for c in &coeffs[1..] {
        sum = sum * t + c;
    }
    sum
}

/// Find solutions for c[0]*x^3 + c[1]*x^2 + c[2]*x + c[3] = 0
/// This function assumes the equation has exactly one solution
fn solve_3rd_degree(c: &[f64;4], mut t: f64) -> f64
{
    let d = &[3.0 * c[0], 2.0 * c[1], c[2]];
    loop {
        let v = polynom_value(c, t);
        if v.abs() < 0.00001 {
            break;
        }
        let k = polynom_value(d, t);
        t -= v / k;
    }
    t
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
fn test_solve_3rd_degree()
{
    let x1 = solve_3rd_degree(&[2.0, 3.0, 3.0, -2.0],0.0);
    assert_relative_eq!(x1, 0.429445, max_relative=0.0001);
}



#[test]
fn test_rotate()
{
    let mut c1=Point{x: 0.0, y: 1.0};
    let mut c2=Point{x: 2.0, y: 3.0};
    let mut p2= Point{x: 3.0, y: 3.0};
    rotate_segment(&mut c1, &mut c2, &mut p2);
    assert_relative_eq!(c1.x, 1.0/2.0f64.sqrt());
    assert_relative_eq!(c1.y, 1.0/2.0f64.sqrt());
    
    assert_relative_eq!(c2.x, 3.0*2.0f64.sqrt() - 1.0/2.0f64.sqrt());
    assert_relative_eq!(c2.y, 1.0/2.0f64.sqrt());
    assert_relative_eq!(p2.x, 3.0*2.0f64.sqrt());
    assert_eq!(p2.y, 0.0);
}

#[test]
fn test_split_segment()
{
    let mut c1=Point{x: 0.0, y: 1.0};
    let mut c2=Point{x: 2.0, y: 3.0};
    let mut p2= Point{x: 3.0, y: 3.0};
    rotate_segment(&mut c1, &mut c2, &mut p2);
    let mut s2 = BezierSegment {
        t_end: 1.0,
        p: p2.x,
        c1: c1,
        c2: c2
    };
    let s1 = split_segment(&mut s2, 0.5, 0.0);
    assert_eq!(s1.t_end, 0.5);
    assert_eq!(s2.t_end, 1.0);
    assert_relative_eq!(s1.p, s2.p, max_relative=0.0001);
    assert_relative_eq!(s1.c1.x, s2.p - s2.c2.x, max_relative=0.0001);
    assert_relative_eq!(s1.c1.y, s2.c2.y, max_relative=0.0001);
    assert_relative_eq!(s1.p - s1.c2.x, s2.c1.x, max_relative=0.0001);
    assert_relative_eq!(s1.c2.y, s2.c1.y, max_relative=0.0001);
}
        

#[test]
fn test_bezier()
{
    let c = Bezier::new(Point{x: 1.0, y: 1.0},
                        Point{x: 2.0, y: -1.0}, 
                        Point{x: 3.0, y: 1.0});
    let len = c.length();
    println!("Length: {}", len);

    for i in 0..=10 {
        let p = f64::from(i) * len / 10.0;
        let (c, _) = c.value(p);
        println!("{} => ({}, {})", p, c.x,c.y);
    }
        
}
#[test]
fn test_bezier_line()
{
    let c = Bezier::new(Point{x: 1.0, y: 1.0},
                        Point{x: 2.0, y: 2.0}, 
                        Point{x: 10.0, y: 10.0});
    let len = c.length();
    assert_relative_eq!(len, 10.0 * 2.0f64.sqrt());
    println!("Length: {}", len);

    for i in 0..=10 {
        let p = f64::from(i) * len / 10.0;
        let (c, _) = c.value(p);
        assert_relative_eq!(p, c.length(),max_relative=0.0001);
        println!("{} => ({}, {})", p, c.x,c.y);
    }
        
}

#[test]
fn test_bezier_linearity()
{
    let p2 = Point{x: 76.444165, y: 45.198655};
    let c1 = Point{x: 34.385054, y: 2.4875446};
    let c2 = Point{ x: -36.419844, y: 10.041719 } + p2;
    let c = Bezier::new(c1,c2,p2);
    println!("{:?}", c);
    let len = c.length();
    let mut p = 0.0;
    let mut prev = Point{x: 0.0, y: 0.0};
    while p <= len - 0.5 {
        p += 0.5;
        let (pos, _) = c.value(p);
        println!("{}: {}", p, pos); 
        assert_relative_eq!((prev - pos).length(), 0.5, max_relative=0.01);
        prev = pos;
    }
}

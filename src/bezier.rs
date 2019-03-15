
pub fn bezier_value(c: &[f64; 4], t: f64) -> f64
{
    let mut p = c[0];
    for k in &c[1..] {
        p = (p * t) + k;
    }
    p
}

pub fn bezier_deriv(d: &[f64; 3], t: f64) -> f64
{
    let mut p = d[0];
    for k in &d[1..] {
        p = (p * t) + k;
    }
    p
}
      
        
pub fn coefficients_from_control(p: &[f64;4]) -> [f64;4]
{
    [
        -p[0]+ 3.0*(p[1] - p[2]) + p[3],
        3.0*(p[0] - 2.0*p[1] + p[2]),
        3.0*(p[1] - p[0]),
        p[0]
    ]
}

pub fn deriv_from_coefficients(c: &[f64;4]) -> [f64;3]
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

fn vlen(x: f64, y: f64) -> f64 
{
    (x * x + y * y).sqrt()
}

fn find_arc_distance(px: &[f64;4], py: &[f64;4], dist: f64, err: f64,
                     t_start:f64, t_end:f64) -> f64
{
    let mut d = dist;
    let mut start = t_start;
    let mut end = t_end;
    let mut x : [f64;4] = px.clone();
    let mut y : [f64;4] = py.clone();
    while d > err {
        // Split curve in half
        let c0x = (x[0] + x[1]) * 0.5;
        let c0y = (y[0] + y[1]) * 0.5;
        let c1x = (x[1] + x[2]) * 0.5;
        let c1y = (y[1] + y[2]) * 0.5;
        let c2x = (x[2] + x[3]) * 0.5;
        let c2y = (y[2] + y[3]) * 0.5;
        
        let d0x = (c0x + c1x) * 0.5;
        let d0y = (c0y + c1y) * 0.5;
        let d1x = (c1x + c2x) * 0.5;
        let d1y = (c1y + c2y) * 0.5;
        
        let e0x = (d0x + d1x) * 0.5;
        let e0y = (d0y + d1y) * 0.5;


        let dlow = vlen(e0x-x[0], e0y- y[0]);
        if d < dlow {
            // Select lower half
            x[1] = c0x;
            y[1] = c0y;
            x[2] = d0x;
            y[2] = d0y;
            x[3] = e0x;
            y[3] = e0y;
            end = (start + end) / 2.0;
        } else {
            x[0] = e0x;
            y[0] = e0y;
            x[1] = d1x;
            y[1] = d1y;
            x[2] = c2x;
            y[2] = c2y;
            start = (start + end) / 2.0;
            d -= dlow;
        }
    }
    return start;
}

pub fn find_arc_distance_vec(segments: &Vec<BezierSegment>, dist: f64, err: f64) -> f64
{
    let mut d = dist;
    for s in segments {
        if s.length <= d {
            return find_arc_distance(&s.px, &s.py, d, err, 
                                     s.t_start, s.t_end);
        }
        d -= s.length;
    }
    return 1.0;
}

pub fn find_equidistant_steps(px: &[f64;4], py: &[f64;4], step: f64, err: f64) 
                              -> (f64, Vec<f64>)
{
    let mut tvec = Vec::<f64>::new();
    let mut segs = Vec::<BezierSegment>::new();
    let l = bezier_subdiv(px,py, 0.0, 1.0, err, &mut segs);
    println!("total length = {}", l);
    let n = (l / step).ceil() as usize;
    let step = l / n as f64;
    let mut seg_iter = segs.iter();
    let mut seg = seg_iter.next().unwrap();
    let mut l = 0.0;
    for i in 1..n {
        let s = i as f64 * step;
        //println!("s={}, l = {}", s, l);
        while l + seg.length < s {
            seg = match seg_iter.next() {
                Some(s) => {
                    l += seg.length;
                    s
                },
                None => break
            }
        }
        let t = find_arc_distance(&seg.px, &seg.py, s - l, err, 
                                  seg.t_start, seg.t_end);
        tvec.push(t);
    }
    tvec.push(1.0); // Last parameter value is always 1, no need to calculate it
    return (step, tvec)
}

/*
Split curve into sufficiently straight (as determined by err) segments
 */

pub fn bezier_subdiv(px: &[f64;4], py: &[f64;4], t_start:f64, t_end:f64, err: f64, segments: &mut Vec<BezierSegment>) -> f64
{
    //println!("{} - {}: {:?} {:?}", t_start, t_end, px, py);
    /* Check if end directions are the same, within 90degrees, as that
    of the line between end points. */

    let dx = px[3]-px[0];
    let dy = py[3]-py[0];

    if dx.abs() < err && dy.abs() < err
        && (px[1]-px[0]).abs() < err && (py[1]-py[0]).abs() < err 
        && (px[2]-px[3]).abs() < err && (py[2]-py[3]).abs() < err
    {
        let l = (dx*dx + dy*dy).sqrt();
        //println!("l = {}", l);
        segments.push(BezierSegment {t_start: t_start, t_end: t_end,
                                     px: px.clone(), py: py.clone(),
                                     length: l});
        return l;
    }
        
    /* Direction at start */
    let dstart = (px[1]-px[0]) * dx + (py[1]-py[0]) * dy;
    /* Direction at start */

    let dend = (px[2]-px[3]) * -dx + (py[2]-py[3]) * -dy;
    if dstart >= 0.0 && dend >= 0.0 {
        let xy = (px[2]-px[3]) * dy - (py[2]-py[3]) * dx;
        let l = (dx*dx + dy*dy).sqrt();
        let d = xy.abs() / l;
        if d <= err {
            //println!("d = {}", d);
            segments.push(BezierSegment {t_start: t_start, t_end: t_end,
                                         px: px.clone(), py: py.clone(),
                                         length: l});
            return l;
        }
    }
    
    // Split curve in half
    let c0x = (px[0] + px[1]) * 0.5;
    let c0y = (py[0] + py[1]) * 0.5;
    let c1x = (px[1] + px[2]) * 0.5;
    let c1y = (py[1] + py[2]) * 0.5;
    let c2x = (px[2] + px[3]) * 0.5;
    let c2y = (py[2] + py[3]) * 0.5;

    let d0x = (c0x + c1x) * 0.5;
    let d0y = (c0y + c1y) * 0.5;
    let d1x = (c1x + c2x) * 0.5;
    let d1y = (c1y + c2y) * 0.5;

    let e0x = (d0x + d1x) * 0.5;
    let e0y = (d0y + d1y) * 0.5;

    let p1x: [f64;4] = [px[0], c0x, d0x, e0x];
    let p1y: [f64;4] = [py[0], c0y, d0y, e0y];
    
    let p2x: [f64;4] = [e0x, d1x, c2x, px[3]];
    let p2y: [f64;4] = [e0y, d1y, c2y, py[3]];


    let t_mid = (t_start + t_end) / 2.0;
    bezier_subdiv(&p1x, &p1y, t_start, t_mid, err, segments) 
        + bezier_subdiv(&p2x, &p2y, t_mid, t_end, err, segments)
}

const ERR_WEIGHT:f64 = 0.2;
fn remove_int_part(a: &mut f64) -> i64
{
    let i = a.round();
    *a -= i;
    i as i64
}

fn limit_acc(a: i64, v:i64, amax: i16, vmax:i32) -> i64
{
    let mut a_adj = a;
    if a_adj > amax as i64 {
        a_adj = amax as i64;
    } else if a_adj < -amax as i64 {
        a_adj = -amax as i64;
    }
    if v + a_adj > vmax as i64 {
        a_adj = vmax as i64 - v;
    } else if  v + a_adj < -vmax as i64 {
        a_adj = -vmax as i64 - v;
    }
    a_adj
}
        
pub fn constant_speed<F>(px: &[f64;4], py: &[f64;4], v:f64, uvec: &[f64], 
                         vx:i32, vy:i32, vmax:i32, amax:i16, cb:&mut F)
                         -> (i64, i64)
    where F: FnMut(i16, i16)
{
    let mut int_x = px[0] as i64;
    let mut int_y = py[0] as i64;
    let px = [px[0]*0.5, px[1]*0.5, px[2]*0.5, px[3]*0.5];
    let py = [py[0]*0.5, py[1]*0.5, py[2]*0.5, py[3]*0.5];
    let cx = coefficients_from_control(&px);
    let cy = coefficients_from_control(&py);
    let cdx = deriv_from_coefficients(&cx);
    let cdy = deriv_from_coefficients(&cy);
    
    let mut int_vx = vx as i64;
    let mut int_vy = vy as i64;
    println!("int_vx: {}, int_vy: {}, v: {}", int_vx, int_vy,v);
    let mut err_x = 0f64;
    let mut err_y = 0f64;
    for &u in uvec {
        let x = bezier_value(&cx, u);
        let y = bezier_value(&cy, u);
        let dx = bezier_deriv(&cdx, u);
        let dy = bezier_deriv(&cdy, u);
        let mut dl = vlen(dx, dy);
        if dl < 0.1 {
            dl = 0.1;
        }
        let v_scale = v / dl;

        let int_vx1 = (dx*v_scale).round() as i64;
        let int_vy1 = (dy*v_scale).round() as i64;
        let int_ax =  limit_acc(int_vx1 - int_vx - remove_int_part(&mut err_x),
                                int_vx, amax, vmax);
        let int_ay = limit_acc(int_vy1 - int_vy - remove_int_part(&mut err_y),
                               int_vy, amax, vmax);
        cb(int_ax as i16, int_ay as i16);
        //println!("int_vx: {}, int_vy: {}", int_vx, int_vy);
        //println!("int_ax: {}, int_ay: {}", int_ax, int_ay);
        int_x += int_ax + 2*int_vx;
        int_y += int_ay + 2*int_vy;
        int_vx += int_ax;
        int_vy += int_ay;
        err_x += (int_x as f64 - 2.0*x) * ERR_WEIGHT;
        err_y += (int_y as f64 - 2.0*y) * ERR_WEIGHT;
        //println!("{:>5.3} ({:>6.2}, {:>6.2}) ({:>6}, {:>6})", u, 2.0*x , 2.0*y, int_x, int_y); 
        //println!("vx: {:>6.2}, vy: {:>6.2}", dx, dy);
    }
    println!("int_vx: {}, int_vy: {}", int_vx, int_vy);
    (int_x, int_y)
}

#[cfg(test)]
    fn assert_close(err:f64, got:f64, expected:f64)
{
    println!("{} == {}", got, expected);
    assert!(((expected-got)/expected).abs() <= err);
}
#[test]
fn test_find_equidistant_steps()
{
    let px = [0.28525391, 1.1596191, 0.8516276, 0.50642903];
    let py = [3.9373211, 2.7115563, 4.8054852, 2.6144045];
    let t = find_equidistant_steps(&px,&py, 0.2, 0.001);
    println!("t= {:?}", t);
}

#[test]
fn test_find_arc_distance()
{
    let px = [1.0, 2.0, 3.0, 4.0];
    let py = [1.0, 2.0, 3.0, 4.0];
    let t = find_arc_distance(&px,&py, 2.0, 0.001, 0.0, 1.0);
    assert_close(0.01, t, 2f64.sqrt() / 3.0);
    let px = [1.0, 3.0, 3.0, 4.0];
    let py = [1.0, 3.1, 3.0, 4.0];
    let t = find_arc_distance(&px,&py, 2.0, 0.001, 0.0, 1.0);
    println!("t={}", t);
    let s = vlen(bezier_value(&coefficients_from_control(&px),t) - px[0],
                 bezier_value(&coefficients_from_control(&py),t) - py[0]);
    assert_close(0.001, s, 2.0);
}

#[test]
fn test_bezier_subdiv()
{
    let mut segs = Vec::<BezierSegment>::new();
    let px = [1.0, 1.0, 2.0, 2.0];
    let py = [1.0, 2.0, 2.0, 1.0];
    let l = bezier_subdiv(&px,&py, 0.0,1.0, 0.01, &mut segs);
    println!("Length: {}\nPath: {}", l, segs.as_svg_path());
    
    let px = [1.0, 0.0, 1.5, 2.0];
    let py = [3.0, 4.0, 5.0, 3.0];
    segs.clear();
    let l = bezier_subdiv(&px,&py, 0.0, 1.0, 0.01, &mut segs);
    println!("Length: {}\nPath: {}", l, segs.as_svg_path());

    let px = [0.28525391, 1.1596191, 0.8516276, 0.50642903];
    let py = [3.9373211, 2.7115563, 4.8054852, 2.6144045];
    segs.clear();
    let l = bezier_subdiv(&px,&py, 0.0, 1.0, 0.01, &mut segs);
    println!("Length: {}\nPath: {}", l, segs.as_svg_path());
}

#[test]
fn test_constant_speed()
{
    let px = [0f64, 30001f64, 18882f64, 30000f64];
    let py = [0f64, 288343f64, 19102f64, -62777f64];
    /*
    let px = [0f64, 10000f64, 20000f64, 30000f64];
    let py = [0f64, -30000f64, -60000f64, -90000f64];
     */
    let v = 3000f64;
    let (v_adj, uvec) = find_equidistant_steps(&px,&py, v, 0.1); 
    constant_speed(&px,&py, v_adj, &uvec, 0,0, 5000, 5000, &mut |x,y| println!("a: ({}, {})", x, y));
}

#[test]
fn test_curve()
{
    let p = [3f64, 2.5, -1f64, 0f64];
    let c = coefficients_from_control(&p);
    let d = deriv_from_coefficients(&c);
    let mut t = 0f64;
    println!("");
    while t < 1f64 - 0.0001 {
        println!("{:>5.2} {:>6.2} {:>6.2}", t, bezier_value(&c, t), bezier_deriv(&d,t));
        t += 0.01;
    }
    t = 1f64;
    println!("{:>5.2} {:>6.2} {:>6.2}", t, bezier_value(&c, t), bezier_deriv(&d,t));
}

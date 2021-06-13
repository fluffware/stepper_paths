use coords::Transform;
use coords::Point;
use std::f64::consts::PI;
use nom::Finish;
use stepper_context::CurveSegment;
use std::error::Error;

mod parser {
    use std::str::FromStr;
    use std::f64::consts::PI;
    use coords::Point;
    use coords::Transform;
    type Input<'a> = &'a str;
    use nom::character::complete::{
        multispace0,
        multispace1,
        digit0,
        digit1,
        self
    };

    named!(wsp_opt<Input, Input >, call!(multispace0));
    named!(comma_wsp<Input, Input >, 
           alt!(recognize!(tuple!(wsp_opt, 
                                  call!(complete::char(',')),
                                  wsp_opt)) 
                | multispace1)); 
    named!(sign_opt<Input, f64>, map!(opt!(alt!(call!(complete::char('+')) 
                                                | call!(complete::char('-')))), 
                                      |s| {
                                          match s {
                                              Some('-') => -1.0,
                                             Some(_) | None => 1.0
                                          }
                                      }));
    named!(fractional_constant<Input, Input>,
           alt!(recognize!(tuple!(digit0, 
                                  call!(complete::char('.')), 
                                  digit1))
                | recognize!(tuple!(digit1, 
                                    call!(complete::char('.'))))
           ));
    named!(floating_point_constant<Input, Input>,
           alt!(recognize!(tuple!(fractional_constant, opt!(exponent))) 
                | recognize!(tuple!(digit1, exponent))));

    named!(exponent<Input, Input>,
           recognize!(tuple!(alt!(call!(complete::char('e')) 
                                  | call!(complete::char('E'))),
                             sign_opt, digit1)));

    named!(floating_point_signed<Input, f64>, 
           map_res!(recognize!(tuple!(sign_opt, floating_point_constant)),
                    |s:Input| {f64::from_str(s)})
           );
    named!(integer<Input, f64>, map_res!(recognize!(tuple!(opt!(alt!(char!('+') | char!('-'))), digit1)), |s:Input| {f64::from_str(s)}));

    named!(number<Input, f64>, alt!(floating_point_signed | integer)); 
    named!(matrix<Input,Transform>,
           delimited!(tuple!(tag!("matrix"),wsp_opt,char!('('),wsp_opt),
                      map!(tuple!(number, comma_wsp, number, comma_wsp, 
                                  number, comma_wsp, number, comma_wsp,
                                  number, comma_wsp, number),
                           |(a,_,b,_,c,_,d,_,e,_,f)| 
                           Transform{matrix:[a,b,c,d,e,f]}),
                      tuple!(wsp_opt,char!(')'))));
    
    named!(translate<Input,Transform>, 
           delimited!(tuple!(tag!("translate"),wsp_opt,char!('('),wsp_opt), 
                      map!(tuple!(number, 
                                  opt!(preceded!(comma_wsp,
                                                 number))),
                                  |(x,y)| {
                                      let y = y.unwrap_or(0.0);
                                      Transform{matrix:[1f64,0f64,
                                                        0f64,1f64,
                                                        x,y]}}),
                      tuple!(wsp_opt,char!(')'))));
    named!(scale<Input,Transform>,  
           delimited!(tuple!(tag!("scale"),wsp_opt,char!('('),wsp_opt),
                      map!(tuple!(number, 
                                  opt!(preceded!(comma_wsp,
                                                 number))),
                           |(x,y)| {
                               let y = match y {
                                   Some(v) => v,
                                   None => x
                               };
                               Transform{matrix:[x,0f64,0f64,y,0.0,0.0]}}),
                      tuple!(wsp_opt,char!(')'))));

    named!(rotate<Input,Transform>,  
           delimited!(tuple!(tag!("rotate"),wsp_opt,char!('('),wsp_opt),
                      map!(tuple!(number, 
                                  opt!(tuple!(preceded!(comma_wsp,
                                                        number),
                                              preceded!(comma_wsp,
                                                        number))
                                  )),
                           |(a,c)| {
                               let a = a * PI / 180.0;
                               match c {
                                   Some((cx,cy)) =>
                                       Transform::rotate_around(a, 
                                                                &Point {x: cx,
                                                                        y: cy}),
                                   None => Transform::rotate(a)
                               } }),
                      tuple!(wsp_opt,char!(')'))));
    
    named!(skewx<Input,Transform>,  
           delimited!(tuple!(tag!("skewX"),wsp_opt,char!('('),wsp_opt),
                      map!(number,
                           |a| {
                               let a = a * PI / 180.0;
                              Transform::skew_x(a)
                           }),
                      tuple!(wsp_opt,char!(')'))));
    
    named!(skewy<Input,Transform>,  
           delimited!(tuple!(tag!("skewY"),wsp_opt,char!('('),wsp_opt),
                      map!(number,
                           |a| {
                               let a = a * PI / 180.0;
                               Transform::skew_y(a)
                           }),
                      tuple!(wsp_opt,char!(')'))));
    named!(transform<Input, Transform>, 
           alt!(matrix 
                | translate 
                | scale 
                | rotate 
                | skewx 
                | skewy));
    
    named!(transforms<Input, Transform>, 
           alt!(map!(complete!(tuple!(transform, comma_wsp, transforms)), 
                     |(a,_,b)| a*b)
                | transform
           ));
    
    named!(pub transform_list<Input, Transform>, preceded!(multispace0,transforms));

    named!(pub path_command<Input, (char, Vec<f64>)>,
               tuple!(preceded!(wsp_opt, alt!(char!('m') | char!('M')
                                              | char!('z') | char!('Z')
                                              | char!('l') | char!('L')
                                              | char!('h') | char!('H')
                                              | char!('v') | char!('V')
                                              | char!('c') | char!('C')
                                              | char!('s') | char!('S')
                                              | char!('q') | char!('Q')
                                              | char!('t') | char!('T')
                                              | char!('a') | char!('A'))),
                      preceded!(wsp_opt,
                                separated_list0!(comma_wsp, number))));


    named!(pub view_box_args<Input, [f64;4]>,
           map!(tuple!(wsp_opt, number, comma_wsp, number, comma_wsp, 
                    number, comma_wsp, number, wsp_opt),
                |(_,a,_,b,_,c,_,d,_)| [a,b,c,d] 
           ));

    // Returns unit length as milimeters
    const INCH:f64 = 25.4;
    const PX:f64 = 1.0;
    const PT:f64 = INCH/72.0;
    
    named!(pub length_unit<Input, f64>,
           alt!(value!(1.0,tag!("mm"))
                | value!(0.01,tag!("cm"))
                | value!(INCH,tag!("in"))
                | value!(PX,tag!("px"))
                | value!(PT,tag!("pt"))
                | value!(PT*12.0,tag!("pc"))
           ));
    named!(pub physical_length<Input, f64>,
        map!(tuple!(wsp_opt, number, wsp_opt, 
                    opt!(terminated!(length_unit,wsp_opt))),
             |(_,value,_,opt_unit)| {
                 match opt_unit {
                     Some(unit) => value * unit,
                     None => value
                 }
             }));
                                
}

pub fn parse_transform(s: &str) -> Result<Transform,nom::error::Error<&str>> {
    
    match parser::transform_list(s).finish() {
        Ok((_,o)) => Ok(o),
        Err(e) => Err(e)
    }
        
}

pub fn parse_view_box(s: &str) -> Result<[f64;4], nom::error::Error<&str> > {
    
    match parser::view_box_args(s).finish() {
        Ok((_,o)) => Ok(o),
        Err(e) => Err(e)
    }
        
}

pub fn parse_length(s: &str) -> Result<f64, nom::error::Error<&str>> {
    
    match parser::physical_length(s).finish() {
        Ok((_,o)) => Ok(o),
        Err(e) => Err(e)
    }
        
}

pub struct PathContext
{
    pos: Point,
    start: Point,
    last_control: Option<Point>
}

impl PathContext
{
    pub fn new() -> PathContext
    {
        PathContext {
            pos: Point {x: 0.0, y: 0.0},
            start: Point{x: 0.0, y: 0.0},
            last_control: None
        }
    }
}

struct PointIterator<I>
where I: Iterator
{
    iter : I
}

impl<'a, I> PointIterator<I> 
    where I: Iterator<Item = &'a f64>
{
    pub fn new(iter: I) -> PointIterator<I> {
        PointIterator{iter}
    }
}

impl<'a, I> Iterator for PointIterator<I>
    where I: Iterator<Item = &'a f64>
{
    type Item = Point;
    fn next(&mut self) -> Option<Self::Item> {
        let a = match self.iter.next() {
            Some(i) => i,
            None => return None
        };
        let b = match self.iter.next() {
            Some(i) => i,
            None => return None
        };
        Some(Point {x:*a,y:*b})
    }
}

fn mul2x2(a: [f64;4], b: [f64;4]) -> [f64;4]
{
    [a[0] * b[0] + a[1] * b[2],
     a[0] * b[1] + a[1] * b[3],
     a[2] * b[0] + a[3] * b[2],
     a[2] * b[1] + a[3] * b[3]]
}

fn invert2x2(m: [f64;4]) -> Option<[f64;4]>
{
    let d = m[0] * m[3] - m[1] * m[2];
    if d == 0.0 {
        None
    } else {
        Some([m[3] / d, -m[1] / d,
              -m[2] / d, m[0] / d])
    }
}

fn transpose2x2(m: [f64;4]) -> [f64;4]
{
    [m[0], m[2],
     m[1], m[3]]
}

#[allow(dead_code)] 
fn ellipse_canonical(mut rx : f64, mut ry: f64, mut rot: f64) -> (f64,f64, f64)
{
    if (rx-ry).abs() < f64::EPSILON {
        (rx,ry, 0.0)
    } else {

        if rot < 0.0 {
            rot += PI;
        }
        if rot >= PI/2.0 {
            std::mem::swap(&mut rx,&mut ry);
            rot -= PI/2.0;
        }
        (rx,ry, rot)
    }
}


// Returns (rx, ry, rot)
#[allow(non_snake_case)] 
fn transform_ellipse(rx: f64, ry: f64, rot: f64, m: [f64;4]) -> (f64,f64, f64)
{
    let (rs, rc) = rot.sin_cos();
    let rs2 = rs *rs;
    let rc2 = rc *rc;
    let rx2 = rx * rx;
    let ry2 = ry * ry;

    // A*x^2 + B*y^2 + 2*C*x*y = 1 */
    let A = rc2/rx2 + rs2/ry2;
    let B = rs2/rx2 + rc2/ry2;
    let C = (1.0/rx2 - 1.0/ry2)*rs*rc;

    let mut q = [A, C,
                 C, B];
    //println!("Before: {}*x^2 + {}*y^2 + {}*x*y = 1", q[0],q[3],2.0*q[2]);
    let invm = match invert2x2(m) {
        Some(m) => m,
        None => panic!("Can't invert transformation matrix for ellipse")
    };

    q = mul2x2(transpose2x2(invm), q);
    q = mul2x2(q, invm);
    //println!("After: {}*x^2 + {}*y^2 + {}*x*y = 1", q[0],q[3],2.0*q[2]);
    let trot = (-q[2]* 2.0).atan2(q[3] - q[0])/2.0;

    let (ts, tc) = trot.sin_cos();
    let ts2 = ts *ts;
    let tc2 = tc *tc;
    let tcs = tc *ts;
    let den = q[0] * tc2 + 2.0 * q[2] * tcs + q[3] * ts2;
    if den <= 0.0 {
        panic!("Failed to compute rx");
    }
    let trx = (1.0 / den).sqrt();

    let den = q[0] * ts2 - 2.0 * q[2] * tcs + q[3] * tc2;
    if den <= 0.0 {
        panic!("Failed to compute ry");
    }
    let try = (1.0 / den).sqrt();
    (trx, try, trot)
    //ellipse_canonical(trx, try, trot)
}    

struct ArcArgs {
    rx: f64,
    ry: f64,
    rot: f64,
    large_arc: bool,
    sweep: bool,
    x: f64,
    y: f64,
}
    
fn get_arc_args<'a, I:Iterator<Item=&'a f64>>(arg_iter: &mut I) -> Option<ArcArgs>
{
    Some(ArcArgs {
        rx: *arg_iter.next()?,
        ry: *arg_iter.next()?,
        rot: arg_iter.next().map(|v| {v * PI / 180.0})?,
        large_arc: arg_iter.next().map(|v| {*v != 0.0})?,
        sweep: arg_iter.next().map(|v| {*v != 0.0})?,
        x: *arg_iter.next()?,
        y: *arg_iter.next()?
    })
}

fn build_path(ctxt: &mut PathContext, cmd: char, args: &[f64],
              segs: &mut Vec<CurveSegment>, tr: &Transform, first_elem: bool) 
              -> Result<(), String>
{
    //println!("{}: {:?}", cmd, args);
    
    let rel = ('a'..= 'z').contains(&cmd);
    match cmd {
        'm' | 'M' => {
            if args.len() < 2 {
                return Err("Move-to command must have at least two arguments".to_string());
            }
            if args.len() % 2 != 0 {
                return Err("Move-to command must have an even number of arguments".to_string());
            }

            let mut coords = PointIterator::new(args.iter());
            if let Some(p) = coords.next() {
                 if rel && !first_elem {
                    ctxt.pos += p;
                } else {
                    ctxt.pos = p;
                 }
                let p1 = *tr*ctxt.pos;
                segs.push(CurveSegment::GoTo(p1));
                ctxt.start = ctxt.pos;
            }
            for p in coords {
                if rel {
                    ctxt.pos += p;
                } else {
                    ctxt.pos = p;
                }
                let p2 = *tr*ctxt.pos;
                segs.push(CurveSegment::LineTo(p2));
            }
            ctxt.last_control = None;
        },

        'l' | 'L' => {
            if args.len() < 2 {
                return Err("Line-to command must have at least two arguments".to_string());
            }
            if args.len() % 2 != 0 {
                return Err("Line-to command must have an even number of arguments".to_string());
            }
            
            for p in PointIterator::new(args.iter()) {
                if rel {
                       ctxt.pos += p;
                } else {
                    ctxt.pos = p;
                }
                let p2 = *tr*ctxt.pos;
                segs.push(CurveSegment::LineTo(p2));
            }
            ctxt.last_control = None;
        },
        'v' | 'V' => {
            if args.is_empty() {
                return Err("Vertical-line-to command must have at least one argument".to_string());
            }
            for &y in args {
                if rel {
                    ctxt.pos.y += y;
                } else {
                    ctxt.pos.y = y;
                }
                let p2 = *tr*ctxt.pos;
                segs.push(CurveSegment::LineTo(p2));
            }
            ctxt.last_control = None;
        },
        'h' | 'H' => {
            if args.is_empty() {
                return Err("Horizontal-line-to command must have at least one argument".to_string());
            }
            for &x in args {
                if rel {
                    ctxt.pos.x += x;
                } else {
                    ctxt.pos.x = x;
                }
                let p2 = *tr*ctxt.pos;
                segs.push(CurveSegment::LineTo(p2));
            }
            ctxt.last_control = None;
        },
        'c' | 'C' => {
            //println!("curve-to");
            if args.len() < 6 {
                return Err("Cubic curve-to command must have at least six arguments".to_string());
            }
            if args.len() % 6 != 0 {
                return Err("The number of arguments to cubic curve-to must be a multiple of six".to_string());
            }
            let mut arg_iter = args.iter();
            
            loop {
                let offset = if rel {ctxt.pos} else {Point{x: 0.0, y:0.0}};
                
                let &x1 = match arg_iter.next() {
                    Some(v) => v,
                    None => break
                };
                let &y1 = match arg_iter.next() {
                    Some(v) => v,
                    None => break
                };
                let &x2 = match arg_iter.next() {
                    Some(v) => v,
                    None => break
                };
                let &y2 = match arg_iter.next() {
                    Some(v) => v,
                    None => break
                };
                let &x = match arg_iter.next() {
                    Some(v) => v,
                    None => break
                };
                let &y = match arg_iter.next() {
                    Some(v) => v,
                    None => break
                };
                let p2 = Point{x, y} + offset;
                let c1 = Point{x: x1, y: y1} + offset;
                let c2 = Point{x: x2, y: y2} + offset;
                
                let c1 = *tr * c1 - *tr * ctxt.pos;
                ctxt.pos = p2;
                ctxt.last_control = Some(c2);
                let p2 = *tr * p2;
                let c2 = *tr * c2 - p2;
                segs.push(CurveSegment::CurveTo(p2, c1, c2));
            }
        },
         'q' | 'Q' => {
            //println!("curve-to");
            if args.len() < 4 {
                return Err("Quadratic curve-to command must have at least four arguments".to_string());
            }
            if args.len() % 4 != 0 {
                return Err("The number of arguments to quadraticcurve-to must be a multiple of four".to_string());
            }
            let mut arg_iter = args.iter();
            
            loop {
                let offset = if rel {ctxt.pos} else {Point{x: 0.0, y:0.0}};
                
                let &x1 = match arg_iter.next() {
                    Some(v) => v,
                    None => break
                };
                let &y1 = match arg_iter.next() {
                    Some(v) => v,
                    None => break
                };
              
                let &x = match arg_iter.next() {
                    Some(v) => v,
                    None => break
                };
                let &y = match arg_iter.next() {
                    Some(v) => v,
                    None => break
                };
                let p2 = Point{x, y} + offset;
                let m1 = Point{x: x1, y: y1} + offset;
                println!("p2: {}, m1: {}", p2, m1);

                ctxt.last_control = Some((m1-p2) * (2.0 / 3.0) + p2);
                println!("last_control: {:?}", ctxt.last_control);
                
                let c1 = (*tr * m1 - *tr * ctxt.pos) * (2.0 / 3.0);
                ctxt.pos = p2;
                let p2 = *tr * p2;
                let c2 = (*tr * m1 - p2) * (2.0 / 3.0);
                

                segs.push(CurveSegment::CurveTo(p2, c1, c2));
            }
        },
        'a' | 'A' => {
            if args.len() < 7 {
                return Err("Arc command must have at least seven arguments".to_string());
            }
            if args.len() % 7 != 0 {
                return Err("The number of arguments to curve-to must be a multiple of seven".to_string());
            }
             let mut arg_iter = args.iter();
            
            while let Some(args) = get_arc_args(&mut arg_iter) {
                let Point{x:dx,y:dy} = if rel {
                    *tr*(Point{x:args.x,y:args.y} + ctxt.pos)
                } else {
                    *tr*Point{x:args.x,y:args.y}
                } - *tr * ctxt.pos;
                
                if dx == 0.0 && dy == 0.0 {
                    return Err("Zero length arc".to_string());
                }
                //println!("rx: {} ry: {}",rx,ry);
                //println!("x: {} y: {}",x,y);

                // Transform to rendering coords
                let (rx, ry, rot) = 
                    transform_ellipse(args.rx,args.ry, args.rot,
                                      [tr.matrix[0], tr.matrix[2],
                                       tr.matrix[1], tr.matrix[3]]);
                                               
                
                let (rs,rc) = args.rot.sin_cos();

                let x1;
                let y1;
                let x2;
                let y2;
                
                // Rotate the intersecting line to match the unrotated
                // ellipse.  This rotates the line in the opposite
                // direction of the ellipse,
                let (dx,dy) = (dx * rc - dy * rs, dx * rs + dy * rc);
                //println!("dx: {} dy: {}",dx,dy);

                if dy.abs() < dx.abs() {
                    let k = dy / dx;
                    let ry2 = ry*ry;
                    let c = ry2 / (rx * rx) + k * k;
                    let m2 = (c * dx *dx / 4.0 - ry2) / (k*k/c - 1.0);
                    //println!("(dy/dx) m2: {}",m2);
                    if m2 < 0.0 {
                        return Err("Ellipse too small for given endpoints of arc".to_string());
                    }
                    let m = m2.sqrt();
                    // Mid-point
                    let mp = -k*m/c;
                    x1 = mp-dx/2.0;
                    x2 = mp+dx/2.0;
                    y1 = k*x1+m;
                    y2 =  k*x2+m;
                } else {
                    let k = dx / dy;
                    let rx2 = rx*rx;
                    let c = rx2 / (ry * ry) + k * k;
                    let m2 = (c * dy *dy / 4.0 - rx2) / (k*k/c - 1.0);
                    //println!("(dx/dy) m2: {}",m2);
                    if m2 < 0.0 {
                        return Err("Ellipse too small for given endpoints of arc".to_string());
                    }
                    let m = m2.sqrt();
                    // Mid-point
                    let mp = -k*m/c;
                    y1 = mp-dy/2.0;
                    y2 = mp+dy/2.0;
                    x1 = k*y1+m;
                    x2 = k*y2+m;
                }
                //println!("({}, {}) ({}, {})",x1,y1, x2,y2);
                let a1 = (y1/ry).atan2(x1/rx);
                let a2 = (y2/ry).atan2(x2/rx);
                
                let a2 = if args.large_arc ^ ((a2 - a1).abs() > PI) {
                    if a1 < a2 {
                        a2 - 2.0 * PI
                    } else {
                        a2 + 2.0 * PI
                    }
                } else {
                    a2
                };
                
                let (a1,a2) = if args.sweep ^ (a2 > a1) {
                    (a2 + PI, a1 + PI) 
                } else {
                    (a1,a2)
                };

                
                /* println!("a1: {} a2: {}",a1,a2);

               
                let xa1 = rx * a1.cos();
                let ya1 = ry * a1.sin();
                let (xa1,ya1) = (xa1*rc - ya1*rs, xa1*rs + ya1*rc);
                
                let xa2 = rx * a2.cos();
                let ya2 = ry * a2.sin();
                let (xa2,ya2) = (xa2*rc - ya2*rs, xa2*rs + ya2*rc);

                println!("Calculated: ({}, {})", xa2 -xa1, ya2-ya1);
                */
                segs.push(CurveSegment::Arc(rx,ry, a1,a2, rot));
                
                if rel {
                    ctxt.pos.x += args.x;
                    ctxt.pos.y += args.y;
                } else {
                    ctxt.pos.x = args.x;
                    ctxt.pos.y = args.y;
                }
            }
            
        },
        'z' | 'Z' => {
            ctxt.pos = ctxt.start;
            let p2 = *tr*ctxt.pos;
            segs.push(CurveSegment::LineTo(p2));
        },
        c => return Err(format!("Unknown command '{}'", c))
    };
    Ok(())
}

pub fn parse_path(path_str: &str, path_ctxt: &mut PathContext,
                  transform: &Transform, segs: &mut Vec<CurveSegment>)
                  -> Result<(), Box<dyn Error + 'static>>
{
    let mut pos = path_str;
    let mut first_elem = true;
    while let Ok((rest, (cmd, args))) = 
        parser::path_command(pos) 
    {
        build_path(path_ctxt,
                   cmd, &args,
                   segs, &transform,
                   first_elem)?;
        //println!("rest: {:?}", rest);
        pos = rest;
        first_elem = false;
    }
    Ok(())
}


#[cfg(test)]
fn matrix_approx_eq(a:&[f64;6], b:&[f64;6]) -> bool
{
    for (a,b) in a.iter().zip(b) {
        if (a-b).abs() > 1e-5 {
            return false;
        }
    }
    return true;
}

#[cfg(test)]
fn assert_matrix_eq(a:&[f64;6], b:&[f64;6])
{
    for (a,b) in a.iter().zip(b) {
        if (a-b).abs() > 1e-5 {
            panic!("{} != {}", a,b);
        }
    }
}

#[cfg(test)]
fn assert_transform_eq(a:&Transform, b:&Transform)
{
    assert_matrix_eq(&a.matrix, &b.matrix);
}

#[cfg(test)]
fn assert_transform_result_eq(a:&nom::IResult<&str,Transform>,
    b:&nom::IResult<&str,Transform>)
{
    if let (&Ok((_,left)), &Ok((_,right))) = (a,b) {
        assert_transform_eq(&left, &right);
    } else {
        panic!("{:?} or {:?} has no value", a, b);
    }
}

#[cfg(test)]
fn assert_parsed_transform_eq(a: &str, b: &str)
{
    let am = match parser::transform_list(a) {
        Ok((_,o)) => o,
        err => panic!("String: '{}' did not parse: {:?}",a, err)
    };
    let bm = match parser::transform_list(b) {
        Ok((_,o)) => o,
        err => panic!("String: '{}' did not parse: {:?}",b, err)
    };
    if !matrix_approx_eq(&am.matrix,&bm.matrix) {
        panic!("{} ({:?}) != {} ({:?})", a, am, b,bm);
    }
}

#[cfg(test)]
fn transform_result(m : &[f64;6]) ->  nom::IResult<&str,Transform>
{
    Ok::<(&str,Transform), _>(("", Transform::new(m)))
}

#[test]
fn test_parser_matrix()
{
    let res = Ok::<(&str,Transform), _>(("", Transform::new(&[-8.0, 0.0, 12.0, 9.3, 1.2e3, 1.2e-2])));
    assert_eq!(parser::transform_list("matrix(-8 0 12 9.3 1.2e3 1.2e-2)"), res);
    assert_eq!(parser::transform_list("translate(-8 0.5)"),
               transform_result(&[1.0, 0.0, 0.0, 1.0, -8.0, 0.5]));
    assert_eq!(parser::transform_list("translate(1e-2)"), 
               transform_result(&[1.0, 0.0, 0.0, 1.0, 0.01, 0.0]));

     assert_transform_result_eq(&parser::transform_list("rotate(30 8.3 1.2)"),
                                &transform_result(&[0.86603, 0.50000,
                                                    -0.50000, 0.86603,
                                                    1.71199,-3.98923]));

    assert_transform_result_eq(&parser::transform_list("rotate(30)"), 
                               &transform_result(&[ 0.86603, 0.50000, 
                                                    -0.50000, 0.86603,
                                                    0.0,0.0]));

    assert_transform_result_eq(&parser::transform_list("skewX(30)"),
                               &transform_result(&[1.0, 0.0, 
                                                   0.57735, 1.0, 
                                                   0.0, 0.0]));
    assert_transform_result_eq(&parser::transform_list("skewY(70.4)"),
                               &transform_result(&[1.0, 2.8083262,
                                                   0.0, 1.0, 
                                                   0.0, 0.0]));
    assert_parsed_transform_eq(
        &"matrix(0.70711 0.70711 -0.70711 0.70711 8.3 1.2)",
        &"translate( 8.3 1.2), rotate(45)");
    assert_parsed_transform_eq(
        &"matrix(0.70711 0.70711 -0.70711 0.70711   7.071067811  7.071067811)",
        &"rotate(45) translate( 10 0)");
    assert_parsed_transform_eq(
        &"matrix(0.70711 0.70711 -0.70711 0.70711  5.0204581464 6.7175144212722)",
        &" rotate(45) translate( 8.3 1.2)");
    
    
    assert_parsed_transform_eq(
        &"rotate(30 8.3 1.2)",
        &"translate(8.3 1.2),rotate(30), translate( -8.3 -1.2)");
    assert_parsed_transform_eq(&"translate(6,5) rotate(23)",
                               "rotate(23,-9.2878926,17.245471)");
    assert_parsed_transform_eq(&"translate(6,5) rotate(23)",
                               &"rotate(23,-9.2878926,17.245471)");
    assert_parsed_transform_eq(&"translate(23,12) skewX(32)",
                               &"matrix(1,0,0.62486935,1,23,12)");
    assert_parsed_transform_eq(&"translate(-23.5,12) skewY(-56)",
                               &"matrix(1,-1.482561,0,1,-23.5,12)");
    
}

#[cfg(test)]
fn assert_ellipse_eq(rx1: f64, ry1: f64, rot1: f64,
                     rx2: f64, ry2: f64, rot2: f64)
{
    const EPSILON: f64 = 1e-6;
    let (rx1,ry1,rot1) = ellipse_canonical(rx1,ry1,rot1);
    let (rx2,ry2,rot2) = ellipse_canonical(rx2,ry2,rot2);
    if (rx1 - rx2).abs() > EPSILON || (ry1 - ry2).abs() > EPSILON 
        || ((rot1 - rot2).abs() > EPSILON && (rx1 - rx2).abs() > EPSILON)
    {
        panic!("Ellipses are not equal: ({}, {}, {}) != ({}, {}, {})",
               rx1,ry1,rot1, rx2,ry2,rot2);
    }
}
#[test]
fn test_ellipse_transform()
{
    let (rx,ry,rot) = transform_ellipse(2.0, 3.0, 0.0, [1.0,0.0,0.0,1.0]);
    assert_ellipse_eq(rx,ry,rot, 2.0, 3.0, 0.0);

    let (rx,ry,rot) = transform_ellipse(2.0, 3.0, 0.0, [3.0,0.0,0.0,3.0]);
    assert_ellipse_eq(rx,ry,rot, 6.0, 9.0, 0.0);

    let a = 0.2f64;
    let (rs, rc) = a.sin_cos();
    let (rx,ry,rot) = transform_ellipse(2.0, 3.0, -0.2, [rc,-rs,rs,rc]);
    assert_ellipse_eq(rx,ry,rot, 2.0, 3.0, 0.0);
    
    let a = 0.2f64;
    let (rs, rc) = a.sin_cos();
    let (rx,ry,rot) = transform_ellipse(2.0, 3.0, 0.0, [rc,-rs,rs,rc]);
    assert_ellipse_eq(rx,ry,rot, 2.0, 3.0, 0.2);

    let a = 0.3f64;
    let (rs, rc) = a.sin_cos();
    let (rx,ry,rot) = transform_ellipse(2.0, 3.0, 0.2, [rc,-rs,rs,rc]);
    assert_ellipse_eq(rx,ry,rot, 2.0, 3.0, 0.5);
    
   
    
    let (rx,ry,rot) = transform_ellipse(2.0, 3.0, 0.0, [1.0, 1.0, 0.0, 1.0]);
    println!("Skewed:  ({}, {}, {})", rx,ry,rot);
    let (rx,ry,rot) = transform_ellipse(rx,ry,rot, [1.0, -1.0, 0.0, 1.0]);
    println!("Unskewed:  ({}, {}, {})", rx,ry,rot);
    assert_ellipse_eq(rx,ry,rot, 2.0, 3.0, 0.0);

    let m = [-1.0, 0.3,  -0.2, 2.0];
    let (rx,ry,rot) = transform_ellipse(2.0, 3.0, 0.0, m);
    println!("Skewed:  ({}, {}, {})", rx,ry,rot);
    let (rx,ry,rot) = transform_ellipse(rx,ry,rot, invert2x2(m).unwrap());
    println!("Unskewed:  ({}, {}, {})", rx,ry,rot);
    assert_ellipse_eq(rx,ry,rot, 2.0, 3.0, 0.0);

    // Circle
    let (rx,ry,rot) = transform_ellipse(3.0, 3.0, 0.0, [3.0,0.0, 0.0,3.0]);
    assert_ellipse_eq(rx,ry,rot, 9.0, 9.0, 0.0);
    
    let m = [-1.0, 0.3,  -0.2, 2.0];
    let (rx,ry,rot) = transform_ellipse(3.0, 3.0, 0.0, m);
    println!("Skewed:  ({}, {}, {})", rx,ry,rot);
    let (rx,ry,rot) = transform_ellipse(rx,ry,rot, invert2x2(m).unwrap());
    println!("Unskewed:  ({}, {}, {})", rx,ry,rot);
    assert_ellipse_eq(rx,ry,rot, 3.0, 3.0, 0.0);

    
}

#[test]
fn test_viewbox_parser()
{
    let s = "0 0 210 297";
    assert_eq!(parser::view_box_args(s), Ok(("", [0.0,0.0,210.0, 297.0])));
}

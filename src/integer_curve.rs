use solve;
use num::Integer;
use std::fmt;
use multi_range::MultiRange;
use std::i64;
use std::convert::TryFrom;
use acc_vector::AccSegment;
use acc_vector::AccVector;



fn min_i64(a:i64, b:i64) ->i64
{
    if a <= b {a} else {b}
}

fn max_i64(a:i64, b:i64) ->i64
{
    if a >= b {a} else {b}
}


//    ___vflat
//   /   \vn
//v0/     

// See DistancePath.asy for a graphical discription
#[derive(Clone)]
pub struct DistancePath {
    pub v0: i64, // Start speed
    pub vn: i64, // End speed
    pub amax:i64, // Maximum acceleration
    pub vflat:i64, // Flat part of speed curve, may be between v0 and vn
    pub t_total:i64, // Total time for path
    pub t_adjust: i64, // Adjust the flat part of the curve so that the speed is increased or decreased by one. The sign decides the direction and the absolute value is the length
}

/* Split acceleration into max acceleration time and adjustment step.
There is always an adjustment step in the range from -max_acc to +max_acc .
 */
pub fn split_acc(v_diff: i64, max_acc: i64) -> (i64, i64) {
    if v_diff > 0 {
            let dr = (v_diff-1).div_rem(&max_acc);
        (dr.0.abs(), dr.1 + 1) 
    } else if v_diff < 0 {
        let dr = (v_diff+1).div_rem(&max_acc);
        (dr.0.abs(), dr.1 - 1)
        } else {
        (0,0)
    }
}

impl DistancePath {
    /*
    pub fn distance(&self) -> i64 {
         let a = self.amax;
        let v0 = self.v0;
        let vn = self.vn;
        let vmax = self.vflat;
        let t1 = ((vmax - v0).abs() + a - 1) / a;
        let t4 = ((vmax - vn).abs() + a - 1) / a;
        let mut s = 2*vmax  * (self.t_total - t1 - t4);
        let v1 = ((vmax - v0) / a) * a + v0;
        let v4 = ((vmax - vn) / a) * a + vn;
        //println!("v1 = {}, v4 = {}", v1,v4);
        if v1 == vmax {
            s += t1*(v0 + v1);
        } else {
            s += (t1 - 1) * (v0 + v1) + (v1 + vmax);;
        }
        if v4 == vmax {
            s += t4*(v4 + vn);
        } else {
            s += (t4 - 1) * (v4 + vn) + (v4 + vmax);;
        }
        s += 2*self.t_adjust;
        s
    }*/

   /*
    pub fn distance(&self) -> i64 {
         let a = self.amax;
        let v0 = self.v0;
        let vn = self.vn;
        let vflat = self.vflat;
        if self.t_total <= 1 {
            if self.t_total == 0 {
                return 0;
            } else { // t_total == 1
                return v0 + vn;
            }
        }
        // t_total >= 2
        let (t1, adj1) = split_acc(vflat - v0,a);
        let (t4,adj4) = split_acc(vflat - vn,a);
        //println!("t1 = {}, adj1 = {}", t1,adj1);
        //println!("t4 = {}, adj4 = {}", t4,adj4);
        assert!(self.t_total >= t1 + t4 + 2);
        
        let v1 = vflat - adj1;
        let v4 = vflat - adj4;
        
        //println!("v1 = {}, v4 = {}", v1,v4);
        let mut s =
        // Flat part of velocity curve
            2*vflat  * (self.t_total - t1 - t4 - 2) 
        // Acceleration before flat part
            + (v0 + v1) * t1 + v1 + vflat
        // Acceleration after flat part
            + (v4 + vn) * t4 + vflat + v4;
        s += 2*self.t_adjust;
        s
    }
     */
    pub fn distance(&self) -> i64 {
        return path_distance(self.v0, self.vn, self.vflat, 
                             self.amax, self.t_total)
            + 2*self.t_adjust;
    }
    pub fn acc_seq(&self) -> Vec<AccSegment>
    {
        let a = self.amax;
        let mut s = Vec::new();
        if self.t_total <= 1 {
            if self.t_total == 0 {
                assert!(self.v0 == self.vn);
            } else {
                assert!((self.v0 - self.vn) <= a);
                s.acc_push(1, i16::try_from(self.vn - self.v0).unwrap());
            }
            return s;
        }
        let mut vflat_start = self.vflat;
        let mut vflat_end = self.vflat;
        let mut t_step = 0;
        if self.t_adjust > 0 {
            if self.vn > self.vflat {
                vflat_end += 1;
                t_step = -self.t_adjust;
                if (vflat_end - self.vn) % a == 0 {
                    t_step -=1;
                }
            } else {
                assert!(self.v0 > self.vflat);
                vflat_start += 1;
                t_step = self.t_adjust;
                if (vflat_start - self.v0) % a == 0 {
                    t_step += 1;
                }
            }
        } else if self.t_adjust < 0 {
            if self.vn < self.vflat {
                vflat_end -= 1;
                t_step = self.t_adjust;
                if (vflat_end  - self.vn) % a == 0 {
                    t_step -= 1;
                }
            } else {
                assert!(self.v0 < self.vflat);
                vflat_start -= 1;
                t_step = -self.t_adjust;
                if (vflat_start - self.v0) % a == 0 {
                    t_step += 1;
                }
            }
        }
        
        
        let (t1, dv1) = div_rem_away(vflat_start - self.v0, a);
        let (t4, dv4) = div_rem_away(vflat_end - self.vn, a);
        let t1 = t1.abs();
        let t4 = t4.abs();
        
        let a1 = if vflat_start >= self.v0 {a} else {-a};
        let a4 = if vflat_end <= self.vn {a} else {-a};
        assert!(self.t_total >= t1 + t4);
        let t_flat = self.t_total - t1 - t4;
        // println!("a1: {}, a4: {}, t_flat: {}",a1,a4, t_flat);
        if dv1 == 0 {
            s.acc_push(u16::try_from(t1).unwrap() , i16::try_from(a1).unwrap());
        } else {
            s.acc_push(u16::try_from(t1 - 1).unwrap(),
                       i16::try_from(a1).unwrap());
            s.acc_push(1, i16::try_from(vflat_start - (a1*(t1-1) + self.v0)).unwrap());
        }

        
        assert!(self.t_adjust <= t_flat);
        if t_step > 0 {
            s.acc_push(u16::try_from(t_step-1).unwrap(), 0);
            s.acc_push(1, i16::try_from(vflat_end - vflat_start).unwrap());
            s.acc_push(u16::try_from(t_flat - t_step).unwrap(), 0);
        } else if t_step < 0 {
            s.acc_push(u16::try_from(t_flat + t_step).unwrap(), 0);
            s.acc_push(1, i16::try_from(vflat_end - vflat_start).unwrap());
            s.acc_push(u16::try_from(-t_step - 1).unwrap(), 0);
        } else {
            s.acc_push(u16::try_from(t_flat).unwrap(), 0);
        }

        if dv4 == 0 {
            s.acc_push(u16::try_from(t4).unwrap(), i16::try_from(a4).unwrap());
        } else {
            s.acc_push(1, i16::try_from((self.vn - a4*(t4-1)) -
                                        vflat_end).unwrap());
            s.acc_push(u16::try_from(t4 - 1).unwrap() ,
                       i16::try_from(a4).unwrap());
        }
        return s;
    }
}

impl fmt::Display for DistancePath {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
       
        write!(f, "v= {} -> {} amax = {}, vflat = {}, t_total = {}, t_adjust = {}, s = {}",
               self.v0, self.vn,
               self.amax, self.vflat, self.t_total, self.t_adjust,
               self.distance()
            )
    }
}

    
// ->(n,a)
pub fn max_acceleration_curve(ds:i64, v0: i64, vn: i64, a: i64) -> MultiRange<i64>
{
    let mut r = MultiRange::<i64>::new();
    let tmin = 
        if v0 != vn {
            ((v0-vn).abs() + a - 1) / a
        } else {
            1
        };
    
    // Constants for the limiting polynom
    // c = cx2 x^2 + cx x;
    let cx2 = a*a;
    let cx = 2*a*(v0+vn);
    let c = (v0-vn)*(v0-vn) + a*a + 2*a * ds;
    
    // Positive acceleration
    //println!("cx2: {}, cx: {}, c: {}", cx2,cx,c);
    let s = c + (v0 + vn)*(v0 + vn);
    if s < 0 { // No solution
        r.or(tmin, i64::MAX);
        return r;
    }
    // Find minima
    let xmin =
    // Round -(v0 + vn) / a to neareast integer
        if v0 + vn >= 0 {
            -((2*(v0 + vn) + a) / (2*a))
        } else {
            (-2*(v0 + vn) + a)  / (2*a)
        };
    
    //println!("xmin: {}", xmin);
    let ymin = solve::sqpoly(cx2,cx,c, xmin);
    if ymin >= 0 {
        // Either a double root or both roots are between adjacent integers
        return MultiRange::range(tmin, i64::MAX);
    }

    let mut limit = 16;
    while solve::sqpoly(cx2,cx,c, xmin + limit) <= 0 {
        limit *= 2;
    }
    limit += 1;
    match (solve::find_root(cx2,cx,c, xmin - limit, xmin), 
           solve::find_root(cx2,cx,c, xmin, xmin + limit)) {
        (Some(min), Some(max)) => {
            r.or(i64::MIN, min);
            r.or(max, i64::MAX);
            r = r.and(tmin, i64::MAX);
        },
        _ => panic!("No solution found, but one exists")
    }
    return r;
}




pub fn speed_limited_curve(ds:i64, v0: i64, vn: i64, a: i64, vmax:i64) -> MultiRange<i64>
{
    let mut r = MultiRange::<i64>::new();
    // Check positive limit
    let (t1, dv1)  = (vmax - v0 + a - 1).div_mod_floor(&a);
    let (t4,dv4) = (vmax - vn + a - 1).div_mod_floor(&a);
    let s0 = (2*v0 + t1*a) * t1 + (dv1 - (a - 1))
        + (2*vn + t4*a) * t4 + (dv4 - (a - 1));
    if ds >= s0 {
        r.or(t1 + t4 + (ds - s0 + 2*vmax  - 1) / (2*vmax), i64::MAX);
    }
    return r;
}

// Divide and round away from zero with remainder.
// If result is (q,r) then b * q + r = a
// b > 0
fn div_rem_away(a:i64,b:i64) -> (i64, i64)
{
    assert!(b > 0);
    let (q1, r1) = a.div_rem(&b);
    if r1 == 0 {(q1,r1)}
    else {
        if r1 > 0 {
            (q1 + 1, r1 - b)
        } else {
            (q1 - 1, r1 + b)
        }
    }
}

fn path_distance(v0: i64, vn:i64, vflat:i64, a:i64, t: i64) -> i64
{
    /*
    println!("path_distance: v0: {}, vn: {}, vflat: {}, a: {}, t: {}",
             v0,vn,vflat,a,t);*/
    let (t1, dv1) = div_rem_away(vflat - v0, a);
    let (t4, dv4) = div_rem_away(vflat - vn, a);
    let t1 = t1.abs();
    let t4 = t4.abs();
    return (vflat-dv1 +v0) * t1  + dv1 
        + 2*vflat* (t - t1 - t4) 
        + dv4 + (vflat-dv4 + vn) * t4;
}

fn adjust_curve_mid(t_total: i64, ds: i64, v0: i64, vn: i64, a: i64) 
                    -> Option<DistancePath>
{
     let (v_low, v_high) = if v0 <= vn {(v0,vn)} else {(vn,v0)};
    // Always keep an adjustment step
    let t_acc = (v_high - v_low - 1) / a;

    let k = 2*(v_low - v_high + a * t_total);
    let m = v_low - v_high + a*t_acc + 2*(v_high-a*t_acc)*t_total 
        + a * t_acc *t_acc; 
    // println!("k= {}, m= {}",k,m);
    if k == 0 {
        /* There is only one possible path, no adjustment can be done. */
        return Some(DistancePath {v0: v0, vn: vn, amax:a, 
                                  vflat:v_low+a, 
                                  t_total:t_total, t_adjust: 0});
    }
    let smax = k * t_acc + m; // Acceleration
    let smin = (2*v_high - t_acc * a) * t_acc // Acceleration
        + v_high - t_acc * a + v_low // Adjustment step
        + 2 * v_low * (t_total - t_acc - 1); // Constant velocity
    // Return unless v_low <= v_flat <= v_high
    if smax < ds || smin > ds {return None;}
    assert!(k > 0);
    assert!(ds- m + k >= 0);
    let t = (ds - m + k) / k;
    let s_high = k * t + m;
    let s_split = (2 * v_low + a * t) * t // Acceleration
        + 2 * (v_low + a* t) * (t_total - t_acc - 1) // Constant velocity
        + v_high - (t_acc - t) * a + v_low + t*a // Adjustment step
        + (2*v_high - (t_acc - t) * a) * (t_acc - t); // Acceleration
    let v_flat;
    let s;
    if ds >= s_split {
        let v_max = v_high - a * (t_acc - t);
        v_flat = v_max - (s_high - ds) / (2*(t_total - t_acc - 1));
        s = s_high - 2 * (v_max - v_flat) * (t_total - t_acc - 1);
    } else {
        let v_min = v_low + a * t;
        v_flat = v_min - (s_split - ds) / (2*(t_total - t_acc));
        s = s_split - 2 * (v_min - v_flat) * (t_total - t_acc);
    }
    // println!("v_flat= {}, s_flat= {}", v_flat, s);

    //println!("t= {}, s_high= {}, s_plit= {}",t, s_high, s_split);
    return Some(DistancePath {v0:v0, vn: vn, amax:a, 
                              vflat:v_flat, 
                              t_total:t_total, t_adjust: (ds - s) / 2});
}

fn adjust_curve(t_total: i64, ds: i64,v0: i64, vn: i64, a: i64) -> DistancePath
{

    /* Path of length 1 is a special case. No adjustment is
     * possible. It either fits or it doesn't.
     */
    if t_total ==  1 {
        let s = v0 + vn;
        assert!((s-ds).abs() <= 1);
        return DistancePath {v0:v0, vn: vn, amax:a, vflat:v0, 
                             t_total:1, t_adjust:0};
    }

    // Solve for when v0 <= vflat <= vn or v0 >= vflat >= vn
    if let Some(path) = adjust_curve_mid(t_total, ds, v0, vn, a) {
        return path;
    }
    

    // Find highest and lowest possible vflat
    // Time from v0 to highest peak
    let tmid_low = (t_total * a - (v0 - vn)) / (2*a);
    // Time from highest peak to vn
    let tmid_high = (t_total * a - (vn - v0)) / (2*a);
    //println!("tmid_low: {}, tmid_high: {}",tmid_low, tmid_high);
    assert!(tmid_low + tmid_high <= t_total);
    assert!(tmid_low + tmid_high >= t_total - 1);
    let mut max = 
        if tmid_high == 0 {
            v0 + tmid_low * a
        } else if tmid_low == 0 {
            vn + tmid_high * a
        } else {
            max_i64(v0 + tmid_low * a, vn + tmid_high * a)
        };
    let mut min = 
        if tmid_high == 0 {
            vn + tmid_low * -a
        } else if tmid_low == 0 {
            v0 + tmid_high * -a
        } else {
            min_i64(v0 + tmid_high * -a, vn + tmid_low * -a)
        };
    assert!(max >= min);
    let mut smax = path_distance(v0,vn,max, a, t_total);
    let mut smin = path_distance(v0,vn,min, a, t_total);
    //println!("({} -> {}) - ({} -> {})", min, smin, max, smax);
    
    assert!(smax >= ds);
    assert!(smin <= ds);
    // Find closest min and max
    loop {
        // println!("({} -> {}) - ({} -> {})", min, smin, max, smax);
        if (max - min) <= 1 {break};
        let mid = (max + min) / 2;
        let s = path_distance(v0,vn,mid, a, t_total);
        if s >= ds {
            max = mid;
            smax = s;
        } else {
            min = mid;
            smin = s;
        }
    }
    let vflat;
    let mut t_adjust = 0;
    if smax == ds {
        vflat = max;
    } else if smin == ds {
        vflat = min;
    } else {
        if v0 > min || vn > min {
            vflat = min;
            t_adjust = (ds - smin) / 2;
        } else  {
            vflat = max;
            t_adjust = (ds - smax) / 2;
        }
    }
    return DistancePath {vflat: vflat, v0: v0, vn:vn, amax: a, 
                         t_adjust: t_adjust, t_total: t_total};
}

pub fn shortest_curve_length(ds:i64, v0: i64, vn: i64, a: i64, vmax:i64) 
                             -> MultiRange<i64>
{
    let mut pos = speed_limited_curve(ds, v0, vn, a, vmax);
    if pos.empty() {
        pos = max_acceleration_curve(ds, v0, vn, a);
    };
    let mut neg = speed_limited_curve(-ds, -v0, -vn, a, vmax);
    if neg.empty() {
        neg = max_acceleration_curve(-ds, -v0, -vn, a);
    };
    return pos.and_range(&neg);
}

pub fn shortest_curve(ds:i64, v0: i64, vn: i64, a: i64, vmax:i64) -> DistancePath
{
    let (t,_) = shortest_curve_length(ds,v0,vn,a,vmax).bounds().unwrap();
    return adjust_curve(t, ds, v0, vn, a);
}

#[derive(Debug)]
pub struct PathLimits
{
    pub ds:i64, 
    pub v0: i32,
    pub vn: i32, 
    pub a: i32, 
    pub vmax:i32
}
    
pub fn shortest_curve_sequences(limits: &[PathLimits]) -> Vec<Vec<AccSegment>>
{
    let mut r = MultiRange::new();
    r.or(i64::MIN, i64::MAX);
    for l in limits {
        r = r.and_range(&shortest_curve_length(l.ds,i64::from(l.v0),
                                               i64::from(l.vn),
                                               i64::from(l.a),
                                               i64::from(l.vmax)));
    }
    let (t,_) = r.bounds().unwrap();
    let mut a: Vec<Vec<AccSegment>> = Vec::new();
    for l in limits {
        let path = adjust_curve(t, l.ds,i64::from(l.v0),i64::from(l.vn),
                                i64::from(l.a));
        a.push(path.acc_seq());
    }
    a
}

#[cfg(test)]
fn check_max_acceleration_curve(ds:i64, v0: i64, vn: i64, a: i64)
{
    println!("\nds= {}, v0= {}, vn= {}, a= {}", ds,v0,vn,a);
    let range = max_acceleration_curve(ds,v0,vn,a);
    println!("=> {:?}", range);
    if let Some((n, _)) = range.bounds() {
        let a2 = a;
        println!("n= {}, a2= {}",n,a2);
        let dt1 = (a*n + vn - v0) / (2*a);
    let dt2 = (a*n + v0 - vn) / (2*a);
        println!("dt1= {}, dt2= {}",dt1, dt2);
        assert!(dt1 >= 0);
        assert!(dt2 >= 0);
        assert!(dt2+ dt1 >= n - 1);
        assert!(dt2+ dt1 <= n);
        let ds2 = (a2*dt1+ 2*v0)*dt1 + (a2*dt2+ 2*vn)*dt2 
            + (n - dt1 -dt2)*(v0+a2*dt1+vn + a2*dt2);
        println!("ds2= {}", ds2); 
        assert!(ds2 >= ds);
    }
}

#[test]
fn test_max_acceleration_curve()
{
    for ds in -20..20 {
        for v0 in -5..5 {
            for vn in -5..5 {
                for a in 1..6 {
                    check_max_acceleration_curve(ds, v0, vn,a);
                }
            }
        }
    }
}

#[test]
fn test_div_rem_away()
{
    for a in -20..20 {
        for b in 1..20 {
            let (q,r) = div_rem_away(a,b);
            assert!((b*q).abs() >= a.abs()); 
            assert_eq!(b*q+r, a);
        }
    }
}

#[test]
fn test_adjust_curve()
{
    let ds = 23;
    let adjusted = adjust_curve(6, ds, -5, 3, 4);
    println!("path = {}, ds = {}", adjusted,ds);
    assert!((adjusted.distance() - ds).abs() <= 1);
}

#[test]
fn test_adjust_curve_mid()
{
    let a = 3;
    for v0 in -5..5 {
        for vn in v0 .. v0 + 10 {
            let n_min = (vn - v0 + a - 1) / a;
            for t_total in n_min .. n_min + 3 {
                let s_min = path_distance(v0, vn, v0, a, t_total);
                let s_max = path_distance(v0, vn, vn, a, t_total);
                for ds in s_min .. s_max {
                    match adjust_curve_mid(t_total, ds, v0, vn, a) {
                        None => {panic!("Outside");}
                        Some(adjusted) => {
                            println!("path = {}, ds = {}", adjusted,ds);
                            assert!((adjusted.distance() - ds).abs() <= 1);
                        }
                    };
                }
            }
        }
    }
}

#[test]
fn test_path_distance()
{
    let v0 = -2;
    let a = 4;
    for dv1 in -7i64..7i64 {
        for dv4 in -7i64..7i64 {
            let vn = v0 + dv1 + dv4;
            let t1 = (dv1.abs() + a - 1) / a;
            let t4 = (dv4.abs() + a - 1) / a;
            let path = DistancePath {v0: v0, vn: vn, amax: a, vflat: v0 + dv1,
                                     t_total: t1 + t4 + 3, t_adjust: 0};
            println!("Path: {}", path);
            assert_eq!(path.distance(), path_distance(v0,vn, path.vflat, a, path.t_total));
        }
    }
}

#[test]
fn test_acc_sequence()
{
    let v0 = -2;
    let a = 4;
    for dv1 in -7i64..7i64 {
        for dv4 in -7i64..7i64 {
            let vn = v0 + dv1 + dv4;
            let t1 = (dv1.abs() + a - 1) / a;
            let t4 = (dv4.abs() + a - 1) / a;
            let path = DistancePath {v0: v0, vn: vn, amax: a, vflat: v0 + dv1,
                                     t_total: t1 + t4 + 3, t_adjust: 0};
            let seq = path.acc_seq();
            let mut s:i64 = 0;
            let mut v = v0;
            for &AccSegment{interval: l,acc: a} in &seq {
                assert!(l > 0);
                let vp = v;
                v += i64::from(l) * i64::from(a);
                s += (vp + v) * i64::from(l);
            }
            println!("{} => {:?} => {}", path, seq,v);
            assert_eq!(s, path_distance(v0,vn, path.vflat, a, path.t_total));
            assert_eq!(s, seq.acc_distance(i32::from(v0)).0);
        }
    }
}


#[cfg(test)]
fn check_shortest_curve(ds:i64, v0: i64, vn: i64, a: i64, vmax: i64)
{
    println!("check_shortest_curve: ds:{}, v0: {}, vn: {}, a: {}, vmax: {}", 
             ds, v0,vn, a, vmax);
    let path = shortest_curve(ds,v0,vn,a, vmax);
    println!("path: {}", path);
    assert!((path.distance() - ds).abs() <= 1);
    assert!(path.vflat.abs() <= vmax);
    assert!((path.vflat + (if path.t_adjust == 0 {
        0
    } else{
        if path.t_adjust > 0 {
            1
        } else {
            -1
        }})).abs() <= vmax);
    let t1 = ((path.vflat - path.v0) + a - 1) / a;
    let t4 = ((path.vflat - path.vn) + a - 1) / a;
    assert!(t1 + t4 <= path.t_total);
    assert!(path.t_total > path.t_adjust);
    assert!(path.t_total - (t1 + t4) >= path.t_adjust);

    if path.t_adjust != 0 {
        // Check that there is room to make adjustments
        let mut t_flat = path.t_total - t1 -t4;
        if t1 == 0 {
            t_flat -= 1;
        }
        if t4 == 0 {
            t_flat -= 1;
        }
        
        assert!(t_flat >=1);
    }
    let seq = path.acc_seq();
    println!("{} => {:?} => {}", path, seq,vn);
    assert_eq!(path.distance(), seq.acc_distance(i32::from(v0)).0);

}

#[test]
fn test_shortest_curve()
{
    for ds in -40..40 {
        for vmax in 1..7 {
            for v0 in -vmax..(vmax+1) {
                for vn in -vmax..(vmax+1) {
                    for a in 1..6 {
                        check_shortest_curve(ds, v0, vn,a,vmax);
                    }
                }
            }
        }
    }
}




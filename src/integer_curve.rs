use solve;
use num::Integer;
use std::fmt;
use std::clone::Clone;
use multi_range::MultiRange;
use std::i64;

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
    v0: i64, // Start speed
    vn: i64, // End speed
    amax:i64, // Maximum acceleration
    vflat:i64, // Flat part of speed curve, may be between v0 and vn
    t_total:i64, // Total time for path
    t_adjust: i64, // Adjust the flat part of the curve so that the speed is increased or decreased by one. The sign decides the direction and the absolute value is the length
}

impl DistancePath {
    fn distance(&self) -> i64 {
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
    // 2 a s_low = cx2 x^2 + cx x + c;
    let cx2 = a*a;
    let cx = 2*a*(v0+vn);
    let c = -(v0-vn)*(v0-vn) - a*a;
    // Positive acceleration
    match solve::find_roots(cx2, cx, -c + 2*a*ds) {
        None => {
            r.or(tmin, i64::MAX);
        },
        Some((min,max)) => {
            r.or(i64::MIN, min);
            r.or(max, i64::MAX);
            r = r.and(tmin, i64::MAX);
        }
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
    println!("pos s0 = {}", s0);
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
    let (t1, dv1) = div_rem_away(vflat - v0, a);
    let (t4, dv4) = div_rem_away(vflat - vn, a);
    let t1 = t1.abs();
    let t4 = t4.abs();
    return (vflat-dv1 +v0) * t1  + dv1 
        + 2*vflat* (t - t1 - t4) 
        + dv4 + (vflat-dv4 + vn) * t4;
}

fn adjust_curve(path: &DistancePath, ds: i64) -> DistancePath
{
    let a = path.amax;
    let v0 = path.v0;
    let vn = path.vn;
    println!("adjust_curve: ds = {}, path = {}",ds, path); 
    // Find highest and lowest possible vflat
    let tmid_low = (path.t_total * a - (v0 - vn)) / (2*a);
    let tmid_high = (path.t_total * a - (vn - v0)) / (2*a);
    let mut max = max_i64(v0 + tmid_low * a, vn + tmid_high * a);
    let mut min = min_i64(v0 + tmid_high * -a, vn + tmid_low * -a);
    assert!(max >= min);
    let mut smax = path_distance(v0,vn,max, a, path.t_total);
    let mut smin = path_distance(v0,vn,min, a, path.t_total);
    println!("{}, {} -> {}, {} -> {}, {}",tmid_low, tmid_high,min,max, smin,smax); 
    assert!(smax >= ds);
    assert!(smin <= ds);
    // Find closest min and max
    loop {
        //println!("({} -> {}) - ({} -> {})", min, smin, max, smax);
        if (max - min) <= 1 {break};
        let mid = (max + min) / 2;
        let s = path_distance(v0,vn,mid, a, path.t_total);
        if s >= ds {
            max = mid;
            smax = s;
        } else {
            min = mid;
            smin = s;
        }
    }
    let mut adjusted = path.clone();
    if smax == ds {
        adjusted.vflat = max;
    } else if smin == ds {
        adjusted.vflat = min;
    } else {
        if smax >= 0 {
            adjusted.vflat = max;
            adjusted.t_adjust = (ds - smax) / 2;
        } else {
            adjusted.vflat = min;
            adjusted.t_adjust = (ds - smin) / 2;
        }
    }
        
    return adjusted;
}

pub fn shortest_curve(ds:i64, v0: i64, vn: i64, a: i64, vmax:i64) -> DistancePath
{
    
    let mut pos = speed_limited_curve(ds, v0, vn, a, vmax);
    if pos.empty() {
        pos = max_acceleration_curve(ds, v0, vn, a);
    };
    let mut neg = speed_limited_curve(-ds, -v0, -vn, a, vmax);
    if neg.empty() {
        neg = max_acceleration_curve(-ds, -v0, -vn, a);
    };
    println!("pos: {:?} neg: {:?}", pos, neg);
    let (t,_) = pos.and_range(&neg).bounds().unwrap();
    let path = DistancePath {v0:v0, vn: vn, amax:a, vflat:0, t_total:t, t_adjust:0};
    return adjust_curve(&path, ds);
}

#[cfg(test)]
fn check_max_acceleration_curve(ds:i64, v0: i64, vn: i64, a: i64)
{
    println!("\nds= {}, v0= {}, vn= {}, a= {}", ds,v0,vn,a);
    let range = max_acceleration_curve(ds,v0,vn,a);
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
    let path = DistancePath {v0: -5, vn: 3, amax: 4, vflat: 0, t_total: 6, t_adjust: 0};
    let ds = 23;
    let adjusted = adjust_curve(&path, ds);
    println!("path = {}, ds = {}", adjusted,ds);
    assert!((adjusted.distance() - ds).abs() <= 1);
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
            assert_eq!(path.distance(), path_distance(v0,vn, path.vflat, a, path.t_total));
        }
    }
}

#[cfg(test)]
fn check_shortest_curve(ds:i64, v0: i64, vn: i64, a: i64, vmax: i64)
{
    println!("check_shortest_curve: ds:{}, v0: {}, vn: {}, a: {}, vmax: {}", 
             ds, v0,vn, a, vmax);
    let path = shortest_curve(ds,v0,vn,a, vmax);
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

}

#[test]
fn test_shortest_curve()
{
    for ds in -20..20 {
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


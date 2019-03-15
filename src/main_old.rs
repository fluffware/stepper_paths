
extern crate libc;
extern crate core;
extern crate paths;
use paths::solve;

use std::env;

macro_rules! permutate {
    ( $f:ident ( $v:expr ) ) => {
        for &q in $v {$f(q)}
    };
    ( $f:ident ( $v:expr, $($r:expr),+ ) ) => {
        for &q in $v {permutate!($f($($r),*;q))} 
    };
    ( $f:ident ( $v:expr, $($r:expr),+ ; $($a:ident),+ )) => {
        for &q in $v {permutate!($f($($r),+;$($a),*,q))} 
    };
    ( $f:ident ( $v:expr ; $($a:expr),+ ) ) => { 
        for &q in $v {$f($($a),+,q)} 
    };
}


#[cfg(test)]
fn eval_acc_seq(s0:i32, v0:i32,  acc: &[i32]) -> (i32, i32)
{
    let mut s = s0;
    let mut v = v0;
    for a in acc {
        s = s + 2*v + a;
        v = v + a;
    }
    (s,v)
}

#[cfg(test)]
fn eval_acc_seq2(s0:i32, v0:i32,  acc: &[i32]) -> (i32, i32)
{
    let mut v = v0;
    let n = acc.len();
    let mut s = s0 + 2 * (n as i32) * v0 ;
    for (k, a) in acc.iter().enumerate()
    {
        v = v + a;
        s = s + (2*((n - k) as i32) - 1) * a;
    }
    (s,v)
    
}

#[allow(dead_code)]
fn eval_acc_time(s0:f64, v0:f64,  acc_time: &[(f64,f64)]) -> (f64, f64)
{
    let mut s = s0;
    let mut v = v0;

    for &(a,dt) in acc_time {
        s = s + 2.0*v * dt + a*dt*dt;
        v = v + a * dt;
    }
    (s,v)
}



#[allow(dead_code)]
fn build_curve(ds: i32, v0: i32, vn:i32, amax:i32, vmax:i32, n:i32)
{
    let dt1 = ((vn - v0) + amax - 1) / amax;
    let ds1h = 2*(dt1-1) * v0 + amax * (dt1*dt1 - dt1) + vn + v0;
    let ds1l = 2*(dt1-1) * vn - amax * (dt1*dt1 - dt1) + vn + v0;

    let dt2 :i32;
   
    println!("ds1h = {} ds1l = {} dt1: {}", ds1h, ds1l, dt1);
}

#[allow(dead_code)]
// ds >= 0, a > 0 , vmax > 0
fn extend_dist(s: f64, v:f64, a: f64, vmax:f64) -> (f64, f64)
{
    // Constant speed time
    let mut dtc = 0.0;
    // Acceleration time
    let mut dta = - v / a + f64::sqrt(s/(2.0 * a) + v*v/(a*a));
    let vtop = v + a*dta;
    println!("v: {} s: {} a: {} vmax: {}", v,s,a,vmax);
    if vtop > vmax {
        dta = (vmax - v)/a;
        dtc = (s - 2.0*(v +vmax)*dta) / (2.0 * vmax);
    } else if vtop < -vmax {
        dta = (vmax - v)/a;
        dtc = (s - 2.0*(v - vmax)*dta) / (2.0 * -vmax);
    }
    (dta, dtc)
}

struct DistancePath {
    amax:i64, // Start acceleration
    t_maxacc: i64, // Maximum acceleration and deceleration time
    /* Acceleration just before and after constant velocity.
    If t_acc_step > 0 then the acceleration before is one greater */
    acc_limited: i64, 
    t_const: i64, // Constant speed time
    t_accstep: i64, // Time offset into constant velocity when velocity decreases by one
}

/*
s = 2 * amax * t_maxacc * t_maxacc + 2*acc_limited + 2*(amax*t_accstep + acc_limited) * t_const + 2*t_accstep

 */
fn build_distance_path(ds:i64,  v0:i64, vmax:i64, amax:i64)
    {
        if ds > 0 {
            match solve::find_roots(amax, 4*v0, 2*ds) {
                Some((min, max)) => {
                    let dt1 = 
                    if min >= 0 {
                            min
                    } else if max >= 0 {
                        max
                    } else {
                        0
                    };
                if 2 * dt1 * v0 > ds {
                    /* Speeds between v0 and vn */
                }
                
            },
            None => {}
        };
        
    } else {
        solve::find_roots(amax, -4*v0, -2*ds);
    }
}

fn build_path(ds:i64,  v0:i64, vn:i64, vmax:i64, amax:i64)
{
    // Try positive acceleration
    match solve::find_roots(amax * amax, 2*amax*(v0+vn), 
                     2*amax*ds + (vn - v0)*(vn - v0)) {
        Some((t1, t2)) => {
            let dshi = (v0+vn) * t2 - (vn - v0)*(vn - v0) 
            println!("ds: {} dt: {}, {}",ds, t1,t2);
        },
        None => {
            println!("No solution");
        }
    }
    // Try negative acceleration
    match solve::find_roots(-amax * amax, 2*amax*(v0+vn), 
                     2*amax*ds - (vn - v0)*(vn - v0)) {
        Some((t1, t2)) => {
            println!("ds: {} dt: {}, {}",ds, t1,t2);
        },
        None => {
            println!("No solution");
        }
    }
}

    

#[allow(dead_code)]
fn build_optimal_path(s0:i32, sn:i32, v0:i32, vn:i32, vmax:i32, amax:i32)
                      ->(f64, [f64;3])
{
    let mut a = (if (vn - v0) >= 0 {amax} else {-amax}) as f64;
    let mut dt = [0.0,0.0,0.0];
    dt[0] = (vn - v0) as f64 / a;
    let s = s0 as f64  + (2 * v0) as f64 * dt[0] + a * dt[0] * dt[0];
    let ds = sn as f64 - s;
    println!("ds: {} a: {}", ds, a);
    if ds < 0.0 {
        if a > 0.0 {
            // start
            let (dta, dtc) = extend_dist(-ds, -v0 as f64, a, vmax as f64);
            dt[2] = dt[0] + dta;
            dt[0] = dta;
            dt[1] = dtc;
            a = -a;
        } else {
            // end
            let (dta, dtc) = extend_dist(-ds, -vn as f64, -a, vmax as f64);
            dt[0] = dt[0] + dta;
            dt[2] = dta;
            dt[1] = dtc;
        }
    } else {
        if a > 0.0 {
            //end
            let (dta, dtc) = extend_dist(ds, vn as f64, a, vmax as f64);
            dt[0] = dt[0] + dta;
            dt[2] = dta;
            dt[1] = dtc;
        } else {
            //start
            let (dta, dtc) = extend_dist(ds, v0 as f64, -a, vmax as f64);
            dt[2] = dt[0] + dta;
            dt[0] = dta;
            dt[1] = dtc;
            a = -a;
        }
    }
    (a, dt)
}

#[allow(dead_code)]
fn build_optimal_int_path(s0:i32, sn:i32, v0:i32, vn:i32, vmax:i32, amax:i32)
                      ->(i32, Vec<i32>)
{
    let (a,dt) = build_optimal_path(s0,sn, v0,vn,vmax, amax);
    let vtop =
        if dt[1] == 0.0 {
            (v0 as f64 + dt[0]*a).floor() as i32
        } else {
            vmax
        };
    let dt: [i32;3] = [dt[0].ceil() as i32, dt[1].ceil() as i32, dt[2].ceil() as i32];
    let a = if a >= 0.0 {amax} else {-amax};
    let mut acc = Vec::with_capacity((dt[0]+dt[1]+dt[2]) as usize);
    for av in &mut acc[0..(dt[0] - 2) as usize] {
        *av = a;
    }
    (0, acc)
}

#[allow(dead_code)]
fn check_path(s:(i32,i32), v:(i32, i32), vmax:i32, amax:i32)
{
    let (s0,sn) = s;
    let (v0,vn) = v;
    println!("s: {} -> {} v: {} -> {}", s0,sn,v0,vn);
    let (a, dt) = build_optimal_path(s0,sn, v0, vn, vmax, amax);
    println!("a: {} dta1: {} dtc: {} dta2: {}", a, dt[0], dt[1], dt[2]);
    let (sn, vn) = eval_acc_time(s0 as f64, v0 as f64,
                                 &[(a,dt[0]), (0.0,dt[1]), (-a,dt[2])]);
    println!("sn: {} vn: {}", sn, vn);
    println!("");
}
    

fn main() {
    let args = env::args();
    let mut args = args.skip(1);
    let ds = args.next().unwrap().parse::<i64>().unwrap();
    let v0 = args.next().unwrap().parse::<i64>().unwrap();
    let vn = args.next().unwrap().parse::<i64>().unwrap();

    build_path(ds,v0,vn,10,3);
}

#[test]
fn seq_eval()
{
    struct TestVec {
        s:i32,
        v: i32,
        acc:Vec<i32>}
    let v = [TestVec{s:3,v:7,acc:vec![3,3,4,3,13]},
             TestVec{s:-1,v:-4,acc:vec![0,-5,1,0,0,2,5,2,1,0,0]},
             TestVec{s:1,v:-4,acc:vec![0,0]},
             TestVec{s:-1,v:4,acc:vec![1,-5,4,0,0,2]}
    ];
    for &TestVec{ref s,ref v,ref acc} in v.iter() {
        let r1 = eval_acc_seq(*s, *v, acc);
        let r2 = eval_acc_seq2(*s, *v, acc);
        assert_eq!(r1,r2);
    }
}

#[cfg(test)]
fn nearly<T>(a:T, b:T, e:T) -> bool
    where T: std::ops::Sub<Output = T> + std::cmp::PartialOrd + std::ops::Neg<Output = T> + std::convert::From<i32>
{
    let d: T = a-b;
    if d >= T::from(0) {d < e} else {-d < e}
}

#[cfg(test)]
    fn check_extend_dist(mut ds:f64, mut v0: f64, a:f64, vmax: f64)
{
    if ds < 0.0 {
        ds = -ds;
        v0 = -v0;
    }
    let (dta, dtc) = extend_dist(ds, v0, a, vmax);
    //println!("dta: {} dtc: {}",dta, dtc);
    let seq = [(a, dta), (0.0,dtc), (-a, dta)];
    let (ds_e, vn_e) = eval_acc_time(0.0, v0, &seq);
    //println!("sn: {} vn: {}",ds_e, vn_e);
    assert!(nearly(ds_e, ds, 1e-6) && nearly(vn_e, v0, 1e-6));
}




#[test]
fn extend()
{
    let a_vec = [0.1, 1.2, 4.0];
    let vmax_vec = [3.0, 1.0 ,0.1];
    let v0_vec = [-21.0, -4.0, 0.0, 2.0, 5.7, 12.5];
    let ds_vec = [-23.0, -4.2, 0.0, 0.2, 6.1,23.0];
  
    permutate!(check_extend_dist(&ds_vec, &v0_vec, &a_vec, &vmax_vec));

}
    

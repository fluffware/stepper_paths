
extern crate libc;
extern crate core;
extern crate paths;
//use paths::solve;
use paths::integer_curve;

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

    

fn main() {
    let args = env::args();
    let mut args = args.skip(1);
    let ds = args.next().unwrap().parse::<i64>().unwrap();
    let v0 = args.next().unwrap().parse::<i64>().unwrap();
    let vn = args.next().unwrap().parse::<i64>().unwrap();
    let vmax = 10;
    let a = 5;
    /*let (n, a2) = integer_curve::max_acceleration_curve(ds, v0, vn, a);
    println!("n= {}  a= {}", n, a2);*/
    let path = integer_curve::shortest_curve(ds,v0,vn,a,vmax);
    println!("Curve: {}", path);
}

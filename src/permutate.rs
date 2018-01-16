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

fn func(x:f64, y:f64,s:&str, z:i32)
{
    println!("{} {} {} {}",x,y,s,z);
}
fn main()
{
    let ds_vec = [-23.0, -4.2, 0.0, 0.2, 6.1,23.0];
    let v0_vec = [-21.0, -4.0, 0.0, 2.0, 5.7, 12.5];
    let s = ["A", "B"];
    let a_vec = [-2, 0, 2];
    permutate!(func(&ds_vec, &v0_vec, &s, &a_vec));
}

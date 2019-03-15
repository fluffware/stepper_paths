
extern crate libc;
extern crate core;
extern crate paths;
use paths::integer_curve;
use paths::acc_vector::AccSegment;
use std::io::Write;
use std::fmt::Write as StrWrite;

extern crate num;
use std::env;

const PREFIX:&'static str = "
void draw_axis(real taxismax, 
	       real vaxismin, 	       
	       real vaxismax) {
  for (real y = vaxismin; y <= vaxismax; ++y) { 
    draw((0, y) -- (taxismax, y),gray +linewidth(0.1pt));
  }
  for (real x = 0; x <= taxismax; ++x) { 
  draw((x, vaxismin) -- (x, vaxismax),gray +linewidth(0.1pt));
  }
  
  draw((0,vaxismin) -- (0,vaxismax), arrow=Arrow);
  draw((0,0) -- (taxismax,0), arrow=Arrow);
}

size(10cm);
real taxismax = 12;
real vaxismax = 7;
real vaxismin = -7;

import curves;

draw_axis(taxismax, vaxismin, vaxismax);

import animation;

animation a;
";

const SUFFIX:&'static str = "
for (int i=0; i < c.length; ++i) { 
  save();
  draw(c[i], red);
  label(Label(format(\"$s=%d$\", dist[i]), align= LeftSide), position=(8,-6), black);
  a.add();
  restore();
}

erase();
a.movie(BBox(0.25cm),delay=100);
";
fn main() {
    let args = env::args();
    let mut args = args.skip(1);
    let v0 = args.next().unwrap().parse::<i64>().unwrap();
    let vn = args.next().unwrap().parse::<i64>().unwrap();
    let vmax = args.next().unwrap().parse::<i64>().unwrap();
    let a = args.next().unwrap().parse::<i64>().unwrap();
    let mut f = std::io::stdout();
    f.write(PREFIX.as_bytes()).unwrap();
    f.write("path c[] = { ".as_bytes()).unwrap();
    let mut dist = String::from("int dist[] = {\n");
    for ds in -70i64 .. 70i64 {
        let path = integer_curve::shortest_curve(ds,v0,vn,a,vmax);
        let seq = path.acc_seq();
        let mut t = 0;
        let mut v = v0 as i32;
        write!(f, "({}, {})", t, v).unwrap();
        for &AccSegment{interval: i, acc: a} in &seq {
            t += i;
            v += a as i32 * i as i32;
            write!(f, " -- ({}, {})", t, v).unwrap();
        }
        f.write(",\n".as_bytes()).unwrap();
        write!(dist, "{},", ds).unwrap();
    }
    dist += "};\n";
    f.write("};\n".as_bytes()).unwrap();
    f.write(dist.as_bytes()).unwrap();
    f.write(SUFFIX.as_bytes()).unwrap();
}

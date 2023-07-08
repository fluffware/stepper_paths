use coords::{Point, Vector};
#[derive(Clone, Debug, PartialEq)]
pub struct BezierSegment {
    pub p: Point, // Start point of this segment, and end point of previous segment
    pub c: Vector, // Control point relative to p
}

impl BezierSegment {
   

}

#[derive(Clone, Debug, PartialEq)]
pub struct BezierSegments(Vec<BezierSegment>);
    
impl BezierSegments{
    pub fn new<C1, C2, P1, P2>(p1: P1, c1: C1, c2: C2, p2: P2) -> BezierSegments
    where
        C1: Into<Vector>,
        C2: Into<Vector>,
        P1: Into<Vector>,
        P2: Into<Vector>,
    {
	BezierSegments(vec!{BezierSegment{p: p1.into(), c: c1.into()},
	     BezierSegment{p: p2.into(), c: -c2.into()}})
     }
    
    pub fn split(&mut self, index: usize, t: f64) {
	let s1 = &self.0[index];
	let s2 = &self.0[index+1];
	let p0 = s1.p;
	let p1 = s1.p + s1.c;
	let p2 = s2.p - s2.c;
	let p3 = s2.p;
	let p00 = p0 + s1.c * t;
        let p01 = p1 + (p2 - p1) * t;
        let p02 = p2 - s2.c * (1.0 - t);
        let p10 = p00 + (p01 - p00) * t;
        let p11 = p01 + (p02 - p01) * t;
        let p20 = p10 + (p11 - p10) * t;
	//assert_eq!(p20.x - p10.x, p11.x - p20.x);
	//assert_eq!(p20.y - p10.y, p11.y - p20.y);

	self.0[index].c = p00 - p0;
	self.0.insert(index + 1,BezierSegment{p: p20, c: p20-p10});
	self.0[index + 2].c =  p2-p20;
    }
}

#[cfg(test)]
use approx::{AbsDiffEq, RelativeEq};

#[test]
fn test_split()
{
    let mut segs = BezierSegments::new((1,0), (0,1), (0,1), (3,0));
    segs.split(0,0.5);
    println!("{segs:?}");
    let mut segs = BezierSegments::new((1,0), (0,1), (0,1), (3,0));
    segs.split(0,0.3);
    println!("{segs:?}");
    /*
    assert_relative_eq!(segs.0[0].p.x, 1.0);
    assert_relative_eq!(segs.0[0].p.y, 0.0);
    assert_relative_eq!(segs.0[0].c.x, 0.0);
    assert_relative_eq!(segs.0[0].c.y, 0.5);
    
    assert_relative_eq!(segs.0[1].p.x, 2.0);
    assert_relative_eq!(segs.0[1].p.y, 0.5);
    assert_relative_eq!(segs.0[1].c.x, 0.0);
    assert_relative_eq!(segs.0[1].c.y, 0.0);
    */
}
    

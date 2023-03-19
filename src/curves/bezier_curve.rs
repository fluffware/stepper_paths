use coords::{Point, Vector};

/// A bezier curve that starts at (0,0) and ends at p2
/// c1 is the control point relative (0,0)
/// c2 is the control point relative p2
#[derive(Clone, Debug, PartialEq)]
pub struct BezierSegment {
    pub c1: Vector,
    pub c2: Vector,
    pub p2: Point,
}

impl BezierSegment {
    pub fn new<C1, C2, P2>(c1: C1, c2: C2, p2: P2) -> BezierSegment
    where
        C1: Into<Vector>,
        C2: Into<Vector>,
        P2: Into<Vector>,
    {
        BezierSegment {
            c1: c1.into(),
            c2: c2.into(),
            p2: p2.into(),
        }
    }

    pub fn point(&self, t: f64) -> Point {
        let t2 = t * t;
        let ti = 1.0 - t;
        let ti2 = ti * ti;
        self.c1 * (3.0 * ti2 * t) + (self.p2 + self.c2) * (3.0 * ti * t2) + self.p2 * (t2 * t)
    }

    pub fn split(&self, t: f64) -> (BezierSegment, BezierSegment) {
        let p00 = self.c1 * t;
        let p01 = self.c1 + (self.p2 + self.c2 - self.c1) * t;
        let p02 = self.p2 + self.c2 * (1.0 - t);
        let p10 = p00 + (p01 - p00) * t;
        let p11 = p01 + (p02 - p01) * t;
        let p20 = p10 + (p11 - p10) * t;
        (
            BezierSegment::new(p00, p10 - p20, p20),
            BezierSegment::new(p11 - p20, p02 - self.p2, self.p2 - p20),
        )
    }

    pub fn to_svg_path(&self) -> String {
        format!(
            "c{},{} {},{} {},{}",
            self.c1.x,
            self.c1.y,
            self.c2.x + self.p2.x,
            self.c2.y + self.p2.y,
            self.p2.x,
            self.p2.y
        )
    }
}

#[test]
fn test_point() {
    let b1 = BezierSegment::new(
        Point { x: 1.0, y: 3.0 },
        Point { x: -2.0, y: -1.0 },
        Point { x: 5.0, y: -2.0 },
    );
    assert_eq!(b1.point(0.0), Point { x: 0.0, y: 0.0 });
    assert_eq!(b1.point(1.0), Point { x: 5.0, y: -2.0 });
    assert_eq!(
        b1.point(0.5),
        Point {
            x: ((1.0 + 3.0) * 3.0 + 5.0) * 0.5 * 0.5 * 0.5,
            y: ((3.0 - 3.0) * 3.0 - 2.0) * 0.5 * 0.5 * 0.5
        }
    );
    assert_eq!(
        b1.point(0.25),
        Point { x: 1.0, y: 3.0 } * 3.0 * 0.75 * 0.75 * 0.25
            + Point { x: 3.0, y: -3.0 } * 3.0 * 0.75 * 0.25 * 0.25
            + Point { x: 5.0, y: -2.0 } * 0.25 * 0.25 * 0.25
    );
}

#[test]
fn test_split() {
    let b1 = BezierSegment::new(Vector::from((1.0, 3.0)), (-2.0, -1.0), (5.0, -2.0));
    let (b11, b12) = b1.split(0.25);

    assert_eq!(
        b11,
        BezierSegment::new((0.25, 0.75), (-0.359375, 0.125), (0.921875, 0.8125))
    );
    assert_eq!(
        b12,
        BezierSegment::new((1.078125, -0.375), (-1.5, -0.75), (4.078125, -2.8125))
    );
    println!("b11: {} b12: {}", b11.to_svg_path(), b12.to_svg_path());
}

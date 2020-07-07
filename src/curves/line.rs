use curve_approx::CurveInfo;
use coords::{Point,Vector};

pub struct Line
{
    dir: Vector,
    length: f64
}

impl Line {
    pub fn new(end: Point) -> Line
    {
        let length = end.length();
        let dir = if length == 0.0 {
            Vector {x: 1.0, y: 0.0} // Arbitrary unit vector
        } else {
            end.unit()
        };
        Line{dir,
             length
        }
    }
}

impl CurveInfo for Line {
    fn length(&self) -> f64
    {
        self.length
    }

    fn value(&self, pos: f64) -> (Point, Vector)
    {
        (self.dir * pos, self.dir)
    }
}

use coords::{Point, Vector};
use curve_approx::CurveInfo;

struct CurveSegment
{
    end_coord: Point,
    end_len: f64,
    curve: Box<dyn CurveInfo>
}

impl std::fmt::Debug for CurveSegment
{
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result
    {
     write!(f," ->({}, {}) @{}", self.end_coord.x, self.end_coord.y, self.end_len)   
    }
}

#[derive(Debug)]
pub struct ConcatCurve
{
    segments: Vec<CurveSegment>
}

impl ConcatCurve
{
    pub fn new() -> ConcatCurve
    {
        ConcatCurve{segments: Vec::new()}
    }

    pub fn add(&mut self, curve: Box<dyn CurveInfo>)
    {
        let clen = curve.length();
        let end_coord = match self.segments.last() {
            None => Point{x: 0.0, y: 0.0},
            Some(seg) => seg.end_coord

        };
        
        let seg = CurveSegment {
            end_coord: curve.value(clen).0 + end_coord,
            end_len: clen + self.length(),
            curve
        };
        self.segments.push(seg);
    }

    pub fn clear(&mut self) {
        self.segments.clear();
    }

    pub fn is_empty(&self) -> bool
    {
        self.segments.is_empty()
    }
        
}    

impl Default for ConcatCurve
{
    fn default() -> ConcatCurve
    {
        ConcatCurve::new()
    }
}

impl CurveInfo for ConcatCurve {
    fn length(&self) -> f64
    {
        if let Some(s) = self.segments.last()
        {
            s.end_len
        } else {0.0}
    }
    
    fn value(&self, pos: f64) -> (Point, Vector)
    {
        match &self.segments.binary_search_by(
            |s| s.end_len.partial_cmp(&pos).unwrap()) {
            Ok(i) | Err(i) => {
                let i = (*i).min(self.segments.len() - 1);
                let seg = &self.segments[i];
                let (start_coord, start_len) =
                    if i  > 0 {
                        let prev = &self.segments[i - 1];
                        (prev.end_coord, prev.end_len)
                    } else {
                        (Point{x: 0.0, y: 0.0}, 0.0)
                    };
                let (point, deriv) = seg.curve.value(pos - start_len);
                (point + start_coord, deriv)
            }
        }
        
    }
}

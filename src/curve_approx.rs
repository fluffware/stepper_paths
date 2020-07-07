use coords::{Vector,Point};

pub trait CurveInfo
{
    /// Returns total length of curve
    fn length(&self) -> f64;
    
    /// Returns coordinates and direction at given position along the curve.
    /// Direction is normalized. If `pos` is 0.0 then the return
    /// cordinates should be (0,0), i.e the curve always starts from
    /// origo.
    ///
    /// # Arguments
    ////
    /// * `pos` - Position on curve in the range 0.0..length(). This is the
    /// distance from the start along the curve
    fn value(&self, pos: f64) -> (Point, Vector);
}

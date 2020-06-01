
pub trait CurveInfo
{
    /// Returns total length of curve
    fn length(&self) -> f64;
    /// Returns coordinates and direction at given relative length.
    /// Direction is normalized. If `pos` is 0.0 then the return
    /// cordinates should be (0,0), i.e the curve always starts from
    /// origo.
    ///
    /// # Arguments
    ////
    /// * `pos` - Position on curve 0.0..1.0. Must be proportional to
    /// curve length.
    fn value(&self, pos: f64) -> ([f64;2], [f64;2]);
}

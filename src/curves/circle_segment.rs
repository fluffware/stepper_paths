use curve_approx::CurveInfo;
use coords::{Point,Vector};
use std::f64::consts::PI;

pub struct CircleSegment
{
    radius: f64,
    start_angle: f64,
    end_angle: f64,
    center: Point,
    angle_scale: f64
}

impl CircleSegment {
    pub fn new(radius: f64,
               start_angle: f64,
               end_angle: f64) -> CircleSegment
    {
        CircleSegment{
            radius,
            start_angle,
            end_angle,
            center: {
                let (ss,sc) = start_angle.sin_cos();
                Point{x: -sc * radius,
                      y: -ss * radius}
            },
            angle_scale: {
                if start_angle <= end_angle {
                    1.0 / radius
                } else {
                    -1.0 / radius
                }
            }
        }
    }

    pub fn new_start_direction(end: Point,
                               start_direction: Vector) -> CircleSegment
    {
        let p0 = end * 0.5;
        let p1 = end.rotate_90_ccw();
        let p2 = start_direction.rotate_90_ccw();
        let b =  (p1.x * p0.y - p0.x * p1.y) / (p1.x * p2.y - p2.x * p1.y);
        let center = p2 * b;
        let radius = center.length();
        let start_angle = (-center).angle();
        let mut end_angle = (end - center).angle();
        if start_direction.diff_sign(&center) < 0.0 {
            // clockwise
            if end_angle < start_angle {
                end_angle += 2.0*PI;
            }
        } else {
            // counter-clockwise
             if end_angle > start_angle {
                end_angle -= 2.0*PI;
            }
        }
            
        CircleSegment{
            radius,
            start_angle, 
            end_angle,
            center,
            angle_scale: {
                if start_angle <= end_angle {
                    1.0 / radius
                } else {
                    -1.0 / radius
                }
            }
        }
    }
    
}

impl CurveInfo for CircleSegment {
    fn length(&self) -> f64
    {
        self.radius * (self.end_angle - self.start_angle).abs()
    }
   
    fn value(&self, pos: f64) -> (Point, Vector)
    {
        let a = pos * self.angle_scale + self.start_angle;
        let (s,c) = a.sin_cos();
        (Point{x: c*self.radius, y: s*self.radius} + self.center,
         if self.start_angle <= self.end_angle {
             Vector{x: -s, y: c}
         } else {
             Vector{x: s, y: -c}
         })
    }
}

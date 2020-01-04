use std;
use std::fmt;

#[derive(Debug)]
pub struct Point
{
    pub x: f64,
    pub y: f64
}
impl std::fmt::Display for Point
{
    fn fmt(self: &Self,f: &mut fmt::Formatter) -> fmt::Result
    {
        write!(f, "({}, {})", self.x, self.y)
    }
}

impl std::cmp::PartialEq for Point
{
    fn eq(&self, t: &Self) ->bool
    {
        self.x == t.x && self.y == t.y
    }
}
impl Clone for Point
{
    fn clone(&self) -> Point {
        Point{x: self.x, y: self.y}
    }
}

impl Copy for Point
{
}
impl std::ops::Add<Point> for Point {
    type Output = Point;
    fn add(self, v: Point) -> Point {
        Point { x: self.x + v.x,
                y: self.y + v.y}
    }
}
impl std::ops::AddAssign<Point> for Point {
    fn add_assign(&mut self, v: Point)  {
        self.x += v.x;
        self.y += v.y;
    }
}

impl std::ops::Sub<Point> for Point {
    type Output = Point;
    fn sub(self, v: Point) -> Point {
        Point { x: self.x - v.x,
                y: self.y - v.y}
    }
}
impl std::ops::SubAssign<Point> for Point {
    fn sub_assign(&mut self, v: Point)  {
        self.x -= v.x;
        self.y -= v.y;
    }
}
impl std::ops::Mul<f64> for Point {
    type Output = Point;
    fn mul(self, s: f64) -> Point {
        Point { x: self.x * s,
                y: self.y * s}
    }
}
#[derive(Debug)]
pub struct Transform {
    pub matrix : [f64;6]
}

fn matrix_mul(a: &[f64; 6], b: &[f64; 6]) -> [f64; 6]
{
    [a[0] * b[0] + a[2]*b[1],
     a[1] * b[0] + a[3]*b[1],
     a[0] * b[2] + a[2]*b[3],
     a[1] * b[2] + a[3]*b[3],
     a[0] * b[4] + a[2]*b[5] + a[4],
     a[1] * b[4] + a[3]*b[5] + a[5]]
}

impl Transform {
    pub fn new(m :&[f64; 6]) -> Transform {
        Transform{matrix: *m}
    }
    pub fn identity() -> Transform {
        Transform{matrix:[1.0, 0.0, 0.0, 1.0, 0.0, 0.0]}
    }
    
    pub fn translate(x: f64, y: f64) -> Transform {
        Transform{matrix:[1.0, 0.0, 0.0, 1.0, x, y]}
    }
    
    pub fn scale(s: f64) -> Transform {
        Transform{matrix:[s, 0.0, 0.0, s, 0.0, 0.0]}
    }
    
    pub fn scale_xy(sx: f64, sy: f64) -> Transform {
        Transform{matrix:[sx, 0.0, 0.0, sy, 0.0, 0.0]}
    }

    pub fn rotate(a: f64) -> Transform {
        let (s,c) = a.sin_cos();
        Transform{matrix:[c, s, -s, c, 0.0, 0.0]}
    }

    pub fn rotate_around(a: f64, pivot: &Point) -> Transform {
        let (s,c) = a.sin_cos();
        Transform{matrix:[c, s, -s, c, 
                          pivot.x*(1.0-c) + pivot.y*s, 
                          -pivot.x*s + pivot.y*(1.0-c)]}
    }
    
    pub fn skew_x(a: f64) -> Transform
    {
        Transform{matrix:[1.0, 0.0, a.tan(), 1.0, 0.0, 0.0]}
    }
    
    pub fn skew_y(a: f64) -> Transform
    {
        Transform{matrix:[1.0, a.tan(), 0.0, 1.0, 0.0, 0.0]}
    }
}
impl std::cmp::PartialEq for Transform
{
    fn eq(&self, t: &Self) ->bool
    {
        self.matrix == t.matrix
    }
}
impl Clone for Transform
{
    fn clone(&self) -> Transform {
        Transform {matrix: self.matrix}
    }
}

impl Copy for Transform
{
}

impl std::ops::Mul for Transform {
    type Output = Transform;
    fn mul(self, t: Self) -> Self {
        Transform  {matrix: matrix_mul(&self.matrix, &t.matrix)}
    }
}

impl std::ops::Mul<Point> for Transform {
    type Output = Point;
    fn mul(self, v: Point) -> Point {
        Point { x: self.matrix[0] * v.x + self.matrix[2] * v.y + self.matrix[4],
                y: self.matrix[1] * v.x + self.matrix[3] * v.y + self.matrix[5]}
    }
}

#[cfg(test)]
fn assert_matrix_eq(a:&[f64;6], b:&[f64;6])
{
    for (a,b) in a.iter().zip(b) {
        if (a-b).abs() > 1e-5 {
            panic!("{} != {}", a,b);
        }
    }
}

#[cfg(test)]
fn assert_transform_eq(a:&Transform, b:&Transform)
{
    assert_matrix_eq(&a.matrix, &b.matrix);
}

#[test]
fn test_transform()
{
    let a = Transform::identity();
    let b = Transform::identity();
    assert_eq!(a*b, Transform::identity());
    let a = Transform::new(&[1.3, 7.1, -23.0, 8.0, 1.45, -12.7]);
    let b = Transform::new(&[4.3, 0.1, 23.0, 18.7, 14.5, 2.7]);
    assert_transform_eq(&(a*b), &Transform::new(&[3.29,31.33,-400.2,312.9,-41.8,111.85]));
    assert_transform_eq(&(b*a), &Transform::new(&[168.89,132.9,85.1,147.3,-271.365,-234.645]));
}

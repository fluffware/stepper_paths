use std;
use std::fmt;

#[derive(Debug)]
pub struct Vector
{
    pub x: f64,
    pub y: f64
}
impl Vector {
    pub fn length(&self) -> f64
    {
        (self.x*self.x + self.y*self.y).sqrt() 
    }

    pub fn scalar_mul(&self, other: Self) -> f64
    {
        self.x * other.x + self.y * other.y
    }

    /// Normalize vector to unit length
    pub fn unit(&self) -> Vector
    {
        let l = self.length();
        Point{x: self.x/l, y: self.y/l}
    }

    /// Angular direction
    pub fn angle(&self) -> f64
    {
        self.y.atan2(self.x)
    }

    /// Rotate vector 90 deg clockwise
    pub fn rotate_90_cw(&self) -> Vector
    {
        Vector{x: self.y, y: -self.x}
    }

    /// Rotate vector 90 deg counter-clockwise
    pub fn rotate_90_ccw(&self) -> Vector
    {
        Vector{x: -self.y, y: self.x}
    }

    /// Get a value that has the same sign as the difference between the
    /// angles
    ///
    /// If self.angle() < other.angle() then the result is negative, otherwise
    /// positive.
    /// The result is zero if the vectors are parallel.
    pub fn diff_sign(&self, other: &Vector) -> f64
    {
        self.y * other.x - self.x * other.y
    }

    /// Returns cosine of the angle between `self` and `other`.
    ///
    /// For parallel vectors the returned value is 1,
    /// for anti-parallel vectors it's -1 and 0 for perpendicular vectors.
    pub fn angle_cos(&self, other: &Vector) -> f64
    {
        (self.x * other.x + self.y * other.y) 
            / ((self.x * self.x + self.y * self.y ) 
               * (other.x * other.x + other.y * other.y )).sqrt()
    }
}

impl std::fmt::Display for Vector
{
    fn fmt(self: &Self,f: &mut fmt::Formatter) -> fmt::Result
    {
        write!(f, "({}, {})", self.x, self.y)
    }
}

impl std::cmp::PartialEq for Vector
{
    fn eq(&self, t: &Self) ->bool
    {
        self.x == t.x && self.y == t.y
    }
}
impl Clone for Vector
{
    fn clone(&self) -> Vector {
        Vector{x: self.x, y: self.y}
    }
}

impl Copy for Vector
{
}

impl std::ops::Index<usize> for Vector
{
    type Output = f64;
    fn index(&self, key: usize) -> &Self::Output
    {
        match key {
            0 => &self.x,
            1 => &self.y,
            _ => panic!("Index out of range")
        }
    }
}
impl std::ops::IndexMut<usize> for Vector
{
    fn index_mut(&mut self, key: usize) -> &mut Self::Output
    {
        match key {
            0 => &mut self.x,
            1 => &mut self.y,
            _ => panic!("Index out of range")
        }
    }
}


impl std::ops::Add<Vector> for Vector {
    type Output = Vector;
    fn add(self, v: Vector) -> Vector {
        Vector { x: self.x + v.x,
                y: self.y + v.y}
    }
}
impl std::ops::AddAssign<Vector> for Vector {
    fn add_assign(&mut self, v: Vector)  {
        self.x += v.x;
        self.y += v.y;
    }
}

impl std::ops::Sub<Vector> for Vector {
    type Output = Vector;
    fn sub(self, v: Vector) -> Vector {
        Vector { x: self.x - v.x,
                y: self.y - v.y}
    }
}
impl std::ops::SubAssign<Vector> for Vector {
    fn sub_assign(&mut self, v: Vector)  {
        self.x -= v.x;
        self.y -= v.y;
    }
}
impl std::ops::Mul<f64> for Vector {
    type Output = Vector;
    fn mul(self, s: f64) -> Vector {
        Vector { x: self.x * s,
                y: self.y * s}
    }
}

impl std::ops::Neg for Vector {
    type Output = Vector;
    fn neg(self) -> Self::Output
    {
        Vector{x: -self.x, y: -self.y}
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

    pub fn rotate_around(a: f64, pivot: &Vector) -> Transform {
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

impl std::ops::Mul<Vector> for Transform {
    type Output = Vector;
    fn mul(self, v: Vector) -> Vector {
        Vector { x: self.matrix[0] * v.x + self.matrix[2] * v.y + self.matrix[4],
                y: self.matrix[1] * v.x + self.matrix[3] * v.y + self.matrix[5]}
    }
}

pub type Point = Vector;

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

#[test]
fn test_diff_sign()
{
    
    assert_eq!(Vector{x: 1.0, y:1.0}.diff_sign(&Vector{x:-4.0, y:-4.0}),
               0.0);
    assert!(Vector{x: 1.0, y:1.0}.diff_sign(&Vector{x:4.0, y:1.0}) > 0.0);
    assert!(Vector{x: 1.0, y:1.0}.diff_sign(&Vector{x:1.0, y:4.0}) < 0.0);
    assert!(Vector{x: 1.0, y:-1.0}.diff_sign(&Vector{x:4.0, y:-1.0}) < 0.0);
    assert!(Vector{x: 1.0, y:-1.0}.diff_sign(&Vector{x:1.0, y:-4.0}) > 0.0);
}

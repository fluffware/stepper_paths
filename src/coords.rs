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
    fn fmt(&self,f: &mut fmt::Formatter) -> fmt::Result
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
    /* Arrange the same way as an SVG matrix
    [0] [2] [4]
    [1] [3] [5]
    0   0   1
     */
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

const MIN_COLUMN_LENGTH: f64 = 1e-6;

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

    pub fn no_translation(&self) -> Transform
    {
	let mut matrix = [0.0; 6];
	matrix[0..4].clone_from_slice(&self.matrix[0..4]);
	Transform{matrix}
    }

    /* Split the transform into (translate, scale, rotate, skew_x, skew_y) */
    pub fn decompose(&self) -> (Vector, Vector, f64, f64, f64)
    {
	/* Algorithms by Frédéric Wang
	(https://frederic-wang.fr/decomposition-of-2d-transform-matrices.html)
	 */
	let m = &self.matrix;
	let rot;
	let scale;
	let mut skew_x = 0.0;
	let mut skew_y = 0.0;
	let det = m[0] * m[3] - m[1] * m[2];
	
	let r2 = m[0] * m[0] + m[1] * m[1];
	let r = r2.sqrt();
	if r >= MIN_COLUMN_LENGTH {
	    scale = Vector{x: r, y: det / r};
	    rot = m[1].atan2(m[0]);
	    skew_x = (m[0]*m[2] + m[1] * m[3]).atan2(r2);
	} else {
	    let s2 = m[2] * m[2] + m[3] * m[3];
	    let s = s2.sqrt();
	    if s >= MIN_COLUMN_LENGTH {
		rot = m[3].atan2(m[4]);
		scale = Vector{x: det / s, y: s};
		skew_y = (m[0]*m[2] + m[1] * m[3]).atan2(s2);
	    } else {
		rot = 0.0;
		    scale = Vector{x: 0.0, y: 0.0};
	    }
	}

	(Vector{x: m[4], y: m[5]}, scale, rot, skew_x, skew_y)
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
use std::f64::consts::PI;

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
#[cfg(test)]
fn check_decompose(tr: Vector, scale: Vector, rot:f64, skew_x: f64, skew_y: f64)
{
    let b = Transform::translate(tr.x,tr.y)
	* Transform::rotate(rot)
	* Transform::scale_xy(scale.x, scale.y)
	* Transform::skew_x(skew_x)
	* Transform::skew_y(skew_y);
    println!("in:  translate({}) scale({}) rotate({}) skewX({}) skewY({})",
	     tr, scale, rot, skew_x, skew_y);
    let (tr, scale, rot, skew_x, skew_y) = b.decompose();

    println!("res: translate({}) scale({}) rotate({}) skewX({}) skewY({})",
	     tr, scale, rot, skew_x, skew_y);
    let c = Transform::translate(tr.x,tr.y)
	* Transform::rotate(rot)
	* Transform::scale_xy(scale.x, scale.y)
	* Transform::skew_x(skew_x)
	* Transform::skew_y(skew_y);
    println!("in: {:?}\nres: {:?}",b.matrix,c.matrix);
    assert_transform_eq(&b, &c);

}

#[test]
fn test_decompose()
{
    let mut a = Transform::identity();
    assert_eq!(a.decompose(), (Vector{x: 0.0, y:0.0}, Vector{x: 1.0, y:1.0},
			       0.0,0.0,0.0));
    check_decompose(Vector{x: 3.0, y:7.0}, Vector{x: 1.0, y:1.0},
		    1.0, 0.0, 0.0);
    check_decompose(Vector{x: 3.0, y:7.0}, Vector{x: -2.0, y:3.0},
		    4.0*PI/5.0, -1.0, 0.0);
    
    let b = Transform{matrix: [0.0, 1.0, 0.0, 2.0, -1.5, 8.0]};
    let (tr, scale, rot, skew_x, skew_y) = b.decompose();
    let c = Transform::translate(tr.x,tr.y)
	* Transform::rotate(rot)
	* Transform::scale_xy(scale.x, scale.y)
	* Transform::skew_x(skew_x)
	* Transform::skew_y(skew_y);
    assert_transform_eq(&b, &c);
}

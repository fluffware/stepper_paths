use coords::{Point, Vector};

fn incenter_point(p1: Point, p2: Point, p3: Point) -> Point {
    let l1 = (p2 - p3).length();
    let l2 = (p1 - p3).length();
    let l3 = (p1 - p2).length();
    (p1 * l1 + p2 * l2 + p3 * l3) / (l1 + l2 + l3)
}

/// Find were the lines (or their extensions) intersect. Return None if they are parallel.
fn intersection(l1: (Point, Point), l2: (Point, Point)) -> Option<Point> {
    let d1 = l1.0 - l1.1;
    let d2 = l2.0 - l2.1;
    let den = d1.cross_mul(d2);
    if den <= f64::EPSILON {
        return None;
    }
    let a1 = l1.0.cross_mul(l1.1);
    let a2 = l2.0.cross_mul(l2.1);
    Some(Point {
        x: (a1 * d2.x - a2 * d1.x) / den,
        y: (a1 * d2.y - a2 * d1.y) / den,
    })
}

#[derive(Debug, PartialEq)]
enum Circle {
    Line, // A circle with infinite radius
    Cw { center: Point, radius: f64 },
    Ccw { center: Point, radius: f64 },
}

/// Find a circle that goes through (0,0), p2, and p3
fn fit_circle(p2: Point, p3: Point) -> Circle {
    let d = p2.cross_mul(p3) * 2.0;
    if d == 0.0 {
        return Circle::Line;
    }
    let d2 = p2.scalar_mul(p2);
    let d3 = p3.scalar_mul(p3);
    let center = Point {
        x: (p3.y * d2 - p2.y * d3) / d,
        y: (p2.x * d3 - p3.x * d2) / d,
    };
    let radius = (p2 - center).length();
    if d > 0.0 {
        Circle::Ccw { center, radius }
    } else {
        Circle::Cw { center, radius }
    }
}


#[cfg(test)]
use approx::{AbsDiffEq, RelativeEq};

#[cfg(test)]
impl AbsDiffEq<Vector> for Vector {
    type Epsilon = f64;
    fn default_epsilon() -> Self::Epsilon {
        f64::EPSILON
    }

    fn abs_diff_eq(&self, other: &Vector, epsilon: Self::Epsilon) -> bool {
        self.x - other.x <= epsilon && self.y - other.y <= epsilon
    }
}

#[cfg(test)]
impl RelativeEq<Vector> for Vector {
    fn default_max_relative() -> Self::Epsilon {
        f64::EPSILON
    }

    fn relative_eq(
        &self,
        other: &Vector,
        epsilon: Self::Epsilon,
        max_relative: Self::Epsilon,
    ) -> bool {
        let max_x = (self.x.abs().max(other.x.abs()) * max_relative).max(epsilon);
        let max_y = (self.y.abs().max(other.y.abs()) * max_relative).max(epsilon);
        (self.x - other.x).abs() <= max_x && (self.y - other.y).abs() <= max_y
    }
}

#[test]
fn test_intersection() {
    assert_relative_eq!(
        intersection(
            (Point { x: -3.0, y: 3.0 }, Point { x: 5.0, y: -5.0 }),
            (Point { x: -3.0, y: -3.0 }, Point { x: 7.0, y: 7.0 })
        )
        .unwrap(),
        Point { x: 0.0, y: 0.0 },
        max_relative = 0.00001
    );

    assert_relative_eq!(
        intersection(
            (Point { x: -3.0, y: 4.0 }, Point { x: 5.0, y: -4.0 }),
            (Point { x: -3.0, y: 3.0 }, Point { x: 6.0, y: -3.0 })
        )
        .unwrap(),
        Point { x: 0.0, y: 1.0 },
        max_relative = 0.00001
    );
    assert_eq!(
        intersection(
            (Point { x: -3.0, y: 4.0 }, Point { x: 5.0, y: -4.0 }),
            (Point { x: -5.0, y: 3.0 }, Point { x: 3.0, y: -5.0 })
        ),
        None
    );
}

#[test]
fn test_fit_circle() {
    assert_eq!(
        fit_circle(Point { x: 3.0, y: 3.0 }, Point { x: 6.0, y: 0.0 }),
        Circle::Cw {
            center: Point { x: 3.0, y: 0.0 },
            radius: 3.0
        }
    );
    assert_eq!(
        fit_circle(Point { x: 3.0, y: -3.0 }, Point { x: 6.0, y: 0.0 }),
        Circle::Ccw {
            center: Point { x: 3.0, y: 0.0 },
            radius: 3.0
        }
    );
}

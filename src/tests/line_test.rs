use super::super::{coords::Point, curve_approx::CurveInfo, curves::line::Line};

#[test]
fn test_line() {
    let end: Point = Point { x: 5.3, y: 1.0 };
    let line = Line::new(end);
    assert_relative_eq!(line.length(), end.length(), max_relative = 0.0001);
    let (p, d) = line.value(0.0);
    assert_eq!(p, Point { x: 0.0, y: 0.0 });
    assert_relative_eq!(d.x, end.x / end.length(), max_relative = 0.0001);
    assert_relative_eq!(d.y, end.y / end.length(), max_relative = 0.0001);

    let (p, d) = line.value(line.length());
    assert_relative_eq!(p.x, end.x, max_relative = 0.0001);
    assert_relative_eq!(p.y, end.y, max_relative = 0.0001);
    assert_relative_eq!(d.x, end.x / end.length(), max_relative = 0.0001);
    assert_relative_eq!(d.y, end.y / end.length(), max_relative = 0.0001);
}

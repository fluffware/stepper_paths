use super::super:: {
    curve_approx::CurveInfo,
    coords::Point,
    curves::circle_segment::CircleSegment
};

#[test]
fn test_circle_segment()
{
    const radius: f64 = 3.0;
    const start: f64 = 1.3;
    const end: f64 = 2.7;
    let cs = CircleSegment::new(radius, start, end);
    assert_relative_eq!(cs.length(), radius * (end - start), 
                        max_relative=0.0001);
    let (p,d) = cs.value(0.0);
    assert_eq!(p, Point{x: 0.0, y: 0.0});
    assert_relative_eq!(d.x, -start.sin(), max_relative=0.0001);
    assert_relative_eq!(d.y, start.cos(), max_relative=0.0001);
    
    let (p,d) = cs.value(cs.length()); 
    assert_relative_eq!(p.x, radius * (end.cos() - start.cos()), 
                        max_relative=0.0001);
    assert_relative_eq!(p.y, radius * (end.sin() - start.sin()), 
                        max_relative=0.0001);
    assert_relative_eq!(d.x, -end.sin(), max_relative=0.0001);
    assert_relative_eq!(d.y, end.cos(), max_relative=0.0001);
}

#[test]
fn test_circle_segment_reverse()
{
    const radius: f64 = 3.0;
    const start: f64 = 1.3;
    const end: f64 = -0.7;
    let cs = CircleSegment::new(radius, start, end);
    assert_relative_eq!(cs.length(), radius * (start - end), 
                        max_relative=0.0001);
    let (p,d) = cs.value(0.0);
    assert_eq!(p, Point{x: 0.0, y: 0.0});
    assert_relative_eq!(d.x, start.sin(), max_relative=0.0001);
    assert_relative_eq!(d.y, -start.cos(), max_relative=0.0001);
    
    let (p,d) = cs.value(cs.length()); 
    assert_relative_eq!(p.x, radius * (end.cos() - start.cos()), 
                        max_relative=0.0001);
    assert_relative_eq!(p.y, radius * (end.sin() - start.sin()), 
                        max_relative=0.0001);
    assert_relative_eq!(d.x, end.sin(), max_relative=0.0001);
    assert_relative_eq!(d.y, -end.cos(), max_relative=0.0001);
}

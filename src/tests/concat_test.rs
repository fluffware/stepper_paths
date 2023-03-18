use super::super::{
    coords::{Point, Vector},
    curve_approx::CurveInfo,
    curves::bezier::Bezier,
    curves::circle_segment::CircleSegment,
    curves::concat_curve::ConcatCurve,
    curves::line::Line,
};
use super::svg_utils;
//use std::f64::consts::PI;
use std::io::Result;

#[test]
fn concat_equidistant() -> Result<()> {
    let mut plot = svg_utils::SvgPlot::new();
    let mut concat = ConcatCurve::new();
    let mut p_start = Point { x: 0.0, y: 0.0 };
    {
        let c1 = Point { x: 4.0, y: 0.0 };
        let c2 = Point { x: 20.0, y: 5.0 };
        let p2 = Point { x: 20.0, y: 10.0 };
        plot.add_bezier(p_start, c1 + p_start, c2 + p_start, p2 + p_start);
        let b = Bezier::new(c1, c2, p2);
        concat.add(Box::new(b));
        p_start += p2;
    }
    {
        let p2 = Point { x: 0.0, y: 5.0 };
        plot.add_line(p_start, p2 + p_start);
        let l = Line::new(p2);
        concat.add(Box::new(l));
        p_start += p2;
    }
    {
        let dir = Vector { x: 0.0, y: 5.0 };
        let end = Point { x: -7.0, y: 0.0 };
        let c = CircleSegment::new_start_direction(end, dir);
        p_start += end;
        concat.add(Box::new(c));
    }
    {
        let dir = Vector { x: 0.0, y: -5.0 };
        let end = Point { x: -7.0, y: 2.0 };
        let c = CircleSegment::new_start_direction(end, dir);
        p_start += end;
        concat.add(Box::new(c));
    }
    {
        let dir = Vector { x: 4.0, y: 5.0 };
        let end = Point { x: 3.0, y: 2.0 };
        let c = CircleSegment::new_start_direction(end, dir);
        p_start += end;
        concat.add(Box::new(c));
    }
    {
        let c1 = Point { x: 17.0, y: 8.0 };
        let c2 = Point { x: 0.0, y: 10.0 } - p_start;
        let p2 = -p_start;
        plot.add_bezier(p_start, c1 + p_start, c2 + p_start, p2 + p_start);
        let b = Bezier::new(c1, c2, p2);
        concat.add(Box::new(b));
        p_start += p2;
    }
    let len = concat.length();
    plot.set_color("blue");
    plot.set_point_radius(0.3);
    let mut pos = 0.0;
    while pos <= len {
        let (p, _dp) = concat.value(pos);
        plot.add_point(p);
        pos += 1.5;
    }
    plot.set_color("red");
    let mut pos = 0.0;
    while pos <= len {
        let (p, dp) = concat.value(pos);
        plot.add_vector(p, dp + p);
        pos += 1.5;
    }
    plot.write("concat_equidistant")?;
    Ok(())
}

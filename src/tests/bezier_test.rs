use super::super::{coords::Point, curve_approx::CurveInfo, curves::bezier::Bezier};
use super::svg_utils;
use std::io::Result;

#[cfg(test)]
fn draw_bezier_equidistant(
    plot: &mut svg_utils::SvgPlot,
    p1: Point,
    c1: Point,
    c2: Point,
    p2: Point,
) {
    let b = Bezier::new(c1 - p1, c2 - p1, p2 - p1);
    let l = b.length();
    plot.set_color("black");
    plot.add_bezier(p1, c1, c2, p2);

    plot.set_point_radius(2.0);
    for i in 0..=10 {
        let (p, _) = b.value(l * f64::from(i) / 10.0);
        plot.add_point(p + p1);
    }
    plot.set_point_radius(1.0);
    plot.set_color("red");
    for i in 0..=10 {
        let (tp, _) = b.bezier_value(f64::from(i) / 10.0);
        plot.add_point(tp + p1);
    }
    plot.set_color("blue");
    for i in 0..=10 {
        let (p, dp) = b.value(l * f64::from(i) / 10.0);
        let p0 = p + p1;
        let p1 = dp * 5.0 + p0;
        plot.add_vector(p0, p1);
    }
}

#[test]
fn test_bezier_equidistant() -> Result<()> {
    let mut plot = svg_utils::SvgPlot::new();
    // Split curve in 10 equal parts
    let p1 = Point { x: 0.0, y: 0.0 };
    let c1 = Point { x: 34.0, y: 72.0 };
    let c2 = Point { x: 119.0, y: 58.0 };
    let p2 = Point { x: 20.0, y: 10.0 };
    draw_bezier_equidistant(&mut plot, p1, c1, c2, p2);

    let p1 = Point { x: 10.0, y: 90.0 };
    let c1 = Point { x: 10.0, y: 90.0 };
    let c2 = Point { x: -5.0, y: -10.0 };
    let p2 = Point { x: 80.0, y: 80.0 };

    draw_bezier_equidistant(&mut plot, p1, c1, c2, p2);

    let p1 = Point { x: 10.0, y: 95.0 };
    let p2 = Point { x: 80.0, y: 85.0 };

    draw_bezier_equidistant(&mut plot, p1, p1, p2, p2);

    let p1 = Point { x: 48.0, y: 12.0 };
    let c1 = Point { x: 87.0, y: 34.0 };
    let c2 = Point { x: 50.0, y: 6.0 };
    let p2 = Point { x: 67.0, y: 4.0 };
    draw_bezier_equidistant(&mut plot, p1, c1, c2, p2);
    plot.write("test_bezier_equidistant")?;
    Ok(())
}

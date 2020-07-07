use super::super:: {
    curve_approx::CurveInfo,
    coords::Point,
    curves::bezier::Bezier
};
use super::svg_utils;
use std::io::Result;

#[test]
fn test_bezier_equidistant() -> Result<()>
{
    let mut plot = svg_utils::SvgPlot::new();
    // Split curve in 10 equal parts
    let p1 = Point{x:0.0, y:0.0};
    let c1 = Point{x: 34.0, y: 72.0};
    let c2 = Point{x:119.0,y:58.0};
    let p2 = Point{x:20.0, y:10.0};
    plot.add_bezier(p1,c1,c2,p2);
    let b = Bezier::new(c1,c2,p2);
    let l = b.length();
    plot.set_point_radius(2.0);
    for i in 0..=10 {
        let (p, _) = b.value(l * f64::from(i) / 10.0);
        plot.add_point(Point{x:p[0], y:p[1]});
    }
    plot.set_point_radius(1.0);
    plot.set_color("red");
    for i in 0..=10 {
        let (tp, _) = b.bezier_value(f64::from(i) / 10.0);
        plot.add_point(Point{x:tp[0], y:tp[1]});
    }
    plot.set_color("blue");
    for i in 0..=10 {
        let (p, dp) = b.value(l * f64::from(i) / 10.0);
        let p0 = Point{x:p[0], y:p[1]};
        let p1 = Point{x:dp[0], y:dp[1]} * 5.0 + p0;
        plot.add_vector(p0,p1);
    }

    plot.write("test_bezier_equidistant")?;
    Ok(())
}

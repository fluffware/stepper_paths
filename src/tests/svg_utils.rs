use coords::{Point, Transform, Vector};
use std::fs::{self, File};
use std::io::Write;
use std::io::{Error, ErrorKind, Result};
use std::path::PathBuf;

pub struct SvgPlot {
    lines: Vec<String>,
    point_raduis: f64,
    color: String,
    line_width: f64,
}

fn find_svg_dir() -> Result<PathBuf> {
    let mut dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    if !dir.is_dir() {
        return Err(Error::new(
            ErrorKind::NotFound,
            "Manifest directory not found",
        ));
    }
    dir.push("target");
    if !dir.is_dir() {
        return Err(Error::new(
            ErrorKind::NotFound,
            "Target directory not found",
        ));
    }
    dir.push("test_files_svg");
    if !dir.is_dir() {
        fs::create_dir(&dir)?;
    }
    Ok(dir)
}

impl SvgPlot {
    pub fn new() -> SvgPlot {
        SvgPlot {
            lines: Vec::new(),
            point_raduis: 1.0,
            color: "black".to_string(),
            line_width: 0.5,
        }
    }

    pub fn set_point_radius(&mut self, radius: f64) {
        self.point_raduis = radius;
    }

    pub fn set_color(&mut self, color: &str) {
        self.color = color.to_string();
    }

    pub fn set_line_width(&mut self, width: f64) {
        self.line_width = width;
    }

    pub fn add_point(&mut self, p: Point) {
        self.lines.push(format!(
            r#"<circle cx="{}" cy="{}" r="{}" style="fill:{}"/>"#,
            p.x, p.y, self.point_raduis, self.color
        ));
    }

    pub fn add_line(&mut self, p1: Point, p2: Point) {
        self.lines.push(format!(
            r#"<line x1="{}" y1="{}" x2="{}" y2="{}"  style="stroke:{};fill:none;stroke-width:{}"/>"#,
            p1.x, p1.y, p2.x, p2.y, self.color, self.line_width
        ));
    }

    pub fn add_bezier(&mut self, p1: Point, c1: Point, c2: Point, p2: Point) {
        self.lines.push(format!(
            r#"<path d="M{},{} C{},{} {},{} {},{}" style="stroke:{};fill:none;stroke-width:{}"/>"#,
            p1.x, p1.y, c1.x, c1.y, c2.x, c2.y, p2.x, p2.y, self.color, self.line_width
        ));
    }

    /// Add a circle segment from start to end with the direction start_dir at the start point
    pub fn add_circle_segment(&mut self, start: Point, end: Point, start_dir: Vector) {
        let rend = end - start;
        let rend_len = rend.length();
        // Calculate direction of radius
        let (dr, sweep) = if rend.cross_mul(start) > 0.0 {
            (
                Point {
                    x: start_dir.y,
                    y: -start_dir.x,
                },
                0,
            )
        } else {
            (
                Point {
                    x: -start_dir.y,
                    y: start_dir.x,
                },
                1,
            )
        };
        let dr_len = dr.length();
        let r = (rend_len * 0.5 * rend_len * dr_len) / dr.scalar_mul(rend);
        let (large_arc, sweep) = if dr.scalar_mul(rend) < 0.0 {
            (1, 1 - sweep)
        } else {
            (0, sweep)
        };
        self.lines.push(format!(
            r#"<path d="M{},{} A{},{} {} {} {} {},{}"
style="stroke:{};fill:none;stroke-width:{}"/>"#,
            start.x,
            start.y,
            r,
            r,
            0.0,
            large_arc,
            sweep,
            end.x,
            end.y,
            self.color,
            self.line_width
        ));
    }

    pub fn add_vector(&mut self, p1: Point, p2: Point) {
        self.lines.push("<g>".to_string());
        self.lines.push(format!(
            r#"<line x1="{}" y1="{}" x2="{}" y2="{}" style="stroke:{};stroke-width:{}"/>"#,
            p1.x, p1.y, p2.x, p2.y, self.color, self.line_width
        ));
        let ul = (p1 - p2).length() / 5.0;
        let a0 = Point { x: -ul, y: -ul };
        let a1 = Point { x: 0.0, y: 0.0 };
        let a2 = Point { x: -ul, y: ul };
        let delta = p2 - p1;
        let l = delta.length();
        let u = delta * (1.0 / l);
        let tm = Transform::new(&[u.x, u.y, -u.y, u.x, p2.x, p2.y]);
        let a0 = tm * a0;
        let a1 = tm * a1;
        let a2 = tm * a2;
        self.lines.push(format!(
            r#"<path d="M{},{} L{},{} L{},{}" style="stroke:{};fill:none;stroke-width:{}"/>"#,
            a0.x, a0.y, a1.x, a1.y, a2.x, a2.y, self.color, self.line_width
        ));
        self.lines.push("</g>".to_string());
    }

    pub fn write(&self, file_name: &str) -> Result<()> {
        let dir = find_svg_dir()?;
        let file_path = dir.join(file_name).with_extension("svg");
        let mut file = File::create(file_path)?;
        file.write(
            br#"<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<svg xmlns="http://www.w3.org/2000/svg"
     width="100px" height="100px">
"#,
        )?;
        for l in &self.lines {
            file.write(b"  ")?;
            file.write(l.as_bytes())?;
            file.write(b"\n")?;
        }
        file.write(b"</svg>\n")?;
        Ok(())
    }
}

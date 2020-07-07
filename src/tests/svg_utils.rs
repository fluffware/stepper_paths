use coords::{Point, Transform};
use std::fs::{File, self};
use std::path::PathBuf;
use std::io::{Error, ErrorKind, Result};
use std::io::Write;

pub struct SvgPlot
{
    lines: Vec<String>,
    point_raduis: f64,
    color: String
}

fn find_svg_dir() -> Result<PathBuf>
{
    let mut dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    if !dir.is_dir() {
        return Err(Error::new(ErrorKind::NotFound, "Manifest directory not found"));
    }
    dir.push("target");
     if !dir.is_dir() {
        return Err(Error::new(ErrorKind::NotFound, "Target directory not found"));
     }
    dir.push("test_files_svg");
    if !dir.is_dir() {
        fs::create_dir(&dir)?;
    }
    Ok(dir)
}

    
impl SvgPlot
{
    pub fn new() -> SvgPlot
    {
        SvgPlot{lines: Vec::new(),
                point_raduis: 1.0,
                color: "black".to_string()
        }
    }
    
    pub fn set_point_radius(&mut self, radius: f64)
    {
        self.point_raduis = radius;
    }
    
    pub fn set_color(&mut self, color: &str)
    {
        self.color = color.to_string();
    }
    
    pub fn add_point(&mut self, p: Point)
    {
        self.lines.push(format!(
            r#"<circle cx="{}" cy="{}" r="{}" style="fill:{}"/>"#,
            p.x, p.y, self.point_raduis, self.color));
    }

    pub fn add_line(&mut self, p1: Point, p2: Point)
    {
        self.lines.push(format!(r#"<line x1="{}" y1="{}" x2="{}" y2="{}"  style="stroke:{};fill:none"/>"#,
                                p1.x, p1.y, p2.x, p2.y, self.color));
    }
    
    pub fn add_bezier(&mut self, p1: Point, c1: Point, c2: Point, p2: Point)
    {
        self.lines.push(format!(r#"<path d="M{},{} C{},{} {},{} {},{}" style="stroke:{};fill:none"/>"#,
                                p1.x, p1.y, c1.x, c1.y, c2.x, c2.y, p2.x, p2.y,
                                self.color));
       
    }

    pub fn add_vector(&mut self, p1: Point, p2: Point)
    {
        self.lines.push("<g>".to_string());
        self.lines.push(format!(
            r#"<line x1="{}" y1="{}" x2="{}" y2="{}" style="stroke:{}"/>"#,
            p1.x, p1.y, p2.x, p2.y, self.color));
        let ul = (p1-p2).length() / 5.0;
        let a0 = Point{x:-ul, y:-ul};
        let a1 = Point{x: 0.0, y: 0.0};
        let a2 = Point{x: -ul, y: ul};
        let delta = p2-p1;
        let l = delta.length();
        let u = delta * (1.0 / l);
        let tm = Transform::new(&[u.x, u.y, -u.y, u.x, p2.x, p2.y]);
        let a0 = tm * a0;
        let a1 = tm * a1;
        let a2 = tm * a2;
        self.lines.push(format!(
            r#"<path d="M{},{} L{},{} L{},{}" style="stroke:{};fill:none"/>"#,
            
            a0.x, a0.y, a1.x, a1.y, a2.x, a2.y,
            self.color));
        self.lines.push("</g>".to_string());
    }
    
    pub fn write(&self, file_name: &str) ->Result<()>
    {
        let dir = find_svg_dir()?;
        let file_path = dir.join(file_name).with_extension("svg");
        let mut file = File::create(file_path)?;
        file.write(
            br#"<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<svg xmlns="http://www.w3.org/2000/svg"
     width="100px" height="100px">
"#);
        for l in &self.lines {
            file.write(b"  ");
            file.write(l.as_bytes());
            file.write(b"\n");
        }
        file.write(b"</svg>\n");
        Ok(())
    }
}

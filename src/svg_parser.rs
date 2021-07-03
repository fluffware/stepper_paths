use std::io::BufReader;
use std::io::Read;

use xml::reader::{EventReader, XmlEvent};
use xml::attribute::OwnedAttribute;
use xml::name::OwnedName;

use stepper_context::CurveSegment;
use coords::Transform;

use svg_path_parser;
use svg_path_parser::PathContext;

const SVG_NS : &str = "http://www.w3.org/2000/svg";

type DynResult<T> = 
    std::result::Result<T, Box<dyn std::error::Error + Send + Sync + 'static>>;

pub fn parse_document<T: Read>(input :T,t0: &Transform,
                                 filter: Box<dyn Fn(&OwnedName,
                                                &Vec<OwnedAttribute>) -> bool>)
                               -> DynResult<Vec<CurveSegment>>
{
    let file = BufReader::new(input);
    let mut trans_stack = Vec::<Transform>::new();
    let mut transform = *t0;
    let parser = EventReader::new(file);
    let mut path = Vec::<CurveSegment>::new();
    let mut path_ctxt = PathContext::new();
    let mut ignore_nest = 0; // Ignore elements when > 0
    for e in parser {
        match e {
            Ok(XmlEvent::StartElement { ref name, ref attributes, .. })
                if ignore_nest > 0 || !filter(name, attributes)=>
            {
                ignore_nest += 1;
            },
            
            Ok(XmlEvent::StartElement { name, attributes, .. }) => {
                if name.namespace == Some(SVG_NS.to_string()) {
                    trans_stack.push(transform);
                    let mut d: Option<String> = None;
                    let mut view_box_rect: Option<[f64;4]> = None;
                    let mut width:Option<String> = None;
                    let mut height:Option<String> = None;
                    // Parse attributes
                    for attr in &attributes {
                        if attr.name.local_name == "transform" {
                            //println!("Transform: '{}'", attr.value);
                            match svg_path_parser::parse_transform(&attr.value) {
                                Ok(t) => {
                                    transform = transform * t;
                                },
                                err => return Err(format!("Failed to parse transform {}: {:?}", attr.value, err).into())
                            }
                        } else if attr.name.local_name == "d" {
                            d = Some(attr.value.clone());
                        } else if attr.name.local_name == "viewBox" {
                            match svg_path_parser::parse_view_box(&attr.value) {
                                Ok(args) => {
                                    view_box_rect = Some(args);
                                },
                                err => return Err(format!("Failed to parse viewBox args {}: {:?}", attr.value, err).into())
                            }
                        } else if attr.name.local_name == "height" {
                            height = Some(attr.value.clone());
                        } else if attr.name.local_name == "width" {
                            width = Some(attr.value.clone());
                        }
                    }

                    if name.local_name == "path"  {
                        match d {
                            None => 
                                return Err("path element has no d attribute".into()),
                            Some(d) => {
                                svg_path_parser::parse_path(
                                    d.as_str(), &mut path_ctxt,
                                    &transform, &mut path)?;
                            }
                        }
                    } else if name.local_name == "svg"  {
                       
                        let (width, height) =
                            match (width, height) {
                                (Some(w), Some(h)) => {
                                    (match svg_path_parser::parse_length(&w) {
                                        Ok(w) => w,
                                        Err(e) => return Err(format!("Failed to parse length {}: {:?}", w, e).into())
                                    },
                                     match svg_path_parser::parse_length(&h) {
                                         Ok(h) => h,
                                         Err(e) => return Err(format!("Failed to parse length {}: {:?}", w, e).into())
                                     })},
                                _ => return Err("svg element needs width and height attributes".into())
                            };
                        let (vbx, vby, vbw, vbh) =
                            match view_box_rect {
                                Some(vb) =>
                                    (vb[0], vb[1], vb[2], vb[3]),
                                None => (0.0,0.0, width, height)
                            };
                        let sc = Transform::scale_xy(width/vbw, height/vbh);
                        let tr = Transform::translate(vbx,vby);
                        transform = transform * sc * tr;
                    }
                }
            }
            Ok(XmlEvent::EndElement { .. }) if ignore_nest > 0 => {
                ignore_nest -= 1;
            },
            
            Ok(XmlEvent::EndElement { name }) => {
                if name.namespace == Some(SVG_NS.to_string()) {
                    transform = trans_stack.pop().unwrap();
                }
            }
            Err(e) => {
                println!("Error: {}", e);
                break;
            }
            _ => {}
        }
    }
    
    Ok(path)
}

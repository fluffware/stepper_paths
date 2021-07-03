extern crate paths;
extern crate serial;
extern crate xml;

use paths::stepper_context::StepperContext;
use paths::stepper_control::play_events;
use paths::stepper_config_parser;

use serial::core::SerialDevice;
use serial::core::SerialPortSettings;
use std::env;

use std::fs::File;

use paths::coords::Transform;
use paths::svg_parser;

use xml::attribute::OwnedAttribute;
use xml::name::OwnedName;

use std::str::FromStr;

extern crate getopts;
use getopts::Options;


fn usage(prg: &str, opts: Options)
{
    let brief = format!("Usage: {} [options] FILE", prg);
    print!("{}", opts.usage(&brief));
}
const INKSCAPE_NS: &str = "http://www.inkscape.org/namespaces/inkscape";
const SVG_NS : &str = "http://www.w3.org/2000/svg";

 

fn main() {
    
    let args: Vec<String> = env::args().collect();
    let program = args[0].clone();
    let mut opts = Options::new();
    opts.optopt("v", "draw-speed", "preferred velocity when drawing", "SPEED");
    opts.optopt("","a-max", "maximum allowed acceleration", "ACC");
    opts.optopt("","v-max", "maximum allowed velocity", "VEL");
    opts.optopt("","intensity", "laser power (%)", "INTENSITY");
    opts.optopt("","repeat", "repeat paths N times", "PATH");
    opts.optopt("","layer", "only draw this Inkscape layer", "LAYER");
    opts.optopt("d", "device", "serial device", "DEV");
    opts.optopt("c", "config", "stepper configuration file", "FILE");
    opts.optflag("h", "help", "print this help menu");
    
    let matches = match opts.parse(&args[1..]) {
        Ok(m) => { m }
        Err(f) => { panic!("{}", f.to_string()) }
    };
    if matches.opt_present("h") {
        usage(&program, opts);
        return;
    }
   

    if matches.free.is_empty() {
        println!("No file name");
        return;
    }

    let config;
    /*
        ticks_per_second: TICKS_PER_SECOND as u32,
        stepper_x: 
        stepper_config::StepperConfig {
            steps_per_millimeter: STEPS_PER_MM * S_SCALE as f64,
            max_velocity: 50.0,
            max_acceleration: 50.0
        },
        stepper_y:
        stepper_config::StepperConfig {
            steps_per_millimeter: STEPS_PER_MM * S_SCALE as f64,
            max_velocity: 50.0,
            max_acceleration: 50.0
        },
        laser:
        stepper_config::LaserConfig {
            max_intensity: 256,
            pwm_period: 256,
            alignment_intensity: 5
        }
    };*/
    if let Some(filename) = matches.opt_str("config")
    {
        config = match stepper_config_parser::read_config(&filename) {
            Err(e) =>
            {
                println!("{}", e);
            return
            },
            Ok(c) => c
        };
    } else {
        println!("No configuration file");
        return;
    }
    
    let port_name = matches.opt_str("device");
    
    let file_name = matches.free[0].clone();

    let vx_scale = config.stepper_x.steps_per_millimeter
        / (2 * config.ticks_per_second) as f64;
    let mut vx_max = (vx_scale.abs() * config.stepper_x.max_velocity).round() as i32;

    let vy_scale = config.stepper_y.steps_per_millimeter
        / (2 * config.ticks_per_second) as f64;
    let mut vy_max = (vy_scale.abs() * config.stepper_y.max_velocity).round() as i32;
    
    if let Some(arg) = matches.opt_str("v-max") {
        match f64::from_str(&arg) {
            Ok(value) => {
                vx_max = (vx_scale.abs() * value) as i32;
                vy_max = (vy_scale.abs() * value) as i32;
            },
            Err(err) => {
                println!("Invalid max velocity: {}", err);
                return
            }
        }
    };
   
    let ax_scale = vx_scale / config.ticks_per_second as f64;
    let ay_scale = vy_scale / config.ticks_per_second as f64;
    
    let mut ax_max =
        (ax_scale * config.stepper_x.max_acceleration).abs().round() as i32;
    let mut ay_max =
        (ay_scale * config.stepper_y.max_acceleration).abs().round() as i32;
    
    if let Some(arg) = matches.opt_str("a-max") {
        match f64::from_str(&arg) {
            Ok(value) => {
                ax_max = (value * ax_scale).abs() as i32;
                ay_max = (value * ay_scale).abs() as i32;
            },
            Err(err) => {
                println!("Invalid max acceleration: {}", err);
                return
            }
        }
    };

    
    let v_draw = match matches.opt_str("draw-speed") {
        Some(arg) => match str::parse(&arg) {
            Ok(value) => value,
            Err(err) => {
                println!("Invalid draw velocity: {}", err);
                return
            }
        },
        None => 20.0
    };

    let n_repeat = match matches.opt_str("repeat") {
        Some(arg) => match str::parse(&arg) {            
            Ok(value) if value >= 1 => value,
            Ok(_)  => {
                println!("Repeat must be >= 1");
                return
            },
            Err(err) => {
                println!("Invalid repeat: {}", err);
                return
            }
        },
        None => 1
    };
    
    let weight = match matches.opt_str("intensity") {
        Some(arg) => match str::parse(&arg) {
            Ok(value) if (0..=100).contains(&value) => value,
            Ok(_)  => {
                println!("Invalid intensity, must be 0 - 100c");
                return
            }
            Err(err) => {
                println!("Invalid intensity: {}", err);
                return
            }
        },
        None => 5
    } * 254 / 100;


    let serport = 
        match port_name {
            Some(port_name) =>
            {
                let mut serport = 
                    match serial::open(&port_name) {
                        Ok(port) => port,
                        Err(e) => panic!("Failed to open {}: {}", port_name, e)
                    };
                match serport.read_settings() {
                    Ok(mut settings) => {
                        settings.set_baud_rate(serial::Baud9600).unwrap();
                        if let Err(e) = serport.write_settings(&settings) {
                            panic!("Failed write serial settings: {}", e)
                    }
                    },
                    Err(e) => panic!("Failed read serial settings: {}", e)
                }
                Some(serport)
            },
            None => None
        };
 
   
    let mut ctxt = StepperContext::new(
        &[ax_max, ay_max],
        &[vx_max, vy_max],
        &[config.stepper_x.steps_per_millimeter,
          config.stepper_y.steps_per_millimeter]);
    ctxt.set_weight(weight);
    
    let file = match File::open(&file_name) {
        Ok(file) => file,
        Err(e) => panic!("Failed to open {}: {}", file_name, e)
    };
    let filter: Box<dyn Fn(&OwnedName,
                       &std::vec::Vec<OwnedAttribute>) -> bool>;
    filter = match matches.opt_str("layer") {
        Some(layer_name) => {
            Box::new(move |name,attrs| {
                match name {
                    OwnedName{local_name ,namespace: Some(namespace),prefix: _} 
                    if local_name == "g" && namespace == SVG_NS => {
                        println!("Group");
                        let mut is_layer = false;
                        let mut layer_name_matches = false;
                        for ref a in attrs {
                            if let OwnedName{
                                ref local_name,
                                namespace: Some(ref namespace),
                                prefix: _} = a.name 
                            {
                                if namespace == INKSCAPE_NS {
                                    if local_name == "groupmode" &&
                                        a.value == "layer" {
                                                is_layer = true;
                                        }
                                    if local_name == "label" &&
                                        a.value == layer_name {
                                            layer_name_matches = true;
                                        }
                                }
                            }
                        }
                        if is_layer && !layer_name_matches {
                            return false;
                        }
                    },
                    _ => {}
                };
                true})
        },
        None => {
            Box::new(|_,_| {
                true})
        }
    };

    let trans = Transform::identity();
    let path = match svg_parser::parse_document(file, &trans, filter) {
        Err(msg) => panic!("Parser failed: {}", msg),
        Ok(path) => {println!("Segments: {:?}",path);path}
    };

    let mut res = Ok([0,0]);
    for _ in 0..n_repeat {
        res = ctxt.draw_curves(&path,v_draw / (config.ticks_per_second as f64));
        match &res {
            Ok(_) => {},
            Err(e) => {
                eprintln!("Failed to draw curve: {}", e);
                return;
            }
        }
    }
     println!("Max errors: {:?}", res.unwrap());
                
    println!("Pos: {:?}",ctxt.position());
    ctxt.step_goto(0,0);
    println!("End time: {} s", ctxt.ticks() / (config.ticks_per_second as u64));
    if let Some(mut serport) = serport {
        let events = ctxt.events();
        
        play_events(&mut serport, events).unwrap();
    }
}


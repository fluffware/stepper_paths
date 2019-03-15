extern crate paths;
extern crate serial;
extern crate xml;

use paths::stepper_context::StepperContext;
use paths::stepper_control::play_events;

use serial::core::SerialDevice;
use serial::core::SerialPortSettings;
use std::env;

use std::fs::File;

use paths::coords::Transform;
use paths::svg_parser;





extern crate getopts;
use getopts::Options;


// const BEZIER_CIRCLE_C: f64 = 0.551915024494;
//const WEIGHT : i32 = 160;
fn usage(prg: &str, opts: Options)
{
    let brief = format!("Usage: {} [options] FILE", prg);
    print!("{}", opts.usage(&brief));
}

const V_SCALE: i32 = 128;
const S_SCALE: i64 = 2*128*128;

const STEPS_PER_MM:f64 = 60.0;

fn main() {
    
    let args: Vec<String> = env::args().collect();
    let program = args[0].clone();
    let mut opts = Options::new();
    opts.optopt("v", "draw-speed", "preferred velocity when drawing", "SPEED");
    opts.optopt("","a-max", "maximum allowed acceleration", "ACC");
    opts.optopt("","v-max", "maximum allowed velocity", "VEL");
    opts.optopt("","intensity", "laser power (%)", "INTENSITY");
    opts.optopt("","repeat", "repeat paths N times", "PATH");
    opts.optopt("d", "device", "serial device", "DEV");
    opts.optflag("h", "help", "print this help menu");
    
    let matches = match opts.parse(&args[1..]) {
        Ok(m) => { m }
        Err(f) => { panic!(f.to_string()) }
    };
    if matches.opt_present("h") {
        usage(&program, opts);
        return;
    }
   

    if matches.free.len() < 1 {
        println!("No file name");
        return;
    }

    let port_name = matches.opt_str("device");
    
    let file_name = matches.free[0].clone();

    let v_max = V_SCALE * match matches.opt_str("v-max") {
        Some(arg) => match i32::from_str_radix(&arg, 10) {
            Ok(value) => value,
            Err(err) => {
                println!("Invalid max velocity: {}", err);
                return
            }
        },
        None => 5000
    };
    
    let a_max = match matches.opt_str("a-max") {
        Some(arg) => match i32::from_str_radix(&arg, 10) {
            Ok(value) => value,
            Err(err) => {
                println!("Invalid max acceleration: {}", err);
                return
            }
        },
        None => 8000
    };
        
    let v_draw = V_SCALE * match matches.opt_str("draw-speed") {
        Some(arg) => match i32::from_str_radix(&arg, 10) {
            Ok(value) => value,
            Err(err) => {
                println!("Invalid draw velocity: {}", err);
                return
            }
        },
        None => 1000
    };

    let n_repeat = match matches.opt_str("repeat") {
        Some(arg) => match u32::from_str_radix(&arg, 10) {            
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
        Some(arg) => match i32::from_str_radix(&arg, 10) {
            Ok(value) if value <= 100 && value >= 0 => value,
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
 
   
    let mut ctxt = StepperContext::new(a_max, v_max);
    ctxt.set_weight(weight);
    
    let file = match File::open(&file_name) {
        Ok(file) => file,
        Err(e) => panic!("Failed to open {}: {}", file_name, e)
    };
    let path = match svg_parser::parse_document(file, &Transform::scale_xy(-STEPS_PER_MM * S_SCALE as f64, -STEPS_PER_MM * S_SCALE as f64)) {
        Err(msg) => panic!("Parser failed: {}", msg),
        Ok(path) => {println!("Segments: {:?}",path);path}
    };
  
    for _ in 0..n_repeat {
        ctxt.draw_curves(&path,v_draw);
    }
                
    println!("Pos: {:?}",ctxt.position());
    ctxt.goto(0,0);
    if let Some(mut serport) = serport {
        let events = ctxt.events();
        
        play_events(&mut serport, events).unwrap();
    }
}


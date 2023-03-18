extern crate image;
extern crate paths;
extern crate serial;

use paths::stepper_config_parser;
use paths::stepper_context::StepperContext;
use paths::stepper_control::play_events;

use serial::core::SerialDevice;
use serial::core::SerialPortSettings;
use std::env;
use std::fs::File;
use std::io::BufReader;

use std::convert::TryFrom;
use std::str::FromStr;

extern crate getopts;
use getopts::Options;

use image::ImageBuffer;
use image::Pixel;

fn usage(prg: &str, opts: Options) {
    let brief = format!("Usage: {} [options] FILE", prg);
    print!("{}", opts.usage(&brief));
}

const DEFAULT_SPOT_DIAMETER: f64 = 0.2;
const DEFAULT_WIDTH: f64 = 20.0;

fn find_row_interval<P, Container, F>(
    image: &ImageBuffer<P, Container>,
    row: u32,
    filter: F,
) -> (u32, u32)
where
    F: Fn(&P) -> bool,
    P: image::Pixel + 'static,
    Container: std::ops::Deref<Target = [P::Subpixel]>,
{
    let (width, _) = image.dimensions();
    let mut x = 0;
    while x < width {
        if filter(image.get_pixel(x, row)) {
            break;
        }
        x += 1;
    }
    if x == width {
        return (0, 0);
    }
    let start = x;
    x = width - 1;
    while x > start {
        if filter(image.get_pixel(x, row)) {
            break;
        }
        x -= 1;
    }
    x += 1;
    (start, x - start)
}

fn main() {
    let args: Vec<String> = env::args().collect();
    let program = args[0].clone();
    let mut opts = Options::new();
    opts.optopt(
        "v",
        "draw-speed",
        "preferred velocity when drawing",
        "SPEED",
    );
    opts.optopt("", "a-max", "maximum allowed acceleration", "ACC");
    opts.optopt("", "v-max", "maximum allowed velocity", "VEL");
    opts.optopt("", "width", "width of image (mm)", "WIDTH");
    opts.optopt("", "height", "height of image (mm)", "HEIGHT");
    opts.optopt("", "spot-diam", "diameter of laser spot (mm)", "DIAMETER");
    opts.optopt("", "intensity", "laser power (%)", "INTENSITY");
    //opts.optopt("","repeat", "repeat paths N times", "PATH");
    opts.optopt("d", "device", "serial device", "DEV");
    opts.optopt("c", "config", "stepper configuration file", "FILE");
    opts.optflag("h", "help", "print this help menu");

    let matches = match opts.parse(&args[1..]) {
        Ok(m) => m,
        Err(f) => {
            panic!("{}", f.to_string())
        }
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
    if let Some(filename) = matches.opt_str("config") {
        config = match stepper_config_parser::read_config(&filename) {
            Err(e) => {
                println!("{}", e);
                return;
            }
            Ok(c) => c,
        };
    } else {
        println!("No configuration file");
        return;
    }

    let port_name = matches.opt_str("device");

    let file_name = matches.free[0].clone();

    let vx_scale = config.stepper_x.steps_per_millimeter / (2 * config.ticks_per_second) as f64;
    let mut vx_max = (vx_scale.abs() * config.stepper_x.max_velocity).round() as i32;

    let vy_scale = config.stepper_y.steps_per_millimeter / (2 * config.ticks_per_second) as f64;
    let mut vy_max = (vy_scale.abs() * config.stepper_y.max_velocity).round() as i32;

    if let Some(arg) = matches.opt_str("v-max") {
        match str::parse::<f64>(&arg) {
            Ok(value) => {
                vx_max = (vx_scale.abs() * value) as i32;
                vy_max = (vy_scale.abs() * value) as i32;
            }
            Err(err) => {
                println!("Invalid max velocity: {}", err);
                return;
            }
        }
    };

    let ax_scale = vx_scale / config.ticks_per_second as f64;
    let ay_scale = vx_scale / config.ticks_per_second as f64;

    let mut ax_max = (ax_scale.abs() * config.stepper_x.max_acceleration).round() as i32;
    let mut ay_max = (ay_scale.abs() * config.stepper_y.max_acceleration).round() as i32;

    if let Some(arg) = matches.opt_str("a-max") {
        match f64::from_str(&arg) {
            Ok(value) => {
                ax_max = (value * ax_scale).abs() as i32;
                ay_max = (value * ay_scale).abs() as i32;
            }
            Err(err) => {
                println!("Invalid max acceleration: {}", err);
                return;
            }
        }
    };

    let v_draw = match matches.opt_str("draw-speed") {
        Some(arg) => match f64::from_str(&arg) {
            Ok(value) => value,
            Err(err) => {
                println!("Invalid draw velocity: {}", err);
                return;
            }
        },
        None => 20.0,
    };

    let weight = match matches.opt_str("intensity") {
        Some(arg) => match str::parse(&arg) {
            Ok(value) if (0..=100).contains(&value) => value,
            Ok(_) => {
                println!("Invalid intensity, must be 0 - 100c");
                return;
            }
            Err(err) => {
                println!("Invalid intensity: {}", err);
                return;
            }
        },
        None => 5,
    } * 254
        / 100;

    let width_opt = match matches.opt_str("width") {
        Some(arg) => match f64::from_str(&arg) {
            Ok(value) => Some(value),
            Err(err) => {
                println!("Invalid width value: {}", err);
                return;
            }
        },
        None => None,
    };

    let height_opt = match matches.opt_str("height") {
        Some(arg) => match f64::from_str(&arg) {
            Ok(value) => Some(value),
            Err(err) => {
                println!("Invalid height value: {}", err);
                return;
            }
        },
        None => None,
    };

    let spot_diam = match matches.opt_str("spot-diam") {
        Some(arg) => match f64::from_str(&arg) {
            Ok(value) => value,
            Err(err) => {
                println!("Invalid height value: {}", err);
                return;
            }
        },
        None => DEFAULT_SPOT_DIAMETER,
    };

    let file = match File::open(&file_name) {
        Ok(file) => file,
        Err(e) => {
            println!("Failed to open {}: {}", file_name, e);
            return;
        }
    };

    let buffer = match image::load(BufReader::new(file), image::ImageFormat::Png) {
        Ok(dyn_buffer) => dyn_buffer.to_luma8(),
        Err(e) => {
            println!("Failed to read image: {}", e);
            return;
        }
    };
    let (x_pixels, y_pixels) = buffer.dimensions();

    println!("Dimensions: {}, {}", x_pixels, y_pixels);

    // Keep aspect ratio if possible
    let (width_mm, height_mm) = match (width_opt, height_opt) {
        (Some(w), Some(h)) => (w, h),
        (Some(w), None) => (w, (w * y_pixels as f64 / x_pixels as f64)),
        (None, Some(h)) => ((h * x_pixels as f64 / y_pixels as f64), h),
        (None, None) => (
            DEFAULT_WIDTH,
            (DEFAULT_WIDTH * y_pixels as f64 / x_pixels as f64),
        ),
    };
    let serport = match port_name {
        Some(port_name) => {
            let mut serport = match serial::open(&port_name) {
                Ok(port) => port,
                Err(e) => panic!("Failed to open {}: {}", port_name, e),
            };
            match serport.read_settings() {
                Ok(mut settings) => {
                    settings.set_baud_rate(serial::Baud9600).unwrap();
                    if let Err(e) = serport.write_settings(&settings) {
                        panic!("Failed write serial settings: {}", e)
                    }
                }
                Err(e) => panic!("Failed read serial settings: {}", e),
            }
            Some(serport)
        }
        None => None,
    };

    let mut ctxt = StepperContext::new(
        &[ax_max, ay_max],
        &[vx_max, vy_max],
        &[
            config.stepper_x.steps_per_millimeter,
            config.stepper_y.steps_per_millimeter,
        ],
    );
    ctxt.set_weight(weight);

    let width = width_mm * config.stepper_x.steps_per_millimeter;
    let res_x = (width * config.ticks_per_second as f64
        / (x_pixels as f64 * (v_draw * config.stepper_x.steps_per_millimeter)))
        .round() as u32;
    let res_x = u32::max(1, res_x);
    let interval = x_pixels * res_x;
    let v_x = (width / (2 * interval) as f64).round() as i32;

    let res_y = (height_mm / (y_pixels as f64 * spot_diam)).round() as u32;
    let step_y = (height_mm * config.stepper_y.steps_per_millimeter / (res_y * y_pixels) as f64)
        .round() as i64;

    // Calculate real width and height
    let width_mm = (v_x as u32 * 2 * interval) as f64 / config.stepper_x.steps_per_millimeter;
    let height_mm =
        (res_y as i64 * step_y * y_pixels as i64) as f64 / config.stepper_y.steps_per_millimeter;
    println!("Width: {} mm, height: {} mm", width_mm, height_mm);

    println!("res_x: {} res_y: {}", res_x, res_y);
    let mut ltr = true;
    for yn in 0..y_pixels {
        let (start_x, length_x) =
            find_row_interval(&buffer, yn, |p| -> bool { p.channels()[0] < 255 });
        if length_x > 0 {
            for repeat_y in 0..res_y {
                let mut prev_weight = 0;

                let line_start;
                let line_v;

                let bitmap_x: Box<dyn Fn(u32) -> u32>;
                if ltr {
                    line_start = start_x;
                    line_v = v_x;
                    bitmap_x = Box::new(|xn| start_x + xn);
                } else {
                    line_start = start_x + length_x;
                    line_v = -v_x;
                    bitmap_x = Box::new(|xn| start_x + length_x - 1 - xn);
                }
                ctxt.step_goto_speed(
                    (res_x * 2 * v_x as u32 * line_start) as i64,
                    (yn * res_y + repeat_y) as i64 * step_y,
                    line_v,
                    0,
                );
                let start = ctxt.ticks();
                for xn in 0..length_x {
                    let w = (255 - buffer.get_pixel(bitmap_x(xn), yn).channels()[0]) as i32
                        * weight
                        / 256;
                    if w != prev_weight {
                        let when = start + u64::from(xn * res_x);
                        let ivl = u32::try_from(when - ctxt.ticks()).unwrap();
                        if ivl != 0 {
                            ctxt.add_acc_interval(0, 0, ivl);
                        }
                        ctxt.add_weight_change(w, when);
                        prev_weight = w;
                    }
                }
                let when = start + u64::from(length_x * res_x);
                let ivl = u32::try_from(when - ctxt.ticks()).unwrap();
                ctxt.add_acc_interval(0, 0, ivl);
                ctxt.add_weight_change(0, when);
                ltr = !ltr;
            }
        }
    }

    ctxt.step_goto(0, 0);
    println!(
        "End time: {} s",
        ctxt.ticks() / (config.ticks_per_second as u64)
    );
    if let Some(mut serport) = serport {
        let events = ctxt.events();

        play_events(&mut serport, events).unwrap();
    }
}

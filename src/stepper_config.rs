use std::io::BufReader;

pub struct StepperConfig
{
    pub steps_per_millimeter: f64,
    pub max_velocity: f64, // mm/s
    pub max_acceleration: f64 // mm/s^2
}


pub struct LaserConfig
{
    pub max_intensity: u32,
    pub pwm_period: u32, 
    pub alignment_intensity: u32 // Just indicate the spot without burning
}
use std::fs::File;

pub struct XYStepperConfig
{
    pub ticks_per_second: u32,
    pub stepper_x: StepperConfig,
    pub stepper_y: StepperConfig,
    pub laser: LaserConfig
}

pub struct Error
{
   descr: String 
}

impl Error
{
    pub fn new(descr: &str) -> Error
    {
        Error {descr: String::from(descr)}
    }
}
            
impl std::error::Error for Error
{
    fn description(&self) -> &str
    {
        return &self.descr;
    }
    
    fn cause(&self) -> Option<&dyn std::error::Error> {None}
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {None}
}

impl std::fmt::Display for Error
{
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{}", self.descr)
    }
}

impl std::fmt::Debug for Error
{
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{}", self.descr)
    }
}

fn set_f64(json: &serde_json::Value, v: &mut f64) -> Result<(), Error> 
{
    match json
    {
        serde_json::Value::Null => Ok(()),
        serde_json::Value::Number(num) => {
            *v = num.as_f64().unwrap();
            Ok(())
        },
        _ => Err(Error::new("Illegal numeric argument"))
    }
}

fn set_u32(json: &serde_json::Value, v: &mut u32) -> Result<(), Error> 
{
    match json
    {
        serde_json::Value::Null => Ok(()),
        serde_json::Value::Number(num) => {
            match num.as_u64()
            {
                Some(i) => {
                    if i > std::u32::MAX as u64 {
                        return Err(Error::new("Integer value to big"))
                    }
                    *v = i as u32
                },
                None => return Err(Error::new("Not an integer value "))
            }
            Ok(())
        },
        _ => Err(Error::new("Illegal numeric argument"))
    }
}
fn update_stepper(stepper: &mut StepperConfig, map: &serde_json::Map<String, serde_json::Value>) -> Result<(), Error>
{
    set_f64(&map["stepsPerMillimeter"],&mut stepper.steps_per_millimeter)?;
    set_f64(&map["maxVelocity"],&mut stepper.max_velocity)?;
    set_f64(&map["maxAcceleration"],&mut stepper.max_acceleration)?;
    Ok(())
}

fn update_laser(laser: &mut LaserConfig, map: &serde_json::Map<String, serde_json::Value>) -> Result<(), Error>
{
    set_u32(&map["maxIntensity"],&mut laser.max_intensity)?;
    set_u32(&map["pwmPeriod"],&mut laser.pwm_period)?;
    set_u32(&map["alignmentIntensity"],&mut laser.alignment_intensity)?;
    Ok(())
}  
pub fn read_config(config: &mut XYStepperConfig, file_name: &str) -> Result<(), Error>
{
    let file = match File::open(&file_name) {
        Ok(file) => file,
        Err(e) =>  return Err(Error::new(&format!("Failed to open file: {}", e)))
    };
    let reader = BufReader::new(file);

    let json: serde_json::Value = match serde_json::from_reader(reader)
    {
        Ok(json) => json,
        Err(e) =>  return Err(Error::new(&format!("Failed to parse config: {}", e)))
    };
    set_u32(&json["ticksPerSecond"], &mut config.ticks_per_second)?;
    match &json["stepperX"]
    {
        serde_json::Value::Object(map) =>
        {
            update_stepper(&mut config.stepper_x, map)?;
        },
        _ => ()
    }
    match &json["stepperY"]
    {
        serde_json::Value::Object(map) =>
        {
            update_stepper(&mut config.stepper_y, map)?;
        },
        _ => ()
    }
    match &json["laser"]
    {
        serde_json::Value::Object(map) =>
        {
            update_laser(&mut config.laser, map)?;
        },
        _ => ()
    }
    Ok(())
}

use std::io::BufReader;
use stepper_config::XyStepperConfig;
use std::fs::File;

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
        &self.descr
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
pub fn read_config(file_name: &str) -> Result<XyStepperConfig, Error>
{
    let file = match File::open(&file_name) {
        Ok(file) => file,
        Err(e) =>  return Err(Error::new(&format!("Failed to open file: {}", e)))
    };
    let reader = BufReader::new(file);

    let config: XyStepperConfig = match serde_json::from_reader(reader)
    {
        Ok(json) => json,
        Err(e) =>  return Err(Error::new(&format!("Failed to parse config: {}", e)))
    };
 
    Ok(config)
}

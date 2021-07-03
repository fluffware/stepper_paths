
use serde::{Serialize, Deserialize};

#[derive(Serialize, Deserialize, Debug)]
#[serde(rename_all = "camelCase")]
pub struct StepperConfig
{
    pub steps_per_millimeter: f64,
    pub max_velocity: f64, // mm/s
    pub max_acceleration: f64 // mm/s^2
}

#[derive(Serialize, Deserialize, Debug)]
#[serde(rename_all = "camelCase")]
pub struct LaserConfig
{
    pub max_intensity: u32,
    pub pwm_period: u32, 
    pub alignment_intensity: u32 // Just indicate the spot without burning
}

#[derive(Serialize, Deserialize, Debug)]
#[serde(rename_all = "camelCase")]
pub struct XyStepperConfig
{
    pub ticks_per_second: u32,
    pub stepper_x: StepperConfig,
    pub stepper_y: StepperConfig,
    pub laser: LaserConfig
}

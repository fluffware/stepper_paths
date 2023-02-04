extern crate num;
#[macro_use]
extern crate nom;
extern crate serde;
extern crate xml;

#[cfg(test)]
#[macro_use]
extern crate approx;

pub mod acc_vector;
pub mod coords;
pub mod curve_approx;
pub mod gauss_quadrature;
pub mod integer_curve;
pub mod multi_range;
pub mod solve;
pub mod stepper_config;
pub mod stepper_config_parser;
pub mod stepper_context;
pub mod stepper_control;
pub mod svg_parser;
pub mod svg_path_parser;

pub mod curves {
    pub mod bezier;
    pub mod circle_segment;
    pub mod concat_curve;
    pub mod line;
    pub mod linear_curve;
}

#[cfg(test)]
mod tests;

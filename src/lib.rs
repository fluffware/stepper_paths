extern crate nom;
extern crate num;
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

pub mod curves;

#[cfg(test)]
mod tests;

extern crate num;
#[macro_use]
extern crate nom;
extern crate xml;

#[cfg(test)]
#[macro_use]
extern crate approx;

pub mod solve;
pub mod integer_curve;
pub mod multi_range;
pub mod acc_vector;
pub mod stepper_context;
pub mod curve_approx;
pub mod stepper_config;
pub mod gauss_quadrature;
pub mod stepper_control;
pub mod coords;
pub mod svg_parser;
pub mod curves {
    pub mod bezier;
    pub mod linear_curve;
    pub mod concat_curve;
    pub mod circle_segment;
    pub mod line;
}

#[cfg(test)]
mod tests;

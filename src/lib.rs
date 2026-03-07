//! # estrannaise-rs
//!
//! Pharmacokinetic (PK) modelling library for estradiol ester formulations,
//! ported 1:1 from the JavaScript source of the *estrannaise* project.
//!
//! ## Feature flags
//! - `wasm` - enables `wasm-bindgen` derives and a thin WASM API surface
//!   (see the `wasm` module).  Requires the `wasm-bindgen` crate.
//!
//! ## Quick start
//! ```rust
//! use estrannaise_rs::models::{Model, pk_function, e2_multidose_3c, RandomMode};
//!
//! // Single dose: 5 mg EV im at t = 0, evaluate at t = 3 days
//! let e2 = pk_function(Model::EvIm, 3.0, 5.0, 1.0, false, 0.0);
//! println!("E2 at day 3: {e2:.1} pg/mL");
//!
//! // Steady state: 5 mg every 7 days
//! let e2_ss = pk_function(Model::EvIm, 3.0, 5.0, 1.0, true, 7.0);
//! println!("E2 (SS) at day 3: {e2_ss:.1} pg/mL");
//!
//! // Multi-dose schedule
//! let e2_multi = e2_multidose_3c(
//!     10.0,
//!     &[5.0, 5.0],
//!     &[0.0, 7.0],
//!     &[Model::EvIm, Model::EvIm],
//!     1.0,
//!     RandomMode::None,
//!     false,
//! );
//! println!("E2 (multi) at day 10: {e2_multi:.1} pg/mL");
//! ```

pub mod menstrual;
pub mod modeldata;
pub mod models;
pub mod spline;
pub mod units;

#[cfg(feature = "wasm")]
pub mod wasm;

// Convenience re-exports
pub use menstrual::{MenstrualPoint, fill_menstrual_cycle_curve};
pub use models::{
    Model, PKQuantities, RandomMode, e2_curve_3c, e2_multidose_3c, e2_patch_3c, e2_ss_average_3c,
    e2_steady_state_3c, e2_steady_state_patch_3c, es_single_dose_3c, fill_curve, fill_target_range,
    get_pk_quantities, get_pk_quantities_3c, golden_section_search, pk_function,
    pk_random_function, random_mcmc_sample, terminal_elimination_time_3c,
};
pub use units::Unit;

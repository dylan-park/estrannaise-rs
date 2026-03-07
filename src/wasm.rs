//! Thin `wasm-bindgen` API surface.
//!
//! Enabled only when compiling with `--features wasm`.
//! Exposes the most commonly-needed functions to JavaScript
//! with the same semantics as the original JS library.

use crate::models::{self, Model, RandomMode};
use crate::units::Unit;
use wasm_bindgen::prelude::*;

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

fn model_from_js(s: &str) -> Result<Model, JsValue> {
    Model::from_str(s).ok_or_else(|| JsValue::from_str(&format!("Unknown model: {s}")))
}

fn unit_cf(unit_str: &str) -> f64 {
    Unit::from_str(unit_str)
        .map(|u| u.conversion_factor())
        .unwrap_or(1.0)
}

// ---------------------------------------------------------------------------
// Core PK functions
// ---------------------------------------------------------------------------

/// Evaluate the median PK curve for a single injection ester.
#[wasm_bindgen]
pub fn pk_function(
    model: &str,
    t: f64,
    dose: f64,
    unit: &str,
    steadystate: bool,
    t_period: f64,
) -> Result<f64, JsValue> {
    let m = model_from_js(model)?;
    let cf = unit_cf(unit);
    Ok(models::pk_function(m, t, dose, cf, steadystate, t_period))
}

/// Same as `pk_function` but draws from the MCMC posterior.
/// Pass `idx = -1` for a random draw; otherwise pins to that index.
#[wasm_bindgen]
pub fn pk_random_function(
    model: &str,
    t: f64,
    dose: f64,
    unit: &str,
    steadystate: bool,
    t_period: f64,
    idx: i32,
) -> Result<f64, JsValue> {
    let m = model_from_js(model)?;
    let cf = unit_cf(unit);
    let idx = if idx < 0 { None } else { Some(idx as usize) };
    Ok(models::pk_random_function(
        m,
        t,
        dose,
        cf,
        steadystate,
        t_period,
        idx,
    ))
}

// ---------------------------------------------------------------------------
// Multi-dose
// ---------------------------------------------------------------------------

/// Compute the combined E2 at time `t` for a multi-dose schedule.
///
/// `random_mode`:
///   - `"none"`    → median PK
///   - `"random"`  → random MCMC sample per dose
///   - `"<N>"`     → MCMC sample at index N
#[wasm_bindgen]
pub fn e2_multidose(
    t: f64,
    doses: &[f64],
    times: &[f64],
    model_strs: Vec<JsValue>,
    unit: &str,
    random_mode: &str,
    intervals: bool,
) -> Result<f64, JsValue> {
    let models: Result<Vec<Model>, _> = model_strs
        .iter()
        .map(|v| model_from_js(v.as_string().as_deref().unwrap_or("")))
        .collect();
    let models = models?;
    let cf = unit_cf(unit);
    let rng = match random_mode {
        "none" => RandomMode::None,
        "random" => RandomMode::Random,
        s => s
            .parse::<usize>()
            .map(RandomMode::Index)
            .map_err(|_| JsValue::from_str("invalid random_mode"))?,
    };
    Ok(models::e2_multidose_3c(
        t, doses, times, &models, cf, rng, intervals,
    ))
}

// ---------------------------------------------------------------------------
// Steady-state average
// ---------------------------------------------------------------------------

#[wasm_bindgen]
pub fn e2_ss_average(model: &str, dose: f64, t_period: f64, unit: &str) -> Result<f64, JsValue> {
    let m = model_from_js(model)?;
    let [d, k1, k2, k3] = crate::modeldata::pk_parameters(m);
    let cf = unit_cf(unit);
    Ok(models::e2_ss_average_3c(cf * dose, t_period, d, k1, k2, k3))
}

// ---------------------------------------------------------------------------
// fill_curve helper
// ---------------------------------------------------------------------------

/// Evaluate the median PK curve at `nb_steps` evenly-spaced points.
/// Returns a flat `[t0, e2_0, t1, e2_1, …]` array for easy JS consumption.
#[wasm_bindgen]
pub fn fill_curve(
    model: &str,
    dose: f64,
    unit: &str,
    steadystate: bool,
    t_period: f64,
    x_min: f64,
    x_max: f64,
    nb_steps: usize,
) -> Result<Vec<f64>, JsValue> {
    let m = model_from_js(model)?;
    let cf = unit_cf(unit);
    let pts = models::fill_curve(
        |t| models::pk_function(m, t, dose, cf, steadystate, t_period),
        x_min,
        x_max,
        nb_steps,
    );
    Ok(pts.iter().flat_map(|&(t, e2)| [t, e2]).collect())
}

// ---------------------------------------------------------------------------
// Menstrual cycle curve
// ---------------------------------------------------------------------------

/// Returns `[time, e2, e2p5, e2p95, …]` interleaved.
#[wasm_bindgen]
pub fn fill_menstrual_cycle_curve(x_min: f64, x_max: f64, nb_steps: usize, unit: &str) -> Vec<f64> {
    let cf = unit_cf(unit);
    let pts = crate::menstrual::fill_menstrual_cycle_curve(x_min, x_max, nb_steps, cf);
    pts.iter()
        .flat_map(|p| [p.time, p.e2, p.e2p5, p.e2p95])
        .collect()
}

// ---------------------------------------------------------------------------
// PK quantities
// ---------------------------------------------------------------------------

/// Exposes PK summary quantities to JavaScript with named fields.
/// `wasm-bindgen` generates JS getters for each `pub` field automatically.
#[wasm_bindgen]
pub struct JsPKQuantities {
    pub t_max: f64,
    pub c_max: f64,
    pub half_life: f64,
    pub half_life_absorption: f64,
}

/// Returns a `JsPKQuantities` object with named fields accessible from JS.
#[wasm_bindgen]
pub fn get_pk_quantities(model: &str) -> Result<JsPKQuantities, JsValue> {
    let m = model_from_js(model)?;
    let q = models::get_pk_quantities(m);
    Ok(JsPKQuantities {
        t_max: q.t_max,
        c_max: q.c_max,
        half_life: q.half_life,
        half_life_absorption: q.half_life_absorption,
    })
}

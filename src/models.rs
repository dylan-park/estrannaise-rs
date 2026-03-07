//! Core pharmacokinetic model functions.
//! Direct port of models.js - logic and variable names are intentionally kept
//! as close to the original JavaScript as possible.

use crate::modeldata;

// ---------------------------------------------------------------------------
// Model enum
// ---------------------------------------------------------------------------

/// All supported estradiol ester / delivery models.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[cfg_attr(feature = "wasm", wasm_bindgen::prelude::wasm_bindgen)]
pub enum Model {
    EvIm,
    EEnIm,
    EcIm,
    EbIm,
    EUnIm,
    EUnCasubq,
    PatchTw,
    PatchOw,
}

impl Model {
    /// Parse from the same string keys used in the JavaScript source.
    pub fn from_str(s: &str) -> Option<Self> {
        match s {
            "EV im" => Some(Self::EvIm),
            "EEn im" => Some(Self::EEnIm),
            "EC im" => Some(Self::EcIm),
            "EB im" => Some(Self::EbIm),
            "EUn im" => Some(Self::EUnIm),
            "EUn casubq" => Some(Self::EUnCasubq),
            "patch tw" => Some(Self::PatchTw),
            "patch ow" => Some(Self::PatchOw),
            _ => None,
        }
    }

    pub fn as_str(self) -> &'static str {
        match self {
            Self::EvIm => "EV im",
            Self::EEnIm => "EEn im",
            Self::EcIm => "EC im",
            Self::EbIm => "EB im",
            Self::EUnIm => "EUn im",
            Self::EUnCasubq => "EUn casubq",
            Self::PatchTw => "patch tw",
            Self::PatchOw => "patch ow",
        }
    }

    /// Width (in days) used internally for transdermal patch models.
    pub fn patch_width(self) -> Option<f64> {
        match self {
            Self::PatchTw => Some(3.5),
            Self::PatchOw => Some(7.0),
            _ => None,
        }
    }

    pub fn is_patch(self) -> bool {
        self.patch_width().is_some()
    }

    pub fn all() -> &'static [Model] {
        &[
            Model::EvIm,
            Model::EEnIm,
            Model::EcIm,
            Model::EbIm,
            Model::EUnIm,
            Model::EUnCasubq,
            Model::PatchTw,
            Model::PatchOw,
        ]
    }
}

// ---------------------------------------------------------------------------
// MCMC random sampling
// ---------------------------------------------------------------------------

/// Return a random MCMC sample `[d, k1, k2, k3]` for `model`.
/// Pass `Some(idx)` to pin to a specific sample (same semantics as
/// the JS `idx` parameter).
pub fn random_mcmc_sample(model: Model, idx: Option<usize>) -> [f64; 4] {
    let samples = modeldata::mcmc_samples(model);
    let i = match idx {
        Some(i) => i % samples.len(),
        None => (pseudo_random_usize()) % samples.len(),
    };
    samples[i]
}

/// Dead-simple LCG - only used when the caller has not provided an external
/// RNG.  For real uncertainty-cloud work you should supply indices from your
/// own RNG (the `idx: Some(i)` path).
fn pseudo_random_usize() -> usize {
    use std::cell::Cell;
    thread_local! {
        static STATE: Cell<u64> = Cell::new(6364136223846793005);
    }
    STATE.with(|s| {
        let v = s
            .get()
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        s.set(v);
        (v >> 33) as usize
    })
}

// ---------------------------------------------------------------------------
// logsubexp helper  (JS: _logsubexp)
// ---------------------------------------------------------------------------

/// Numerically stable `x + ln(1 - exp(y - x))`.
/// Panics when `y > x` (mirrors the JS throw).
#[inline]
fn logsubexp(x: f64, y: f64) -> f64 {
    assert!(y <= x, "logsubexp: y must be <= x");
    if x == y {
        f64::NEG_INFINITY
    } else {
        x + (1.0 - (y - x).exp()).ln()
    }
}

// ---------------------------------------------------------------------------
// Single-dose 3-compartment - injections  (JS: e2Curve3C)
// ---------------------------------------------------------------------------

/// Single-dose E2 curve for injection esters (3-compartment model).
///
/// # Parameters
/// - `t`            – time (days)
/// - `dose`         – dose (mg), already scaled by the caller
/// - `d, k1, k2, k3` – PK parameters
/// - `ds`           – initial condition for the second compartment (default 0)
/// - `d2`           – initial condition for the third compartment  (default 0)
/// - `steadystate`  – if true, delegates to `e2_steady_state_3c`
/// - `t_period`     – dosing interval `T` used in steady-state mode
pub fn e2_curve_3c(
    t: f64,
    dose: f64,
    d: f64,
    k1: f64,
    k2: f64,
    k3: f64,
    ds: f64,
    d2: f64,
    steadystate: bool,
    t_period: f64,
) -> f64 {
    if steadystate {
        return e2_steady_state_3c(t, dose, t_period, d, k1, k2, k3);
    }

    if t < 0.0 {
        return 0.0;
    }

    let mut ret = 0.0;

    if d2 > 0.0 {
        ret += d2 * (-k3 * t).exp();
    }

    if ds > 0.0 {
        if k2 == k3 {
            ret += ds * k2 * t * (-k2 * t).exp();
        } else {
            ret += ds * k2 / (k2 - k3) * ((-k3 * t).exp() - (-k2 * t).exp());
        }
    }

    if dose > 0.0 && d > 0.0 {
        if k1 == k2 && k2 == k3 {
            ret += dose * d * k1 * k1 * t * t * (-k1 * t).exp() / 2.0;
        } else if k1 == k2 && k2 != k3 {
            ret += dose * d * k1 * k1 * ((-k3 * t).exp() - (-k1 * t).exp() * (1.0 + (k1 - k3) * t))
                / (k1 - k3).powi(2);
        } else if k1 != k2 && k1 == k3 {
            ret += dose * d * k1 * k2 * ((-k2 * t).exp() - (-k1 * t).exp() * (1.0 + (k1 - k2) * t))
                / (k1 - k2).powi(2);
        } else if k1 != k2 && k2 == k3 {
            ret += dose * d * k1 * k2 * ((-k1 * t).exp() - (-k2 * t).exp() * (1.0 - (k1 - k2) * t))
                / (k1 - k2).powi(2);
        } else {
            ret += dose
                * d
                * k1
                * k2
                * ((-k1 * t).exp() / (k1 - k2) / (k1 - k3)
                    - (-k2 * t).exp() / (k1 - k2) / (k2 - k3)
                    + (-k3 * t).exp() / (k1 - k3) / (k2 - k3));
        }
    }

    if ret.is_nan() { 0.0 } else { ret }
}

// ---------------------------------------------------------------------------
// Steady-state - injections  (JS: e2SteadyState3C)
// ---------------------------------------------------------------------------

pub fn e2_steady_state_3c(
    t: f64,
    dose: f64,
    t_period: f64,
    d: f64,
    k1: f64,
    k2: f64,
    k3: f64,
) -> f64 {
    let phase = t - t_period * (t / t_period).floor();
    dose * d
        * k1
        * k2
        * ((-k1 * phase).exp() / (1.0 - (-k1 * t_period).exp()) / (k1 - k2) / (k1 - k3)
            - (-k2 * phase).exp() / (1.0 - (-k2 * t_period).exp()) / (k1 - k2) / (k2 - k3)
            + (-k3 * phase).exp() / (1.0 - (-k3 * t_period).exp()) / (k1 - k3) / (k2 - k3))
}

// ---------------------------------------------------------------------------
// Single-dose intermediate (Es) compartment  (JS: esSingleDose3C)
// ---------------------------------------------------------------------------

pub fn es_single_dose_3c(t: f64, dose: f64, d: f64, k1: f64, k2: f64, ds: f64) -> f64 {
    if t < 0.0 {
        return 0.0;
    }
    let mut ret = 0.0;
    if ds > 0.0 {
        ret += ds * (-k2 * t).exp();
    }
    if dose > 0.0 && d > 0.0 {
        if k1 == k2 {
            ret += dose * d * k1 * t * (-k1 * t).exp();
        } else {
            ret += dose * d * k1 / (k1 - k2) * ((-k2 * t).exp() - (-k1 * t).exp());
        }
    }
    ret
}

// ---------------------------------------------------------------------------
// Transdermal patch - single dose  (JS: e2Patch3C)
// ---------------------------------------------------------------------------

pub fn e2_patch_3c(
    t: f64,
    dose: f64,
    d: f64,
    k1: f64,
    k2: f64,
    k3: f64,
    w: f64, // patch wear duration (days): 3.5 or 7.0
    steadystate: bool,
    t_period: f64,
) -> f64 {
    if steadystate {
        return e2_steady_state_patch_3c(t, dose, t_period, d, k1, k2, k3, w);
    }
    if t < 0.0 {
        return 0.0;
    } else if t <= w {
        e2_curve_3c(t, dose, d, k1, k2, k3, 0.0, 0.0, false, 0.0)
    } else {
        let es_w = es_single_dose_3c(w, dose, d, k1, k2, 0.0);
        let e2_w = e2_curve_3c(w, dose, d, k1, k2, k3, 0.0, 0.0, false, 0.0);
        e2_curve_3c(t - w, 0.0, 0.0, k1, k2, k3, es_w, e2_w, false, 0.0)
    }
}

// ---------------------------------------------------------------------------
// Transdermal patch - steady state  (JS: e2SteadyStatePatch3C)
// ---------------------------------------------------------------------------

pub fn e2_steady_state_patch_3c(
    t: f64,
    dose: f64,
    t_period: f64,
    d: f64,
    k1: f64,
    k2: f64,
    k3: f64,
    w: f64,
) -> f64 {
    let es_w = es_single_dose_3c(w, dose, d, k1, k2, 0.0);
    let e2_w = e2_curve_3c(w, dose, d, k1, k2, k3, 0.0, 0.0, false, 0.0);

    let phase_t = t - t_period * (t / t_period).floor();
    let phase_t_w = t - t_period * ((t - w) / t_period).floor();

    let mut ret = dose
        * d
        * k1
        * k2
        * (logsubexp(-k1 * phase_t, -k1 * phase_t_w).exp()
            / (1.0 - (-k1 * t_period).exp())
            / (k1 - k2)
            / (k1 - k3)
            - logsubexp(-k2 * phase_t, -k2 * phase_t_w).exp()
                / (1.0 - (-k2 * t_period).exp())
                / (k1 - k2)
                / (k2 - k3)
            + logsubexp(-k3 * phase_t, -k3 * phase_t_w).exp()
                / (1.0 - (-k3 * t_period).exp())
                / (k1 - k3)
                / (k2 - k3));

    let tw_phase = t - w - t_period * ((t - w) / t_period).floor();
    ret += es_w * k2 / (k2 - k3)
        * ((-k3 * tw_phase).exp() / (1.0 - (-k3 * t_period).exp())
            - (-k2 * tw_phase).exp() / (1.0 - (-k2 * t_period).exp()));

    ret += e2_w * (-k3 * tw_phase).exp() / (1.0 - (-k3 * t_period).exp());

    ret
}

// ---------------------------------------------------------------------------
// Steady-state average (injection only)  (JS: e2ssAverage3C)
// ---------------------------------------------------------------------------

/// `_k1` and `_k2` are accepted but unused (kept for API symmetry with
/// splatting `PKParameters` the same way the JS does).
pub fn e2_ss_average_3c(dose: f64, t_period: f64, d: f64, _k1: f64, _k2: f64, k3: f64) -> f64 {
    dose * d / k3 / t_period
}

// ---------------------------------------------------------------------------
// Terminal elimination time estimate  (JS: terminalEliminationTime3C)
// ---------------------------------------------------------------------------

pub fn terminal_elimination_time_3c(k1: f64, k2: f64, k3: f64, nb_half_lives: f64) -> f64 {
    nb_half_lives * 2_f64.ln() * (1.0 / k1 + 1.0 / k2 + 1.0 / k3)
}

// ---------------------------------------------------------------------------
// High-level PK function dispatch  (JS: PKFunctions / PKRandomFunctions)
// ---------------------------------------------------------------------------

/// Evaluate the median PK curve for `model` at time `t`.
///
/// - `dose`        – dose in mg (or mcg/day for patches), **not yet scaled**
/// - `cf`          – unit conversion factor (1.0 → pg/mL)
/// - `steadystate` – use steady-state formula
/// - `t_period`    – dosing interval T (days), used only in steady-state mode
pub fn pk_function(
    model: Model,
    t: f64,
    dose: f64,
    cf: f64,
    steadystate: bool,
    t_period: f64,
) -> f64 {
    let [d, k1, k2, k3] = modeldata::pk_parameters(model);
    let scaled_dose = cf * dose;
    if let Some(w) = model.patch_width() {
        e2_patch_3c(t, scaled_dose, d, k1, k2, k3, w, steadystate, t_period)
    } else {
        e2_curve_3c(
            t,
            scaled_dose,
            d,
            k1,
            k2,
            k3,
            0.0,
            0.0,
            steadystate,
            t_period,
        )
    }
}

/// Same as `pk_function` but draws PK parameters from the MCMC posterior.
///
/// `idx` pins the sample index; `None` draws randomly (see `random_mcmc_sample`).
pub fn pk_random_function(
    model: Model,
    t: f64,
    dose: f64,
    cf: f64,
    steadystate: bool,
    t_period: f64,
    idx: Option<usize>,
) -> f64 {
    let [d, k1, k2, k3] = random_mcmc_sample(model, idx);
    let scaled_dose = cf * dose;
    if let Some(w) = model.patch_width() {
        e2_patch_3c(t, scaled_dose, d, k1, k2, k3, w, steadystate, t_period)
    } else {
        e2_curve_3c(
            t,
            scaled_dose,
            d,
            k1,
            k2,
            k3,
            0.0,
            0.0,
            steadystate,
            t_period,
        )
    }
}

// ---------------------------------------------------------------------------
// Multi-dose summation  (JS: e2multidose3C)
// ---------------------------------------------------------------------------

/// Compute the combined E2 level at time `t` for an arbitrary schedule of
/// doses with (possibly) different esters.
///
/// # Parameters
/// - `doses`            – dose amounts (mg or mcg/day)
/// - `times`            – absolute or interval times (days), depending on `intervals`
/// - `models`           – ester model for each dose entry
/// - `cf`               – unit conversion factor
/// - `random`           – sampling mode (see [`RandomMode`])
/// - `intervals`        – when `true`, `times` are inter-dose intervals;
///                        the first entry is the initial offset and is skipped
///                        (mirrors the JS cumsum logic)
pub fn e2_multidose_3c(
    t: f64,
    doses: &[f64],
    times: &[f64],
    models: &[Model],
    cf: f64,
    random: RandomMode,
    intervals: bool,
) -> f64 {
    assert_eq!(doses.len(), times.len());
    assert_eq!(doses.len(), models.len());

    // Convert intervals → absolute times  (JS cumsum logic)
    let abs_times: Vec<f64> = if intervals {
        let mut acc = 0.0;
        times
            .iter()
            .enumerate()
            .map(|(i, &v)| {
                if i == 0 {
                    0.0
                } else {
                    acc += v;
                    acc
                }
            })
            .collect()
    } else {
        times.to_vec()
    };

    let mut sum = 0.0;
    for i in 0..doses.len() {
        sum += match random {
            RandomMode::None => pk_function(models[i], t - abs_times[i], doses[i], cf, false, 0.0),
            RandomMode::Random => {
                pk_random_function(models[i], t - abs_times[i], doses[i], cf, false, 0.0, None)
            }
            RandomMode::Index(idx) => pk_random_function(
                models[i],
                t - abs_times[i],
                doses[i],
                cf,
                false,
                0.0,
                Some(idx),
            ),
        };
    }
    sum
}

/// Controls how `e2_multidose_3c` samples PK parameters.
///
/// Mirrors the JS `random` parameter which accepted `false`, `true`, or a
/// non-negative integer index.
#[derive(Debug, Clone, Copy)]
pub enum RandomMode {
    /// Use median PK parameters.
    None,
    /// Draw a random MCMC sample for each dose.
    Random,
    /// Use a specific MCMC sample index.
    Index(usize),
}

// ---------------------------------------------------------------------------
// Curve helpers  (JS: fillCurve, fillTargetRange)
// ---------------------------------------------------------------------------

/// Evaluate `f` at `nb_steps` evenly-spaced points in `[x_min, x_max]`.
pub fn fill_curve<F>(f: F, x_min: f64, x_max: f64, nb_steps: usize) -> Vec<(f64, f64)>
where
    F: Fn(f64) -> f64,
{
    (0..nb_steps)
        .map(|i| {
            let t = x_min + i as f64 * (x_max - x_min) / (nb_steps - 1) as f64;
            (t, f(t))
        })
        .collect()
}

/// Target E2 range (100–200 pg/mL) scaled by `cf`.
pub fn fill_target_range(x_min: f64, x_max: f64, cf: f64) -> [(f64, f64, f64); 2] {
    [
        (x_min, cf * 100.0, cf * 200.0),
        (x_max, cf * 100.0, cf * 200.0),
    ]
}

// ---------------------------------------------------------------------------
// Golden-section search  (JS: goldenSectionSearch)
// ---------------------------------------------------------------------------

pub fn golden_section_search<F>(
    f: F,
    mut a: f64,
    mut b: f64,
    tolerance: f64,
    max_iterations: usize,
) -> f64
where
    F: Fn(f64) -> f64,
{
    let phi = (1.0 + 5_f64.sqrt()) / 2.0;
    let mut c = b - (b - a) / phi;
    let mut d = a + (b - a) / phi;

    for _ in 0..max_iterations {
        if f(c) < f(d) {
            b = d;
        } else {
            a = c;
        }
        c = b - (b - a) / phi;
        d = a + (b - a) / phi;
        if (b - a).abs() < tolerance {
            break;
        }
    }
    (b + a) / 2.0
}

// ---------------------------------------------------------------------------
// PK summary quantities  (JS: getPKQuantities3C / getPKQuantities)
// ---------------------------------------------------------------------------

#[derive(Debug, Clone)]
pub struct PKQuantities {
    pub t_max: f64,
    pub c_max: f64,
    pub half_life: f64,
    pub half_life_absorption: f64,
}

pub fn get_pk_quantities_3c(d: f64, k1: f64, k2: f64, k3: f64) -> PKQuantities {
    let terminal = terminal_elimination_time_3c(k1, k2, k3, 5.0);
    let t_max = golden_section_search(
        |t| -e2_curve_3c(t, 1.0, 1.0, k1, k2, k3, 0.0, 0.0, false, 0.0),
        0.0,
        terminal,
        1e-5,
        100,
    );
    let c_max = e2_curve_3c(t_max, 1.0, 1.0, k1, k2, k3, 0.0, 0.0, false, 0.0);
    let c_half = c_max / 2.0;
    let t_half = golden_section_search(
        |t| (e2_curve_3c(t, 1.0, 1.0, k1, k2, k3, 0.0, 0.0, false, 0.0) - c_half).powi(2),
        t_max,
        terminal,
        1e-5,
        100,
    );
    let t_half_absorption = golden_section_search(
        |t| (e2_curve_3c(t, 1.0, 1.0, k1, k2, k3, 0.0, 0.0, false, 0.0) - c_half).powi(2),
        0.0,
        t_max,
        1e-5,
        100,
    );
    PKQuantities {
        t_max,
        c_max: d * c_max,
        half_life: t_half - t_max,
        half_life_absorption: t_half_absorption,
    }
}

pub fn get_pk_quantities(model: Model) -> PKQuantities {
    let [d, k1, k2, k3] = modeldata::pk_parameters(model);
    get_pk_quantities_3c(d, k1, k2, k3)
}

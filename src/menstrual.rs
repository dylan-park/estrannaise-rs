//! Menstrual cycle E2 reference curve.
//!
//! Port of `fillMenstrualCycleCurve` from models.js.

use crate::modeldata::{
    MENSTRUAL_CYCLE_E2, MENSTRUAL_CYCLE_E2P5, MENSTRUAL_CYCLE_E2P95, MENSTRUAL_CYCLE_T,
};
use crate::spline::CubicSpline;

use std::sync::OnceLock;

static SPLINE_E2: OnceLock<CubicSpline> = OnceLock::new();
static SPLINE_P05: OnceLock<CubicSpline> = OnceLock::new();
static SPLINE_P95: OnceLock<CubicSpline> = OnceLock::new();

fn spline_e2() -> &'static CubicSpline {
    SPLINE_E2.get_or_init(|| CubicSpline::new(MENSTRUAL_CYCLE_T, MENSTRUAL_CYCLE_E2))
}
fn spline_p05() -> &'static CubicSpline {
    SPLINE_P05.get_or_init(|| CubicSpline::new(MENSTRUAL_CYCLE_T, MENSTRUAL_CYCLE_E2P5))
}
fn spline_p95() -> &'static CubicSpline {
    SPLINE_P95.get_or_init(|| CubicSpline::new(MENSTRUAL_CYCLE_T, MENSTRUAL_CYCLE_E2P95))
}

/// One point in the menstrual cycle curve.
#[derive(Debug, Clone, Copy)]
pub struct MenstrualPoint {
    pub time: f64,
    pub e2: f64,
    pub e2p5: f64,
    pub e2p95: f64,
}

/// Generate the menstrual cycle reference curve (median + 5th/95th percentile).
///
/// The 28-day cycle repeats via `t mod 28`.
///
/// # Parameters
/// - `x_min`, `x_max` – time range (days)
/// - `nb_steps`       – number of evaluation points
/// - `cf`             – unit conversion factor (1.0 → pg/mL)
pub fn fill_menstrual_cycle_curve(
    x_min: f64,
    x_max: f64,
    nb_steps: usize,
    cf: f64,
) -> Vec<MenstrualPoint> {
    (0..nb_steps)
        .map(|i| {
            let t = x_min + i as f64 * (x_max - x_min) / (nb_steps - 1) as f64;
            // ((t % 28) + 28) % 28  - same as JS for negative t
            let phase = ((t % 28.0) + 28.0) % 28.0;
            MenstrualPoint {
                time: t,
                e2: cf * spline_e2().at(phase),
                e2p5: cf * spline_p05().at(phase),
                e2p95: cf * spline_p95().at(phase),
            }
        })
        .collect()
}

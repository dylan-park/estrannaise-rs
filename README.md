# estrannaise-rs

A Rust library that is a 1:1 port of the pharmacokinetic (PK) model from the [estrannaise](https://git.gay/Estra/estrannaise) JavaScript project. It models serum estradiol levels for various estradiol ester formulations used in transfeminine HRT.

Compiles to native Rust **and** WebAssembly (via `wasm-bindgen`).

---

## Crate layout

```
src/
  lib.rs          - crate root, re-exports
  models.rs       - core PK math (direct port of models.js)
  modeldata.rs    - auto-generated PK parameters + MCMC samples (from modeldata.js)
  spline.rs       - cubic hermite spline (direct port of cubic-spline.js)
  menstrual.rs    - menstrual cycle reference curve (fillMenstrualCycleCurve)
  units.rs        - unit conversion (pg/mL, pmol/L, ng/L, firkin/furlong³)
  wasm.rs         - thin wasm-bindgen API (feature = "wasm")
tests/
  parity.rs       - cross-validation tests against JS reference values
```

---

## Supported models

| Rust `Model` variant | JS key           | Description                                     |
|----------------------|------------------|-------------------------------------------------|
| `Model::EvIm`        | `"EV im"`        | Estradiol Valerate, intramuscular               |
| `Model::EEnIm`       | `"EEn im"`       | Estradiol Enanthate, intramuscular              |
| `Model::EcIm`        | `"EC im"`        | Estradiol Cypionate, intramuscular              |
| `Model::EbIm`        | `"EB im"`        | Estradiol Benzoate, intramuscular               |
| `Model::EUnIm`       | `"EUn im"`       | Estradiol Undecylate, intramuscular             |
| `Model::EUnCasubq`   | `"EUn casubq"`   | Estradiol Undecylate in castor oil, subcutaneous|
| `Model::PatchTw`     | `"patch tw"`     | Transdermal patch, twice-weekly                 |
| `Model::PatchOw`     | `"patch ow"`     | Transdermal patch, once-weekly                  |

---

## Usage (native)

```toml
[dependencies]
estrannaise = { path = "." }
```

```rust
use estrannaise::{Model, pk_function, e2_multidose_3c, RandomMode};

// Single dose: 5 mg EV im, evaluate at t = 3 days
let e2 = pk_function(Model::EvIm, 3.0, 5.0, 1.0, false, 0.0);

// Steady state: 5 mg EV every 7 days, evaluate at t = 3 days
let e2_ss = pk_function(Model::EvIm, 3.0, 5.0, 1.0, true, 7.0);

// Multi-dose schedule (absolute times)
let e2_multi = e2_multidose_3c(
    10.0,
    &[5.0, 5.0],
    &[0.0, 7.0],
    &[Model::EvIm, Model::EvIm],
    1.0,         // pg/mL
    RandomMode::None,
    false,       // absolute times, not intervals
);

// With uncertainty (random MCMC sample)
let e2_rand = pk_function_random(Model::EvIm, 3.0, 5.0, 1.0, false, 0.0, None);

// Fill a curve (900 points)
let curve = fill_curve(
    |t| pk_function(Model::EvIm, t, 5.0, 1.0, false, 0.0),
    0.0, 28.0, 900,
);
```

---

## Unit conversion

```rust
use estrannaise::Unit;
let cf = Unit::PmolPerL.conversion_factor(); // 3.6713
let e2_pmol = pk_function(Model::EvIm, 3.0, 5.0, cf, false, 0.0);
```

---

## Building for WebAssembly

```bash
# Install wasm-pack once
cargo install wasm-pack

# Build (outputs to pkg/)
wasm-pack build --release --target web --features wasm

# Or for Node.js
wasm-pack build --release --target nodejs --features wasm
```

The WASM API mirrors the Rust API. All functions that accept a `model` string use the same keys as the original JS (`"EV im"`, `"patch tw"`, etc.).

```js
import init, { pk_function, e2_multidose, fill_curve } from './pkg/estrannaise.js';
await init();

const e2 = pk_function("EV im", 3.0, 5.0, "pg/mL", false, 0.0);
const curve = fill_curve("EV im", 5.0, "pg/mL", true, 7.0, 0.0, 28.0, 900);
// curve is Float64Array: [t0, e2_0, t1, e2_1, ...]
```

---

## Running tests

```bash
cargo test
```

---

## Mathematical model

All models use a 3-compartment open pharmacokinetic system:

- **Depot → Absorption → Plasma** (injection esters)
- **Patch-on → Absorption → Plasma** (transdermal patches, with a wear-time cutoff W)

Steady-state solutions are derived analytically. Uncertainty is represented via pre-computed MCMC posterior samples (313 samples per model), drawn at evaluation time.

The model parameters and MCMC samples are the original work of the estrannaise project author and differ from those used in other simulators.

---

## Correspondence table (JS → Rust)

| JavaScript function             | Rust equivalent                                |
|---------------------------------|------------------------------------------------|
| `e2Curve3C`                     | `e2_curve_3c`                                  |
| `e2SteadyState3C`               | `e2_steady_state_3c`                           |
| `esSingleDose3C`                | `es_single_dose_3c`                            |
| `e2Patch3C`                     | `e2_patch_3c`                                  |
| `e2SteadyStatePatch3C`          | `e2_steady_state_patch_3c`                     |
| `e2ssAverage3C`                 | `e2_ss_average_3c`                             |
| `e2multidose3C`                 | `e2_multidose_3c`                              |
| `PKFunctions(cf)[model](…)`     | `pk_function(model, t, dose, cf, …)`           |
| `PKRandomFunctions(cf)[m](…)`   | `pk_random_function(model, t, dose, cf, …)`    |
| `randomMCMCSample`              | `random_mcmc_sample`                           |
| `fillCurve`                     | `fill_curve`                                   |
| `fillMenstrualCycleCurve`       | `fill_menstrual_cycle_curve`                   |
| `fillTargetRange`               | `fill_target_range`                            |
| `terminalEliminationTime3C`     | `terminal_elimination_time_3c`                 |
| `goldenSectionSearch`           | `golden_section_search`                        |
| `getPKQuantities3C`             | `get_pk_quantities_3c`                         |
| `getPKQuantities`               | `get_pk_quantities`                            |

## Future Work

- [ ] Modify parity tests to directly use the values from the js script instead of hard coding

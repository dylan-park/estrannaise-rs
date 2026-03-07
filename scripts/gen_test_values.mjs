// gen_test_values.mjs
// Run with: node gen_test_values.mjs
// Place this file alongside models.js and modeldata.js

import {
  e2Curve3C,
  e2SteadyState3C,
  esSingleDose3C,
  e2Patch3C,
  e2SteadyStatePatch3C,
  e2ssAverage3C,
  e2multidose3C,
  PKFunctions,
  PKParameters,
  fillMenstrualCycleCurve,
} from "./models.js";

function section(title) {
  console.log(`\n${"=".repeat(60)}`);
  console.log(`// ${title}`);
  console.log("=".repeat(60));
}

function row(label, value) {
  console.log(`${label}: ${value}`);
}

// -----------------------------------------------------------------------
// e2Curve3C - all injection models, single dose 5mg
// -----------------------------------------------------------------------
section("e2Curve3C - single dose 5mg, all injection models");
const injectionModels = [
  "EV im",
  "EEn im",
  "EC im",
  "EB im",
  "EUn im",
  "EUn casubq",
];
for (const model of injectionModels) {
  console.log(`\n// ${model}`);
  const [d, k1, k2, k3] = PKParameters[model];
  for (const t of [0, 1, 3, 7, 14, 21]) {
    const v = e2Curve3C(t, 5, d, k1, k2, k3, 0, 0, false, 0);
    row(`  t=${t}`, v);
  }
}

// -----------------------------------------------------------------------
// e2Curve3C - initial conditions (Ds and D2 paths)
// -----------------------------------------------------------------------
section("e2Curve3C - non-zero initial conditions (EV im params)");
{
  const [d, k1, k2, k3] = PKParameters["EV im"];
  console.log("\n// Ds=50, D2=0 (second compartment IC)");
  for (const t of [0, 1, 3, 7]) {
    const v = e2Curve3C(t, 0, d, k1, k2, k3, 50, 0, false, 0);
    row(`  t=${t}`, v);
  }
  console.log("\n// Ds=0, D2=100 (third compartment IC)");
  for (const t of [0, 1, 3, 7]) {
    const v = e2Curve3C(t, 0, d, k1, k2, k3, 0, 100, false, 0);
    row(`  t=${t}`, v);
  }
  console.log("\n// Ds=50, D2=100, dose=5 (all terms active)");
  for (const t of [0, 1, 3, 7]) {
    const v = e2Curve3C(t, 5, d, k1, k2, k3, 50, 100, false, 0);
    row(`  t=${t}`, v);
  }
}

// -----------------------------------------------------------------------
// e2Curve3C - equal-rate edge cases
// -----------------------------------------------------------------------
section("e2Curve3C - equal-rate edge cases (dose=1, d=1)");
{
  const k = 0.5;
  console.log("\n// k1==k2==k3");
  for (const t of [1, 3, 7]) {
    const v = e2Curve3C(t, 1, 1, k, k, k, 0, 0, false, 0);
    row(`  t=${t}`, v);
  }
  console.log("\n// k1==k2, k2!=k3");
  for (const t of [1, 3, 7]) {
    const v = e2Curve3C(t, 1, 1, k, k, 1.5, 0, 0, false, 0);
    row(`  t=${t}`, v);
  }
  console.log("\n// k1!=k2, k1==k3");
  for (const t of [1, 3, 7]) {
    const v = e2Curve3C(t, 1, 1, k, 1.5, k, 0, 0, false, 0);
    row(`  t=${t}`, v);
  }
  console.log("\n// k1!=k2, k2==k3");
  for (const t of [1, 3, 7]) {
    const v = e2Curve3C(t, 1, 1, 1.5, k, k, 0, 0, false, 0);
    row(`  t=${t}`, v);
  }
}

// -----------------------------------------------------------------------
// esSingleDose3C - EV im params
// -----------------------------------------------------------------------
section("esSingleDose3C - EV im, dose=5");
{
  const [d, k1, k2, k3] = PKParameters["EV im"];
  for (const t of [0, 1, 3, 7]) {
    const v = esSingleDose3C(t, 5, d, k1, k2, k3, 0);
    row(`  t=${t}`, v);
  }
  console.log("\n// with Ds=50 IC");
  for (const t of [0, 1, 3, 7]) {
    const v = esSingleDose3C(t, 5, d, k1, k2, k3, 50);
    row(`  t=${t}`, v);
  }
}

// -----------------------------------------------------------------------
// esSingleDose3C - k1==k2 degenerate branch
// -----------------------------------------------------------------------
section("esSingleDose3C - k1==k2 degenerate branch (dose=1, d=1)");
{
  const k = 0.5;
  console.log("\n// k1==k2, no IC");
  for (const t of [1, 3, 7]) {
    const v = esSingleDose3C(t, 1, 1, k, k, k, 0);
    row(`  t=${t}`, v);
  }
  console.log("\n// k1==k2, Ds=50 IC");
  for (const t of [0, 1, 3, 7]) {
    const v = esSingleDose3C(t, 1, 1, k, k, k, 50);
    row(`  t=${t}`, v);
  }
}

// -----------------------------------------------------------------------
// e2SteadyState3C - all injection models
// -----------------------------------------------------------------------
section("e2SteadyState3C - 5mg every 7 days, all injection models");
for (const model of injectionModels) {
  console.log(`\n// ${model}`);
  const [d, k1, k2, k3] = PKParameters[model];
  for (const t of [0.5, 2, 4, 6.9]) {
    const v = e2SteadyState3C(t, 5, 7, d, k1, k2, k3);
    row(`  t=${t}`, v);
  }
}

// -----------------------------------------------------------------------
// e2Patch3C - single dose, both patch models
// -----------------------------------------------------------------------
section("e2Patch3C - single dose 100mcg/day, both patch models");
for (const [model, W] of [
  ["patch tw", 3.5],
  ["patch ow", 7.0],
]) {
  console.log(`\n// ${model} (W=${W})`);
  const [d, k1, k2, k3] = PKParameters[model];
  for (const t of [0, 1, W - 0.1, W, W + 0.1, W + 2, W * 2]) {
    const v = e2Patch3C(t, 100, d, k1, k2, k3, W, false, 0);
    row(`  t=${t.toFixed(1)}`, v);
  }
}

// -----------------------------------------------------------------------
// e2SteadyStatePatch3C - both patch models
// -----------------------------------------------------------------------
section("e2SteadyStatePatch3C - 100mcg/day steady state");
for (const [model, W, T] of [
  ["patch tw", 3.5, 3.5],
  ["patch ow", 7.0, 7.0],
]) {
  console.log(`\n// ${model} (W=${W}, T=${T})`);
  const [d, k1, k2, k3] = PKParameters[model];
  for (const t of [0.5, 1.5, 2.5, T - 0.1, T, T + 0.5]) {
    const v = e2SteadyStatePatch3C(t, 100, T, d, k1, k2, k3, W);
    row(`  t=${t.toFixed(1)}`, v);
  }
}

// -----------------------------------------------------------------------
// PKFunctions dispatch - all models, single dose and steady state
// -----------------------------------------------------------------------
section("PKFunctions dispatch - single dose and steady state, all models");
const pf = PKFunctions(1.0);
for (const model of [...injectionModels, "patch tw", "patch ow"]) {
  console.log(`\n// ${model}`);
  for (const t of [1, 3, 7]) {
    const v = pf[model](t, 5, false, 0);
    row(`  single t=${t}`, v);
  }
  for (const t of [1, 3, 6]) {
    const v = pf[model](t, 5, true, 7);
    row(`  ss     t=${t}`, v);
  }
}

// -----------------------------------------------------------------------
// e2ssAverage3C - all injection models
// -----------------------------------------------------------------------
section("e2ssAverage3C - 5mg every 7 days, all injection models");
for (const model of injectionModels) {
  const v = e2ssAverage3C(5, 7, ...PKParameters[model]);
  row(`  ${model}`, v);
}

// -----------------------------------------------------------------------
// e2multidose3C - mixed models, intervals mode
// -----------------------------------------------------------------------
section("e2multidose3C - mixed models");
{
  console.log("\n// EV im + EC im, absolute times [0, 5]");
  for (const t of [3, 5, 8, 14]) {
    const v = e2multidose3C(
      t,
      [5, 5],
      [0, 5],
      ["EV im", "EC im"],
      1.0,
      false,
      false,
    );
    row(`  t=${t}`, v);
  }
  console.log(
    "\n// EV im + EV im, intervals mode [0, 7] (same as absolute [0,7])",
  );
  for (const t of [3, 7, 10, 14]) {
    const v = e2multidose3C(
      t,
      [5, 5],
      [0, 7],
      ["EV im", "EV im"],
      1.0,
      false,
      true,
    );
    row(`  t=${t}`, v);
  }
  console.log(
    "\n// three doses: EV im 5mg at 0, EEn im 3mg at 7, EC im 4mg at 14",
  );
  for (const t of [5, 10, 16, 21]) {
    const v = e2multidose3C(
      t,
      [5, 3, 4],
      [0, 7, 14],
      ["EV im", "EEn im", "EC im"],
      1.0,
      false,
      false,
    );
    row(`  t=${t}`, v);
  }
}

// -----------------------------------------------------------------------
// Unit conversion - pmol/L (cf=3.6713)
// -----------------------------------------------------------------------
section("Unit conversion - EV im 5mg single dose in pmol/L (cf=3.6713)");
{
  const cf = 3.6713;
  const [d, k1, k2, k3] = PKParameters["EV im"];
  for (const t of [1, 3, 7, 14]) {
    const v = e2Curve3C(t, cf * 5, d, k1, k2, k3, 0, 0, false, 0);
    row(`  t=${t}`, v);
  }
}

// -----------------------------------------------------------------------
// fillMenstrualCycleCurve - knot values (integer days)
// -----------------------------------------------------------------------
section(
  "fillMenstrualCycleCurve - knot values (cf=1.0, 29 steps over 28 days)",
);
{
  const curve = fillMenstrualCycleCurve(0, 28, 29, 1.0);
  for (const pt of curve.filter((p) => Number.isInteger(p.Time))) {
    console.log(
      `  t=${pt.Time}: e2=${pt.E2.toFixed(6)}, p5=${pt.E2p5.toFixed(6)}, p95=${pt.E2p95.toFixed(6)}`,
    );
  }
}

// -----------------------------------------------------------------------
// fillMenstrualCycleCurve - mid-interval interpolated values
// This is the critical spline test: at knot points any implementation
// returns the data exactly. Only between knots does the algorithm matter.
// -----------------------------------------------------------------------
section("fillMenstrualCycleCurve - mid-interval interpolated values (cf=1.0)");
{
  // Generate 10000 points so chosen targets land very close to exact values
  const N = 10000;
  const curve = fillMenstrualCycleCurve(0, 28, N, 1.0);
  const step = 28 / (N - 1);

  const midpoints = [
    0.5,
    1.5,
    2.5,
    3.5,
    4.5,
    5.5,
    6.5, // early follicular
    7.5,
    8.5,
    9.5,
    10.5, // late follicular
    11.5,
    12.5,
    13.5, // pre-ovulation peak (sharpest curve)
    14.5,
    15.5,
    16.5, // post-ovulation drop
    17.5,
    18.5,
    19.5,
    20.5, // luteal rise
    21.5,
    22.5,
    23.5,
    24.5,
    25.5,
    26.5,
    27.5, // luteal decline
  ];

  for (const target of midpoints) {
    const idx = Math.round(target / step);
    const pt = curve[idx];
    console.log(
      `  t=${pt.Time.toFixed(4)}: e2=${pt.E2.toFixed(8)}, p5=${pt.E2p5.toFixed(8)}, p95=${pt.E2p95.toFixed(8)}`,
    );
  }
}

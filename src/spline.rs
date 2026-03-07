//! Hermite cubic spline - exact port of `cubic-spline.js`
//!
//! The package solves for first derivatives (`ks`) at each knot using a
//! tridiagonal system with Gaussian elimination + partial pivoting, then
//! evaluates using the Hermite basis form.  This is NOT a natural spline
//! (zero second derivative at endpoints); the boundary conditions are:
//!
//!   row 0:   2/dx_0 * k_0  +  1/dx_0 * k_1  =  3*dy_0/dx_0²
//!   row n:   1/dx_n * k_{n-1}  +  2/dx_n * k_n  =  3*dy_n/dx_n²

/// A Hermite cubic spline whose algorithm matches `cubic-spline.js`
#[derive(Debug, Clone)]
pub struct CubicSpline {
    xs: Vec<f64>,
    ys: Vec<f64>,
    /// First derivatives at each knot, solved by `get_natural_ks`
    ks: Vec<f64>,
}

impl CubicSpline {
    pub fn new(xs: &[f64], ys: &[f64]) -> Self {
        assert_eq!(xs.len(), ys.len());
        assert!(xs.len() >= 2);
        let ks = get_natural_ks(xs, ys);
        Self {
            xs: xs.to_vec(),
            ys: ys.to_vec(),
            ks,
        }
    }

    /// Evaluate the spline at `x`.
    pub fn at(&self, x: f64) -> f64 {
        let i = self.get_index_before(x);
        let dx = self.xs[i] - self.xs[i - 1];
        let dy = self.ys[i] - self.ys[i - 1];
        let t = (x - self.xs[i - 1]) / dx;
        let a = self.ks[i - 1] * dx - dy;
        let b = -self.ks[i] * dx + dy;
        (1.0 - t) * self.ys[i - 1] + t * self.ys[i] + t * (1.0 - t) * (a * (1.0 - t) + b * t)
    }

    /// Binary search: returns the index `s` such that `xs[s-1] <= x < xs[s]`.
    /// Matches the JS `getIndexBefore` exactly.
    fn get_index_before(&self, x: f64) -> usize {
        let mut lo = 0usize;
        let mut hi = self.xs.len();
        loop {
            let mid = (lo + hi) / 2;
            if self.xs[mid] < x && mid != lo {
                lo = mid;
            } else if self.xs[mid] >= x && mid != hi {
                hi = mid;
            } else {
                break;
            }
        }
        // clamp to valid interval
        (lo + 1).min(self.xs.len() - 1)
    }
}

/// Build the `ks` (first-derivative) vector.
/// Direct port of `getNaturalKs` + `solveLowerTriangular` (`q`) from the JS.
fn get_natural_ks(xs: &[f64], ys: &[f64]) -> Vec<f64> {
    let n = xs.len() - 1; // number of intervals
    let sz = n + 1; // number of knots
    let aug = sz + 1; // augmented matrix width

    // Augmented matrix: sz rows × (sz+1) cols, row-major
    let mut mat = vec![vec![0.0_f64; aug]; sz];

    // Interior rows
    for i in 1..n {
        let dx_l = xs[i] - xs[i - 1];
        let dx_r = xs[i + 1] - xs[i];
        let dy_l = ys[i] - ys[i - 1];
        let dy_r = ys[i + 1] - ys[i];
        mat[i][i - 1] = 1.0 / dx_l;
        mat[i][i] = 2.0 * (1.0 / dx_l + 1.0 / dx_r);
        mat[i][i + 1] = 1.0 / dx_r;
        mat[i][sz] = 3.0 * (dy_l / (dx_l * dx_l) + dy_r / (dx_r * dx_r));
    }

    // First row (boundary)
    let dx0 = xs[1] - xs[0];
    let dy0 = ys[1] - ys[0];
    mat[0][0] = 2.0 / dx0;
    mat[0][1] = 1.0 / dx0;
    mat[0][sz] = 3.0 * dy0 / (dx0 * dx0);

    // Last row (boundary)
    let dxn = xs[n] - xs[n - 1];
    let dyn_ = ys[n] - ys[n - 1];
    mat[n][n - 1] = 1.0 / dxn;
    mat[n][n] = 2.0 / dxn;
    mat[n][sz] = 3.0 * dyn_ / (dxn * dxn);

    solve_gauss(&mut mat, sz)
}

/// Gaussian elimination with partial pivoting.
/// Port of the `q()` function in cubic-spline.js.
fn solve_gauss(mat: &mut Vec<Vec<f64>>, sz: usize) -> Vec<f64> {
    let aug = sz + 1;
    let mut row = 0usize;
    let mut col = 0usize;

    while row < sz && col <= sz {
        // Find pivot
        let mut max_val = f64::NEG_INFINITY;
        let mut max_row = row;
        for r in row..sz {
            let v = mat[r][col].abs();
            if v > max_val {
                max_val = v;
                max_row = r;
            }
        }

        if mat[max_row][col] == 0.0 {
            col += 1;
            continue;
        }

        mat.swap(row, max_row);

        for r in (row + 1)..sz {
            let factor = mat[r][col] / mat[row][col];
            mat[r][col] = 0.0;
            for c in (col + 1)..aug {
                let sub = mat[row][c] * factor;
                mat[r][c] -= sub;
            }
        }
        row += 1;
        col += 1;
    }

    // Back-substitution
    let mut ks = vec![0.0_f64; sz];
    for i in (0..sz).rev() {
        let f = if mat[i][i] != 0.0 {
            mat[i][sz] / mat[i][i]
        } else {
            0.0
        };
        ks[i] = f;
        for l in (0..i).rev() {
            mat[l][sz] -= mat[l][i] * f;
            mat[l][i] = 0.0;
        }
    }
    ks
}

#[cfg(test)]
mod tests {
    use super::*;

    // Knot values are reproduced exactly
    #[test]
    fn interpolates_knots_exactly() {
        let xs = vec![0.0, 1.0, 2.0, 3.0, 4.0];
        let ys = vec![0.0, 1.0, 0.0, 1.0, 0.0];
        let s = CubicSpline::new(&xs, &ys);
        for (&x, &y) in xs.iter().zip(ys.iter()) {
            let v = s.at(x);
            assert!((v - y).abs() < 1e-10, "at x={x}: got {v}, want {y}");
        }
    }

    // Linear data is reproduced exactly between knots
    #[test]
    fn linear_data_is_exact() {
        let xs = vec![0.0, 1.0, 2.0, 3.0];
        let ys = vec![0.0, 2.0, 4.0, 6.0];
        let s = CubicSpline::new(&xs, &ys);
        assert!((s.at(1.5) - 3.0).abs() < 1e-10);
    }

    // Mid-interval value matches JS reference at t=0.5013
    // (full menstrual cycle knots, JS value = 39.86835433)
    #[test]
    fn menstrual_cycle_midpoint_matches_js() {
        let xs: Vec<f64> = (0..30).map(|i| i as f64).collect();
        let ys = vec![
            37.99, 40.59, 37.49, 34.99, 35.49, 39.54, 41.99, 44.34, 53.43, 58.58, 71.43, 98.92,
            132.31, 177.35, 255.88, 182.80, 85.23, 70.98, 87.97, 109.92, 122.77, 132.56, 150.30,
            133.81, 137.16, 134.96, 92.73, 85.68, 46.34, 41.19,
        ];
        let s = CubicSpline::new(&xs, &ys);
        let got = s.at(0.5013);
        let js = 39.86835433_f64;
        let err = (got - js).abs() / js;
        assert!(err < 1e-4, "got {got:.8}, js {js:.8}, err {err:.2e}");
    }
}

// JIBES implemented in Rust, this code is a mirror of the python code in cellranger.analysis.jibes
// with the regression and model fitting steps optimized to use significantly less memory and run
// a lot faster.  Any changes here should be reflected in the python code and vice-versa.

use ndarray::prelude::*;
use ndarray::{Array, Array1, Array2};
use ndarray_stats::QuantileExt;
use numpy::{IntoPyArray, PyArray1, PyArray2, PyReadonlyArray1};
use pyo3::prelude::*;
use pyo3::types::PyList;
//use pyo3::{wrap_pyfunction};
//use statrs::function::gamma::ln_gamma;
use statrs::distribution::{Continuous, Normal};
use statrs::function::factorial::ln_factorial;
use std::cmp::max;
use std::collections::HashMap;
use std::iter::zip;

static _MAX_K_LETS_TO_CONSIDER: i32 = 3;
static MIN_FOREGROUND_DELTA: f64 = 0.01;

fn lfactorial(x: f64) -> f64 {
    // Log factorial function
    // Via gamma
    // let gam_version = ln_gamma(x + 1.0);

    // Via cached recursition results

    //assert_eq!(gam_version, fac_version);
    ln_factorial(x as u64)
}

fn lmultinomial_comb(cnts: &Array1<f64>) -> f64 {
    // Log version of multinomial coefficient
    let numerator = lfactorial(cnts.sum());
    let denominator: f64 = cnts.iter().map(|x| lfactorial(*x)).sum();
    numerator - denominator
}
#[pyclass]
#[derive(Debug, Clone)]
struct JibesModelO3 {
    #[pyo3(get)]
    pub num_tags: usize,
    #[pyo3(get)]
    pub blank_freq: f64,
    #[pyo3(get)]
    pub background: Vec<f64>,
    #[pyo3(get)]
    pub foreground: Vec<f64>,
    #[pyo3(get)]
    pub std_devs: Vec<f64>,
    #[pyo3(get)]
    pub frequencies: Vec<f64>,
    #[pyo3(get)]
    pub n_gems: i32,
}

#[pymethods]
impl JibesModelO3 {
    #[new]
    fn new(
        backgrounds: &PyArray1<f64>,
        foregrounds: &PyArray1<f64>,
        std_devs: &PyArray1<f64>,
        frequencies: Option<&PyArray1<f64>>,
        blank_prob: Option<f64>,
        n_gems: Option<i32>,
    ) -> Self {
        let num_tags = backgrounds.len();
        assert!(
            num_tags == std_devs.len(),
            "Backgrounds does not match std_devs"
        );
        assert!(
            num_tags == foregrounds.len(),
            "Tags does not match std_devs"
        );
        // In case we ever want to enable frequencies
        let freqs = match frequencies {
            Some(x) => {
                assert!(num_tags == x.len(), "Tags does not match std_devs");
                x.to_vec().unwrap()
            }
            None => Array::from_elem(num_tags, 1.0 / (num_tags as f64)).to_vec(),
        };
        let blank_freq = blank_prob.unwrap_or(0.04);
        let background = backgrounds.to_vec().unwrap();
        let foreground = foregrounds.to_vec().unwrap();
        for &z in &foreground {
            assert!(
                z >= MIN_FOREGROUND_DELTA,
                "Cannot initialize with foreground value {z} less than {MIN_FOREGROUND_DELTA}"
            );
        }
        let std_dev = std_devs.to_vec().unwrap();
        JibesModelO3 {
            num_tags,
            blank_freq,
            background,
            foreground,
            std_devs: std_dev,
            frequencies: freqs,
            n_gems: n_gems.unwrap_or(95000),
        }
    }
}

#[pyclass]
struct JibesDataO3 {
    #[pyo3(get)]
    pub counts: PyObject,
    #[pyo3(get)]
    pub column_names: Py<PyList>,
    #[pyo3(get)]
    pub barcodes: Py<PyList>,
}

#[pyclass]
#[allow(non_snake_case)]
struct JibesEMO3 {
    #[pyo3(get)]
    pub n_gems: i32,
    #[pyo3(get)]
    pub max_k_let_setting: i16,
    #[pyo3(get)]
    pub optimize_pop_frequencies: bool,
    #[pyo3(get)]
    pub k_let_limited: bool,
    #[pyo3(get)]
    pub max_modeled_k_let: i16,
    // This has no "intercept" term for the regression
    pub latent_states: Array2<f64>,
    // Rust copy of the data held in
    // `data` on the python side
    pub counts: Array2<f64>,
    #[pyo3(get)]
    pub data: Option<PyObject>,
    #[pyo3(get)]
    pub estimated_cells: f64,
    #[pyo3(get)]
    pub model: JibesModelO3,
    #[pyo3(get)]
    pub LL: f64,
    #[pyo3(get)]
    pub converged: bool,
    #[pyo3(get)]
    pub iterations: i32,
    // Barcodes x Latent State Matrix
    pub posterior: Array2<f64>,
}

impl JibesEMO3 {
    fn calculate_latent_state_weights(&self) -> Array1<f64> {
        let cnts_vec: Vec<f64> = get_multiplet_counts_unrounded(self.estimated_cells, self.n_gems)
            .iter()
            .take(self.max_modeled_k_let as usize)
            .copied()
            .collect();
        let cnts = Array1::from(cnts_vec);
        let total = cnts.sum();

        let mut p_k_let: Array1<f64> = cnts / total;
        if self.k_let_limited {
            // We mush all the probability for the higher states into the last state that we want
            // to account for all that data
            let array_end_index = p_k_let.len() - 1;
            p_k_let[array_end_index] = p_k_let.iter().skip(self.max_k_let_setting as usize).sum();
        }
        let klet: Array<usize, Ix1> = self.latent_states.sum_axis(Axis(1)).map(|x| *x as usize);
        let state: Vec<f64> = klet.iter().skip(1).map(|k| p_k_let[k - 1].ln()).collect();
        let pis: Vec<f64> = self.model.frequencies.iter().map(|x| x.ln()).collect();
        let z = self.latent_states.dim().0;
        let mut comb_probs = Array1::zeros(z);
        comb_probs[0] = self.model.blank_freq.ln();
        let log_not_blank = (1.0 - self.model.blank_freq).ln();
        for zi in 1..z {
            let latent_state = self.latent_states.slice(s![zi, ..]);
            let relevant_ps: Vec<usize> = latent_state
                .iter()
                .enumerate()
                .filter(|x| *x.1 != 0.0)
                .map(|z| z.0)
                .collect();
            let cnts: Array1<f64> = relevant_ps.iter().map(|i| latent_state[*i]).collect();
            let ps: Array1<f64> = relevant_ps.iter().map(|i| pis[*i]).collect();
            let multiplier = lmultinomial_comb(&cnts);
            let p: f64 = zip(cnts, ps).map(|(cnt, p)| cnt * p).sum();
            comb_probs[zi] = p + multiplier + state[zi - 1] + log_not_blank;
        }
        comb_probs
    }
    fn maximize_parameters(&mut self) {
        // The maximization step in this algorithm is equivalent to a weighted linear
        // regression with weights given by the latent state probabilities, and covariates
        // given by the cells in each latent state.  Since each tag is independent, this
        // becomes a 2 parameter weighted least squares regression, which we calculate
        // more efficiently (though further optimizations possible) by reusing the same
        // covariate matrix for each data point instead of tiling it like on the python side.
        // Further details are in the IDF for the JIBES model.

        // Solution to find coefficents in weighted least squares model is (X^T W X)^-1 X^T W Y
        for k in 0..(self.k() as usize) {
            // solve for each tag independently

            // First calculate X^T(transpose) W Y, a 2 x 1 matrix
            let mut first_elem: f64 = 0.0; // This should just be the sum of the tag counts
                                           // if normalization is such that sum(w) = 1 for all data points, we still calculate
                                           // it below to avoid enforcing that assumption and to enable easier testing with
                                           // different weights.
            let mut second_elem: f64 = 0.0;
            for i in 0..(self.n() as usize) {
                let tag_cnt: f64 = self.counts[[i, k]];
                first_elem += tag_cnt * self.posterior.row(i).sum();
                second_elem += tag_cnt
                    * zip(self.posterior.row(i), self.latent_states.column(k))
                        .map(|x| *x.0 * *x.1)
                        .sum::<f64>();
            }
            #[allow(non_snake_case)]
            let XTWY = array![first_elem, second_elem];

            //Next calculate (X^T W X), a 2x2 symmetrical matrix
            // typically just self.n() if weights sum to one for each point
            let elem_11 = self.posterior.sum();
            let mut elem_12 = 0.0;
            let mut elem_22 = 0.0;
            for i in 0..(self.n() as usize) {
                // Calculate the off diagonal and diagonal element for the 2x2 matrix
                let prod = zip(self.posterior.row(i), self.latent_states.column(k))
                    .map(|x| (*x.0 * *x.1, *x.1))
                    .map(|z| (z.0, z.0 * z.1))
                    .fold((0.0, 0.0), |acc, val| (acc.0 + val.0, acc.1 + val.1));
                elem_12 += prod.0;
                elem_22 += prod.1;
            }
            // let XTWY = array![[elem_11, elem_12], [elem_12, elem_22]];
            // Get inverse of this 2 x 2
            let determinant = 1.0 / (elem_11 * elem_22 - elem_12 * elem_12);
            let inverse = determinant * array![[elem_22, -elem_12], [-elem_12, elem_11]];

            // Multiply inverse by other matrix to get regression coefficients
            let mut coefficients = inverse.dot(&XTWY);

            // Check if the coefficient is negative and correct it if so
            if coefficients[1] < MIN_FOREGROUND_DELTA {
                // Force the coefficient to be positive
                coefficients[1] = MIN_FOREGROUND_DELTA;
                // Re-estimate intercept term conditioned on the slope value.
                // which is the classic B_0 = mean(Y_w) - B_1 * mean(X_w)
                // where mean(Z_w) = sum(w*z) / sum(w)
                let x_bar = elem_12 / elem_11;
                let y_bar = first_elem / elem_11;
                coefficients[0] = y_bar - MIN_FOREGROUND_DELTA * x_bar;
            }

            self.model.background[k] = coefficients[0];
            self.model.foreground[k] = coefficients[1];

            // Calculate residuals to estimate variance
            let predicted = self
                .latent_states
                .column(k)
                .map(|x| coefficients[1] * x + coefficients[0]);
            let mut residuals = 0.0;
            for i in 0..(self.n() as usize) {
                let weights = self.posterior.row(i);
                let tag_cnt = array![self.counts[[i, k]]];
                let cur_res = &predicted - &tag_cnt;
                residuals += zip(cur_res, weights).map(|x| x.0 * x.0 * x.1).sum::<f64>();
            }

            // Note we are using the biased MLE estimate here, R uses the unbiased estimate
            // which divides by N-2 instead of N so consider that when comparing results.
            // In practice for a reasonable number of barcodes (>1000)
            // the difference will be trivial.
            let std_dev_estimate = (residuals / elem_11).sqrt();
            // Never allow the variance to go to 0, which would only happen if
            // only one data point is in this component, but would cause the LL to go to infinity
            self.model.std_devs[k] = std_dev_estimate.max(0.0001);
            assert!(
                !self.optimize_pop_frequencies,
                "Frequency optimization not added yet."
            );
        }
    }
}

#[pymethods]
impl JibesEMO3 {
    #[getter]
    #[allow(non_snake_case)]
    fn get_X<'py>(&self, py: Python<'py>) -> &'py PyArray2<f64> {
        // Get the regression matrix present in the python code (need to add intercept term to latent states)

        let mut array = Array2::<f64>::zeros((self.z() as usize, (self.k() + 1) as usize));
        //  = PyArray2::<f64>::zeros(gil.python(), [self.z() as usize, (self.k() + 1) as usize], false );
        // add intercept term
        for i in 0..array.nrows() {
            array[[i, 0]] = 1.0;
            for k in 0..self.latent_states.ncols() {
                array[[i, k + 1]] = self.latent_states[[i, k]];
            }
        }
        //let gil = Python::acquire_gil();
        array.into_pyarray(py)
    }

    #[getter]
    fn get_posterior<'py>(&self, py: Python<'py>) -> &'py PyArray2<f64> {
        // Get the regression matrix present in the python code (need to add intercept term to latent states)
        PyArray2::from_array(py, &self.posterior)
    }

    #[new]
    #[pyo3(signature=(data, model, optimize_pop_freqs=false, max_k_lets=3, n_gems=95000))]
    fn new(
        py: Python<'_>,
        data: PyObject,
        model: &JibesModelO3,
        optimize_pop_freqs: bool,
        max_k_lets: usize,
        n_gems: i32,
    ) -> Self {
        // Copy this to rust to make it easier to deal with
        let cnt_data =
            pyarray2_to_array2_f64(data.getattr(py, "counts").unwrap().extract(py).unwrap());
        // What's the largest k-let we're likely to see?
        let n = cnt_data.dim().0;
        let estimated_cells = calculate_expected_total_cells(n as i32, n_gems)
            .expect("calculate_expected_total_cells");
        let exp_cnts =
            get_multiplet_counts(py, estimated_cells, n_gems).expect("get_multiplet_counts");

        let max_found_multiplets = exp_cnts
            .as_array()
            .iter()
            .enumerate()
            .filter(|x| *x.1 > 0.0)
            .map(|e| e.0)
            .max();
        let mut max_multiplets = match max_found_multiplets {
            Some(x) => max(x + 1, 2),
            None => {
                println!(
                    "Search for prob < 1 failed, crazy high cell load or zero cells has occurred."
                );
                4
            }
        };
        let mut k_let_limited = false;
        if max_multiplets > max_k_lets {
            k_let_limited = true;
            max_multiplets = max_k_lets;
            println!("Limiting Latent States to K-lets of {max_k_lets} with a blank state");
        }
        let latent_states =
            generate_all_multiplets(model.num_tags as i32, max_multiplets as i32, k_let_limited);
        let max_modeled_k_let = max(model.num_tags, max_multiplets);
        let z = latent_states.dim().0;
        println!("Total Latent States = {z}");
        println!("Observed Barcodes = {n}");
        println!("Inferred Cells = {estimated_cells}");
        JibesEMO3 {
            n_gems,
            max_k_let_setting: max_k_lets as i16,
            optimize_pop_frequencies: optimize_pop_freqs,
            k_let_limited,
            max_modeled_k_let: max_modeled_k_let as i16,
            latent_states,
            estimated_cells,
            counts: cnt_data,
            data: Some(data),
            model: model.clone(),
            LL: f64::NEG_INFINITY,
            converged: false,
            iterations: 0,
            posterior: Array::zeros((n, z)),
        }
    }

    fn _calculate_latent_state_weights(&self) -> Vec<f64> {
        let weights = self.calculate_latent_state_weights();
        weights.to_vec()
    }

    fn _calculate_posterior_by_state(&mut self) {
        // Make a barcodes x latent state matrix with the the log
        // likelihood of coming from any state

        // Clear matrix just for safety
        self.posterior.fill(0.0);

        // Calculate the prior and fill all the latent states with it
        let z_prior = self.calculate_latent_state_weights();
        for z in 0..z_prior.len() {
            for bc in 0..self.posterior.dim().0 {
                self.posterior[[bc, z]] = z_prior[z];
            }
        }

        // Iterate through and calculate the posterior for each latent state
        for state in 0..(self.z() as usize) {
            let cur_state = self.latent_states.slice(s![state, ..]);
            for tag in 0..(self.k() as usize) {
                let mu = self.model.background[tag] + self.model.foreground[tag] * cur_state[tag];
                let std_dev = self.model.std_devs[tag];
                let dist = Normal::new(mu, std_dev).unwrap();
                for bc in 0..(self.n() as usize) {
                    // Slowest portion of the code right now, a simple way to
                    // speed this up would be to cache results per barcode as
                    // most latent states have similar marginal distributions for a given
                    // tag.
                    self.posterior[[bc, state]] += dist.ln_pdf(self.counts[[bc, tag]]);
                }
            }
        }
        // Now to normalize, log-sum-exp operation
        // 1- subtract max - matches python
        let ll_max: Array1<f64> = self
            .posterior
            .axis_iter(Axis(0))
            .map(|x| *x.max().unwrap())
            .collect();
        for bc in 0..self.posterior.nrows() {
            //subtract max
            self.posterior.row_mut(bc).map_inplace(|x| *x -= ll_max[bc]);
        }
        // 2 - calc exp
        self.posterior.mapv_inplace(f64::exp); // okay here
                                               // 3 - normalize
        let mut marginal: Array1<f64> =
            self.posterior.axis_iter(Axis(0)).map(|x| x.sum()).collect();
        for (mut x, marg) in zip(self.posterior.axis_iter_mut(Axis(0)), &marginal) {
            x.map_inplace(|c| *c /= marg);
        }
        // 4 - Calculate LL
        marginal.mapv_inplace(f64::ln);
        self.LL = marginal.sum() + ll_max.sum();
    }

    #[allow(non_snake_case)]
    fn one_EM_step(&mut self) -> f64 {
        // Perform one step of the EM algorithm.
        // :returns The current LL of the model.
        if self.LL == f64::NEG_INFINITY {
            // initialize
            self._calculate_posterior_by_state();
        }
        println!("LL = {}", self.LL);
        self.maximize_parameters();
        self._calculate_posterior_by_state();
        self.iterations += 1;
        self.LL
    }
    #[allow(non_snake_case)]
    #[pyo3(signature=(max_reps=50000, abs_tol=1e-2, rel_tol=1e-7))]
    fn perform_EM(&mut self, max_reps: i32, abs_tol: f64, rel_tol: f64) -> PyResult<f64> {
        //  Run the EM algorithm until either max_reps is exceeded or the change in  the LL after
        // a run is <= tol or 1 - new/old <= rel_tol
        // :param max_reps: Maximum number of repetitions
        // :param abs_tol: Absolute difference in LL which would terminate the algorithm
        // :param rel_tol: Relative difference in LL which would terminate the algorithm
        // :return: The LL of the model after fitting
        let mut last_ll = self.LL;
        let mut rep = 0;
        if self.counts.dim().0 == 0 {
            // Don't perform EM if there's no cells, just set to zero.
            self.LL = 0.0;
            self.converged = true;
            for k in 0..(self.k() as usize) {
                self.model.foreground[k] = 0.0;
            }
        } else {
            loop {
                self.one_EM_step();
                rep += 1;
                let rel_change = 1.0 - self.LL / last_ll;
                let abs_change = self.LL - last_ll;
                if rep > max_reps {
                    break;
                }
                if !(last_ll == f64::NEG_INFINITY)
                    & ((abs_change < abs_tol) | (rel_change < rel_tol))
                {
                    println!("EM algorithm has converged.");
                    self.converged = true;
                    break;
                }
                last_ll = self.LL;
            }
        }
        Ok(self.LL)
    }

    fn get_snr_dictionary(&self, py: Python<'_>) -> HashMap<String, f64> {
        assert!(
            self.converged,
            "Model must converge before SNRs can be reported."
        );
        let snrs: Vec<f64> = zip(&self.model.foreground, &self.model.std_devs)
            .map(|(signal, noise)| signal / noise)
            .collect();
        let column_names: Vec<Vec<u8>> = self
            .data
            .as_ref()
            .unwrap()
            .getattr(py, "column_names")
            .unwrap()
            .extract(py)
            .unwrap();
        let mut snr_mapping: HashMap<String, f64> = HashMap::new();
        for tag in column_names.iter().enumerate() {
            let key = format!("snr_{}", std::str::from_utf8(tag.1).unwrap());
            snr_mapping.insert(key, snrs[tag.0]);
        }
        snr_mapping
    }

    #[getter]
    fn z(&self) -> i32 {
        self.latent_states.nrows() as i32
    }
    #[getter]
    fn k(&self) -> i32 {
        self.counts.ncols() as i32
    }
    #[getter]
    fn n(&self) -> i32 {
        self.counts.nrows() as i32
    }
}

// Rather convoluted way to call the original method, perhaps there's something better?
#[pyfunction]
fn calculate_expected_total_cells(n_bcs: i32, n_gems: i32) -> PyResult<f64> {
    Python::with_gil(|py| {
        let cr_fa = PyModule::import(py, "cellranger.feature.feature_assigner")?;
        let total = cr_fa
            .getattr("calculate_expected_total_cells")?
            .call1((n_bcs, n_gems))?
            .extract()?;
        Ok(total)
    })
}

fn generate_all_multiplets(num_tags: i32, max_multiplets: i32, k_let_limited: bool) -> Array2<f64> {
    Python::with_gil(|py| {
        let combinatorics = PyModule::import(py, "cellranger.analysis.combinatorics").unwrap();
        let total: &PyArray2<i32> = combinatorics
            .getattr("generate_all_multiplets_as_array")
            .unwrap()
            .call1((num_tags, max_multiplets, k_let_limited))
            .unwrap()
            .downcast::<PyArray2<i32>>()
            .unwrap();
        pyarray2_to_array2_i32(total)
    })
}

fn get_multiplet_counts_unrounded(estimated_cells: f64, n_gems: i32) -> Vec<f64> {
    Python::with_gil(|py| {
        let cr_fa = PyModule::import(py, "cellranger.feature.feature_assigner").unwrap();
        let probs: &PyArray1<f64> = cr_fa
            .getattr("get_multiplet_counts_unrounded")
            .unwrap()
            .call1((estimated_cells, n_gems))
            .unwrap()
            .extract()
            .unwrap();
        //println!("Probs = {}", probs);
        probs.to_vec().unwrap()
    })
}

fn get_multiplet_counts<'py>(
    py: Python<'py>,
    estimated_cells: f64,
    n_gems: i32,
) -> PyResult<PyReadonlyArray1<'py, f64>> {
    let cr_fa = PyModule::import(py, "cellranger.feature.feature_assigner")?;
    let probs: PyReadonlyArray1<'py, f64> = cr_fa
        .getattr("get_multiplet_counts")?
        .call1((estimated_cells, n_gems))?
        .extract()?;
    Ok(probs)
}

// These should be generic but ran into some small issue with PyArray<T> when I tried,
// just replicated two functions for now.
fn pyarray2_to_array2_i32(src: &PyArray2<i32>) -> Array2<f64> {
    let dims = src.dims();
    let rows = dims[0];
    let cols = dims[1];
    let mut result = Array2::<f64>::zeros(dims);
    for row in 0..rows {
        for col in 0..cols {
            unsafe {
                result[[row, col]] = *src.get([row, col]).unwrap() as f64;
            }
        }
    }
    result
}

fn pyarray2_to_array2_f64(src: &PyArray2<f64>) -> Array2<f64> {
    let dims = src.dims();
    let rows = dims[0];
    let cols = dims[1];
    let mut result = Array2::<f64>::zeros(dims);
    for row in 0..rows {
        for col in 0..cols {
            unsafe {
                result[[row, col]] = *src.get([row, col]).unwrap();
            }
        }
    }
    result
}

/// The jibes python classes implemented in Rust
#[pymodule]
fn jibes_o3(_py: Python<'_>, m: &PyModule) -> PyResult<()> {
    m.add_class::<JibesModelO3>()?;
    m.add_class::<JibesEMO3>()?;
    // This function is was added for one unit test
    m.add_function(wrap_pyfunction!(calculate_expected_total_cells, m)?)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    macro_rules! assert_delta {
        ($x:expr, $y:expr, $d:expr) => {
            if !($x - $y < $d || $y - $x < $d) {
                panic!();
            }
        };
    }
    use super::*;
    // use pyo3::types::IntoPyDict;
    #[test]
    fn test_wls_regression() {
        // Cheesy tests with equal "weights" to the wls regression, quick sanity check
        let latent_states = array![[0.0, 0.0], [0.0, 1.0], [1.0, 0.0], [1.0, 1.0]];
        let counts = array![[10.0, 400.0], [1000.0, 50.0]];
        let posterior = Array2::ones((2, 4));
        let model = JibesModelO3 {
            n_gems: 95000,
            std_devs: vec![1.0, 2.0],
            num_tags: 2,
            foreground: vec![5.0, 5.0],
            background: vec![1.0, 1.0],
            blank_freq: 0.04,
            frequencies: vec![0.5, 0.5],
        };
        // make a fake pyobject
        //let gil = Python::acquire_gil();
        //let py = gil.python();
        //let dummy = (0..3 as u64).map(|i| (i, i * 2)).into_py_dict(py).into();
        let mut fitter = JibesEMO3 {
            n_gems: 95000,
            max_k_let_setting: 2,
            optimize_pop_frequencies: false,
            k_let_limited: false,
            max_modeled_k_let: 2,
            latent_states,
            counts,
            data: None,
            estimated_cells: -2.0,
            model,
            LL: f64::NEG_INFINITY,
            converged: false,
            iterations: 0,
            posterior,
        };
        fitter.maximize_parameters();
        assert_delta!(fitter.model.background[0], 505.0, 0.00001);
        assert_delta!(fitter.model.foreground[0], 0.0, 0.00001);
    }
}

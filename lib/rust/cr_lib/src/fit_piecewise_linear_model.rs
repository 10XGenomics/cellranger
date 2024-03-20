use itertools::Itertools;
use ndarray::prelude::*;
use ndarray::{concatenate, Array1, NdFloat, Zip};
use num::FromPrimitive;

/// Local helper function that computes the dot product between
/// two one dimensional vectors (ndarray::Array1).
/// Using this rather than ndarray::dot as the latter calls intel-mkl
/// which we do not currently have in cr_lib
fn dot_product<F>(a: &Array1<F>, b: &Array1<F>) -> F
where
    F: NdFloat,
{
    Zip::from(a)
        .and(b)
        .fold(F::zero(), |acc: F, w, z| acc + *w * *z)
}

/// The solution of a two dimensional linear regression
/// Using offset and slope. The model is
/// y = x  * offset + slope.
struct TwoDimLinearRegressionParameter<F: FromPrimitive + NdFloat> {
    pub offset: F,
    pub slope: F,
}

/// Local helper function that solves a two dimensional linear regression
/// The model being solved is
/// [1 x] *[offset, slope]^T = y
/// TODO: Solve things using ndarray-linalg. The bottleneck there is that
/// we can not import ndarray-linalg to cr_lib.
fn solve_two_dim_linear_regression<F>(
    x: &Array1<F>,
    y: &Array1<F>,
) -> Option<TwoDimLinearRegressionParameter<F>>
where
    F: FromPrimitive + NdFloat,
{
    let dim = x.len();
    let ones_vec = Array1::ones(dim);
    // We compute A = [1 x]^T [1 x] . And B = [1 x]^T  y
    // where A = [ a b ]
    //           [ b c ],
    // And B = [y1 y2]^T
    // We see that a = 1^T 1; b = 1^T x; c = x^T x.
    // And that y1 = 1^T y, y2 = x^T y.
    // Thus the solution is
    // 1/det * [  c  -b ] * [y1 y2]^T.
    //         [ -b   a ]
    let a = dot_product(&ones_vec, &ones_vec);
    let b = dot_product(&ones_vec, x);
    let c = dot_product(x, x);
    let y1 = dot_product(&ones_vec, y);
    let y2 = dot_product(x, y);
    let det = a * c - b * b;
    // If the determinant is zero that means that
    // x is zero or some multiple of 1.
    // hence we fit a 1D linear regression of y against 1.
    // which ends up being slope of zero and offset being
    // the parameter of 1D linear regression.
    // If a is zero (which should not happen) we return None.
    match (det.is_zero(), a.is_zero()) {
        (true, true) => None,
        (true, false) => Some(TwoDimLinearRegressionParameter {
            offset: y1 / a,
            slope: F::zero(),
        }),
        (false, _) => Some(TwoDimLinearRegressionParameter {
            offset: (c * y1 - b * y2) / det,
            slope: (-b * y1 + a * y2) / det,
        }),
    }
}

/// Struct describing the parameters of a piecewise linear model
/// The model is the following
/// f(x) =  self.constant, if x < self.critical_point;
///         self.constant + self.slope*(x-self.critical_point), if x >= self.critical_point.
#[derive(Copy, Clone)]
pub struct PiecewiseLinearModel<F: FromPrimitive + NdFloat> {
    pub constant: F,
    pub slope: F,
    pub critical_point: F,
}

///Struct describing the parameters  of an estimated model.
/// This holds a PieceWiseModel plus a residual error.
#[derive(Copy, Clone)]
pub struct EstimatedModel<F: FromPrimitive + NdFloat> {
    pub model: PiecewiseLinearModel<F>,
    pub rss: F,
}

/// Struct holding the data input into the model
#[derive(Clone)]
pub struct PiecewiseLinearData<F: FromPrimitive + NdFloat> {
    /// Assumes that there are at least 3 such genes
    /// Unspliced counts of genes with both spliced and unspliced probes
    /// log transformed. Fit assumes that these are sorted in increasing
    /// order.
    spliced_counts: Array1<F>,

    /// Spliced counts of genes with both spliced and unspliced probes
    /// log transformed.
    unspliced_counts: Array1<F>,
}

impl<F> PiecewiseLinearData<F>
where
    F: FromPrimitive + NdFloat,
{
    /// Given vectors of spliced and unspliced counts, returns a new PiecewiseLinearData
    /// Makes sure that we have at least 3 data points, and that the spliced
    /// and unspliced vectors are of the same length.
    pub fn new(spliced_counts: Vec<F>, unspliced_counts: Vec<F>) -> Self {
        assert_eq!(spliced_counts.len(), unspliced_counts.len());
        assert!(
            spliced_counts.len() >= 3,
            "Vector lengths are {}. Need to be at least 3.",
            spliced_counts.len()
        );
        let sorted_unspliced_counts = unspliced_counts
            .into_iter()
            .enumerate()
            .sorted_by(|(a, _), &(b, _)| {
                spliced_counts[*a].partial_cmp(&spliced_counts[b]).unwrap()
            })
            .map(|x| x.1);
        let sorted_spliced_counts = spliced_counts
            .into_iter()
            .sorted_by(|a, b| a.partial_cmp(b).unwrap());
        Self {
            spliced_counts: Array1::from_iter(sorted_spliced_counts),
            unspliced_counts: Array1::from_iter(sorted_unspliced_counts),
        }
    }

    /// Returns spliced count as a vector
    pub fn get_spliced_counts(&self) -> Vec<F> {
        self.spliced_counts.clone().to_vec()
    }

    /// Returns unspliced counts as a vector
    pub fn get_unspliced_counts(&self) -> Vec<F> {
        self.unspliced_counts.clone().to_vec()
    }

    /// Fit a piecewise linear model using self.spliced_counts[pivot_index] as knee point
    /// Will panic if pivot_index is larger than length of self.spliced_counts
    /// Crucially depends upon self.spliced_counts being sorted (else
    /// we'd have to sort it every fit doing unnecessary work)
    /// This reduces to solve a 2D linear regression as
    /// [1 X_covariate] * [offset slope]^T = self.unspliced_counts
    /// where X_covariate = [ 0, 0, ... 0, self.spliced_counts[pivot_index+1]-z, ..,  self.spliced_counts[dimension-1]-z]^T
    /// where z = self.spliced_counts[pivot_index].
    /// Even if model is ill defined, the covariates have rank 1. Thus
    /// solves the 1D linear regression.
    fn fit_at_pivot(&self, pivot_index: usize) -> EstimatedModel<F> {
        let x_covariate = concatenate!(
            Axis(0),
            Array::zeros(pivot_index),
            &self.spliced_counts.slice(s![pivot_index..]) - self.spliced_counts[pivot_index]
        );
        let soln = solve_two_dim_linear_regression(&x_covariate, &self.unspliced_counts).unwrap();
        let residuals = &self.unspliced_counts - x_covariate * soln.slope - soln.offset;
        let rss = dot_product(&residuals, &residuals);
        EstimatedModel {
            model: PiecewiseLinearModel {
                critical_point: self.spliced_counts[pivot_index],
                constant: soln.offset,
                slope: soln.slope,
            },
            rss,
        }
    }

    /// Fit a piecewise linear model to the data
    /// Crucially depends upon self.spliced_counts being sorted
    /// Uses the fact that only points in self.spliced_counts are possible
    /// knee points that we care about
    pub fn fit(&self) -> EstimatedModel<F> {
        let n_samples = self.spliced_counts.len();
        (1..n_samples - 1)
            .map(|x| self.fit_at_pivot(x))
            .min_by(|a, b| (a.rss).partial_cmp(&(b.rss)).unwrap())
            .unwrap()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use assert_approx_eq::assert_approx_eq;
    use ndarray::Array;
    use ndarray_rand::rand_distr::{StandardNormal, Uniform};
    use ndarray_rand::RandomExt;
    use rand::rngs::StdRng;
    use rand::SeedableRng;

    impl<F> PiecewiseLinearModel<F>
    where
        F: FromPrimitive + NdFloat,
    {
        /// Not completely necessary, but having fun return a closure rather than
        /// just a value
        pub fn get_predictor(self) -> impl Fn(F) -> F {
            move |b| match b.lt(&self.critical_point) {
                true => self.constant,
                false => self.constant + (b - self.critical_point) * self.slope,
            }
        }
    }

    // Generate Test data
    fn get_data(
        n_samples: usize,
        predictor: &impl Fn(f64) -> f64,
        state: u64,
    ) -> (Array1<f64>, Array1<f64>) {
        let mut rng = StdRng::seed_from_u64(state);
        let x = Array::random_using(n_samples, Uniform::new(0., 10.), &mut rng);
        let noise: Array1<f64> = Array::random_using(n_samples, StandardNormal, &mut rng);
        (x.clone(), x.map(|&z| predictor(z)) + 0.25 * noise)
    }

    #[test]
    fn test_piecewise_linear_with_random_sample() {
        let model_gt = PiecewiseLinearModel {
            constant: 4.0,
            slope: 1.0,
            critical_point: 4.0,
        };
        let gt_predictor = model_gt.get_predictor();
        let data_in = {
            let (x, y) = get_data(5000, &gt_predictor, 77);
            PiecewiseLinearData::new(x.to_vec(), y.to_vec())
        };
        let fit_model = data_in.fit();

        // Using assert_approx_eq! rather than assert_eq! as we're comparing
        // floating point numbers. Values returned in linux and macOS are
        // slightly different
        assert_approx_eq!(fit_model.rss, 312.47644505361575);
        assert_approx_eq!(fit_model.model.constant, 3.9995689041095446);
        assert_approx_eq!(fit_model.model.slope, 0.9991637367907095);
        assert_approx_eq!(fit_model.model.critical_point, 3.9993728477509705);
    }
}

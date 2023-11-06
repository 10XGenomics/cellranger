use itertools::Itertools;
use ndarray::{Array2, ShapeBuilder};
use num_traits::Num;
use sprs::{CsMatI, SpIndex};
use std::ops::AddAssign;

/// Get representatives for each cluster in a clustering.
/// The current implementation returns the centroid of the features in all
/// barcodes in a cluster as the cluster representative.
/// Index to cluster map is a an array of cluster values in order of
/// barcodes in the feature bc matrix
/// Cluster numbers are 1..max_clusters
pub fn get_cluster_representatives<N, I>(
    index_to_cluster_map: &[i64],
    feature_bc_matrix: &CsMatI<N, I>,
) -> Array2<f64>
where
    N: Num + Clone,
    f64: std::convert::From<N>,
    I: SpIndex,
{
    let size_factors = index_to_cluster_map.iter().counts();

    let mut cluster_reps = Array2::<f64>::zeros((feature_bc_matrix.rows(), size_factors.len()).f());
    for (index, bc_expressions) in feature_bc_matrix.outer_iterator().enumerate() {
        let cluster_number = index_to_cluster_map[index];
        let cluster_index = (cluster_number - 1) as usize;
        cluster_reps.column_mut(cluster_index).add_assign(
            &(bc_expressions.to_dense().mapv(|x| f64::from(x))
                / size_factors[&cluster_number] as f64),
        );
    }
    cluster_reps
}

#[cfg(test)]
mod tests {
    use super::*;
    use ndarray::array;
    use sprs::CsMat;

    #[test]
    fn test_get_cluster_reps() {
        // # Python code to reconstruct this test
        // import numpy as np
        // from scipy.cluster import hierarchy
        // arr = np.array([
        //             [4, 5, 10, 4, 3, 11, 14, 6, 10, 12],
        //             [21, 19, 24, 17, 16, 25, 24, 22, 21, 21],
        //             [13, 10, 42, 7, 1, 17, 14, 20, 11, 9]
        //         ])
        // label_list = np.array([3, 2, 1, 3, 3, 3, 1, 3, 4, 4])
        // unique_labels = sorted(set(label_list))
        // cluster_reps = np.zeros((arr.shape[0], len(unique_labels)))
        // for ind, clstr_name in enumerate(unique_labels):
        //     cluster_reps[:, ind] = arr[:,label_list==clstr_name].mean(axis=1)
        // print(cluster_reps)
        // >> array([[12. ,  5. ,  5.6, 11. ],
        //          [24. , 19. , 20.2, 21. ],
        //          [28. , 10. , 11.6, 10. ]])
        //
        // To get sparse reps for the input array note that
        // import scipy.sparse as sp_sparse
        // sp_arr = sp_sparse.csc_matrix(arr)
        // sp_arr.indptr
        // >> array([ 0,  3,  6,  9, 12, 15, 18, 21, 24, 27, 30], dtype=int32)
        // sp_arr.indices
        // >> array([0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0,
        //              1, 2, 0, 1, 2, 0, 1, 2], dtype=int32)
        // sp_arr.data
        // >> array([ 4, 21, 13,  5, 19, 10, 10, 24, 42,  4, 17,  7,  3, 16,  1, 11, 25,
        //              17, 14, 24, 14,  6, 22, 20, 10, 21, 11, 12, 21,  9])

        let mat = CsMat::new_csc(
            (3, 10),
            vec![0, 3, 6, 9, 12, 15, 18, 21, 24, 27, 30],
            vec![
                0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0,
                1, 2,
            ],
            vec![
                4., 21., 13., 5., 19., 10., 10., 24., 42., 4., 17., 7., 3., 16., 1., 11., 25., 17.,
                14., 24., 14., 6., 22., 20., 10., 21., 11., 12., 21., 9.,
            ],
        );
        let label_list = vec![3, 2, 1, 3, 3, 3, 1, 3, 4, 4];
        let expected_cluster_reps_out = array![
            [12., 5., 5.6, 11.],
            [24., 19., 20.2, 21.],
            [28., 10., 11.6, 10.]
        ];

        let cluster_reps = get_cluster_representatives(&label_list, &mat);
        assert!(expected_cluster_reps_out.abs_diff_eq(&cluster_reps, 1e-6));
    }
}

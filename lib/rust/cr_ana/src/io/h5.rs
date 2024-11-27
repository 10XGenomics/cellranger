//! I/O HDF5 helper functions

use crate::types::{
    clustering_key, ClusteringKey, ClusteringResult, ClusteringType, EmbeddingResult,
    EmbeddingType, H5File, PcaResult,
};
use anyhow::{Context, Result};
use cr_types::reference::feature_reference::FeatureType;
use itertools::FoldWhile::{Continue, Done};
use itertools::Itertools;
use ndarray::{s, Array2};
use std::collections::BTreeMap;

const VERSION_DS: &str = "version";
const VERSION: i64 = 2;

pub(crate) type FA = hdf5::types::FixedAscii<256>;
pub(crate) type FA32 = hdf5::types::FixedAscii<32>;
pub(crate) type FA16 = hdf5::types::FixedAscii<16>;

pub(crate) fn make_fixed_ascii(s: &str) -> Result<FA> {
    Ok(FA::from_ascii(s)?)
}

pub(crate) fn make_fixed_ascii32(s: &str) -> Result<FA32> {
    Ok(FA32::from_ascii(s)?)
}

pub(crate) mod matrix {
    pub(crate) const GROUP: &str = "matrix";
    pub(crate) const SHAPE: &str = "shape";
    pub(crate) const FEATURE_GROUP: &str = "features";
    pub(crate) const FEATURE_TYPE: &str = "feature_type";
}

/// Return (num_features, num_barcodes, nnz).
pub(crate) fn matrix_shape(matrix: &H5File) -> Result<(isize, isize, isize)> {
    let matrix = hdf5::File::open(matrix)?.group(matrix::GROUP)?;
    let shape = matrix.dataset(matrix::SHAPE)?.read_1d::<isize>()?;
    let nnz = matrix.dataset("data")?.size() as isize;
    Ok((shape[0], shape[1], nnz))
}

pub(crate) fn matrix_feature_types(matrix: &H5File) -> Result<BTreeMap<FeatureType, usize>> {
    let dataset = hdf5::File::open(matrix)?
        .group(matrix::GROUP)?
        .group(matrix::FEATURE_GROUP)?
        .dataset(matrix::FEATURE_TYPE)?;
    let feature_types = dataset
        .read_1d::<FA>()?
        .iter()
        .map(|x| x.parse::<FeatureType>())
        .fold_while(Ok(BTreeMap::default()), |mut acc, x| match x {
            Ok(x) => {
                *acc.as_mut().unwrap().entry(x).or_insert(0) += 1;
                Continue(acc)
            }
            Err(e) => Done(Err(e)),
        })
        .into_inner()?;
    Ok(feature_types)
}

/// Estimate the memory (GiB) required to load a matrix from the number of non-zero entries (NNZ).
pub(crate) fn estimate_mem_gib_from_nnz(matrix: &H5File) -> Result<f64> {
    let (_num_features, _num_cells, nnz) = matrix_shape(matrix)?;
    Ok(20e-9 * nnz as f64)
}

/// Estimate the memory (GiB) required to store the PCA outputs.
pub(crate) fn estimate_mem_gib_for_pca_outs(matrix: &H5File, num_pcs: usize) -> Result<f64> {
    let (num_features, num_cells, _nnz) = matrix_shape(matrix)?;
    Ok(20e-9 * ((num_features + num_cells) * num_pcs as isize) as f64)
}

pub(crate) mod pca {
    pub(crate) const GROUP: &str = "pca";
    pub(crate) const COMPONENTS: &str = "components";
    pub(crate) const DISPERSION: &str = "dispersion";
    pub(crate) const FEATURES: &str = "features_selected";
    pub(crate) const MATRIX: &str = "transformed_pca_matrix";
    pub(crate) const VARIANCE: &str = "variance_explained";
}

fn feature_group(group: hdf5::Group, feature_type: FeatureType) -> Result<hdf5::Group> {
    let member = group
        .member_names()?
        .into_iter()
        .find(|m| m.contains(feature_type.as_snake_case()));
    if member.is_none() {
        anyhow::bail!(
            "Could not find {:?} feature in group {:?}",
            feature_type,
            group
        )
    }
    Ok(group.group(&member.unwrap())?)
}

pub(crate) struct PcaMatrixShape {
    pub num_bcs: usize,
    pub num_components: usize,
}

pub(crate) fn load_transformed_pca(
    pca_h5: &H5File,
    feature_type: FeatureType,
    input_pcs: Option<usize>,
) -> Result<(hdf5::Group, PcaMatrixShape)> {
    let file = hdf5::File::open(pca_h5)?;
    let group = file.group(pca::GROUP)?;
    let n_components = feature_group(group, feature_type)?;
    let sh = n_components.dataset(pca::MATRIX)?.shape();
    let mut shape = PcaMatrixShape {
        num_bcs: sh[0],
        num_components: sh[1],
    };
    if let Some(input_pcs) = input_pcs {
        if shape.num_components > input_pcs {
            shape.num_components = input_pcs;
        }
    }
    Ok((n_components, shape))
}

pub(crate) fn load_transformed_pca_matrix(
    pca_h5: &H5File,
    feature_type: FeatureType,
    input_pcs: Option<usize>,
) -> Result<Array2<f64>> {
    let (n_components, shape) = load_transformed_pca(pca_h5, feature_type, input_pcs)?;
    let proj = if input_pcs.is_some() {
        // Get the minimum of the number of available and requested components.
        n_components
            .dataset(pca::MATRIX)?
            .read_slice_2d::<f64, _>(s![.., 0..shape.num_components])?
    } else {
        n_components.dataset(pca::MATRIX)?.read_2d::<f64>()?
    };
    Ok(proj)
}

pub(crate) fn save_pca(pca_h5: &H5File, result: &PcaResult<'_>) -> Result<()> {
    let PcaResult {
        components,
        dispersion,
        features_selected,
        transformed_pca_matrix,
        variance_explained,
        ..
    } = result;
    let file = hdf5::File::create(pca_h5)?;
    file.new_dataset::<i64>()
        .create(VERSION_DS)?
        .write_scalar(&VERSION)?;
    let group = file
        .create_group(pca::GROUP)?
        .create_group(&result.h5_key())?;
    group
        .new_dataset::<f64>()
        .shape(components.dim())
        .create(pca::COMPONENTS)?
        .write(components)?;
    group
        .new_dataset::<f64>()
        .shape(dispersion.dim())
        .create(pca::DISPERSION)?
        .write(dispersion)?;
    let features_selected: Vec<_> = features_selected
        .iter()
        .map(|f| make_fixed_ascii(f))
        .try_collect()?;
    group
        .new_dataset::<FA>()
        .shape((features_selected.len(),))
        .create(pca::FEATURES)?
        .write(&features_selected)?;
    group
        .new_dataset::<f64>()
        .shape(transformed_pca_matrix.dim())
        .create(pca::MATRIX)?
        .write(transformed_pca_matrix)?;
    group
        .new_dataset::<f64>()
        .shape(variance_explained.dim())
        .create(pca::VARIANCE)?
        .write(variance_explained)?;
    Ok(())
}

pub(crate) fn combine_pcas<'a>(
    path: &H5File,
    parts: impl Iterator<Item = &'a H5File>,
) -> Result<()> {
    let output = hdf5::File::create(path)?;
    output
        .new_dataset::<i64>()
        .create(VERSION_DS)?
        .write_scalar(&VERSION)?;
    let out_group = output.create_group(pca::GROUP)?;
    for part in parts {
        let in_group = hdf5::File::open(part)?.group(pca::GROUP)?;
        for num_pcs in in_group.member_names()? {
            let in_group = in_group.group(&num_pcs)?;
            let components_staging_area = in_group.dataset(pca::COMPONENTS)?;
            let mut components = ndarray::Array2::<f64>::zeros((1, 1));
            let mut dispersion = ndarray::Array1::<f64>::zeros(1);
            let mut features_selected = ndarray::Array1::<FA>::from_elem(1, FA::new());
            let mut variance_explained = ndarray::Array1::<f64>::zeros(1);
            // if components is 2d, pca information is present above transformed_pca_matrix
            if components_staging_area.shape()[0] == 2 {
                components = components_staging_area.read_2d::<f64>()?;
                dispersion = in_group.dataset(pca::DISPERSION)?.read_1d::<f64>()?;
                features_selected = in_group.dataset(pca::FEATURES)?.read_1d::<FA>()?;
                variance_explained = in_group.dataset(pca::VARIANCE)?.read_1d::<f64>()?;
            }
            let transformed_pca_matrix = in_group.dataset(pca::MATRIX)?.read_2d::<f64>()?;
            let out_group = out_group.create_group(&num_pcs)?;
            out_group
                .new_dataset::<f64>()
                .shape(components.dim())
                .create(pca::COMPONENTS)?
                .write(&components)?;
            out_group
                .new_dataset::<f64>()
                .shape(dispersion.dim())
                .create(pca::DISPERSION)?
                .write(&dispersion)?;
            out_group
                .new_dataset::<FA>()
                .shape((features_selected.len(),))
                .create(pca::FEATURES)?
                .write(&features_selected)?;
            out_group
                .new_dataset::<f64>()
                .shape(transformed_pca_matrix.dim())
                .create(pca::MATRIX)?
                .write(&transformed_pca_matrix)?;
            out_group
                .new_dataset::<f64>()
                .shape(variance_explained.dim())
                .create(pca::VARIANCE)?
                .write(&variance_explained)?;
        }
    }
    Ok(())
}

pub(crate) mod clustering {
    pub(crate) const GROUP: &str = "clustering";
    pub(crate) const SCORE: &str = "cluster_score";
    pub(crate) const TYPE: &str = "clustering_type";
    pub(crate) const CLUSTERS: &str = "clusters";
    pub(crate) const DESCRIPTION: &str = "description";
    pub(crate) const SORT_KEY: &str = "global_sort_key";
    pub(crate) const NUM: &str = "num_clusters";
}

pub(crate) fn load_clustering(
    path: &H5File,
    clustering_type: ClusteringType,
    feature_type: FeatureType,
) -> Result<ClusteringResult> {
    let file = hdf5::File::open(path)?;
    let group = file.group(clustering::GROUP)?;

    let key = clustering_key(clustering_type, feature_type);
    let query_key = format!("_{key}");
    let labels = group
        .group(&query_key)?
        .dataset(clustering::CLUSTERS)?
        .read_1d::<i64>()?
        .to_vec();
    Ok(ClusteringResult {
        clustering_type,
        feature_type,
        labels,
        key,
    })
}

pub(crate) fn load_clusterings(
    path: &H5File,
    is_antibody_only: bool,
) -> Result<Vec<ClusteringResult>> {
    let file = hdf5::File::open(path)?;
    let group = file.group(clustering::GROUP)?;
    let members = group.member_names()?;
    let mut result = Vec::with_capacity(members.len());
    for member in members {
        let ClusteringKey {
            clustering_type,
            feature_type,
        } = member
            .parse::<ClusteringKey>()
            .with_context(|| format!("failed to parse '{member}' as clustering key"))?;

        let feature_type = if is_antibody_only {
            FeatureType::Barcode(cr_types::FeatureBarcodeType::Antibody)
        } else {
            feature_type
        };
        let key = clustering_key(clustering_type, feature_type);
        let labels = group
            .group(&member)?
            .dataset(clustering::CLUSTERS)?
            .read_1d::<i64>()?
            .to_vec();
        result.push(ClusteringResult {
            clustering_type,
            feature_type,
            labels,
            key,
        });
    }
    Ok(result)
}

pub(crate) fn save_clustering(path: &H5File, result: &ClusteringResult) -> Result<()> {
    let ClusteringResult {
        labels,
        clustering_type,
        ..
    } = result;
    let file = hdf5::File::create(path)?;
    file.new_dataset::<i64>()
        .create(VERSION_DS)?
        .write_scalar(&VERSION)?;
    let clustering_key = &result.key;
    let group = file
        .create_group(clustering::GROUP)?
        .create_group(&format!("_{clustering_key}"))
        .with_context(|| format!("unable to create group : {clustering_key}"))?;
    group
        .new_dataset::<f64>()
        .create(clustering::SCORE)?
        .write_scalar(&0.0)?;
    group
        .new_dataset::<FA32>()
        .create(clustering::TYPE)?
        .write_scalar(&make_fixed_ascii32(&clustering_type.lc()).with_context(|| {
            format!(
                "unable to create convert dataset name to fixed string size : {}",
                &clustering_type.lc()
            )
        })?)?;
    group
        .new_dataset::<i64>()
        .shape((labels.len(),))
        .create(clustering::CLUSTERS)?
        .write(&labels[..])?;
    group
        .new_dataset::<FA32>()
        .create(clustering::DESCRIPTION)?
        .write_scalar(&make_fixed_ascii32(&result.desc()).with_context(|| {
            format!(
                "unable to create convert description to fixed string size : {}",
                &result.desc()
            )
        })?)?;
    group
        .new_dataset::<f64>()
        .create(clustering::SORT_KEY)?
        .write_scalar(&f64::NEG_INFINITY)?;
    group
        .new_dataset::<i64>()
        .create(clustering::NUM)?
        .write_scalar(labels.iter().max().unwrap_or(&0))?;
    Ok(())
}

pub(crate) fn combine_clusterings<'a>(
    path: &H5File,
    parts: impl IntoIterator<Item = &'a H5File>,
) -> Result<()> {
    let output = hdf5::File::create(path)?;
    output
        .new_dataset::<i64>()
        .create(VERSION_DS)?
        .write_scalar(&VERSION)?;
    let out_group = output.create_group(clustering::GROUP)?;
    for part in parts {
        let in_group = hdf5::File::open(part)?.group(clustering::GROUP)?;
        for clustering_key in in_group.member_names()? {
            let in_group = in_group.group(&clustering_key)?;
            let cluster_score = in_group.dataset(clustering::SCORE)?.read_scalar::<f64>()?;
            let clustering_type = in_group.dataset(clustering::TYPE)?.read_scalar::<FA16>()?;
            let clusters = in_group.dataset(clustering::CLUSTERS)?.read_1d::<i64>()?;
            let description = in_group
                .dataset(clustering::DESCRIPTION)?
                .read_scalar::<FA32>()?;
            let global_sort_key = in_group
                .dataset(clustering::SORT_KEY)?
                .read_scalar::<f64>()?;
            let num_clusters = in_group.dataset(clustering::NUM)?.read_scalar::<i64>()?;
            let out_group = out_group
                .create_group(&clustering_key)
                .with_context(|| format!("could not create group {clustering_key}"))?;
            out_group
                .new_dataset::<f64>()
                .create(clustering::SCORE)?
                .write_scalar(&cluster_score)?;
            out_group
                .new_dataset::<FA16>()
                .create(clustering::TYPE)?
                .write_scalar(&clustering_type)?;
            out_group
                .new_dataset::<i64>()
                .shape((clusters.len(),))
                .create(clustering::CLUSTERS)?
                .write(&clusters)?;
            out_group
                .new_dataset::<FA32>()
                .create(clustering::DESCRIPTION)?
                .write_scalar(&description)?;
            out_group
                .new_dataset::<f64>()
                .create(clustering::SORT_KEY)?
                .write_scalar(&global_sort_key)?;
            out_group
                .new_dataset::<i64>()
                .create(clustering::NUM)?
                .write_scalar(&num_clusters)?;
        }
    }
    Ok(())
}

pub(crate) mod embedding {
    pub(crate) const KEY: &str = "key";
    pub(crate) const NAME: &str = "name";
}

pub(crate) fn save_embedding(path: &H5File, result: &EmbeddingResult<'_>) -> Result<()> {
    let EmbeddingResult {
        embedding,
        embedding_type,
        feature_type,
        key,
        ..
    } = result;
    let dims = embedding.dim().1;
    let file = hdf5::File::create(path)?;
    let group = file
        .create_group(embedding_type.lc())?
        .create_group(&format!("_{key}"))?; // prefixed with a _ for some reason..
    group
        .new_dataset::<f64>()
        .shape(embedding.dim())
        .create(format!("transformed_{}_matrix", embedding_type.lc()).as_ref())?
        .write(embedding)?;
    group
        .new_dataset::<FA32>()
        .create(embedding::KEY)?
        .write_scalar(&make_fixed_ascii32(key)?)?;
    group
        .new_dataset::<FA>()
        .create(embedding::NAME)?
        .write_scalar(&make_fixed_ascii(&format!(
            "{}_{}-d",
            feature_type.as_snake_case(),
            dims,
        ))?)?;
    Ok(())
}

pub(crate) fn combine_embeddings<'a>(
    path: &H5File,
    parts: impl Iterator<Item = &'a H5File>,
    embedding_type: EmbeddingType,
) -> Result<()> {
    let h5 = hdf5::File::create(path)?.create_group(embedding_type.lc())?;
    for part in parts {
        let group = hdf5::File::open(part)?.group(embedding_type.lc())?;
        for embedding_key in group.member_names()? {
            let group = group.group(&embedding_key)?;
            let matrix = group
                .dataset(&format!("transformed_{}_matrix", embedding_type.lc()))?
                .read_2d::<f64>()?;
            let key = group.dataset(embedding::KEY)?.read_scalar::<FA32>()?;
            let name = group.dataset(embedding::NAME)?.read_scalar::<FA>()?;
            let group = h5.create_group(&embedding_key)?;
            group
                .new_dataset::<f64>()
                .shape(matrix.dim())
                .create(format!("transformed_{}_matrix", embedding_type.lc()).as_ref())?
                .write(&matrix)?;
            group
                .new_dataset::<FA32>()
                .create(embedding::KEY)?
                .write_scalar(&key)?;
            group
                .new_dataset::<FA>()
                .create(embedding::NAME)?
                .write_scalar(&name)?;
        }
    }
    Ok(())
}

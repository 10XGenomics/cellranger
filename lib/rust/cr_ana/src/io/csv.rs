//! I/O CSV helper functions

use crate::types::{clustering_key, ClusteringResult, ClusteringType, EmbeddingResult, PcaResult};
use anyhow::{anyhow, Result};
use cr_types::reference::feature_reference::FeatureType;
use ndarray::Array2;
use std::fs::{copy, create_dir_all, File};
use std::io::{BufWriter, Write};
use std::path::Path;

pub(crate) fn save_pca(
    pca_csv: &Path,
    result: &PcaResult<'_>,
    barcodes: &[String],
    feature_ids: &[String],
) -> Result<()> {
    let PcaResult {
        components,
        dispersion,
        features_selected,
        transformed_pca_matrix,
        variance_explained,
        ..
    } = result;
    let (num_bcs, num_pcs) = transformed_pca_matrix.dim();
    let subdir = result.csv_key();
    let component_dir = pca_csv.join(subdir);
    create_dir_all(&component_dir)?;
    {
        let path = component_dir.join("components.csv");
        let mut file = BufWriter::new(File::create(path)?);
        write!(file, "PC")?;
        for feature_id in feature_ids {
            write!(file, ",{feature_id}")?;
        }
        writeln!(file)?;
        for i in 0..num_pcs {
            write!(file, "{}", i + 1)?;
            for j in 0..feature_ids.len() {
                write!(file, ",{}", components[[i, j]])?;
            }
            writeln!(file)?;
        }
    }
    {
        let path = component_dir.join("dispersion.csv");
        let mut file = BufWriter::new(File::create(path)?);
        writeln!(file, "Feature,Normalized.Dispersion")?;
        for (i, feature_id) in feature_ids.iter().enumerate() {
            writeln!(file, "{feature_id},{}", dispersion[i])?;
        }
    }
    {
        // these are weirdly reversed b/c of how the python does it:
        //   np.argsort(...)[-pca_features:]
        let path = component_dir.join("features_selected.csv");
        let mut file = BufWriter::new(File::create(path)?);
        writeln!(file, "Feature")?;
        for (i, feature_id) in features_selected.iter().enumerate() {
            writeln!(file, "{},{feature_id}", i + 1)?;
        }
    }
    {
        let path = component_dir.join("projection.csv");
        let mut file = BufWriter::new(File::create(path)?);
        write!(file, "Barcode")?;
        for i in 1..=num_pcs {
            write!(file, ",PC-{i}")?;
        }
        writeln!(file)?;
        for i in 0..num_bcs {
            write!(file, "{}", barcodes[i])?;
            for j in 0..num_pcs {
                write!(file, ",{}", transformed_pca_matrix[[i, j]])?;
            }
            writeln!(file)?;
        }
    }
    {
        let path = component_dir.join("variance.csv");
        let mut file = BufWriter::new(File::create(path)?);
        writeln!(file, "PC,Proportion.Variance.Explained")?;
        for i in 0..num_pcs {
            writeln!(file, "{},{}", i + 1, variance_explained[i])?;
        }
    }
    Ok(())
}

pub(crate) fn combine_pcas(
    pca_csv: &Path,
    parts: impl Iterator<Item = impl AsRef<Path>>,
) -> Result<()> {
    for part in parts {
        let mut dirents = std::fs::read_dir(&part)?;
        let subdir = loop {
            match dirents.next() {
                Some(path) => {
                    let subdir = path?.file_name();
                    if subdir.to_str().is_some_and(|s| s.ends_with("_components")) {
                        break subdir;
                    }
                }
                None => {
                    panic!("failed to find any components!")
                }
            }
        };
        let in_component_dir = part.as_ref().join(&subdir);
        let out_component_dir = pca_csv.join(subdir);
        create_dir_all(&out_component_dir)?;
        for file in &[
            "components.csv",
            "dispersion.csv",
            "features_selected.csv",
            "projection.csv",
            "variance.csv",
        ] {
            if Path::exists(in_component_dir.join(file).as_ref()) {
                std::fs::copy(in_component_dir.join(file), out_component_dir.join(file))?;
            }
        }
    }
    Ok(())
}

pub(crate) fn save_clustering(
    path: &Path,
    result: &ClusteringResult,
    barcodes: &[String],
) -> Result<()> {
    let ClusteringResult { labels, .. } = result;
    let clustering_dir = path.join(&result.key);
    create_dir_all(&clustering_dir)?;
    {
        let path = clustering_dir.join("clusters.csv");
        let mut file = BufWriter::new(File::create(path)?);
        writeln!(file, "Barcode,Cluster")?;
        for i in 0..labels.len() {
            writeln!(file, "{},{}", barcodes[i], labels[i])?;
        }
    }
    Ok(())
}

pub(crate) fn combine_clusterings(
    path: &Path,
    parts: impl IntoIterator<Item = impl AsRef<Path>>,
) -> Result<()> {
    for part in parts {
        for dirent in std::fs::read_dir(&part)? {
            let subdir = dirent?.file_name();
            create_dir_all(path.join(&subdir))?;
            let src = part.as_ref().join(&subdir).join("clusters.csv");
            if src.exists() {
                let dst = path.join(subdir).join("clusters.csv");
                std::fs::copy(src, dst)?;
            }
        }
    }
    Ok(())
}

pub(crate) fn save_differential_expression(
    path: &Path,
    clustering_type: ClusteringType,
    feature_type: FeatureType,
    feature_ids: &[String],
    feature_names: &[String],
    data: &Array2<f64>,
) -> Result<()> {
    let key = clustering_key(clustering_type, feature_type);
    let clustering_dir = path.join(key);
    create_dir_all(&clustering_dir)?;
    {
        let (num_bcs, num_clusters) = data.dim();
        let num_clusters = num_clusters / 3;
        let path = clustering_dir.join("differential_expression.csv");
        let mut file = BufWriter::new(File::create(path)?);
        write!(file, "Feature ID,Feature Name")?;
        for j in 1..=num_clusters {
            write!(
                file,
                ",Cluster {j} Mean Counts,Cluster {j} Log2 fold change,Cluster {j} Adjusted p value"
            )?;
        }
        writeln!(file)?;
        for i in 0..num_bcs {
            write!(file, "{},{}", feature_ids[i], feature_names[i])?;
            for j in 0..num_clusters {
                write!(
                    file,
                    ",{},{},{}",
                    data[[i, 3 * j]],
                    data[[i, 3 * j + 1]],
                    data[[i, 3 * j + 2]]
                )?;
            }
            writeln!(file)?;
        }
    }
    Ok(())
}

pub(crate) fn save_embedding(path: &Path, result: &EmbeddingResult<'_>) -> Result<()> {
    let EmbeddingResult {
        barcodes,
        embedding,
        embedding_type,
        key,
        ..
    } = result;
    let (num_bcs, num_components) = embedding.dim();
    let key = format!("{key}_components");
    let component_dir = path.join(key);
    create_dir_all(&component_dir)?;
    {
        let path = component_dir.join("projection.csv");
        let mut file = BufWriter::new(File::create(path)?);
        write!(file, "Barcode")?;
        for i in 1..=num_components {
            write!(file, ",{}-{i}", embedding_type.uc())?;
        }
        writeln!(file)?;
        for i in 0..num_bcs {
            write!(file, "{}", barcodes[i])?;
            for j in 0..num_components {
                write!(file, ",{}", embedding[[i, j]])?;
            }
            writeln!(file)?;
        }
    }
    Ok(())
}

pub(crate) fn combine_embeddings(
    path: &Path,
    parts: impl Iterator<Item = impl AsRef<Path>>,
) -> Result<()> {
    for part in parts {
        let part = part.as_ref();
        let key = part
            .read_dir()?
            .filter_map(|d| d.ok().and_then(|d| d.file_name().into_string().ok()))
            .find(|d| d.ends_with("_components"))
            .ok_or_else(|| anyhow!("unable to find projection for {}", part.display()))?;
        let src = part.join(&key).join("projection.csv");
        let dst = path.join(&key).join("projection.csv");
        create_dir_all(path.join(&key))?;
        copy(src, dst)?;
    }
    Ok(())
}

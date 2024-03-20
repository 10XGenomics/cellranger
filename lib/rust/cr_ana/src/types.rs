use anyhow::{anyhow, bail, Result};
use cr_types::reference::feature_reference::FeatureType;
use cr_types::FeatureBarcodeType;
use martian_derive::{martian_filetype, MartianStruct};
use ndarray::{Array1, Array2};
use serde::{Deserialize, Serialize};
use std::borrow::Cow;
use std::convert::TryFrom;
use std::str::FromStr;

martian_filetype!(H5File, "h5");

#[derive(Debug, Serialize, Deserialize, MartianStruct)]
pub struct Unit {}

pub(crate) struct PcaResult<'a> {
    pub components: Array2<f64>,
    pub dispersion: Array1<f64>,
    pub features_selected: Vec<&'a str>,
    pub transformed_pca_matrix: Array2<f64>,
    pub variance_explained: Array1<f64>,
    pub key: String,
}

impl<'a> PcaResult<'a> {
    pub(crate) fn new(
        components: Array2<f64>,
        dispersion: Array1<f64>,
        feature_type: FeatureType,
        features_selected: Vec<&'a str>,
        transformed_pca_matrix: Array2<f64>,
        variance_explained: Array1<f64>,
    ) -> Self {
        let num_pcs = transformed_pca_matrix.dim().1;
        let key = format!("{}_{num_pcs}", feature_type.as_snake_case());
        PcaResult {
            components,
            dispersion,
            features_selected,
            transformed_pca_matrix,
            variance_explained,
            key,
        }
    }
    pub(crate) fn csv_key(&self) -> String {
        format!("{}_components", self.key)
    }
    pub(crate) fn h5_key(&self) -> String {
        format!("_{}", self.key)
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord, Serialize, Deserialize)]
#[serde(try_from = "&str", into = "Cow<'static, str>")]
pub enum ClusteringType {
    KMeans(usize),
    Louvain,
    Hierarchical(usize),
}

impl ClusteringType {
    pub(crate) fn lc(&self) -> Cow<'static, str> {
        use ClusteringType::{Hierarchical, KMeans, Louvain};
        match self {
            KMeans(k) => Cow::Owned(format!("kmeans_{k}_clusters")),
            Louvain => Cow::Borrowed("graphclust"),
            Hierarchical(num_clusters) => Cow::Owned(format!("hcluster_{num_clusters}_clusters")),
        }
    }
    pub(crate) fn desc(&self) -> Cow<'static, str> {
        use ClusteringType::{Hierarchical, KMeans, Louvain};
        match self {
            KMeans(k) => Cow::Owned(format!("K-means (K={k})")),
            Louvain => Cow::Borrowed("Graph-based"),
            Hierarchical(num_clusters) => Cow::Owned(format!("Hclust (K={num_clusters})")),
        }
    }
}

impl FromStr for ClusteringType {
    type Err = anyhow::Error;

    fn from_str(s: &str) -> Result<Self> {
        Ok(match s.strip_prefix('_').unwrap_or(s) {
            "graphclust" => ClusteringType::Louvain,
            t if t.starts_with("kmeans") => {
                let k = t
                    .split('_')
                    .nth(1)
                    .ok_or_else(|| anyhow!("bad kmeans clustering type: {s}"))
                    .and_then(|k| Ok(k.parse::<usize>()?))?;
                ClusteringType::KMeans(k)
            }
            t if t.starts_with("hcluster") => {
                let k = t
                    .split('_')
                    .rev()
                    .nth(1)
                    .ok_or_else(|| anyhow!("bad hierarchical clustering type: {s}"))
                    .and_then(|k| Ok(k.parse::<usize>()?))?;
                ClusteringType::Hierarchical(k)
            }
            _ => bail!("unknown clustering type: {s}"),
        })
    }
}

impl<'a> TryFrom<&'a str> for ClusteringType {
    type Error = anyhow::Error;

    fn try_from(s: &str) -> Result<Self> {
        s.parse::<ClusteringType>()
    }
}

impl From<ClusteringType> for Cow<'static, str> {
    fn from(val: ClusteringType) -> Cow<'static, str> {
        val.lc()
    }
}

pub(crate) struct ClusteringKey {
    pub clustering_type: ClusteringType,
    pub feature_type: FeatureType,
}

impl FromStr for ClusteringKey {
    type Err = anyhow::Error;

    fn from_str(s: &str) -> Result<Self> {
        let t = s.strip_prefix('_').unwrap_or(s);
        let parts = t.split('_').collect::<Vec<_>>();
        let clustering_key = match parts
            .iter()
            .position(|&s| s == "graphclust" || s == "kmeans" || s == "hcluster")
        {
            None => {
                bail!("invalid clustering key: {s}");
            }
            Some(0) => ClusteringKey {
                clustering_type: t.parse::<ClusteringType>()?,
                feature_type: FeatureType::Gene,
            },
            Some(i) => {
                let feature_type = FeatureType::from_snake_case(&parts[0..i].join("_"))?;
                let clustering_type = parts[i..].join("_").parse::<ClusteringType>()?;
                ClusteringKey {
                    clustering_type,
                    feature_type,
                }
            }
        };
        Ok(clustering_key)
    }
}

pub(crate) struct ClusteringResult {
    pub clustering_type: ClusteringType,
    pub feature_type: FeatureType,
    pub labels: Vec<i64>,
    pub key: String,
}

pub(crate) fn clustering_key(clustering_type: ClusteringType, feature_type: FeatureType) -> String {
    format!("{}_{}", feature_type.as_snake_case(), clustering_type.lc())
}

impl ClusteringResult {
    pub(crate) fn new(
        clustering_type: ClusteringType,
        feature_type: FeatureType,
        labels: Vec<i64>,
    ) -> Self {
        let key = clustering_key(clustering_type, feature_type);
        ClusteringResult {
            clustering_type,
            feature_type,
            labels,
            key,
        }
    }

    pub(crate) fn num_clusters(&self) -> i64 {
        self.labels.iter().copied().max().unwrap_or(0)
    }
    pub(crate) fn desc(&self) -> Cow<'static, str> {
        match self.feature_type {
            FeatureType::Gene => self.clustering_type.desc(),
            FeatureType::Barcode(FeatureBarcodeType::Antibody) => {
                Cow::Owned(format!("Antibody {}", self.clustering_type.desc()))
            }
            FeatureType::Barcode(_) => unimplemented!(),
        }
    }
}

pub(crate) enum EmbeddingType {
    Tsne,
    Umap,
}

impl EmbeddingType {
    pub(crate) fn lc(&self) -> &'static str {
        use EmbeddingType::{Tsne, Umap};
        match self {
            Tsne => "tsne",
            Umap => "umap",
        }
    }
    pub(crate) fn uc(&self) -> &'static str {
        use EmbeddingType::{Tsne, Umap};
        match self {
            Tsne => "TSNE",
            Umap => "UMAP",
        }
    }
}

pub(crate) struct EmbeddingResult<'a> {
    pub barcodes: Vec<Cow<'a, String>>,
    pub embedding: Array2<f64>,
    pub embedding_type: EmbeddingType,
    pub feature_type: FeatureType,
    pub key: String,
}

impl<'a> EmbeddingResult<'a> {
    pub(crate) fn new(
        barcodes: &'a [String],
        embedding: Array2<f64>,
        embedding_type: EmbeddingType,
        feature_type: FeatureType,
    ) -> Self {
        let dims = embedding.dim().1;
        let key = format!("{}_{dims}", feature_type.as_snake_case());
        EmbeddingResult {
            barcodes: barcodes.iter().map(Cow::Borrowed).collect::<Vec<_>>(),
            embedding,
            embedding_type,
            feature_type,
            key,
        }
    }
}

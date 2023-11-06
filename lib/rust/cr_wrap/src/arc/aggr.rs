use crate::arc::types::{validate_distance, MAX_CLUSTERS_RANGE};
use crate::mrp_args::MrpArgs;
use crate::utils::{validate_id, CliPath};
use anyhow::{anyhow, bail, Context, Result};
use clap::{self, value_parser, Parser};
use csv::StringRecord;
use ordered_float::NotNan;
use serde::{self, Deserialize, Serialize};
use std::collections::HashMap;
use std::iter::zip;

#[derive(Parser, Debug, Clone)]
pub struct AggrArgs {
    /// A unique run id and output folder name [a-zA-Z0-9_-]+ of maximum length
    /// 64 characters
    #[clap(long, value_name = "ID", required = true, value_parser = validate_id)]
    pub id: String,

    /// Sample description to embed in output files
    #[clap(long = "description", default_value = "", value_name = "TEXT")]
    sample_desc: String,

    /// Path to folder containing cellranger-arc-compatible reference. Reference
    /// packages can be downloaded from support.10xgenomics.com or constructed using
    /// the `cellranger-arc mkref` command. Note: this reference must match the
    /// reference used for the initial `cellranger-arc count` run.
    #[clap(long, value_name = "PATH", required = true)]
    reference: CliPath,

    /// Path to CSV file enumerating 'cellranger-arc count' outputs required for aggregation.
    ///
    /// For example, a CSV for aggregating two samples would look as follows (blank lines are
    /// ignored):
    ///
    /// library_id,atac_fragments,per_barcode_metrics,gex_molecule_info
    ///
    /// L1,/data/L1/outs/atac_fragments.tsv.gz,/data/L1/outs/per_barcode_metrics.csv,/data/L1/outs/gex_molecule_info.h5
    ///
    /// L2,/data/L2/outs/atac_fragments.tsv.gz,/data/L2/outs/per_barcode_metrics.csv,/data/L2/outs/gex_molecule_info.h5
    ///
    /// Optionally, metadata associated with these libraries can
    /// be specified using additional columns. This information is not used by the pipeline but will
    /// be available in the Loupe file for visualization.
    #[clap(long, value_name = "CSV", required = true)]
    csv: CliPath,

    /// Override peak caller: specify peaks to use in downstream analyses from
    /// supplied 3-column BED file.
    /// The supplied peaks file must be sorted by position and not contain overlapping peaks;
    /// comment lines beginning with `#` are allowed.
    #[clap(long, value_name = "BED")]
    peaks: Option<CliPath>,

    /// Library depth normalization mode.
    #[clap(
        long = "normalize",
        default_value = "depth",
        value_name = "MODE",
        value_parser = ["none", "depth"],
    )]
    normalization: String,

    /// Disable secondary analysis, e.g. clustering
    #[clap(long = "nosecondary")]
    no_secondary_analysis: bool,

    /// Hidden option: change feature linkage window size
    #[clap(long, hide = true, value_parser = validate_distance)]
    feature_linkage_max_dist_mb: Option<NotNan<f64>>,

    /// Hidden option: change the max k of kmeans clustering
    #[clap(long, hide = true, value_parser = value_parser!(u64).range(MAX_CLUSTERS_RANGE))]
    k_means_max_clusters: Option<u64>,

    /// Do not execute the pipeline.
    /// Generate a pipeline invocation (.mro) file and stop.
    #[clap(long)]
    pub dry: bool,

    #[clap(flatten)]
    pub mrp: MrpArgs,
}

impl AggrArgs {
    pub fn to_mro_args(&self) -> Result<AggrMro> {
        let args = self.clone();

        Ok(AggrMro {
            sample_id: args.id,
            sample_desc: args.sample_desc,
            reference_path: args.reference,
            aggr_defs: AggrDefs::from_csv(&args.csv)?,
            no_secondary_analysis: args.no_secondary_analysis,
            normalization: args.normalization,
            k_means_max_clusters: args.k_means_max_clusters,
            feature_linkage_max_dist_mb: args.feature_linkage_max_dist_mb,
            custom_peaks: args.peaks,
        })
    }
}

#[derive(Debug, Serialize, Deserialize)]
#[serde(transparent)]
pub struct AggrDefs(Vec<AggrDef>);

impl AggrDefs {
    pub fn from_csv(path: &CliPath) -> Result<AggrDefs> {
        let mut aggr_defs = Vec::new();
        let mut reader = csv::ReaderBuilder::new()
            .has_headers(false)
            .from_path(path)
            .with_context(|| format!("Unable to read provided csv: {path:?}"))?;
        let mut header = StringRecord::new();
        let mut metadata_names = Vec::new();
        let mut metadata_indices = Vec::new();
        let non_metadata_cols = [
            "library_id",
            "atac_fragments",
            "gex_molecule_info",
            "per_barcode_metrics",
        ];
        for (ri, rec) in reader.records().enumerate() {
            let rec = rec.context("Unable to read aggr CSV")?;
            if ri == 0 {
                header = rec.clone();
                for (i, h) in rec.iter().enumerate() {
                    if non_metadata_cols.iter().any(|&x| x == h) {
                        continue;
                    }
                    metadata_indices.push(i);
                    metadata_names.push(h.to_string());
                }
                continue;
            }

            let aggr_def_base: AggrDefBase = rec
                .deserialize(Some(&header))
                .with_context(|| format!("Invalid aggr CSV row = `{}`", rec.as_slice()))?;
            let mut metadata = HashMap::new();
            for (name, i) in zip(&metadata_names, &metadata_indices) {
                let entry = rec.get(*i).ok_or_else(|| {
                    anyhow!(
                        "Unable to read metadata column `{name}` for CSV row `{}`",
                        rec.as_slice()
                    )
                })?;
                metadata.insert(name.to_string(), entry.to_string());
            }
            let row = AggrDef {
                aggr_def_base,
                metadata,
            };
            aggr_defs.push(row);
        }
        if aggr_defs.is_empty() {
            bail!("Unable to read provided csv: {:?}\nNo records found", path,);
        }
        for col in &metadata_names {
            if aggr_defs.iter().all(|adef| adef.metadata[col].is_empty()) {
                bail!(
                    "Error processing aggr CSV: metadata column `{}` has empty entries in every row.",
                    col
                );
            }
        }
        Ok(AggrDefs(aggr_defs))
    }
}

#[derive(Debug, Serialize, Deserialize, Clone)]
struct AggrDefBase {
    library_id: String,
    atac_fragments: CliPath,
    per_barcode_metrics: CliPath,
    gex_molecule_info: CliPath,
}

#[derive(Debug, Serialize, Deserialize)]
struct AggrDef {
    #[serde(flatten)]
    aggr_def_base: AggrDefBase,
    metadata: HashMap<String, String>,
}

#[derive(Serialize)]
pub struct AggrMro {
    sample_id: String,
    sample_desc: String,
    reference_path: CliPath,
    aggr_defs: AggrDefs,
    no_secondary_analysis: bool,
    normalization: String,
    k_means_max_clusters: Option<u64>,
    feature_linkage_max_dist_mb: Option<NotNan<f64>>,
    custom_peaks: Option<CliPath>,
}

#[cfg(test)]
pub mod tests {
    use super::*;
    use std::io::Write;
    use tempfile;

    #[test]
    pub fn test_metadata() {
        let mut input = tempfile::NamedTempFile::new().unwrap();
        input.write_all(b"\n").unwrap();
        input.flush().unwrap();
        let path = input.path().to_str().unwrap();

        // Verify that we correctly read metadata
        {
            let mut tmp = tempfile::NamedTempFile::new().unwrap();
            tmp.write_all(
                b"library_id,atac_fragments,per_barcode_metrics,gex_molecule_info,batch\n",
            )
            .unwrap();
            tmp.write_all(format!("L1,{path},{path},{path},1\n").as_bytes())
                .unwrap();
            tmp.write_all(format!("L2,{path},{path},{path},2\n").as_bytes())
                .unwrap();
            tmp.flush().unwrap();

            let aggr_defs = AggrDefs::from_csv(&CliPath::from(tmp.path())).unwrap().0;
            let batch = aggr_defs
                .iter()
                .map(|def| def.metadata["batch"].as_str())
                .collect::<Vec<_>>();
            assert_eq!(batch, ["1", "2"]);
        }

        // Exit if metadata is empty
        {
            let mut tmp = tempfile::NamedTempFile::new().unwrap();
            tmp.write_all(
                b"library_id,atac_fragments,per_barcode_metrics,gex_molecule_info,batch\n",
            )
            .unwrap();
            tmp.write_all(format!("L1,{path},{path},{path},\n").as_bytes())
                .unwrap();
            tmp.write_all(format!("L2,{path},{path},{path},\n").as_bytes())
                .unwrap();
            tmp.flush().unwrap();

            if let Err(err) = AggrDefs::from_csv(&CliPath::from(tmp.path())) {
                assert_eq!(
                    err.to_string(),
                    "Error processing aggr CSV: metadata column `batch` has empty entries in \
                    every row."
                        .to_string()
                );
            } else {
                panic!("empty metadata test failed");
            }
        }
    }

    #[test]
    pub fn test_empty_csv() {
        let mut tmp = tempfile::NamedTempFile::new().unwrap();
        tmp.flush().unwrap();
        assert!(AggrDefs::from_csv(&CliPath::from(tmp.path())).is_err());
    }
}

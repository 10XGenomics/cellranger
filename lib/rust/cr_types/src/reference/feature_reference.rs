#![expect(missing_docs)]
use super::probe_set_reference::TargetSetFile;
use crate::probe_set::{ProbeSetReferenceMetadata, is_deprecated_probe, strip_deprecated_prefix};
use crate::reference::reference_info::ReferenceInfo;
use crate::{FeatureBarcodeType, FeatureID, FeatureName, GenomeName, LibraryType};
use anyhow::{Context, Result, anyhow, bail, ensure};
use csv::{self, StringRecord};
use fastq_set::read_pair::WhichRead;
use itertools::{Itertools, chain};
use martian::{AsMartianPrimaryType, MartianPrimaryType};
use martian_derive::{MartianStruct, MartianType, martian_filetype};
use martian_filetypes::LazyFileTypeIO;
use martian_filetypes::tabular_file::CsvFileNoHeader;
use metric::{AsMetricPrefix, TxHashMap, TxHashSet};
use regex::Regex;
use serde::{Deserialize, Serialize};
use serde_with::DeserializeFromStr;
use std::collections::hash_map::Entry;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufReader, Read, Write};
use std::path::Path;
use std::str::{self, FromStr};
use std::string::String;
use strum::Display;
use strum::IntoEnumIterator;
use transcriptome::{Gene, Transcriptome};

pub const REQUIRED_FEATURE_TAGS: &[&str] = &[
    "id",
    "name",
    "genome",
    "feature_type",
    "read",
    "pattern",
    "sequence",
];

pub const REQUIRED_FEATURE_REF_COLS: &[&str] =
    &["id", "name", "read", "pattern", "sequence", "feature_type"];

pub const TARGETING_ANTIGEN: &str = "targeting_antigen";
pub const FUNCTIONAL_NAME: &str = "functional_name";

#[derive(PartialEq, Eq, Debug, Clone, Deserialize, Serialize, MartianStruct)]
pub struct SpecificityControls {
    pub control_for_allele: HashMap<String, String>,
    pub has_mhc_allele_column: bool,
}

// information related to beam extracted from multi config CSV
#[derive(PartialEq, Eq, Debug, Clone, Deserialize, Serialize, MartianStruct)]
pub struct FeatureConfig {
    pub beam_mode: Option<BeamMode>,
    pub specificity_controls: Option<SpecificityControls>,
    pub functional_map: Option<HashMap<String, String>>,
    pub hashtag_ids: Option<Vec<String>>,
}

pub const MHC_ALLELE: &str = "mhc_allele";
pub const NO_ALLELE: &str = "no_allele";
pub const HASHTAG: &str = "hashtag";

#[derive(PartialEq, Eq, Hash, Debug, Clone, Copy, Serialize, Deserialize, MartianType, Display)]
#[serde(rename_all = "lowercase")]
pub enum BeamMode {
    #[serde(rename = "beam_ab")]
    #[strum(to_string = "beam_ab")]
    BeamAB,
    #[serde(rename = "beam_t")]
    #[strum(to_string = "beam_t")]
    BeamT,
}

impl FromStr for BeamMode {
    type Err = anyhow::Error;

    fn from_str(s: &str) -> Result<BeamMode> {
        Ok(match s {
            "beam_ab" => BeamMode::BeamAB,
            "beam_t" => BeamMode::BeamT,
            _ => bail!(
                "Unknown variant '{s}' for beam mode. Supported variants are: [beam_ab, beam_t]"
            ),
        })
    }
}

pub type TargetGeneIndicesFile = CsvFileNoHeader<u32>;

/// A target set, which is composed of a name and a vector of on-target feature indices.
#[derive(PartialEq, Eq, Serialize, Deserialize, Clone, Debug)]
pub struct TargetSet {
    /// The name of the feature reference.
    name: String,

    /// A set of the on-target feature indices.
    feature_indices: TxHashSet<u32>,
}

impl TargetSet {
    /// Create a TargetSet from a name and a set of on-target feature indices.
    pub fn from_indices(name: &str, feature_indices: TxHashSet<u32>) -> TargetSet {
        TargetSet {
            name: name.to_string(),
            feature_indices,
        }
    }

    /// Create a TargetSet from a name and a slice of on-target bools.
    pub fn from_bools(name: &str, on_target: &[bool]) -> TargetSet {
        TargetSet::from_indices(
            name,
            on_target
                .iter()
                .enumerate()
                .filter_map(|(i, boolean)| boolean.then_some(i as u32))
                .collect(),
        )
    }

    /// Create a TargetSet from a name and a feature reference path.
    pub fn load(name: &str, path: &TargetGeneIndicesFile) -> Result<TargetSet> {
        Ok(TargetSet::from_indices(name, path.read_all()?))
    }

    // Return the name of the feature reference.
    pub fn name(&self) -> &str {
        &self.name
    }

    // Return a set of on-target feature indices.
    pub fn feature_indices(&self) -> &TxHashSet<u32> {
        &self.feature_indices
    }

    /// Return whether the specified feature index is on target.
    pub fn is_on_target(&self, feature_index: u32) -> bool {
        self.feature_indices.contains(&feature_index)
    }

    // Return a sorted vector of on-target feature indices.
    pub fn to_feature_indices_vec(&self) -> Vec<u32> {
        self.feature_indices.iter().copied().sorted().collect()
    }
}

#[derive(
    Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord, Serialize, DeserializeFromStr,
)]
#[serde(into = "&str")]
pub enum FeatureType {
    Gene,
    Barcode(FeatureBarcodeType),
    ProteinExpression, // used by Xenium for in situ protein expression
    Peaks,
}

impl AsMartianPrimaryType for FeatureType {
    fn as_martian_primary_type() -> MartianPrimaryType {
        MartianPrimaryType::Str
    }
}

impl FeatureType {
    pub fn is_multiplexing(&self) -> bool {
        *self == Self::Barcode(FeatureBarcodeType::Multiplexing)
    }

    const GENE_EXP: &'static str = "Gene Expression";
    const PROTEIN_EXP: &'static str = "Protein Expression";
    const PEAKS: &'static str = "Peaks";

    /// Return the string representation of this feature type.
    pub fn as_str(&self) -> &'static str {
        match self {
            Self::Gene => Self::GENE_EXP,
            Self::Barcode(x) => x.as_str(),
            Self::ProteinExpression => Self::PROTEIN_EXP,
            Self::Peaks => Self::PEAKS,
        }
    }

    const GENE_EXP_SNAKE_CASE: &'static str = "gene_expression";
    const PROTEIN_EXP_SNAKE_CASE: &'static str = "protein_expression";
    const PEAKS_SNAKE_CASE: &'static str = "peaks";

    /// Return a snake_case representation of this feature type.
    pub fn as_snake_case(&self) -> &'static str {
        match self {
            Self::Gene => Self::GENE_EXP_SNAKE_CASE,
            Self::Barcode(bc) => bc.as_snake_case(),
            Self::ProteinExpression => Self::PROTEIN_EXP_SNAKE_CASE,
            Self::Peaks => Self::PEAKS_SNAKE_CASE,
        }
    }

    /// Parse a snake-case representation of this feature type.
    pub fn from_snake_case(s: &str) -> Result<Self> {
        Ok(match s {
            Self::GENE_EXP_SNAKE_CASE => Self::Gene,
            Self::PROTEIN_EXP_SNAKE_CASE => Self::ProteinExpression,
            Self::PEAKS_SNAKE_CASE => Self::Peaks,
            _ => FeatureType::Barcode(FeatureBarcodeType::from_snake_case(s)?),
        })
    }
}

impl From<FeatureType> for &'static str {
    fn from(value: FeatureType) -> &'static str {
        value.as_str()
    }
}

impl FromStr for FeatureType {
    type Err = <FeatureBarcodeType as FromStr>::Err;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Ok(match s {
            Self::GENE_EXP => Self::Gene,
            Self::PROTEIN_EXP => Self::ProteinExpression,
            Self::PEAKS => Self::Peaks,
            _ => Self::Barcode(FeatureBarcodeType::from_str(s)?),
        })
    }
}

impl From<FeatureType> for LibraryType {
    fn from(value: FeatureType) -> Self {
        match value {
            FeatureType::Gene => LibraryType::Gex,
            FeatureType::Barcode(ft) => LibraryType::FeatureBarcodes(ft),
            FeatureType::Peaks => LibraryType::Atac,
            FeatureType::ProteinExpression => {
                unimplemented!(
                    "Protein Expression can't be converted to a LibraryType; the concept of library is not valid here."
                )
            }
        }
    }
}

impl AsMetricPrefix for FeatureType {
    /// Return the metric prefix for this feature type.
    fn as_metric_prefix(&self) -> Option<&str> {
        match self {
            Self::Gene => None,
            Self::Peaks => None,
            Self::Barcode(ft) => ft.as_metric_prefix(),
            Self::ProteinExpression => {
                unimplemented!(
                    "Protein Expression can't be converted to a metric prefix as there should be no metrics computed based on it here."
                )
            }
        }
    }
}

/// Feature reference entry
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct FeatureDef {
    pub index: usize,
    pub id: String,
    pub name: String,
    pub genome: GenomeName,
    pub sequence: String,
    pub pattern: String,
    pub read: WhichRead,
    pub feature_type: FeatureType,
    pub tags: HashMap<String, String>,
}

// a mask for which ascii characters are valid, equivalent to python
// > set(string.printable) - set(string.whitespace) - set("/,'\"\\`")
const VALID_ID_CHARS: &[bool] = &[
    false, false, false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, true, false, true, true, true, true, false,
    true, true, true, true, false, true, true, false, true, true, true, true, true, true, true,
    true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true,
    true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true,
    true, true, true, true, true, false, true, true, true, false, true, true, true, true, true,
    true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true,
    true, true, true, true, true, true, true, true, true, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false,
];

impl FeatureDef {
    pub fn from_string_record(
        data: &StringRecord,
        headers: &StringRecord,
        index: usize,
    ) -> Result<FeatureDef> {
        let get_col = |col_name: &str| {
            for i in 0..headers.len() {
                if &headers[i] == col_name {
                    return &data[i];
                }
            }

            panic!("column not found: {col_name}");
        };

        let mut tags = HashMap::new();

        for i in 0..headers.len() {
            if REQUIRED_FEATURE_REF_COLS.contains(&&headers[i]) {
                continue;
            }
            // skip extra CSV columns with an empty header name: fixes CR-4396
            if !headers[i].is_empty() {
                tags.insert(headers[i].to_string(), data[i].to_string());
            }
        }

        let feature_type_str = get_col("feature_type");
        let feature_type = FeatureBarcodeType::from_str(feature_type_str).map_err(|_| {
            anyhow!(
                "Unknown feature_type '{feature_type_str}' in the feature reference must be one of {}",
                FeatureBarcodeType::iter().format(", ")
            )
        })?;

        let read = match WhichRead::from_str(get_col("read")) {
            Ok(WhichRead::R1) => WhichRead::R1,
            Ok(WhichRead::R2) => WhichRead::R2,
            _ => bail!(
                "The feature definition file contains a read type value '{}' which is not one of the allowed read types 'R1' or 'R2'.",
                get_col("read")
            ),
        };

        let id = get_col("id").to_string();
        for (i, c) in id.chars().enumerate() {
            ensure!(
                !c.is_ascii_whitespace(),
                "Feature id field cannot contain whitespace: '{c}'",
            );
            let c_index = c as usize;
            ensure!(
                c_index < VALID_ID_CHARS.len() && VALID_ID_CHARS[c_index],
                "Feature id field contains illegal character at position {}: '{id}'. \
                 Feature ids may only ASCII characters, and must not use whitespace, \
                 slash, quote, or comma characters.",
                i + 1
            );
        }

        Ok(FeatureDef {
            index,
            id: get_col("id").to_string(),
            name: get_col("name").to_string(),
            genome: GenomeName::default(),
            sequence: get_col("sequence").to_string(),
            pattern: get_col("pattern").to_string(),
            read,
            feature_type: FeatureType::Barcode(feature_type),
            tags,
        })
    }

    pub fn add_tag(&mut self, tag_key: String, tag_value: String) -> Result<()> {
        // Insert returns the existing value if present
        match self.tags.insert(tag_key.clone(), tag_value) {
            Some(value) => {
                bail!("Duplicate value in feature reference: {tag_key} with value {value}")
            }
            None => Ok(()),
        }
    }

    fn len(&self) -> usize {
        // We previously only supported the regex ^/$ markers for start and end,
        // but in order to be more human-friendly, we now also support 5p and 3p too
        // replace them here to make a valid regex
        let re5p = Regex::new("^5[Pp]?[-_]?").unwrap();
        let pat = re5p.replace(&self.pattern, "");

        let re3p = Regex::new("[-_]?3[Pp]?$").unwrap();
        let pat = re3p.replace(&pat, "");

        let bc_regex = Regex::new(r"\(BC\)").unwrap();
        let pat = bc_regex.replace(&pat, &self.sequence);

        pat.chars()
            .fold(0usize, |len, c| match c.to_ascii_uppercase() {
                'A' | 'C' | 'G' | 'T' | 'N' => len + 1,
                _ => len,
            })
    }
}

fn validate_headers(record: &StringRecord) -> Result<()> {
    let hdrs: TxHashSet<_> = record.iter().collect();
    let (found, missing) = REQUIRED_FEATURE_REF_COLS
        .iter()
        .partition::<Vec<&str>, _>(|&req| hdrs.contains(req));
    if !missing.is_empty() {
        bail!(
            r#"The feature reference file header does not contain one or more required comma-separated fields: "{}".
The following fields were found: "{}".
Please check that your file is in CSV format and has the required field names."#,
            missing.join(", "),
            found.join(", ")
        );
    }
    Ok(())
}

martian_filetype! { FeatureReferenceFile, "csv" }

impl FeatureReferenceFile {
    /// Read the feature reference file, using an optional feature config.
    pub fn read(&self, feature_config: Option<&FeatureConfig>) -> Result<FeatureReference> {
        FeatureReference::from_csv(self, feature_config)
    }
}

/// Return the genome name of this feature ID.
pub fn get_genome_from_feature_id<'a>(id: &str, genomes: &'a [GenomeName]) -> &'a GenomeName {
    match genomes {
        [] => panic!("Trying to get genome from feature ID without any genomes"),
        [genome] => genome,
        _ => {
            let stripped_id = GenomeName::from(strip_deprecated_prefix(id));
            genomes
                .iter()
                .find(|&g| stripped_id.starts_with(g.as_str()))
                .with_context(|| {
                    format!(
                        "Genome for feature {id} not in genome list {}",
                        genomes.iter().format(", ")
                    )
                })
                .unwrap()
        }
    }
}

/// Contains all the 'features' in a GEX analysis.  A feature
/// is either a gene from the standard genome reference, or
/// a Feature Barcode declared in the feature reference.
/// Each feature is given an id, which is an index into the
/// `features_defs` field.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq, Default)]
pub struct FeatureReference {
    pub feature_defs: Vec<FeatureDef>,
    /// For each feature barcode type, a hashmap from the pattern sequence to the feature_id
    pub feature_maps: HashMap<FeatureBarcodeType, HashMap<String, Vec<usize>>>,
    // Lookup a gene index given a gene
    pub gene_to_index: HashMap<Gene, usize>,
    // Optional set of on-target features
    pub target_set: Option<TargetSet>,
}

impl FeatureReference {
    pub fn estimate_mem_bytes_for_feature_count(num_features: usize) -> usize {
        // The feature reference bincode file typically comes under 10MB for ~30k
        // features. Assuming 1kB per feature (slight over-estimate)
        1_024 * num_features
    }

    /// Return the number of features.
    pub fn num_features(&self) -> usize {
        self.feature_defs.len()
    }

    /// Create a new feature reference from a CSV file.
    fn from_csv(
        path: &FeatureReferenceFile,
        feature_config: Option<&FeatureConfig>,
    ) -> Result<Self> {
        let rdr = BufReader::new(File::open(path)?);
        let fref = FeatureReference::new(&[], None, Some(rdr), None, None, None, feature_config)?;
        Ok(fref)
    }

    /// Read a feature reference from a transcriptome reference, probe set, and feature reference CSV.
    pub fn from_paths(
        transcriptome_reference_path: Option<&Path>,
        feature_reference: Option<&FeatureReferenceFile>,
        target_set_name: Option<&str>,
        target_set_path: Option<&TargetSetFile>,
        target_features_path: Option<&TargetGeneIndicesFile>,
        feature_config: Option<&FeatureConfig>,
    ) -> Result<Self> {
        assert!(
            transcriptome_reference_path.is_some()
                || target_set_path.is_some()
                || feature_reference.is_some()
        );

        let transcriptome_genomes = transcriptome_reference_path
            .map(ReferenceInfo::from_reference_path)
            .transpose()?
            .map(|x| x.genomes);
        let probe_set_genomes = target_set_path
            .map(ProbeSetReferenceMetadata::load_from)
            .transpose()?
            .map(|x| x.reference_genomes());
        let genomes = transcriptome_genomes
            .or(probe_set_genomes)
            .unwrap_or_default();

        let txome = transcriptome_reference_path
            .map(Transcriptome::from_reference_path)
            .transpose()?;

        let target_set = target_set_name
            .zip(target_features_path)
            .map(|(name, path)| TargetSet::load(name, path))
            .transpose()?;

        let feature_ref_stream = feature_reference
            .map(|x| File::open(x).with_context(|| x.display().to_string()))
            .transpose()?;

        let target_genes_and_included = target_set_path
            .map(|x| x.read_genes_and_included(transcriptome_reference_path))
            .transpose()?;

        FeatureReference::new(
            &genomes,
            txome.as_ref(),
            feature_ref_stream,
            target_set_name,
            target_set,
            target_genes_and_included.as_deref(),
            feature_config,
        )
    }

    /// Create a new feature reference. See `from_paths` for details.
    #[allow(clippy::too_many_arguments)]
    pub fn new<R: Read>(
        genomes: &[GenomeName],
        txome: Option<&Transcriptome>,
        csv_stream: Option<R>,
        target_set_name: Option<&str>,
        target_set: Option<TargetSet>,
        target_genes_and_included: Option<&[(FeatureID, FeatureName, bool)]>,
        feature_config: Option<&FeatureConfig>,
    ) -> Result<FeatureReference> {
        // Create gene features
        let mut fdefs: Vec<_> = txome
            .map(|x| x.genes.as_slice())
            .unwrap_or_default()
            .iter()
            .enumerate()
            .map(|(index, gene)| FeatureDef {
                index,
                id: gene.id.clone(),
                name: gene.name.clone(),
                genome: get_genome_from_feature_id(&gene.id, genomes).clone(),
                sequence: String::default(),
                pattern: String::default(),
                read: WhichRead::R2,
                feature_type: FeatureType::Gene,
                tags: HashMap::new(),
            })
            .collect();

        // Create non-gene probe features.
        let probe_features = target_genes_and_included
            .unwrap_or_default()
            .iter()
            .filter(|(id, _, _)| !txome.is_some_and(|x| x.gene_id_to_idx.contains_key(id)))
            .zip(fdefs.len()..)
            .map(|((id, name, _included), index)| FeatureDef {
                index,
                id: id.clone(),
                name: name.clone(),
                genome: get_genome_from_feature_id(id, genomes).clone(),
                sequence: String::default(),
                pattern: String::default(),
                read: WhichRead::R2,
                feature_type: FeatureType::Gene,
                tags: HashMap::new(),
            });
        fdefs.extend(probe_features);

        // Construct the map of genes to indices.
        let gene_to_index: HashMap<Gene, usize> = fdefs
            .iter()
            .map(|x| {
                (
                    Gene {
                        id: x.id.clone(),
                        name: x.name.clone(),
                    },
                    x.index,
                )
            })
            .collect();

        // Construct the set of on-target feature indices.
        let target_set = if let Some(target_genes_and_included) = target_genes_and_included {
            // target_set contains all genes targeted by probes ignoring included when filter_probes is false.
            let target_set_indices_orig = target_set.map(|x| x.feature_indices).unwrap_or_default();
            let target_set_indices_probe = target_genes_and_included
                .iter()
                .filter(|&&(ref id, _, included)| included && !is_deprecated_probe(id))
                .cloned()
                .map(|(id, name, _)| gene_to_index[&Gene { id, name }] as u32);
            let target_set_indices =
                chain(target_set_indices_orig, target_set_indices_probe).collect();
            Some(TargetSet::from_indices(
                target_set_name.unwrap(),
                target_set_indices,
            ))
        } else {
            assert!(target_set.is_none());
            None
        };

        // Validate contents of specificty_controls in feature_config
        let mut seen_this_control_id = HashMap::new();
        let mut seen_this_mhc_allele = HashMap::new();
        if let Some(feature_config) = feature_config
            && let Some(specificty_controls) = &feature_config.specificity_controls
        {
            let mhc_alleles: Vec<String> = specificty_controls
                .control_for_allele
                .keys()
                .cloned()
                .collect();
            let control_ids: Vec<String> = specificty_controls
                .control_for_allele
                .values()
                .cloned()
                .collect();
            for mhc_allele in mhc_alleles {
                seen_this_mhc_allele.insert(mhc_allele, false);
            }
            for control_id in control_ids {
                seen_this_control_id.insert(control_id, false);
            }
            match (
                specificty_controls.has_mhc_allele_column,
                feature_config.beam_mode,
            ) {
                (true, Some(BeamMode::BeamAB)) => bail!(
                    "Error parsing Multi config CSV: The `mhc_allele` column in \
                     [antigen-specificity] section is invalid for BCR Antigen Capture."
                ),
                (false, Some(BeamMode::BeamT)) => bail!(
                    "Error parsing Multi config CSV: The `mhc_allele` column in \
                     [antigen-specificity] section is required for TCR Antigen Capture."
                ),
                _ => (),
            };
        }

        let mut fmaps = HashMap::new();
        let mut all_feature_ids = HashSet::new();
        if let Some(csv) = csv_stream {
            let reader = BufReader::new(csv);
            let mut csv_reader = csv::ReaderBuilder::new()
                .trim(csv::Trim::All)
                .from_reader(reader);

            let header = csv_reader.headers()?.clone();
            validate_headers(&header)?;

            let hashtag_ids = feature_config.and_then(|f| f.hashtag_ids.as_ref());

            let mut feature_types_seen: HashSet<FeatureBarcodeType> = HashSet::new();
            let mut last_feature_type: Option<FeatureBarcodeType> = None;
            let mut cmo_ids = HashMap::new();
            // Create feature barcode (fBC) features
            // TODO - needs a bunch of details brought over from Python code
            // see lib/python/cellranger/rna/feature_ref.py:parse_feature_def_file
            for record in csv_reader.records() {
                let record = record?;
                let feat_idx = fdefs.len();
                let mut fdef = FeatureDef::from_string_record(&record, &header, feat_idx)?;
                all_feature_ids.insert(fdef.id.clone());

                let FeatureType::Barcode(feature_barcode_type) = fdef.feature_type else {
                    unreachable!();
                };

                if feature_barcode_type == FeatureBarcodeType::Multiplexing {
                    let line = record.position().unwrap().line();
                    match cmo_ids.entry(fdef.id.clone()) {
                        Entry::Occupied(o) => {
                            bail!(
                                "Multiplexing Capture features sharing id \"{}\" found on lines {} and {}",
                                fdef.id,
                                o.get(),
                                line
                            )
                        }
                        Entry::Vacant(v) => {
                            v.insert(line);
                        }
                    }
                }
                if feature_barcode_type == FeatureBarcodeType::Antibody
                    && hashtag_ids.is_some_and(|x| x.contains(&fdef.id))
                {
                    fdef.add_tag(HASHTAG.to_string(), "True".to_string())?;
                }
                if feature_barcode_type == FeatureBarcodeType::Antigen {
                    if let Some(feature_config) = feature_config {
                        if let Some(specificity_controls) = &feature_config.specificity_controls {
                            // Add tags corresponding to control id
                            // First, find allele (if beam-ab, set to None)
                            let allele = match fdef.tags.get(MHC_ALLELE) {
                                Some(allelle) => {
                                    if feature_config.beam_mode == Some(BeamMode::BeamAB) {
                                        bail!(
                                            "Error parsing feature reference: The `mhc_allele` column is invalid for BCR Antigen Capture."
                                        );
                                    }
                                    // Throw an error if feature reference has MHC allele but not multi config csv
                                    if !specificity_controls.has_mhc_allele_column {
                                        bail!(
                                            "Feature reference CSV contains `{MHC_ALLELE}` column but `{MHC_ALLELE}` is not present in the \
                                    [antigen-specificity] section of the multi config CSV."
                                        );
                                    }
                                    allelle.as_str()
                                }
                                None => NO_ALLELE,
                            };
                            // Check if multi config CSV has correct allele and control feature id pairing
                            if allele != NO_ALLELE {
                                for (al, ctrl) in &specificity_controls.control_for_allele {
                                    if ctrl == &fdef.id && al != allele {
                                        bail!(
                                            "Error parsing feature reference: Feature id {ctrl} has MHC allele \
                                        parameter {allele} which is incompatible with the MHC allele {al} in [antigen-specificity] section of \
                                        Multi config CSV."
                                        );
                                    }
                                }
                            }

                            if let Some(seen) = seen_this_control_id.get_mut(&fdef.id) {
                                *seen = true;
                            }
                            if let Some(seen) = seen_this_mhc_allele.get_mut(allele) {
                                *seen = true;
                            }

                            let control_id = specificity_controls.control_for_allele.get(allele);

                            if let Some(control_id) = control_id {
                                // If the control_id for the current feature is equel to the current feature id,
                                // The current feature is non-targeting.
                                let is_targeting = if control_id.eq(&fdef.id) {
                                    "False"
                                } else {
                                    "True"
                                };

                                fdef.add_tag(
                                    TARGETING_ANTIGEN.to_string(),
                                    is_targeting.to_string(),
                                )?;
                            } else {
                                // if allele is not present in the specificity_controls, then this allele does not have
                                // a control.
                                fdef.add_tag(TARGETING_ANTIGEN.to_string(), "Null".to_string())?;
                            }
                        } else {
                            fdef.add_tag(TARGETING_ANTIGEN.to_string(), "Null".to_string())?;
                        }

                        // Only add tag to feature ref if [feature-functional-map] section
                        if let Some(functional_map) = &feature_config.functional_map {
                            if let Some(functional_name) = functional_map.get(&fdef.id) {
                                fdef.add_tag(FUNCTIONAL_NAME.to_string(), functional_name.clone())?;
                            } else {
                                // If [feature-functional-map] exists but it doesn't include this feature, add Null
                                fdef.add_tag(FUNCTIONAL_NAME.to_string(), "Null".to_string())?;
                            }
                        }
                    } else {
                        fdef.add_tag(TARGETING_ANTIGEN.to_string(), "Null".to_string())?;
                    }
                }
                match last_feature_type {
                    None => {
                        feature_types_seen.insert(feature_barcode_type);
                    }
                    Some(lft) => {
                        if feature_barcode_type != lft
                            && feature_types_seen.contains(&feature_barcode_type)
                        {
                            // If we have switched from the last feature type, but the new feature type is already in the feature types seen,
                            // it means we have had a block of the new feature type before the block of last feature type:
                            // Feature ref example:
                            // Feature type A
                            // Feature type B
                            // Feature type B -> last feature type
                            // Feature type A -> new feature type
                            bail!(
                                "Features of the same type must be continous in the feature reference file."
                            );
                        }
                    }
                };
                last_feature_type = Some(feature_barcode_type);
                fmaps
                    .entry(feature_barcode_type)
                    .or_insert_with(HashMap::new)
                    .entry(fdef.sequence.clone())
                    .or_insert_with(Vec::new)
                    .push(feat_idx);

                fdefs.push(fdef);
            }

            let duplicate_ids: Vec<_> = fdefs.iter().map(|x| &x.id).duplicates().collect();
            ensure!(
                duplicate_ids.is_empty(),
                "{} duplicate feature ID(s) found in feature reference: {}",
                duplicate_ids.len(),
                duplicate_ids.into_iter().format(", ")
            );
        }

        // Validate the content of feature-functional-map in feature_config
        let mut feature_ids_in_functional_map = HashSet::new();
        if let Some(feature_config) = feature_config
            && let Some(functional_map) = &feature_config.functional_map
        {
            for feature_id in functional_map.keys() {
                feature_ids_in_functional_map.insert(feature_id.to_string());
            }
        }
        if !feature_ids_in_functional_map.is_subset(&all_feature_ids) {
            bail!(
                "Error parsing feature reference: Feature ids {:?} are in [feature-functional-map] section \
                of Multi config CSV but are not in the feature reference.",
                &feature_ids_in_functional_map - &all_feature_ids
            );
        }

        let mut missing_control_ids = Vec::new();
        for (control_id, seen) in seen_this_control_id {
            if !seen {
                missing_control_ids.push(control_id);
            }
        }
        if !missing_control_ids.is_empty() {
            missing_control_ids.sort();
            bail!(
                "Antigen Capture feature id(s) {missing_control_ids:?} provided in \
                 [antigen-specificity] section of multi config CSV not found in feature reference.",
            );
        }
        let mut missing_mhc_alleles = Vec::new();
        for (mhc_allele, seen) in seen_this_mhc_allele {
            if !seen {
                missing_mhc_alleles.push(mhc_allele);
            }
        }
        if !missing_mhc_alleles.is_empty() {
            missing_mhc_alleles.sort();
            bail!(
                "Antigen Capture MHC allele(s) {:?} provided in [antigen-specificity] section of multi config CSV not found in feature reference.",
                &missing_mhc_alleles
            );
        }

        Ok(FeatureReference {
            feature_defs: fdefs,
            feature_maps: fmaps,
            gene_to_index,
            target_set,
        })
    }

    /// This function is only safe to append Multiplexing Capture features,
    /// other feature invariants may not be preserved
    pub fn append_feature_defs(&mut self, fdefs: &[FeatureDef]) {
        // this does not need to validate for existing CMO ids because CMOs are only allowed from
        // only one of:
        //   1) the cmo-set
        //   2) the feature-ref
        //   3) builtins
        for fdef in fdefs {
            use FeatureBarcodeType::Multiplexing;
            assert_eq!(fdef.feature_type, FeatureType::Barcode(Multiplexing));
            let index = self.feature_defs.len();
            self.feature_maps
                .entry(Multiplexing)
                .or_default()
                .entry(fdef.sequence.clone())
                .or_default()
                .push(index);
            self.feature_defs.push(FeatureDef {
                index,
                ..fdef.clone()
            });
        }
    }

    pub fn gene_index(&self, gene: &Gene) -> usize {
        self.gene_to_index[gene]
    }

    pub fn has_target_features(&self) -> bool {
        self.target_set.is_some()
    }

    /// Get the set of Genes corresponding to self.target_set, if not None
    /// TODO: TargetSet should backed by this HashSet of Genes, not u32
    pub fn target_genes(&self) -> Option<HashSet<Gene>> {
        let target_set = self.target_set.as_ref()?;
        Some(
            self.gene_to_index
                .iter()
                .filter_map(|(gene, &index)| {
                    target_set
                        .is_on_target(index as u32)
                        .then_some(gene.clone())
                })
                .collect(),
        )
    }

    /// Get the set of Multiplexing Capture feature ids
    pub fn multiplexing_ids(&self) -> TxHashSet<String> {
        self.feature_defs
            .iter()
            .filter(|x| x.feature_type == FeatureType::Barcode(FeatureBarcodeType::Multiplexing))
            .map(|x| x.id.clone())
            .collect()
    }

    /// Get the set of Antibody feature ids
    pub fn antibody_ids(&self) -> TxHashSet<String> {
        self.feature_defs
            .iter()
            .filter(|x| x.feature_type == FeatureType::Barcode(FeatureBarcodeType::Antibody))
            .map(|x| x.id.clone())
            .collect()
    }

    /// Get the minimum read lengths per feature type
    pub fn min_feature_read_lengths(&self) -> TxHashMap<FeatureType, TxHashMap<WhichRead, usize>> {
        self.feature_defs
            .iter()
            .fold(TxHashMap::default(), |mut set, x| {
                let read = x.read;
                let len = x.len();
                match set.entry(x.feature_type) {
                    Entry::Occupied(mut e) => match e.get_mut().entry(read) {
                        Entry::Occupied(mut e) if len < *e.get() => {
                            e.insert(len);
                        }
                        Entry::Vacant(e) => {
                            e.insert(len);
                        }
                        Entry::Occupied(_) => {}
                    },
                    Entry::Vacant(e) => {
                        let mut v = TxHashMap::default();
                        v.insert(read, len);
                        e.insert(v);
                    }
                }
                set
            })
    }

    /// Write the feature reference to a CSV file.
    pub fn to_csv(&self, w: &mut impl Write) -> Result<()> {
        let custom_tags: Vec<_> = self
            .feature_defs
            .iter()
            .flat_map(|d| d.tags.keys())
            .unique()
            .sorted()
            .cloned()
            .collect();
        write!(w, "{}", REQUIRED_FEATURE_REF_COLS.join(","))?;
        for ctag in &custom_tags {
            write!(w, ",{ctag}")?;
        }
        writeln!(w)?;
        for fdef in &self.feature_defs {
            // we're writing out the feature reference, skip gene features
            if fdef.feature_type == FeatureType::Gene {
                continue;
            }
            let read = match fdef.read {
                WhichRead::R1 => "R1",
                WhichRead::R2 => "R2",
                WhichRead::I1 => "I1",
                WhichRead::I2 => "I2",
            };
            // "id", "name", "read", "pattern", "sequence", "feature_type"
            write!(
                w,
                "{},{},{read},{},{},{}",
                fdef.id,
                fdef.name,
                fdef.pattern,
                fdef.sequence,
                fdef.feature_type.as_str(),
            )?;
            for ctag in &custom_tags {
                let s = fdef.tags.get(ctag).map_or("", String::as_str);
                write!(w, ",{s}")?;
            }
            writeln!(w)?;
        }
        Ok(())
    }

    /// Write the feature reference to a TSV file. Creates the `features.tsv.gz` file in the output folder.
    pub fn to_tsv(&self, w: &mut impl Write) -> Result<()> {
        for fd in &self.feature_defs {
            write!(w, "{}\t{}\t{}", fd.id, fd.name, fd.feature_type.as_str())?;
        }
        Ok(())
    }

    /// Validate Beam feature reference (make sure at most one non-targeting antigen_type per beam-ab feature ref
    /// and one non-targeting per allele for beam-t)
    pub fn validate_beam_feature_ref(&self, beam_mode: BeamMode) -> Result<()> {
        let antigen_type_per_allele = self
            .feature_defs
            .iter()
            .filter(|fdef| fdef.feature_type == FeatureType::Barcode(FeatureBarcodeType::Antigen))
            .map(|fdef| {
                let allele = fdef.tags.get(MHC_ALLELE).map_or("", String::as_str);
                let antigen_type = fdef
                    .tags
                    .get(TARGETING_ANTIGEN)
                    .map_or("", String::as_str)
                    .to_lowercase();
                (allele, antigen_type)
            })
            .into_group_map();

        for (allele, antigen_types) in antigen_type_per_allele {
            match beam_mode {
                BeamMode::BeamAB => ensure!(
                    allele.is_empty(),
                    "Error parsing feature reference: Value '{allele}' provided for mhc_allele. \
                     The `mhc_allele` column is invalid for BCR Antigen Capture.",
                ),
                BeamMode::BeamT => ensure!(
                    !allele.is_empty(),
                    "Error parsing feature reference: Missing value for mhc_allele. \
                     The `mhc_allele` column is required for TCR Antigen capture"
                ),
            }

            let num_empty = antigen_types.iter().filter(|&f| f.is_empty()).count();
            ensure!(
                num_empty == 0 || num_empty == antigen_types.len(),
                if allele.is_empty() {
                    anyhow!(
                        "Error parsing feature reference: Empty values for antigen_type detected"
                    )
                } else {
                    anyhow!(
                        "Error parsing feature reference: Empty values for antigen_type detected \
                         for mhc_allele: '{allele}'"
                    )
                }
            );

            let num_non_targeting = antigen_types.into_iter().filter(|f| f == "False").count();
            ensure!(
                num_non_targeting < 2,
                if allele.is_empty() {
                    anyhow!(
                        "Error parsing feature reference: More than one non-targeting \
                         antigen_type detected"
                    )
                } else {
                    anyhow!(
                        "Error parsing feature reference: More than one non-targeting \
                         antigen_type detected for mhc_allele: '{allele}'"
                    )
                }
            );
        }

        Ok(())
    }

    pub fn iter_feature_defs<F>(&self, feature_type_filter: F) -> impl Iterator<Item = &FeatureDef>
    where
        F: Fn(FeatureType) -> bool,
    {
        self.feature_defs
            .iter()
            .filter(move |fdef| feature_type_filter(fdef.feature_type))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::{BufWriter, Seek};
    use std::slice;

    // initialize insta test harness
    #[ctor::ctor]
    fn init() {
        // this ensures insta knows where to find its snap tests
        let cwd = std::env::current_dir().unwrap();
        let workspace_root = cwd.parent().unwrap();
        unsafe { std::env::set_var("INSTA_WORKSPACE_ROOT", workspace_root) }
    }

    #[test]
    fn test_serialize_feature_type() -> Result<()> {
        let gex = "\"Gene Expression\"";
        assert_eq!(serde_json::to_string(&FeatureType::Gene)?, gex);
        assert_eq!(serde_json::from_str::<FeatureType>(gex)?, FeatureType::Gene);

        let mut f = tempfile::tempfile()?;
        bincode::serialize_into(BufWriter::new(&f), &FeatureType::Gene)?;
        f.rewind()?;
        let gex: FeatureType = bincode::deserialize_from(BufReader::new(f))?;
        assert_eq!(gex, FeatureType::Gene);

        Ok(())
    }

    #[test]
    fn test_get_genome_from_feature_id() -> Result<()> {
        let x = GenomeName::from("X");
        let ab = GenomeName::from("A_B");
        assert_eq!(get_genome_from_feature_id("A_B_C", slice::from_ref(&x)), &x);
        assert_eq!(
            get_genome_from_feature_id("A_B_C", &[x.clone(), ab.clone()]),
            &ab
        );
        assert_eq!(
            get_genome_from_feature_id("DEPRECATED_A_B_C", slice::from_ref(&x)),
            &x
        );
        assert_eq!(
            get_genome_from_feature_id("DEPRECATED_A_B_C", &[x.clone(), ab.clone()]),
            &ab
        );
        Ok(())
    }

    #[test]
    #[should_panic]
    fn test_get_genome_from_feature_id_genome_not_found() {
        let x = GenomeName::from("X");
        let y = GenomeName::from("Y");
        get_genome_from_feature_id("DEPRECATED_C", &[x.clone(), y.clone()]);
    }

    #[test]
    #[should_panic]
    fn test_get_genome_from_deprecated_feature_id_genome_not_found() {
        let x = GenomeName::from("X");
        let y = GenomeName::from("Y");
        get_genome_from_feature_id("DEPRECATED_C", &[x.clone(), y.clone()]);
    }

    #[test]
    #[should_panic]
    fn test_get_genome_from_feature_id_empty_genomes() {
        get_genome_from_feature_id("A_B_C", &[]);
    }
}

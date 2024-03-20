//! Martian stage WRITE_MOLECULE_INFO

use crate::probe_barcode_matrix::ProbeCounts;
use crate::stages::collate_metrics::SampleMetrics;
use crate::stages::parse_multi_config::CellCallingParam;
use crate::types::FeatureReferenceFormat;
use crate::BcUmiInfoShardFile;
use anyhow::Result;
use barcode::{Barcode, BarcodeContent};
use cr_h5::molecule_info::MoleculeInfoWriter;
use cr_types::barcode_index::BarcodeIndex;
use cr_types::chemistry::{ChemistryDefs, ChemistryDefsExt};
use cr_types::filtered_barcodes::{FilteredBarcodesCsv, FilteredBarcodesCsvRow};
use cr_types::probe_set::Probe;
use cr_types::reference::feature_reference::FeatureReference;
use cr_types::reference::reference_info::ReferenceInfo;
use cr_types::rna_read::{LibraryInfo, RnaChunk};
use cr_types::target_panel_summary::{TargetPanelSummary, TargetPanelSummaryFormat};
use cr_types::{
    AlignerParam, BarcodeIndexFormat, BarcodeToSample, BcUmiInfo, GemWell, GenomeName, H5File,
    MetricsFile, SampleAssignment, SampleBarcodes, SampleBarcodesFile,
};
use itertools::Itertools;
use json_report_derive::JsonReport;
use martian::{MartianRover, MartianStage, MartianVoid, Resource, StageDef};
use martian_derive::{make_mro, MartianStruct};
use martian_filetypes::json_file::JsonFile;
use martian_filetypes::tabular_file::CsvFile;
use martian_filetypes::FileTypeRead;
use metric::{join_metric_name, CountMetric, JsonReport, PercentMetric, TxHashMap, TxHashSet};
use serde::{Deserialize, Serialize};
use serde_json::{json, Map, Value};
use shardio::ShardReader;
use std::cmp::max;
use std::collections::HashMap;
use std::path::PathBuf;
use strum_macros::IntoStaticStr;

/// Describes whether a molecule info is raw (uber), sample, or count.
/// Stored in molecule info metrics and used to set aggr preflights.
#[derive(Clone, Copy, IntoStaticStr)]
#[strum(serialize_all = "snake_case")]
pub enum MoleculeInfoType {
    Raw,
    PerSample,
    Count,
}

pub struct WriteMoleculeInfo;

#[derive(Clone, Deserialize, MartianStruct)]
pub struct WriteMoleculeInfoStageInputs {
    pub chemistry_defs: ChemistryDefs,
    pub gem_well: GemWell,
    pub counts_bc_order: Vec<BcUmiInfoShardFile>,
    pub reference_path: PathBuf,
    pub read_chunks: Vec<RnaChunk>,
    pub feature_reference: FeatureReferenceFormat,
    pub filtered_barcodes: FilteredBarcodesCsv,
    pub per_probe_metrics: Option<CsvFile<ProbeCounts>>,
    pub target_panel_summary: Option<TargetPanelSummaryFormat>,
    pub matrix_computer_summary: MetricsFile,
    pub recovered_cells: Option<CellCallingParam>,
    pub force_cells: Option<CellCallingParam>,
    pub include_introns: bool,
    pub disable_ab_aggregate_detection: bool,
    pub disable_high_occupancy_gem_detection: bool,
    pub filter_probes: Option<bool>,
    pub multi_config_sha: Option<String>,
    pub sample_barcodes: Option<SampleBarcodesFile>,
    pub per_sample_metrics: Option<Vec<SampleMetrics>>,
    pub barcode_index: BarcodeIndexFormat,
    pub slide_serial_capture_area: Option<String>,
}

#[derive(Serialize, MartianStruct)]
pub struct SampleMoleculeInfo {
    pub sample: SampleAssignment,
    pub h5_file: H5File,
    pub summary: JsonFile<()>, // molecule info metrics summary
}

#[derive(Serialize, MartianStruct)]
pub struct WriteMoleculeInfoStageOutputs {
    pub single_mol_info: Option<SampleMoleculeInfo>,
    pub multi_mol_info: Option<Vec<SampleMoleculeInfo>>,
}

/// A JSON object, the inner type of serde_json::Value::Object
type JMap = Map<String, Value>;

/// Convert value to a JSON object, and panic if serialization fails or does not produce an Object.
fn to_object(value: impl Serialize) -> JMap {
    if let Value::Object(map) = serde_json::to_value(value).unwrap() {
        map
    } else {
        panic!("Argument of to_object is not a JSON object");
    }
}

/// Calculate the metric frac_usable_confidently_mapped_reads_on_target.
/// Return None if there is no targeted gene expression library.
fn calculate_frac_usable_confidently_mapped_reads_on_target(metrics: &Value) -> Option<f64> {
    metrics.as_object().unwrap()[LIBRARIES_METRIC]
        .as_object()
        .unwrap()
        .iter()
        .find_map(|(_library_index, library_metrics)| {
            library_metrics.get(ON_TARGET_USABLE_READS_METRIC).map(|x| {
                let on_target_usable_read_pairs = x.as_i64().unwrap();
                let raw_read_pairs = library_metrics[TOTAL_READS_METRIC].as_i64().unwrap();
                PercentMetric::from_parts(on_target_usable_read_pairs, raw_read_pairs)
                    .fraction()
                    .unwrap_or(f64::NAN)
            })
        })
}

/// Return the number of valid barcodes from a metrics JSON.
/// An antibody-only analysis has no metric barcodes_detected,
/// so return the maximum of barcodes_detected and ANTIBODY_barcodes_detected.
fn get_valid_barcodes(metrics: &HashMap<String, Value>) -> u64 {
    let gex_barcodes_detected = metrics.get("barcodes_detected").and_then(Value::as_u64);
    let antibody_barcodes_detected = metrics
        .get("ANTIBODY_barcodes_detected")
        .and_then(Value::as_u64);
    max(gex_barcodes_detected, antibody_barcodes_detected).unwrap()
}

#[make_mro(volatile = strict)]
impl MartianStage for WriteMoleculeInfo {
    type StageInputs = WriteMoleculeInfoStageInputs;
    type ChunkInputs = MartianVoid;
    type ChunkOutputs = MartianVoid;
    type StageOutputs = WriteMoleculeInfoStageOutputs;

    fn split(
        &self,
        args: Self::StageInputs,
        _rover: MartianRover,
    ) -> Result<StageDef<Self::ChunkInputs>> {
        // Multi uses more memory for sample_barcodes and other data.
        let sample_barcodes_mem_gib = if args.sample_barcodes.is_some() {
            // For MULTI_WRITE_PER_SAMPLE_MOLECULE_INFO.
            /// Loading sample_barcodes.json uses this many bytes per valid barcode.
            const MEM_BYTES_PER_VALID_BARCODE: u64 = 384;
            let mem_bytes = MEM_BYTES_PER_VALID_BARCODE
                * get_valid_barcodes(&args.matrix_computer_summary.read()?);
            1 + mem_bytes / 1024 / 1024 / 1024
        } else {
            // For WRITE_MOLECULE_INFO.
            0
        };
        let mem_gib = 8 + isize::try_from(sample_barcodes_mem_gib).unwrap();
        Ok(StageDef::with_join_resource(Resource::with_mem_gb(mem_gib)))
    }

    fn main(
        &self,
        _args: Self::StageInputs,
        _chunk_args: MartianVoid,
        _rover: MartianRover,
    ) -> Result<Self::ChunkOutputs> {
        unreachable!()
    }

    fn join(
        &self,
        args: Self::StageInputs,
        _chunk_defs: Vec<MartianVoid>,
        _chunk_outs: Vec<MartianVoid>,
        rover: MartianRover,
    ) -> Result<Self::StageOutputs> {
        let target_panel_summary = args
            .target_panel_summary
            .as_ref()
            .map(TargetPanelSummaryFormat::read)
            .transpose()?;
        let target_panel_summary = target_panel_summary.as_ref();

        let library_info = cr_types::rna_read::make_library_info(
            &args.read_chunks,
            target_panel_summary.map(|x| x.target_panel_name.as_str()),
        );

        let feature_reference = args.feature_reference.read()?;
        let targets_per_lib = load_targets_per_lib(&feature_reference, &library_info)?;

        let (probes, filtered_probes) =
            if let Some(per_probe_metrics_path) = &args.per_probe_metrics {
                let per_probe_metrics = per_probe_metrics_path.read()?;
                let filtered_probes: Vec<bool> =
                    per_probe_metrics.iter().map(|x| x.pass_filter).collect();
                let probes: Vec<Probe> = per_probe_metrics.into_iter().map_into().collect();
                (Some(probes), Some(filtered_probes))
            } else {
                (None, None)
            };

        let sample_to_metrics: TxHashMap<&SampleAssignment, &MetricsFile> =
            if let Some(sample_metrics_lst) = &args.per_sample_metrics {
                sample_metrics_lst
                    .iter()
                    .map(|x| (&x.sample, &x.summary))
                    .collect()
            } else {
                TxHashMap::default()
            };

        let ref_info = ReferenceInfo::from_reference_path(&args.reference_path)?;

        let library_filtered_genome_barcodes = args.filtered_barcodes.read()?;
        let library_filtered_barcodes: TxHashSet<&Barcode> = library_filtered_genome_barcodes
            .iter()
            .map(|x| &x.barcode as &Barcode)
            .collect();

        // This function call take a long time and can consume a lot of memory
        // TODO: Look for optimizations
        let reader = ShardReader::<BcUmiInfo, BcUmiInfo>::open_set(&args.counts_bc_order)?;

        let barcode_index = args.barcode_index.read()?;
        let sample_barcodes = SampleBarcodes::read_from_json(args.sample_barcodes.as_ref())?;
        let barcode_to_sample = BarcodeToSample::construct(&sample_barcodes);

        let molecule_info_type = if sample_barcodes.is_multiplexed() {
            // "multiplexed" in this context describes any multi run
            // where we are writing a per-sample molecule info rather than a count or raw mol info
            // even if it's "non-multiplexed" multi run, that single sample is treated like a single muxed sample internally
            MoleculeInfoType::PerSample
        } else {
            // if not multiplexed, it's either a count molecule info or an raw (uber) molecule info from multi.
            // if it's a multi run then the multi_config_sha will be something but it'll be None for Count runs
            if args.multi_config_sha.is_some() {
                MoleculeInfoType::Raw
            } else {
                MoleculeInfoType::Count
            }
        };
        // for each sample in the multiplex experiment, create a mol info file, visitor, and writer
        let mut multi_tmp = HashMap::new();
        let mut multi_visitor = HashMap::new();
        let mut multi_writer = HashMap::new();
        let library_info: Vec<_> = library_info
            .into_iter()
            .map(cr_types::LibraryInfo::Count)
            .collect();
        let valid_barcodes: Vec<BarcodeContent> = barcode_index
            .sorted_barcodes()
            .iter()
            .map(|barcode| *barcode.content())
            .collect();

        for sample in sample_barcodes.get_samples() {
            let h5_path: H5File = rover.make_path(format!("molecule_info_{sample}"));
            let mut writer = MoleculeInfoWriter::new(
                &h5_path,
                &feature_reference,
                probes.as_deref(),
                filtered_probes.as_deref(),
                &valid_barcodes,
                &library_info,
            )?;

            let sample_filtered_barcodes = match sample_barcodes.get_barcodes(sample) {
                Some(barcodes) => barcodes
                    .iter()
                    .filter(|barcode| library_filtered_barcodes.contains(barcode))
                    .collect(),
                None => library_filtered_barcodes.clone(),
            };

            // pass_filter table of molecule_info
            let genomes: Vec<_> = ref_info.genomes.iter().sorted().cloned().collect();
            let barcode_info_pass_filter = build_barcode_info(
                &library_filtered_genome_barcodes,
                &library_info,
                &barcode_index,
                &sample_filtered_barcodes,
                &genomes,
            );
            writer.write_barcode_info(&barcode_info_pass_filter, &genomes)?;
            multi_tmp.insert(sample.clone(), h5_path);
            multi_writer.insert(sample.clone(), writer);
            multi_visitor.insert(
                sample,
                MoleculeInfoVisitor {
                    sample_filtered_barcodes,
                    targets_per_lib: targets_per_lib.as_deref(),
                    metrics: MoleculeInfoMetrics::default(),
                },
            );
        }
        // At this point the memory for this stage seems to be near its maximum
        // for each barcode, write to the appropriate sample molecule info
        for bc_umi_info in reader.iter()? {
            let bc_umi_info = bc_umi_info?;
            let sample = barcode_to_sample.get_sample(&bc_umi_info.barcode);
            if sample != &SampleAssignment::Unassigned {
                multi_visitor.get_mut(sample).unwrap().visit(&bc_umi_info);
                multi_writer
                    .get_mut(sample)
                    .unwrap()
                    .fill(&barcode_index, &bc_umi_info)?;
            }
        }

        // for each sample, write metrics summary to file and molecule info h5
        let multi_mol_info_iter = sample_barcodes.get_samples().into_iter().map(|sample| {
            // safe to unwrap, sample_assignments is the set of keys used to populate these
            let visitor = multi_visitor.get_mut(sample).unwrap();
            let writer = multi_writer.get_mut(sample).unwrap();
            let h5_file = multi_tmp[sample].clone();
            writer.flush()?;

            // various metrics dumped into molecule_info.h5 to support aggr
            // if this is non-multiplexed, should use the raw metrics
            // otherwise use sliced metrics
            let per_sample_metrics_json = if sample_barcodes.is_multiplexed() {
                sample_to_metrics
                    .get(sample)
                    .expect("Sample did not have per-sample metrics JSON.")
            } else {
                &args.matrix_computer_summary
            };

            let metrics = get_molecule_info_metrics(
                &args,
                per_sample_metrics_json,
                &visitor.metrics,
                molecule_info_type,
                &ref_info,
                target_panel_summary,
            )?;
            writer.write_metrics(&serde_json::to_string(&metrics)?)?;

            // Write the SAMPLE_summary.json.
            let summary = rover.make_path(format!("{sample}_summary"));
            let mut reporter = visitor.metrics.to_json_reporter();
            if let Some(x) = calculate_frac_usable_confidently_mapped_reads_on_target(&metrics) {
                reporter.insert("frac_usable_confidently_mapped_reads_on_target", x);
            }
            reporter.report(&summary)?;

            anyhow::Ok(SampleMoleculeInfo {
                sample: sample.clone(),
                h5_file,
                summary,
            })
        });
        let multi_mol_info = multi_mol_info_iter.try_collect()?;

        Ok(if sample_barcodes.is_multiplexed() {
            Self::StageOutputs {
                single_mol_info: None,
                multi_mol_info: Some(multi_mol_info),
            }
        } else {
            // legacy output, for convenience with martian until it's not needed
            Self::StageOutputs {
                single_mol_info: Some(multi_mol_info.into_iter().next().unwrap()),
                multi_mol_info: None,
            }
        })
    }
}

const REFERENCE_METRICS: &[&str] = &[
    "reference_mkref_version",
    "reference_fasta_hash",
    "reference_gtf_hash",
    "reference_gtf_hash.gz",
];

fn get_library_metrics(
    read_chunks: &[RnaChunk],
    summary: &JMap,
    metrics: &MoleculeInfoMetrics,
    target_set_name: Option<&str>,
) -> JMap {
    // Get the library info for unique libraries
    let library_info = cr_types::rna_read::make_library_info(read_chunks, target_set_name);

    // Ensure that there is at most one gene expression library.
    let number_of_gene_expression_libraries = library_info
        .iter()
        .filter(|x| x.library_type.is_gex())
        .count();
    assert!(number_of_gene_expression_libraries <= 1);

    // Get the total reads per library metric for each library
    let mut lib_metrics = JMap::new();
    for (lib_idx, lib) in library_info.iter().enumerate() {
        let lib_idx = &(lib_idx as u16);

        let mut lm = JMap::new();
        lm.insert(
            TOTAL_READS_METRIC.to_string(),
            summary[&join_metric_name(lib.library_type, "total_read_pairs")].clone(),
        );

        // total_read_pairs_in_filtered_barcodes (same as above but restricted to molecules associated with cell barcodes)
        // not always present in metrics, so only put into molecule info file if present
        let total_read_pairs_in_filtered_barcodes = summary.get(&join_metric_name(
            lib.library_type,
            "total_read_pairs_in_filtered_barcodes",
        ));
        if let Some(value) = total_read_pairs_in_filtered_barcodes {
            lm.insert(
                TOTAL_READS_IN_FILTERED_BARCODES_METRIC.to_string(),
                value.clone(),
            );
        }

        // usable reads -- computed in this stage.
        lm.insert(
            USABLE_READS_METRIC.to_string(),
            metrics
                .usable_read_pairs
                .get(lib_idx)
                .map_or(0, |x| x.count())
                .into(),
        );

        // feature reads -- also computed in this stage
        lm.insert(
            FEATURE_READS_METRIC.to_string(),
            metrics
                .feature_read_pairs
                .get(lib_idx)
                .map_or(0, |x| x.count())
                .into(),
        );

        if target_set_name.is_some() {
            lm.insert(
                ON_TARGET_USABLE_READS_METRIC.to_string(),
                metrics
                    .on_target_usable_read_pairs
                    .get(lib_idx)
                    .map_or(0, |x| x.count())
                    .into(),
            );
        }

        lib_metrics.insert(lib.library_id.to_string(), Value::Object(lm));
    }

    lib_metrics
}

const GEM_GROUPS_METRIC: &str = "gem_groups";
const LIBRARIES_METRIC: &str = "libraries";
const GG_RECOVERED_CELLS_METRIC: &str = "recovered_cells";
const GG_FORCE_CELLS_METRIC: &str = "force_cells";
const SLIDE_SERIAL_CAPTURE_AREA: &str = "slide_serial_capture_area";

const TOTAL_READS_METRIC: &str = "raw_read_pairs";
const TOTAL_READS_IN_FILTERED_BARCODES_METRIC: &str = "raw_read_pairs_in_filtered_barcodes";
const USABLE_READS_METRIC: &str = "usable_read_pairs";
const FEATURE_READS_METRIC: &str = "feature_read_pairs";
const ON_TARGET_USABLE_READS_METRIC: &str = "on_target_usable_read_pairs";

const ANALYSIS_PARAMETERS_METRIC: &str = "analysis_parameters";

#[derive(Serialize)]
pub struct AnalysisParameters {
    include_introns: Option<bool>,
    filter_probes: Option<bool>,
    filter_aggregates: bool,
    filter_high_occupancy_gems: bool,
}

// assemble a collection of metrics that go into the molecule_info.h5 to support aggr.
fn get_gem_group_metrics(
    gem_groups: &[u16],
    force_cells: &Option<CellCallingParam>,
    recovered_cells: &Option<CellCallingParam>,
) -> JMap {
    let mut gg_metrics = JMap::new();
    for gg in gem_groups {
        let recovered_cells = recovered_cells
            .as_ref()
            .and_then(|rc| rc.sum().map(|x| x / gem_groups.len() as f64));
        let force_cells = force_cells
            .as_ref()
            .and_then(|fc| fc.sum().map(|x| x / gem_groups.len() as f64));

        let gm = json!({
            GG_RECOVERED_CELLS_METRIC: recovered_cells,
            GG_FORCE_CELLS_METRIC: force_cells,
        });

        gg_metrics.insert(gg.to_string(), gm);
    }

    gg_metrics
}

// get_molecule_info_metrics needs both the sliced and unsliced
// metrics because the sliced metrics don't have cellranger_version
// and alignment_aligner
fn get_molecule_info_metrics(
    args: &WriteMoleculeInfoStageInputs,
    sample_metrics_json: &MetricsFile,
    metrics: &MoleculeInfoMetrics,
    molecule_info_type: MoleculeInfoType,
    ref_info: &ReferenceInfo,
    target_panel_summary: Option<&TargetPanelSummary>,
) -> Result<Value> {
    let gem_well_summary_map: JMap = args.matrix_computer_summary.read()?;

    let aligner: AlignerParam = gem_well_summary_map["alignment_aligner"]
        .as_str()
        .unwrap()
        .parse()?;

    let mut final_metrics: JMap = ref_info
        .clone()
        .into_json_report()
        .into_iter()
        .filter(|(k, _)| REFERENCE_METRICS.contains(&k.as_str()))
        .collect();

    final_metrics.insert(
        "cellranger_version".to_string(),
        gem_well_summary_map["cellranger_version"].clone(),
    );

    final_metrics.extend(args.chemistry_defs.primary().to_json_reporter());
    final_metrics.insert(
        "chemistry_defs".to_string(),
        serde_json::to_value(&args.chemistry_defs).unwrap(),
    );

    // insert the molecule info type into the metrics
    final_metrics.insert(
        "molecule_info_type".to_string(),
        <&str>::from(molecule_info_type).into(),
    );

    // Include information about the the target set.
    if let Some(target_panel_summary) = target_panel_summary {
        final_metrics.extend(to_object(target_panel_summary));
    }

    let sample_summary_map: JMap = sample_metrics_json.read()?;
    let lib_metrics = get_library_metrics(
        &args.read_chunks,
        &sample_summary_map,
        metrics,
        target_panel_summary.map(|x| x.target_panel_name.as_str()),
    );
    final_metrics.insert(LIBRARIES_METRIC.to_string(), Value::Object(lib_metrics));

    let gem_groups = [args.gem_well.inner()];
    let gg_metrics = get_gem_group_metrics(&gem_groups, &args.force_cells, &args.recovered_cells);
    final_metrics.insert(GEM_GROUPS_METRIC.to_string(), Value::Object(gg_metrics));

    // analysis parameters metrics could come from either sample or gem-well
    final_metrics.insert(
        ANALYSIS_PARAMETERS_METRIC.to_string(),
        serde_json::to_value(AnalysisParameters {
            include_introns: (aligner == AlignerParam::Star).then_some(args.include_introns),
            filter_probes: args.filter_probes,
            filter_aggregates: !args.disable_ab_aggregate_detection,
            filter_high_occupancy_gems: !args.disable_high_occupancy_gem_detection,
        })
        .unwrap(),
    );

    if let Some(slide_serial_capture_area) = &args.slide_serial_capture_area {
        final_metrics.insert(
            SLIDE_SERIAL_CAPTURE_AREA.to_string(),
            Value::String(slide_serial_capture_area.to_string()),
        );
    }

    if let Some(sha) = args.multi_config_sha.as_deref() {
        final_metrics.insert("multi_config_sha".to_string(), sha.into());
    }

    Ok(Value::Object(final_metrics))
}

fn load_targets_per_lib(
    feature_reference: &FeatureReference,
    library_info: &[LibraryInfo],
) -> Result<Option<Vec<TxHashSet<u32>>>> {
    if feature_reference.target_set.is_none() {
        return Ok(None);
    }
    let target_indices = feature_reference
        .target_set
        .as_ref()
        .unwrap()
        .feature_indices();
    let mut targets_per_lib = Vec::new();
    for l in library_info {
        let targets = match &l.target_set_name {
            Some(_) => target_indices.clone(),
            None => TxHashSet::default(),
        };
        targets_per_lib.push(targets);
    }

    Ok(Some(targets_per_lib))
}

pub struct MoleculeInfoVisitor<'a> {
    sample_filtered_barcodes: TxHashSet<&'a Barcode>,
    targets_per_lib: Option<&'a [TxHashSet<u32>]>,
    metrics: MoleculeInfoMetrics,
}

#[derive(Default, JsonReport)]
pub struct MoleculeInfoMetrics {
    usable_read_pairs: TxHashMap<u16, CountMetric>,
    feature_read_pairs: TxHashMap<u16, CountMetric>,
    on_target_usable_read_pairs: TxHashMap<u16, CountMetric>,
}

impl<'a> MoleculeInfoVisitor<'a> {
    fn visit(&mut self, bc_umi_info: &BcUmiInfo) {
        let is_cell_barcode = self.sample_filtered_barcodes.contains(&bc_umi_info.barcode);

        for c in &bc_umi_info.umi_counts {
            // everything is "feature"
            self.metrics
                .feature_read_pairs
                .entry(c.library_idx)
                .or_default()
                .increment_by(c.read_count);

            if is_cell_barcode {
                self.metrics
                    .usable_read_pairs
                    .entry(c.library_idx)
                    .or_default()
                    .increment_by(c.read_count);

                if let Some(targets) = self.targets_per_lib {
                    let lib_targets = &targets[c.library_idx as usize];

                    if lib_targets.contains(&c.feature_idx) {
                        self.metrics
                            .on_target_usable_read_pairs
                            .entry(c.library_idx)
                            .or_default()
                            .increment_by(c.read_count);
                    }
                }
            }
        }
    }
}

/// Return the /barcode_info/pass_filter and /barcode_info/genomes datasets.
/// pass_filter has three columns: barcode, library, genome.
fn build_barcode_info(
    library_filtered_genome_barcodes: &[FilteredBarcodesCsvRow],
    library_info: &[cr_types::LibraryInfo],
    barcode_index: &BarcodeIndex,
    sample_filtered_barcodes: &TxHashSet<&Barcode>,
    genomes: &[GenomeName],
) -> ndarray::Array2<u64> {
    let genome_to_idx: HashMap<&GenomeName, usize> = genomes
        .iter()
        .enumerate()
        .map(|(genome_idx, genome)| (genome, genome_idx))
        .collect();

    let gem_group_to_libraries: HashMap<u16, Vec<usize>> = library_info
        .iter()
        .enumerate()
        .map(|(lib_idx, library)| match library {
            cr_types::LibraryInfo::Aggr(_) => unimplemented!(),
            cr_types::LibraryInfo::Count(library) => (library.gem_group, lib_idx),
        })
        .into_group_map();

    let barcode_library_genomes: Vec<[u64; 3]> = library_filtered_genome_barcodes
        .iter()
        .filter(|x| sample_filtered_barcodes.contains(&x.barcode as &Barcode))
        .flat_map(|barcode_genome| {
            let barcode: &Barcode = &barcode_genome.barcode;
            let barcode_idx = barcode_index.get_index(barcode);
            let genome_idx = genome_to_idx[&barcode_genome.genome];
            gem_group_to_libraries[&barcode.gem_group()]
                .iter()
                .map(move |&lib_idx| [barcode_idx as u64, lib_idx as u64, genome_idx as u64])
        })
        .sorted()
        .collect();
    ndarray::Array2::from(barcode_library_genomes)
}

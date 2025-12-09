#![deny(missing_docs)]
use crate::fit_piecewise_linear_model::{EstimatedModel, PiecewiseLinearData};
use crate::probe_barcode_matrix::ProbeCounts;
use anyhow::Result;
use barcode::Barcode;
use cr_h5::molecule_info::MoleculeInfoIterator;
use cr_types::probe_set::{ProbeRegion, ProbeSetReference};
use cr_types::reference::probe_set_reference::TargetSetFile;
use cr_types::types::PROBE_IDX_SENTINEL_VALUE;
use cr_types::utils::calculate_median_of_sorted;
use cr_types::{CountShardFile, ProbeBarcodeCount};
use itertools::Itertools;
use metric::{MeanMetric, TxHashMap, TxHashSet};
use ndarray::Array2;
use ndarray::prelude::*;
use serde::Serialize;
use shardio::ShardReader;
use std::collections::{BTreeMap, HashMap};
use std::path::Path;

/// Calculate the median of an unsorted iterator.
/// Return 0 if the iterator is empty.
fn calculate_median_of_unsorted(xs: impl Iterator<Item = usize>) -> usize {
    calculate_median_of_sorted(&xs.sorted().collect::<Vec<_>>()).unwrap_or(0)
}

#[derive(Serialize)]
pub struct GdnaCorrectedMetrics {
    mean_gdna_corrected_genes: f64,
    mean_gdna_corrected_umis: f64,
    median_gdna_corrected_genes: usize,
    median_gdna_corrected_umis: usize,
}

/// GDNA_GENE_THRESHOLD is the minimum number of genes with both spliced and
/// unspliced probes necessary to run gDNA analysis.
/// Also used in mro/rna/stages/targeted/disable_targeted_stages/__init__.py
const GDNA_GENE_THRESHOLD: usize = 10;
/// Struct returned to the stage. Has estimates and things needed to plot
pub struct MetricsComputed {
    pub estimated_gdna_per_probe: f64,
    pub estimated_percentage_of_gdna_umi: f64,
    pub unspliced_counts: Vec<f64>,
    pub spliced_counts: Vec<f64>,
    pub estimated_model: EstimatedModel<f64>,
}

/// Function computing gDNA metrics and returning the metrics currently
/// used.
pub fn compute_gdna_metrics(
    mol_info_path: &Path,
    probe_set_path: &TargetSetFile,
) -> Result<MetricsComputed> {
    // Reading the probe set reference
    let probe_set_reference: ProbeSetReference =
        ProbeSetReference::from_path(probe_set_path, None, 1)?;
    let probes = &probe_set_reference.sorted_probes();

    // Dictionary with number of spliced probes that every gene has
    let mut num_spliced_probes_per_gene = TxHashMap::default();
    // Dictionary with number of unspliced probes that every gene has
    let mut num_unspliced_probes_per_gene = TxHashMap::default();
    let mut number_of_unspliced_probes = 0.0;
    let mut total_number_of_umis = 0;

    // Going over the list of probes and figuring out the number of spliced
    // and unspliced probes each gene has
    for x in probes {
        if x.is_excluded_probe() {
            continue;
        }
        match &x.region {
            Some(ProbeRegion::Spliced) => {
                *num_spliced_probes_per_gene
                    .entry(x.gene.id.as_str())
                    .or_insert(0) += 1;
            }
            Some(ProbeRegion::Unspliced) => {
                number_of_unspliced_probes += 1.0;
                *num_unspliced_probes_per_gene
                    .entry(x.gene.id.as_str())
                    .or_insert(0) += 1;
            }
            _ => (),
        }
    }

    // Computing the set of genes which have both spliced and unspliced probes
    let set_of_genes_with_spliced_probes: TxHashSet<_> =
        num_spliced_probes_per_gene.keys().copied().collect();
    let set_of_genes_with_unspliced_probes: TxHashSet<_> =
        num_unspliced_probes_per_gene.keys().copied().collect();
    let set_of_genes_with_spliced_unspliced_probes: TxHashSet<_> = set_of_genes_with_spliced_probes
        .intersection(&set_of_genes_with_unspliced_probes)
        .copied()
        .collect();

    // Finally reading in the mol info
    let mol_info_iter = MoleculeInfoIterator::new(mol_info_path)
        .unwrap()
        .cell_barcodes_only(true)
        .unwrap();
    // Number of unspliced umis seen in a gene - averaged by the number of unspliced probes the gene has
    // Initialising it here. Only computing on genes with both spliced and unspliced probes.
    // Need to initialise as need things not seen as zero counts.
    let mut num_unspliced_umis_per_gene: TxHashMap<_, _> =
        set_of_genes_with_spliced_unspliced_probes
            .iter()
            .map(|&x| (x, 0.0))
            .collect();

    // Number of spliced umis seen in a gene - averaged by the number of unspliced probes the gene has
    // Initialising it here. Only computing on genes with both spliced and unspliced probes.
    // Need to initialise as need things not seen as zero counts.
    let mut num_spliced_umis_per_gene: TxHashMap<_, _> = set_of_genes_with_spliced_unspliced_probes
        .iter()
        .map(|&x| (x, 0.0))
        .collect();

    // Iterating over the mol info
    for full_umi_record in mol_info_iter {
        let Some(probe_idx) = full_umi_record.umi_data.probe_idx else {
            continue;
        };
        if probe_idx == PROBE_IDX_SENTINEL_VALUE {
            continue;
        }
        assert!(probe_idx >= 0);
        let probe = probes[probe_idx as usize];
        if probe.is_excluded_probe() {
            continue;
        }

        total_number_of_umis += 1;
        let gene_id = probe.gene.id.as_str();
        if !set_of_genes_with_spliced_unspliced_probes.contains(gene_id) {
            continue;
        }

        match &probes[probe_idx as usize].region {
            Some(ProbeRegion::Spliced) => {
                *num_spliced_umis_per_gene.get_mut(gene_id).unwrap() +=
                    1.0 / num_spliced_probes_per_gene[gene_id] as f64;
            }
            Some(ProbeRegion::Unspliced) => {
                *num_unspliced_umis_per_gene.get_mut(gene_id).unwrap() +=
                    1.0 / (num_unspliced_probes_per_gene[gene_id] as f64);
            }
            _ => (),
        }
    }

    // Sorting the genes which have both spliced and unspliced counts by name
    let mut sorted_control_genes: Vec<_> = set_of_genes_with_spliced_unspliced_probes
        .into_iter()
        .collect();
    sorted_control_genes.sort();

    // Getting the spliced and unspliced counts from the dicts in the same order
    // and log transforming them
    let control_log_spliced_umis: Vec<_> = sorted_control_genes
        .iter()
        .map(|&x| num_spliced_umis_per_gene[x].ln_1p())
        .collect();

    let control_log_unspliced_umis: Vec<_> = sorted_control_genes
        .iter()
        .map(|&x| num_unspliced_umis_per_gene[x].ln_1p())
        .collect();

    // fit the piecewise linear model
    let data_to_fit =
        PiecewiseLinearData::new(control_log_spliced_umis, control_log_unspliced_umis);
    let estimated_model = data_to_fit.fit();

    // Compute the two metrics needed
    let estimated_gdna_per_probe = estimated_model.model.constant.exp() - 1.0;
    let estimated_percentage_of_gdna_umi = total_number_of_umis
        .min((estimated_gdna_per_probe * number_of_unspliced_probes).round() as i64)
        as f64
        / total_number_of_umis as f64;

    Ok(MetricsComputed {
        unspliced_counts: data_to_fit.get_unspliced_counts(),
        spliced_counts: data_to_fit.get_spliced_counts(),
        estimated_gdna_per_probe,
        estimated_percentage_of_gdna_umi,
        estimated_model,
    })
}

/// Struct to store spliced and unspliced counts corresponding
/// to a gene. We use floats as the counts are averaged across
/// number of spliced and unspliced probes
#[derive(Default, Clone, Copy)]
struct SplicedUnsplicedCounts {
    spliced: f64,
    unspliced: f64,
}

/// Take in CSF files corresponding to UMI counts of every probe in every barcode
/// and a list of filtered barcodes. Then collates pseudobulk probe counts
/// over all filtered barcodes.
/// Will work for all probe sets.
fn collate_probe_metrics(
    probe_barcode_path: &[CountShardFile],
    filtered_barcodes: &TxHashSet<Barcode>,
    probe_set_path: &TargetSetFile,
    reference_path: Option<&Path>,
) -> Result<Vec<ProbeCounts>> {
    // Read in the probe set reference and get them in sorted order - the order of the probe_idx
    let psr: ProbeSetReference = ProbeSetReference::from_path(probe_set_path, reference_path, 1)
        .expect("ProbeSetReference could not be made from path");
    let probes: Vec<_> = psr.sorted_probes();
    let num_probes = probes.len();

    // Vector where we collate probe counts
    let mut umis_in_filtered_barcodes = vec![0; num_probes];
    let mut umis_in_all_barcodes = vec![0; num_probes];

    // Read in probe barcode shards
    let probe_barcode_reader: ShardReader<ProbeBarcodeCount> =
        ShardReader::open_set(probe_barcode_path)?;

    // Iterate over the probe-barcode records
    for count in probe_barcode_reader.iter()? {
        // Check if the barcode is a filtered barcode
        let count_unwrapped = count.as_ref().expect("Error reading probe barcode count");
        let probe_index = usize::try_from(count_unwrapped.probe_idx)?;
        // Count UMIs in all barcodes
        umis_in_all_barcodes[probe_index] += usize::try_from(count_unwrapped.umi_count).unwrap();
        if filtered_barcodes.contains(&count_unwrapped.barcode) {
            // Increment the probe count by the UMI count of the record
            umis_in_filtered_barcodes[probe_index] +=
                usize::try_from(count_unwrapped.umi_count).unwrap();
        }
    }

    // Returns a vector of probes and their counts.
    // Set `pass_filter` to `false` when it is an excluded probe.
    Ok(probes
        .into_iter()
        .enumerate()
        .map(|(probe_idx, probe)| ProbeCounts {
            umis_in_filtered_barcodes: umis_in_filtered_barcodes[probe_idx],
            umis_in_all_barcodes: umis_in_all_barcodes[probe_idx],
            probe_id: probe.probe_id.clone(),
            region: probe.region.clone(),
            included: probe.included,
            gene_id: probe.gene.id.clone(),
            gene_name: probe.gene.name.clone(),
            pass_filter: !probe.is_excluded_probe(),
            probe_idx,
            probe_type: probe.probe_type,
            ref_sequence_name: probe.ref_sequence_name.clone(),
            ref_sequence_pos: probe.ref_sequence_pos,
            cigar_string: probe.cigar_string.clone(),
            genome: probe.genome.clone(),
        })
        .collect())
}

/// Determine which probes are excluded by the gDNA filter and set `pass_filter` to `false`.
/// `pass_filter` is `true` when the probe
/// - is not deprecated and
/// - passes the gDNA filter
fn compute_filtered_probes(per_probe_metrics: Vec<ProbeCounts>) -> Result<Vec<ProbeCounts>> {
    // Dictionary with number of spliced probes that every gene has
    let mut num_spliced_probes_per_gene = TxHashMap::default();
    // Dictionary with number of unspliced probes that every gene has
    let mut num_unspliced_probes_per_gene = TxHashMap::default();

    // Going over the list of probes and figuring out the number of spliced
    // and unspliced probes each gene has
    for probe_record in &per_probe_metrics {
        if probe_record.is_excluded() {
            continue;
        }
        match &probe_record.region {
            Some(ProbeRegion::Spliced) => {
                *num_spliced_probes_per_gene
                    .entry(probe_record.gene_id.as_str())
                    .or_insert(0) += 1;
            }
            Some(ProbeRegion::Unspliced) => {
                *num_unspliced_probes_per_gene
                    .entry(probe_record.gene_id.as_str())
                    .or_insert(0) += 1;
            }
            _ => (),
        }
    }

    // Computing the set of genes which have both spliced and unspliced probes
    let set_of_genes_with_spliced_probes: TxHashSet<_> =
        num_spliced_probes_per_gene.keys().copied().collect();
    let set_of_genes_with_unspliced_probes: TxHashSet<_> =
        num_unspliced_probes_per_gene.keys().copied().collect();

    // Number of unspliced umis and spliced umis seen in a gene -
    // averaged by the number of unspliced and spliced probes the gene has
    // Initialising it here. Only computing on genes with both spliced and unspliced probes.
    // Need to initialise as need things not seen as zero counts.
    // Use a BTreeMap instead of a HashMap as it preserves the order in every iteration
    // for a given set of keys.
    let mut num_spliced_and_unspliced_umis_per_gene: BTreeMap<_, _> =
        set_of_genes_with_spliced_probes
            .intersection(&set_of_genes_with_unspliced_probes)
            .map(|&x| (x, SplicedUnsplicedCounts::default()))
            .collect();

    // Iterating over the probe records
    for probe_record in &per_probe_metrics {
        // Probe IDX cant be -ve as we have set the dtype to usize
        // PROBE_SENTINEL_VALUE is handled upstream
        if probe_record.is_excluded() {
            continue; // ignoring probes with prefixes we ignore
        }
        let gene_id = probe_record.gene_id.as_str();
        if let Some(num_spliced_and_unspliced) =
            num_spliced_and_unspliced_umis_per_gene.get_mut(gene_id)
        {
            // only considering genes with both spliced and unspliced probes
            match &probe_record.region {
                Some(ProbeRegion::Spliced) => {
                    num_spliced_and_unspliced.spliced += (probe_record.umis_in_filtered_barcodes
                        as f64)
                        / (num_spliced_probes_per_gene[gene_id] as f64);
                }
                Some(ProbeRegion::Unspliced) => {
                    num_spliced_and_unspliced.unspliced += (probe_record.umis_in_filtered_barcodes
                        as f64)
                        / (num_unspliced_probes_per_gene[gene_id] as f64);
                }
                _ => (),
            }
        }
    }

    // Getting the spliced and unspliced counts from the dicts in the same order
    // and log transforming them. We can do this over two different iterations of
    // num_spliced_and_unspliced_umis_per_gene as it is a BTreeMap (rather than a HashMap)
    // and thus preserves the order in which iteration is done
    let control_log_spliced_unspliced_umis: Vec<[f64; 2]> = num_spliced_and_unspliced_umis_per_gene
        .into_values()
        .map(|x| [x.spliced.ln_1p(), x.unspliced.ln_1p()])
        .collect();
    let control_log_spliced_unspliced_umis_to_pass: Array2<_> =
        control_log_spliced_unspliced_umis.into();

    // fit the piecewise linear model
    let data_to_fit = PiecewiseLinearData::new(
        control_log_spliced_unspliced_umis_to_pass
            .slice(s![.., 0])
            .to_vec(),
        control_log_spliced_unspliced_umis_to_pass
            .slice(s![.., 1])
            .to_vec(),
    );

    // Compute the two metrics needed
    let estimated_model = data_to_fit.fit();
    let estimated_gdna_per_probe = (estimated_model.model.constant.exp() - 1.0).ceil() as usize;

    // Determine which probes are excluded by the gDNA filter. Deprecated probes are also excluded.
    Ok(per_probe_metrics
        .into_iter()
        .map(|x| ProbeCounts {
            pass_filter: x.pass_filter
                && !(x.region == Some(ProbeRegion::Unspliced)
                    && x.umis_in_filtered_barcodes <= estimated_gdna_per_probe),
            ..x
        })
        .collect())
}

/// Takes in probe barcode counts, filtered barcodes, the probe sets and the refernce
/// Decides if gDNA analysis should be run
/// Collates probe metrics.
/// If gDNA analysis should be run, it is and per probe metrics are generated
pub fn get_filtered_per_probe_metrics(
    probe_barcode_counts: &[CountShardFile],
    filtered_barcodes: &TxHashSet<Barcode>,
    probe_set: &TargetSetFile,
    reference_path: Option<&Path>,
) -> Result<Vec<ProbeCounts>> {
    let raw_per_probe_metrics = collate_probe_metrics(
        probe_barcode_counts,
        filtered_barcodes,
        probe_set,
        reference_path,
    )?;

    // IDs of genes with spliced probes
    let genes_with_spliced_probes: TxHashSet<_> = raw_per_probe_metrics
        .iter()
        .filter(|&x| x.is_spliced())
        .map(|x| x.gene_id.as_str())
        .collect();
    // IDs of genes with unspliced probes
    let genes_with_unspliced_probes: TxHashSet<_> = raw_per_probe_metrics
        .iter()
        .filter(|&x| x.is_unspliced())
        .map(|x| x.gene_id.as_str())
        .collect();
    // Number of genes with spliced and unspliced probes
    let num_genes_with_spliced_and_unspliced_probes = genes_with_unspliced_probes
        .intersection(&genes_with_spliced_probes)
        .count();

    // We run the gDNA analysis only if there are at least GDNA_GENE_THRESHOLD genes
    // with spliced and unspliced probes
    if num_genes_with_spliced_and_unspliced_probes >= GDNA_GENE_THRESHOLD {
        compute_filtered_probes(raw_per_probe_metrics)
    } else {
        Ok(raw_per_probe_metrics)
    }
}

/// Given the count shards of probe-barcode counts, filtered probes and
/// filtered barcodes; go over the count shards and compute the median genes
/// in filtered barcodes that come only from filtered probes
pub fn compute_gdna_corrected_median_genes_per_spot(
    probe_barcode_counts: &[CountShardFile],
    filtered_per_probe_metrics: &[ProbeCounts],
    filtered_barcodes: &TxHashSet<Barcode>,
) -> Result<GdnaCorrectedMetrics> {
    // Get a map from filtered probes to the genes they correspond to
    let filtered_probe_to_gene_id: TxHashMap<_, _> = filtered_per_probe_metrics
        .iter()
        .filter(|x| x.pass_filter)
        .map(|x| (x.probe_idx, x.gene_id.as_str()))
        .collect();

    // Maintain two hashmaps to keep track of gene and UMI counts.
    // Initialize both to zero as we may have barcodes that are filtered but have
    // zero counts that we may not see in the count-shards.
    let mut bc_gene_counts: HashMap<&Barcode, usize> =
        filtered_barcodes.iter().map(|x| (x, 0)).collect();
    let mut bc_umi_counts = bc_gene_counts.clone();

    // Iterate through the count shards
    let probe_barcode_reader: ShardReader<ProbeBarcodeCount> =
        ShardReader::open_set(probe_barcode_counts)?;
    probe_barcode_reader.iter()?.process_results(|iter| {
        // group by barcodes. Note that the count-shards are expected to be sorted by barcode
        for (barcode, counts) in &iter.chunk_by(|x| x.barcode) {
            // If barcode is not a filtered barcode, ignore it.
            if !filtered_barcodes.contains(&barcode) {
                continue;
            }

            // Get a vector of counts corresponding to filtered probes in a barcode
            let filtered_counts: Vec<_> = counts
                .filter(|x| filtered_probe_to_gene_id.contains_key(&(x.probe_idx as usize)))
                .collect();
            // Compute the number of gDNA corrected genes in the barcode.
            let number_of_gdna_corrected_genes_in_bc = filtered_counts
                .iter()
                .map(|x| filtered_probe_to_gene_id[&(x.probe_idx as usize)])
                .unique()
                .count();
            // Compute the total UMIs in the barcode in filtered probes. These are
            // the gDNA corrected UMIs in the barcode
            let number_of_gdna_corrected_umis_in_bc: u32 =
                filtered_counts.into_iter().map(|x| x.umi_count).sum();
            // Update the gDNA corrected gene and UMI counts for the barcode
            *bc_gene_counts.get_mut(&barcode).unwrap() = number_of_gdna_corrected_genes_in_bc;
            *bc_umi_counts.get_mut(&barcode).unwrap() =
                number_of_gdna_corrected_umis_in_bc as usize;
        }
    })?;

    let mean_gdna_corrected_genes: MeanMetric =
        bc_gene_counts.values().map(|&x| x as i64).collect();
    let mean_gdna_corrected_umis: MeanMetric = bc_umi_counts.values().map(|&x| x as i64).collect();

    Ok(GdnaCorrectedMetrics {
        mean_gdna_corrected_genes: mean_gdna_corrected_genes.mean(),
        mean_gdna_corrected_umis: mean_gdna_corrected_umis.mean(),
        median_gdna_corrected_genes: calculate_median_of_unsorted(bc_gene_counts.into_values()),
        median_gdna_corrected_umis: calculate_median_of_unsorted(bc_umi_counts.into_values()),
    })
}

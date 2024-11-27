//! Probe barcode matrix I/O

use anyhow::{bail, Context, Result};
use barcode::Barcode;
use cr_h5::count_matrix::{write_barcodes_column, MAT_H5_BUF_SZ};
use cr_h5::feature_reference_io::{make_fixed_ascii, write_target_set_group, FA_LEN};
use cr_h5::{extend_dataset, make_column_ds, write_column_ds};
use cr_types::probe_set::{is_deprecated_probe, Probe, ProbeRegion, ProbeSetReference, ProbeType};
use cr_types::reference::feature_reference::TargetSet;
use cr_types::reference::probe_set_reference::TargetSetFile;
use cr_types::{BarcodeIndex, CountShardFile, ProbeBarcodeCount};
use hdf5::types::FixedAscii;
use hdf5::{File, Group};
use itertools::{process_results, Itertools};
use martian_filetypes::json_file::JsonFile;
use martian_filetypes::FileTypeRead;
use metric::TxHashSet;
use serde::{Deserialize, Serialize};
use shardio::ShardReader;
use std::collections::HashMap;
use std::path::Path;
use std::str::FromStr;
use transcriptome::Gene;

const MATRIX_GROUP: &str = "matrix";

/// HDF5 group name of the probe set reference
const PROBE_SET_REF_GROUP: &str = "features";
const PROBE_SET_FEATURE_TYPE_NAME: &str = "Gene Expression";

#[derive(Serialize, Deserialize, Clone)]
/// Struct of a probe_idx, probe and probe_count to write out in
/// per_probe_metrics.csv
pub struct ProbeCounts {
    pub probe_idx: usize,
    pub probe_id: String,
    pub gene_name: String,
    pub gene_id: String,
    pub included: bool,
    pub region: Option<ProbeRegion>,
    pub umis_in_all_barcodes: usize,
    pub umis_in_filtered_barcodes: usize,
    pub pass_filter: bool,
    pub probe_type: ProbeType,
    pub ref_sequence_name: String,
    pub ref_sequence_pos: Option<usize>,
    pub cigar_string: String,
}

impl ProbeCounts {
    /// Return true if this probe ID is excluded based on its probe ID or
    /// if the probe the included field is false.
    pub fn is_excluded(&self) -> bool {
        is_deprecated_probe(&self.probe_id) || !self.included
    }

    /// Return true if the Probe is spliced
    pub fn is_spliced(&self) -> bool {
        self.region == Some(ProbeRegion::Spliced)
    }

    /// Return true if the Probe is unspliced
    pub fn is_unspliced(&self) -> bool {
        self.region == Some(ProbeRegion::Unspliced)
    }
}

impl From<ProbeCounts> for Probe {
    /// Convert a ProbeCounts object to a Probe.
    fn from(probe: ProbeCounts) -> Probe {
        Probe {
            probe_id: probe.probe_id,
            gene: Gene {
                name: probe.gene_name,
                id: probe.gene_id,
            },
            included: probe.included,
            region: probe.region,
            probe_type: probe.probe_type,
            ref_sequence_name: probe.ref_sequence_name,
            ref_sequence_pos: probe.ref_sequence_pos,
            cigar_string: probe.cigar_string,
        }
    }
}

/// Write out a Probe Set Reference as a group into the h5 matrix
/// filtered_probes_set is a Set of Probe IDs which are filtered
/// The written out group is as follows
///     /probes
///         /gene_id         (ID of gene that the probe measures)
///         /gene_name       (Name of gene that the probe measures)
///         /genome          (Name of genome)
///         /id              (ID of the Probe)
///         /name            (Name of the Probe)
///         /probe_region    (Region of the probe: 'spliced' or 'unspliced'. In addition has 'other' or 'null')
///         /filtered_probes (boolean column which is TRUE if the probe is filtered)
fn write_probe_columns(
    group: &mut Group,
    psr: &ProbeSetReference,
    filtered_probe_set: &TxHashSet<String>,
    probe_set_name: &str,
) -> Result<()> {
    // Obtain probes from the probe set reference
    let probes: Vec<_> = psr.sorted_probes();

    // Reserve vectors that we will write out into the h5
    let mut gene_id_strings = Vec::with_capacity(probes.len());
    let mut gene_name_strings = Vec::with_capacity(probes.len());
    let mut probe_id_strings = Vec::with_capacity(probes.len());
    let mut probe_name_strings = Vec::with_capacity(probes.len());
    let mut probe_region_strings = Vec::with_capacity(probes.len());
    let mut feature_type_strings = Vec::with_capacity(probes.len());
    let mut probes_pass_filter = Vec::with_capacity(probes.len());
    let mut probes_included = Vec::with_capacity(probes.len());
    let mut genome_strings = Vec::with_capacity(probes.len());

    // Iterate over all the probes
    for probe in probes {
        // Get gene IDs
        gene_id_strings.push(make_fixed_ascii(&probe.gene.id)?);

        // Get Gene Names
        gene_name_strings.push(make_fixed_ascii(&probe.gene.name)?);

        // Feature Type in a string
        feature_type_strings.push(make_fixed_ascii(PROBE_SET_FEATURE_TYPE_NAME)?);

        // Get probe IDs and probe name
        let probe_id = probe.probe_id.as_str();
        probe_id_strings.push(make_fixed_ascii(probe_id)?);
        // The probe ID is GENE_ID|GENE_NAME|HASH. The probe name is GENE_NAME|HASH.
        let probe_name = probe_id
            .split_once('|')
            .map_or(probe_id, |(_gene_id, probe_name)| probe_name);
        probe_name_strings.push(make_fixed_ascii(probe_name)?);

        // See if the probe IDs passed the filter by checking if they are in the
        // filtered_probe_set
        probes_pass_filter.push(filtered_probe_set.contains(probe_id));

        // Vector of included field of probes
        probes_included.push(!probe.is_excluded_probe());

        // Get regions of the probe
        // 'other' and 'none' are handled separately
        // 'spliced' and 'unspliced' are the usual options
        let probe_region_string = match &probe.region {
            Some(x) => x.to_string(),
            None => String::from("none"),
        };
        probe_region_strings.push(make_fixed_ascii(&probe_region_string)?);

        // Put in the genome as another dataset
        genome_strings.push(make_fixed_ascii(&psr.genome)?);
    }

    // Write vectors to the H5 matrix
    write_column_ds::<FixedAscii<FA_LEN>>(group, "genome", &genome_strings)?;
    write_column_ds::<FixedAscii<FA_LEN>>(group, "gene_id", &gene_id_strings)?;
    write_column_ds::<FixedAscii<FA_LEN>>(group, "gene_name", &gene_name_strings)?;
    write_column_ds::<FixedAscii<FA_LEN>>(group, "id", &probe_id_strings)?;
    write_column_ds::<FixedAscii<FA_LEN>>(group, "name", &probe_name_strings)?;
    write_column_ds::<FixedAscii<FA_LEN>>(group, "feature_type", &feature_type_strings)?;
    write_column_ds::<FixedAscii<FA_LEN>>(group, "probe_region", &probe_region_strings)?;
    write_column_ds::<bool>(group, "filtered_probes", &probes_pass_filter)?;
    write_target_set_group(
        group,
        &TargetSet::from_bools(probe_set_name, &probes_included),
    )?;

    Ok(())
}

/// Save per probe UMI count data to the Cell Ranger HDF5 sparse matrix format.
/// `group` must be the h5 group named 'matrix' where the data will be deposited.
/// /matrix
///     /data (int32 count values)
///     /indices (int32 feature ids)
///     /indptr (int64 indexeinto data for each barcode)
///     /barcodes (barcode info)
///     /filtered_barcodes (boolean vec of barcodes above that are filtered)
///     /probes (probe info)
fn write_probe_matrix_h5_helper(
    group: &mut Group,
    reader: &ShardReader<ProbeBarcodeCount>,
    psr: &ProbeSetReference,
    barcode_index: &BarcodeIndex,
    filtered_probe_set: &TxHashSet<String>,
    filtered_barcodes: &[bool],
    probe_set_name: &str,
) -> Result<()> {
    // Avoid failure in writing a hdf5 dataset of length 0 with gzip chunk size 1
    if barcode_index.is_empty() {
        bail!(
            "No 10x barcodes were observed in the experiment. This is likely the consequence of a \
            sample mixup or very poor sequencing quality on the barcode bases. Further execution \
            is halted."
        );
    }

    // Write the Barcode column in to the h5
    write_barcodes_column(group, "barcodes", barcode_index.sorted_barcodes())?;
    write_column_ds::<bool>(group, "filtered_barcodes", filtered_barcodes)?;

    // Create a group of probes and write the probe info into the h5
    let mut probe_group = group.create_group(PROBE_SET_REF_GROUP)?;
    write_probe_columns(&mut probe_group, psr, filtered_probe_set, probe_set_name)?;

    // Create columns in the H5
    let data = make_column_ds::<i32>(group, "data")?;
    let indices = make_column_ds::<i64>(group, "indices")?;

    // Create vectors for the sparse matrix representations
    // These will be written out into the h5
    let mut data_buf = Vec::with_capacity(MAT_H5_BUF_SZ);
    let mut indices_buf = Vec::with_capacity(MAT_H5_BUF_SZ);
    // Barcode_counts is used to construct the indptr of the sparse matrix
    let mut barcode_counts = vec![0; barcode_index.len()];

    process_results(reader.iter()?, |iter| {
        for (barcode, counts) in &iter.group_by(|x| x.barcode) {
            // If not in the filtered barcodes, we dont process the the UMI
            if !barcode_index.contains_barcode(&barcode) {
                continue;
            }

            // Count all counts corresponding to a barcode
            // probe_idx is the index of the probe. Thus indices_buf maintains
            // the sequence of probe ids that are non-zero as we need in a sparse matrix
            // representation. data_buf maintains the data
            // n counts the number of nonzero probes we have for a barcode
            let mut n = 0;
            for count in counts {
                data_buf.push(count.umi_count as i32);
                indices_buf.push(count.probe_idx as i32);
                n += 1;
            }

            // barcode_counts maintain difference of adjacent indptr
            let barcode_index = barcode_index.get_index(&barcode);
            // Just making sure that the shards are indeed sorted by barcodes
            assert_eq!(barcode_counts[barcode_index], 0);
            barcode_counts[barcode_index] = n;

            if data_buf.len() > MAT_H5_BUF_SZ * 3 / 4 {
                extend_dataset::<i32>(&data, &data_buf)?;
                data_buf.clear();

                extend_dataset(&indices, &indices_buf)?;
                indices_buf.clear();
            }
        }
        // Writing the data and indices to the h5
        extend_dataset(&data, &data_buf)?;
        extend_dataset(&indices, &indices_buf)?;

        // Constructing the indptr from the barcode_counts
        let mut total = 0i64;
        let mut indptr = vec![0i64];
        indptr.extend(barcode_counts.iter().map(|x| {
            total += x;
            total
        }));
        // Writing the indptr to the h5
        write_column_ds(group, "indptr", &indptr)?;

        // Also putting in the shape of the matrix into the h5
        write_column_ds(
            group,
            "shape",
            &[psr.number_of_probes() as i32, barcode_counts.len() as i32],
        )?;
        anyhow::Ok(())
    })??;
    Ok(())
}

#[allow(clippy::too_many_arguments)]
/// Writing a probe barcode matrix in h5 format
/// All probes are written into the matrix, with filtered_probes being annotated
/// filtered_probes_set is a Set of the probe_ids of probes that are filtered
/// Uses a barcode index of barcodes to include in the H5 matrix
pub fn write_probe_bc_matrix(
    probe_set_path: &TargetSetFile,
    reference_path: Option<&Path>,
    probe_barcode_path: &[CountShardFile],
    h5_path: &Path,
    bc_index: &BarcodeIndex,
    filtered_probe_set: &TxHashSet<String>,
    filtered_barcodes: &[bool],
    probe_set_name: &str,
) -> Result<()> {
    let psr: ProbeSetReference = ProbeSetReference::from_path(probe_set_path, reference_path, 1)?;
    let probe_barcode_reader: ShardReader<ProbeBarcodeCount> =
        ShardReader::open_set(probe_barcode_path)?;
    let f = File::create(h5_path).with_context(|| h5_path.display().to_string())?;
    let mut group = f.create_group(MATRIX_GROUP)?;
    write_probe_matrix_h5_helper(
        &mut group,
        &probe_barcode_reader,
        &psr,
        bc_index,
        filtered_probe_set,
        filtered_barcodes,
        probe_set_name,
    )?;
    Ok(())
}

/// Function to read a sample_barcode json.
/// Expect to have a Hashmap of {sample: vec of barcodes}
pub fn read_bc_json(
    json_path: &JsonFile<HashMap<String, Vec<String>>>,
) -> Result<HashMap<String, Vec<Barcode>>> {
    let tmp_hashmap = json_path.read()?;
    let sample_barcodes_hashmap: HashMap<String, Vec<Barcode>> = tmp_hashmap
        .into_iter()
        .map(|(name, str_vc)| {
            (
                name,
                str_vc
                    .into_iter()
                    .map(|x| Barcode::from_str(&x).unwrap())
                    .collect(),
            )
        })
        .collect();
    Ok(sample_barcodes_hashmap)
}

#[cfg(test)]
mod tests {
    use super::*;
    use barcode::BcSeq;
    use std::fs::write;
    use tempfile::tempdir;

    fn barcode(seq: &[u8]) -> Barcode {
        Barcode::with_seq(1, BcSeq::from_bytes(seq), true)
    }

    #[test]
    /// Test function fo read_bc_json
    fn test_read_bc_json() {
        let dir = tempdir().expect("Can not create directory");
        let file_path = dir.path().join("tmp_json.json");
        let json_string = r#"
        {
            "Sample1":[
                "AAACTGGGTAAAGGCCATCCCAAC-1",
                "AAAGCGAAGACCATGGATCCCAAC-1"
            ],
            "Sample2":[
                "TTGAGCGAGTCAACTTAACGCCGA-1",
                "TTGAGGTAGGGATGAAAACGCCGA-1",
                "TTGAGGAGTTAGTGAGAACGCCGA-1"
            ]
        }"#;
        write(&file_path, json_string).expect("Unable to write file");

        let sample_barcode_expected = HashMap::from([
            (
                String::from("Sample1"),
                vec![
                    barcode(b"AAACTGGGTAAAGGCCATCCCAAC"),
                    barcode(b"AAAGCGAAGACCATGGATCCCAAC"),
                ],
            ),
            (
                String::from("Sample2"),
                vec![
                    barcode(b"TTGAGCGAGTCAACTTAACGCCGA"),
                    barcode(b"TTGAGGTAGGGATGAAAACGCCGA"),
                    barcode(b"TTGAGGAGTTAGTGAGAACGCCGA"),
                ],
            ),
        ]);

        let sample_barcodes_hashmap_read = read_bc_json(&JsonFile::from(file_path)).unwrap();
        dir.close().expect("Could not delete directory");
        assert_eq!(sample_barcodes_hashmap_read, sample_barcode_expected);
    }
}

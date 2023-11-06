use crate::barcode_counter::FilteredBarcodeFilter;
use cr_h5::molecule_info::MoleculeInfoIterator;
use cr_types::probe_set::ProbeSetReference;
use cr_types::target_panel_summary::{TargetPanelSummary, TargetPanelSummaryFormat};
use cr_types::types::PROBE_IDX_SENTINEL_VALUE;
use cr_types::TargetingMethod;
use itertools::Itertools;
use martian_filetypes::FileTypeRead;
use pyo3::prelude::*;
use std::path::PathBuf;

#[allow(clippy::type_complexity)]
#[pyfunction]
pub(crate) fn count_umis_per_probe(
    _py: Python<'_>,
    mol_info_path: PathBuf,
    target_panel_summary_path: PathBuf,
    reference_path: PathBuf,
) -> PyResult<Option<(Vec<String>, Vec<String>, Vec<String>, Vec<i32>)>> {
    let tps: TargetPanelSummary = TargetPanelSummaryFormat::from(&target_panel_summary_path)
        .read()
        .unwrap();
    if tps.targeting_method != TargetingMethod::TemplatedLigation {
        return Ok(None);
    }

    let mut valid_bc_filter = FilteredBarcodeFilter::new(&mol_info_path);

    let psr: ProbeSetReference =
        ProbeSetReference::from_path(&tps.target_panel_path, &reference_path, 1)
            .expect("ProbeSetReference could not be made from path");
    let probes: Vec<_> = psr.sorted_probes();

    if probes.is_empty() | (probes[0].region.is_none()) {
        return Ok(None);
    }

    // iterate over mol info records and increment probe counts for each library,
    // provided the molecule passes certain filters (if given)
    let mut result = vec![0; probes.len()];
    for x in MoleculeInfoIterator::new(&mol_info_path).unwrap() {
        if let Some(probe_idx) = x.umi_data.probe_idx {
            let bc_lib = (x.barcode_idx, x.umi_data.library_idx);
            if valid_bc_filter.is_cell_or_spot(bc_lib) & (probe_idx != PROBE_IDX_SENTINEL_VALUE) {
                result[probe_idx as usize] += 1;
            }
        }
    }

    let (probe_ids, gene_names, regions): (Vec<String>, Vec<String>, Vec<String>) = probes
        .into_iter()
        .cloned()
        .map(|x| {
            (
                x.probe_id,
                x.gene.id,
                x.region.map_or("None", |x| x.as_str()).to_owned(),
            )
        })
        .multiunzip();
    Ok(Some((probe_ids, gene_names, regions, result)))
}

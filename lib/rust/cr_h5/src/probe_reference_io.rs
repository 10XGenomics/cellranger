use crate::feature_reference_io::encode_ascii_xml;
use anyhow::{Context, Result};
use cr_types::probe_set::{Probe, ProbeRegion};
use hdf5::types::FixedAscii;
use hdf5::Group;
use itertools::izip;
use std::boxed::Box;
use std::iter::Iterator;
use transcriptome::Gene;

pub const PROBE_DATA_LEN: usize = 64;
const PROBE_ID: &str = "probe_id";
const PROBE_FEATURE_NAME: &str = "feature_name";
const PROBE_FEATURE_ID: &str = "feature_id";
const PROBE_REGION: &str = "region";
const FILTERED_PROBES_NAME: &str = "filtered_probes";
const INCLUDED_PROBES_NAME: &str = "included_probes";

pub fn to_h5(probes: &[Probe], filtered_probes: &[bool], group: &mut Group) -> Result<()> {
    let item = probes
        .first()
        .expect("Probes vector was empty")
        .region
        .as_ref();
    if item.is_some() {
        //if one probe has a region annotation, they should all have annotations.
        let regions: Vec<_> = probes
            .iter()
            .map(|probe| probe.region.as_ref().unwrap().to_string())
            .collect();
        mk_string_col(group, PROBE_REGION, regions)?;
    };

    let (ids, gene_names, gene_ids, includeds): (Vec<_>, Vec<_>, Vec<_>, Vec<_>) =
        itertools::multiunzip(probes.iter().map(|probe| {
            (
                probe.probe_id.clone(),
                probe.gene.name.clone(),
                probe.gene.id.clone(),
                probe.included,
            )
        }));
    mk_string_col(group, PROBE_ID, ids)?;
    mk_string_col(group, PROBE_FEATURE_NAME, gene_names)?;
    mk_string_col(group, PROBE_FEATURE_ID, gene_ids)?;
    mk_bool_col(group, FILTERED_PROBES_NAME, filtered_probes)?;
    mk_bool_col(group, INCLUDED_PROBES_NAME, &includeds)?;
    Ok(())
}

pub fn from_h5(group: &Group) -> Result<Vec<Probe>> {
    //region is only optionally present
    let ds = group.dataset(PROBE_REGION);
    let regions: Box<dyn Iterator<Item = Option<ProbeRegion>>> = if ds.is_ok() {
        let regions_arr = ds.unwrap().read_1d::<FixedAscii<PROBE_DATA_LEN>>()?;
        Box::new(regions_arr.into_iter().map(|x| Some(ProbeRegion::new(&x))))
    } else {
        Box::new(std::iter::repeat(None))
    };
    // Included is only optionally present.
    let includeds: Box<dyn Iterator<Item = bool>> =
        if let Ok(dataset) = group.dataset(INCLUDED_PROBES_NAME) {
            Box::new(dataset.read_1d::<bool>()?.into_iter())
        } else {
            Box::new(std::iter::repeat(true))
        };
    let ds = group.dataset(PROBE_ID)?;
    let ids = ds.read_1d::<FixedAscii<PROBE_DATA_LEN>>()?;
    let ids: Vec<String> = ids.into_iter().map(|x| x.to_string()).collect();

    let ds = group.dataset(PROBE_FEATURE_NAME)?;
    let gene_names = ds.read_1d::<FixedAscii<PROBE_DATA_LEN>>()?;
    let gene_names: Vec<String> = gene_names.into_iter().map(|x| x.to_string()).collect();

    let ds = group.dataset(PROBE_FEATURE_ID)?;
    let gene_ids = ds.read_1d::<FixedAscii<PROBE_DATA_LEN>>()?;
    let gene_ids = gene_ids.into_iter().map(|x| x.to_string());

    let output: Vec<_> = izip!(regions, includeds, ids, gene_names, gene_ids)
        .map(|(region, included, probe_id, gene_name, gene_id)| {
            let gene = Gene {
                id: gene_id,
                name: gene_name,
            };
            Probe {
                probe_id,
                gene,
                included,
                region,
            }
        })
        .collect();

    Ok(output)
}

/// Writes a bool column into the h5
fn mk_bool_col(group: &mut Group, name: &str, data: &[bool]) -> Result<()> {
    let ds = group
        .new_dataset::<bool>()
        .deflate(1)
        .shape((data.len(),))
        .create(name)?;
    ds.as_writer().write(&data)?;
    Ok(())
}

fn mk_string_col(group: &mut Group, name: &str, data: Vec<String>) -> Result<()> {
    let mut strings = Vec::with_capacity(data.len());

    for datapoint in data {
        let v = make_fixed_ascii(&datapoint)?;
        strings.push(v);
    }

    let ds = group
        .new_dataset::<FixedAscii<PROBE_DATA_LEN>>()
        .deflate(1)
        .shape((strings.len(),))
        .create(name)?;

    ds.as_writer().write(&strings)?;

    Ok(())
}

#[inline]
fn make_fixed_ascii(s: &str) -> Result<FixedAscii<PROBE_DATA_LEN>> {
    let s = encode_ascii_xml(s);
    let probe_err = format!(
        "Error saving probe information. All probe names and ids must be less than {PROBE_DATA_LEN} characters"
    );
    let res = FixedAscii::<PROBE_DATA_LEN>::from_ascii(s.as_bytes()).with_context(|| probe_err)?;
    Ok(res)
}

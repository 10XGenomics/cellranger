use crate::feature_reference_io::encode_ascii_xml;
use anyhow::{Context, Result};
use cr_types::probe_set::{Probe, ProbeRegion, ProbeType};
use hdf5::types::FixedAscii;
use hdf5::Group;
use std::boxed::Box;
use std::convert::AsRef;
use std::iter::Iterator;
use std::str::FromStr;
use transcriptome::Gene;

pub const PROBE_DATA_LEN: usize = 64;
const PROBE_ID: &str = "probe_id";
const PROBE_FEATURE_NAME: &str = "feature_name";
const PROBE_FEATURE_ID: &str = "feature_id";
const PROBE_REGION: &str = "region";
const FILTERED_PROBES_NAME: &str = "filtered_probes";
const INCLUDED_PROBES_NAME: &str = "included_probes";
const PROBE_TYPE: &str = "probe_type";
const PROBE_REF_NAME: &str = "ref_name";
const PROBE_REF_POS: &str = "ref_pos";
const PROBE_CIGAR: &str = "cigar";

pub fn to_h5(probes: &[Probe], filtered_probes: &[bool], group: &mut Group) -> Result<()> {
    let item = probes
        .first()
        .expect("Probes vector was empty")
        .region
        .as_ref();
    if item.is_some() {
        //if one probe has a region annotation, they should all have annotations.
        mk_string_col(group, PROBE_REGION, probes, |probe| {
            probe.region.as_ref().unwrap().to_string()
        })?;
    };

    mk_string_col(group, PROBE_ID, probes, |probe| probe.probe_id.clone())?;
    mk_string_col(group, PROBE_FEATURE_NAME, probes, |probe| {
        probe.gene.name.clone()
    })?;
    mk_string_col(group, PROBE_FEATURE_ID, probes, |probe| {
        probe.gene.id.clone()
    })?;
    mk_string_col(group, PROBE_TYPE, probes, |probe| {
        probe.probe_type.as_ref().to_owned()
    })?;
    mk_string_col(group, PROBE_REF_NAME, probes, |probe| {
        probe.ref_sequence_name.clone()
    })?;
    mk_string_col(group, PROBE_REF_POS, probes, |probe| {
        probe
            .ref_sequence_pos
            .map(|x| x.to_string())
            .unwrap_or_default()
    })?;
    mk_string_col(group, PROBE_CIGAR, probes, |probe| {
        probe.cigar_string.clone()
    })?;

    mk_bool_col(group, FILTERED_PROBES_NAME, filtered_probes)?;
    let includes: Vec<bool> = probes.iter().map(|probe| probe.included).collect();
    mk_bool_col(group, INCLUDED_PROBES_NAME, &includes)?;
    Ok(())
}

pub fn from_h5(group: &Group) -> Result<Vec<Probe>> {
    //region is only optionally present
    let ds = group.dataset(PROBE_REGION);
    let mut regions: Box<dyn Iterator<Item = Option<ProbeRegion>>> = if ds.is_ok() {
        let regions_arr = ds.unwrap().read_1d::<FixedAscii<PROBE_DATA_LEN>>()?;
        Box::new(regions_arr.into_iter().map(|x| Some(ProbeRegion::new(&x))))
    } else {
        Box::new(std::iter::repeat(None))
    };

    // Included is only optionally present.
    let mut includeds: Box<dyn Iterator<Item = bool>> =
        if let Ok(dataset) = group.dataset(INCLUDED_PROBES_NAME) {
            Box::new(dataset.read_1d::<bool>()?.into_iter())
        } else {
            Box::new(std::iter::repeat(true))
        };

    let ids = read_str_col(group, PROBE_ID)?;
    let mut gene_names = read_str_col(group, PROBE_FEATURE_NAME)?;
    let mut gene_ids = read_str_col(group, PROBE_FEATURE_ID)?;

    let mut probe_types = read_str_col(group, PROBE_TYPE)?
        .map(|x| ProbeType::from_str(x.as_ref()).unwrap_or_default());
    let mut ref_sequence_names = read_str_col(group, PROBE_REF_NAME)?;
    let mut ref_sequence_positions = read_str_col(group, PROBE_REF_POS)?;
    let mut cigar_strings = read_str_col(group, PROBE_CIGAR)?;

    let mut output = Vec::new();
    for probe_id in ids {
        let probe = Probe {
            probe_id,
            gene: Gene {
                id: gene_ids
                    .next()
                    .expect("Error while loading probeset from h5"),
                name: gene_names
                    .next()
                    .expect("Error while loading probeset from h5"),
            },
            included: includeds.next().unwrap(),
            region: regions.next().unwrap(),
            probe_type: probe_types.next().unwrap(),
            ref_sequence_name: ref_sequence_names.next().unwrap(),
            ref_sequence_pos: ref_sequence_positions
                .next()
                .map(|num| num.parse::<usize>().ok())
                .unwrap(),
            cigar_string: cigar_strings.next().unwrap(),
        };
        output.push(probe);
    }

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

fn mk_string_col(
    group: &mut Group,
    name: &str,
    probes: &[Probe],
    func: impl Fn(&Probe) -> String,
) -> Result<()> {
    let data: Vec<String> = probes.iter().map(func).collect();
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

fn read_str_col(group: &Group, name: &str) -> Result<Box<dyn Iterator<Item = String>>> {
    if let Ok(ds) = group.dataset(name) {
        Ok(Box::new(
            ds.read_1d::<FixedAscii<PROBE_DATA_LEN>>()?
                .map(std::string::ToString::to_string)
                .into_iter(),
        ))
    } else {
        Ok(Box::new(std::iter::repeat(String::new())))
    }
}

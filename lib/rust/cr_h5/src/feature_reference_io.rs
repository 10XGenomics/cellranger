use anyhow::{bail, Context, Result};
use cr_types::reference::feature_reference::{
    FeatureDef, FeatureReference, FeatureType, TargetSet, REQUIRED_FEATURE_TAGS,
};
use fastq_set::read_pair::WhichRead;
use hdf5::types::FixedAscii;
use hdf5::{Group, H5Type};
use itertools::Itertools;
use std::borrow::Cow;
use std::collections::HashMap;
use std::fmt::Write;
use std::str::FromStr;
use transcriptome::Gene;
/// Write the feature reference to an HDF5 group. For use in the matrix.h5 and molecule_info.h5
pub fn to_h5(featref: &FeatureReference, group: &mut Group) -> Result<()> {
    let mut all_tag_keys = vec!["genome"];

    let have_feature_bc_cols = featref.feature_defs.iter().any(|x| !x.sequence.is_empty());

    // what are the 'all_tag_keys'?
    if have_feature_bc_cols {
        all_tag_keys.push("read");
        all_tag_keys.push("pattern");
        all_tag_keys.push("sequence");
    }

    // write the custom tags that we don't have built-in support for
    let custom_tags: Vec<_> = featref
        .feature_defs
        .iter()
        .flat_map(|d| d.tags.keys())
        .map(String::as_str)
        .unique()
        .sorted()
        .collect();
    all_tag_keys.extend(&custom_tags);

    for tag in custom_tags {
        mk_string_col(featref, group, tag, |feat| {
            feat.tags.get(tag).map_or("", String::as_str)
        })?;
    }

    // write the list of tags that are present in this feature reference
    let all_tag_keys_h5: Vec<_> = all_tag_keys
        .iter()
        .map(|x| make_fixed_ascii(x))
        .try_collect()?;
    let feature_types: HashMap<_, _> = featref
        .feature_maps
        .keys()
        .map(|&x| (x, x.to_string()))
        .collect();

    group
        .new_dataset::<FixedAscii<FA_LEN>>()
        .shape((all_tag_keys.len(),))
        .create("_all_tag_keys")?
        .as_writer()
        .write(&all_tag_keys_h5)?;

    // write the 'always on' tag
    mk_string_col(featref, group, "feature_type", |feat| {
        &feature_types[&feat.feature_type]
    })?;
    mk_string_col(featref, group, "genome", |feat| &feat.genome)?;
    mk_string_col(featref, group, "id", |feat| &feat.id)?;
    mk_string_col(featref, group, "name", |feat| &feat.name)?;

    // write the option feature bc columns
    if have_feature_bc_cols {
        mk_byte_col(featref, group, "sequence", |feat| feat.sequence.as_bytes())?;
        mk_string_col(featref, group, "pattern", |feat| &feat.pattern)?;
        mk_string_col(featref, group, "read", |feat| {
            if feat.feature_type == FeatureType::Gene {
                ""
            } else {
                match feat.read {
                    WhichRead::R1 => "R1",
                    WhichRead::R2 => "R2",
                    WhichRead::I1 => "I1",
                    WhichRead::I2 => "I2",
                }
            }
        })?;
    }

    // write target set if we have one
    if let Some(ref ts) = featref.target_set {
        write_target_set_group(group, ts)?;
    }

    Ok(())
}

/// Write a TargetSet to a HDF5 group.
/// Creates a new group target_sets and puts a dataset with the
/// TargetSet name in it. The data set contains the ids of features kept
pub fn write_target_set_group(group: &mut Group, target_set: &TargetSet) -> Result<()> {
    let feature_id_array = ndarray::Array1::from(target_set.to_feature_indices_vec());
    let target_set_group = group.create_group("target_sets")?;
    target_set_group
        .new_dataset::<u32>()
        .shuffle()
        .deflate(1)
        .shape(feature_id_array.shape())
        .create(target_set.name())?
        .as_writer()
        .write(feature_id_array.view())?;
    Ok(())
}

#[allow(dead_code)]
fn mk_col<T: H5Type, F: (Fn(&FeatureDef) -> T)>(
    featref: &FeatureReference,
    group: &mut Group,
    name: &str,
    extract_field: F,
) -> Result<()> {
    let mut v = Vec::with_capacity(featref.feature_defs.len());

    for f in &featref.feature_defs {
        v.push(extract_field(f));
    }

    let ds = group.new_dataset::<T>().shape((v.len(),)).create(name)?;
    ds.as_writer().write(&v)?;
    Ok(())
}

fn mk_string_col<'a, F: (Fn(&'a FeatureDef) -> &'a str)>(
    featref: &'a FeatureReference,
    group: &mut Group,
    name: &str,
    extract_field: F,
) -> Result<()> {
    let mut strings = Vec::with_capacity(featref.feature_defs.len());

    for f in &featref.feature_defs {
        let s = extract_field(f);
        let v = make_fixed_ascii(s)?;
        strings.push(v);
    }

    let ds = group
        .new_dataset::<FixedAscii<FA_LEN>>()
        .shuffle()
        .deflate(1)
        .shape((strings.len(),))
        .create(name)?;
    ds.as_writer().write(&strings)?;
    Ok(())
}

fn mk_byte_col<'a, F: (Fn(&'a FeatureDef) -> &[u8])>(
    featref: &'a FeatureReference,
    group: &mut Group,
    name: &str,
    extract_field: F,
) -> Result<()> {
    let mut strings = Vec::with_capacity(featref.feature_defs.len());

    for f in &featref.feature_defs {
        let s = extract_field(f);
        let v = make_fixed_ascii_u8(s)?;
        strings.push(v);
    }

    let ds = group
        .new_dataset::<FixedAscii<FA_LEN>>()
        .shape((strings.len(),))
        .create(name)?;
    ds.as_writer().write(&strings)?;
    Ok(())
}

pub fn from_h5(group: &Group) -> Result<FeatureReference> {
    let id = require_string_dataset(group, "id")?;
    if id.is_empty() {
        bail!("Couldn't find 'id' dataset in molecule_info.h5 file")
    }

    let name = require_string_dataset(group, "name")?;
    let feature_types: Vec<FeatureType> =
        try_dataset(group, "feature_type", FeatureType::Gene, id.len())?;
    let genome: Vec<String> = try_string_dataset(group, "genome", "", id.len())?;

    // get feature-barcoding elements -- these are empty strings for
    let sequence: Vec<String> = try_string_dataset(group, "sequence", "", id.len())?;
    let read: Vec<WhichRead> = try_string_dataset(group, "read", "", id.len())?
        .into_iter()
        .map(|x| {
            if x.is_empty() {
                Ok(WhichRead::R2)
            } else {
                x.parse::<WhichRead>()
            }
        })
        .try_collect()?;
    let pattern: Vec<String> = try_string_dataset(group, "pattern", "", id.len())?;

    let all_tag_keys = require_string_dataset(group, "_all_tag_keys")?;

    // tags not in the list of 'known' tags
    let extra_tags = get_extra_tags(&all_tag_keys);

    let mut tag_datasets = HashMap::new();
    for c in &extra_tags {
        let ds = require_string_dataset(group, c)?;
        tag_datasets.insert(c, ds);
    }

    let mut features = Vec::new();

    for i in 0..id.len() {
        let mut tags = HashMap::new();
        for c in &extra_tags {
            let val = &tag_datasets[c][i];
            if !val.is_empty() {
                tags.insert(c.clone(), val.clone());
            }
        }

        let f = FeatureDef {
            index: i,
            id: id[i].clone(),
            name: name[i].clone(),
            genome: genome[i].as_str().into(),
            sequence: sequence[i].clone(),
            pattern: pattern[i].clone(),
            read: read[i],
            feature_type: feature_types[i],
            tags,
        };

        features.push(f);
    }

    let mut feature_maps = HashMap::new();
    feature_maps.insert(FeatureType::Gene, HashMap::new());
    let mut gene_to_index = HashMap::new();

    for def in &features {
        if def.feature_type != FeatureType::Gene {
            feature_maps
                .entry(def.feature_type)
                .or_insert_with(HashMap::new)
                .entry(def.sequence.clone())
                .or_insert_with(Vec::new)
                .push(def.index);
        }

        if def.feature_type == FeatureType::Gene {
            let gene = Gene {
                id: def.id.clone(),
                name: def.name.clone(),
            };
            gene_to_index.insert(gene, def.index);
        }
    }

    let target_set = group
        .group("target_sets")
        .ok()
        .map(|target_set_group| -> Result<_> {
            let sets = target_set_group.member_names()?;
            assert_eq!(sets.len(), 1, "SLFE only support a single target set!");
            let name = &sets[0];
            let feature_indices = {
                if let Ok(fids) = target_set_group.dataset(name)?.as_reader().read_1d::<u32>() {
                    fids.into_iter().collect()
                } else {
                    // fallback for older aggrs that generate these as strings
                    // max size of u32 in ascii is 10
                    target_set_group
                        .dataset(name)?
                        .as_reader()
                        .read_1d::<FixedAscii<10>>()?
                        .iter()
                        .map(|v| v.parse::<u32>())
                        .try_collect()?
                }
            };
            Ok(TargetSet::from_indices(name, feature_indices))
        })
        .transpose()?;

    Ok(FeatureReference {
        feature_defs: features,
        feature_maps,
        gene_to_index,
        target_set,
    })
}

pub fn encode_ascii_xml(s: &str) -> Cow<'_, str> {
    if s.chars().all(|c| c.is_ascii()) {
        return Cow::Borrowed(s);
    }
    let mut t = String::new();
    for c in s.chars() {
        if c.is_ascii() {
            t.push(c);
        } else {
            write!(&mut t, "&#{};", c as u32).unwrap();
        }
    }
    Cow::Owned(t)
}

pub const FA_LEN: usize = 256;
const FA_ERR: &str = "Error saving gene / feature information. All gene and features names and ids must be less than 256 characters";
#[inline]
pub fn make_fixed_ascii(s: &str) -> Result<FixedAscii<FA_LEN>> {
    let s = encode_ascii_xml(s);
    let res = FixedAscii::<FA_LEN>::from_ascii(s.as_bytes()).with_context(|| FA_ERR)?;
    Ok(res)
}

#[inline]
fn make_fixed_ascii_u8(s: &[u8]) -> Result<FixedAscii<FA_LEN>> {
    let s = encode_ascii_xml(std::str::from_utf8(s).with_context(|| FA_ERR)?);
    let res = FixedAscii::<FA_LEN>::from_ascii(s.as_bytes()).with_context(|| FA_ERR)?;
    Ok(res)
}

fn try_string_dataset(group: &Group, name: &str, default: &str, len: usize) -> Result<Vec<String>> {
    Ok(if group.link_exists(name) {
        group
            .dataset(name)?
            .read_1d::<FixedAscii<FA_LEN>>()?
            .into_iter()
            .map(|x| x.to_string())
            .collect()
    } else {
        vec![default.into(); len]
    })
}

fn try_dataset<T, D>(group: &Group, name: &str, default: D, len: usize) -> Result<Vec<T>>
where
    T: FromStr + Clone,
    D: Into<T>,
    <T as FromStr>::Err: std::error::Error + Send + Sync + 'static,
{
    Ok(if group.link_exists(name) {
        group
            .dataset(name)?
            .read_1d::<FixedAscii<FA_LEN>>()?
            .into_iter()
            .map(|x| x.parse::<T>())
            .try_collect()?
    } else {
        vec![default.into(); len]
    })
}

fn require_string_dataset(group: &Group, name: &str) -> Result<Vec<String>> {
    Ok(if group.link_exists(name) {
        group
            .dataset(name)?
            .read_1d::<FixedAscii<FA_LEN>>()?
            .into_iter()
            .map(|x| x.to_string())
            .collect()
    } else {
        bail!("Couldn't find '{name}' dataset in molecule_info.h5 file")
    })
}

fn get_extra_tags(all_tag_keys: &[String]) -> Vec<String> {
    all_tag_keys
        .iter()
        .filter(|&x| !REQUIRED_FEATURE_TAGS.contains(&x.as_str()))
        .cloned()
        .collect()
}

#[cfg(test)]
#[allow(dead_code)]
mod tests {
    use super::*;
    use cr_types::reference::feature_reference::FeatureReferenceFile;
    use hdf5::File;
    use pretty_assertions::assert_eq;
    use regex::Regex;
    use std::convert::TryFrom;
    use std::path::Path;
    use tempfile::NamedTempFile;
    use test_refdata::{refdata_available, refdata_path, showroom_available, showroom_path};

    /// Compare feature refs
    fn compare_feature_refs(f1: &FeatureReference, f2: &FeatureReference) {
        assert_eq!(f1.feature_defs.len(), f2.feature_defs.len(), "lens differ");
        for i in 0..f1.feature_defs.len() {
            assert_eq!(f1.feature_defs[i], f2.feature_defs[i]);
        }

        assert_eq!(f1.feature_maps, f2.feature_maps);
        assert_eq!(f1.gene_to_index.len(), f2.gene_to_index.len());
        assert_eq!(f1, f2);
    }

    #[test]
    fn gex_only_features() -> Result<()> {
        if !refdata_available() {
            return Ok(());
        }

        // read from reference path
        let fr_new = FeatureReference::from_paths(
            &refdata_path("GRCh38-3.0.0/"),
            None,
            None,
            None,
            None,
            None,
            None,
        )?;

        // read from legacy molelcule_info.h5 & compare (if possible)
        let mat_path = showroom_path("pbmc_1k_v2/filtered_feature_bc_matrix.h5");

        if showroom_available() {
            let file = File::open(mat_path)?;
            let fr_mat = from_h5(&file.group("matrix/features")?)?;
            compare_feature_refs(&fr_new, &fr_mat);
        }

        // write to h5
        let tmp = NamedTempFile::new_in(".")?;
        let out_file = File::create(tmp.path())?;
        let mut group = out_file.create_group("features")?;
        to_h5(&fr_new, &mut group)?;
        // read it back in
        let fr_roundtrip = from_h5(&group)?;

        compare_feature_refs(&fr_roundtrip, &fr_new);
        Ok(())
    }

    #[test]
    fn gex_abs_features() -> Result<()> {
        test_feature_io(
            "test/feature/pbmc_1k_protein_v3.csv",
            refdata_path("GRCh38-3.0.0/"),
            Some(
                showroom_path("pbmc_1k_protein_v3/filtered_feature_bc_matrix.h5")
                    .to_str()
                    .unwrap(),
            ),
        )
    }

    #[test]
    fn test_extra_cols_csv() -> Result<()> {
        test_feature_io(
            "test/feature/extra_commas.csv",
            refdata_path("GRCh38-3.0.0/"),
            None,
        )
    }

    fn test_feature_io(
        feature_csv: impl AsRef<Path>,
        reference: impl AsRef<Path>,
        legacy_mat_h5: Option<&str>,
    ) -> Result<()> {
        if !refdata_available() {
            return Ok(());
        }

        // read from reference path
        let fr_new = FeatureReference::from_paths(
            reference.as_ref(),
            Some(&FeatureReferenceFile::from(feature_csv.as_ref())),
            None,
            None,
            None,
            None,
            None,
        )?;

        // read from legacy molelcule_info.h5 & compare, if we have access to this file
        if let Some(mat_fn) = legacy_mat_h5 {
            let mat_path = Path::new(mat_fn);

            if mat_path.exists() {
                let file = File::open(mat_path)?;
                let fr_mat = from_h5(&file.group("matrix/features")?)?;
                compare_feature_refs(&fr_new, &fr_mat);
            }
        }

        // write to h5
        let tmp = NamedTempFile::new_in(".")?;
        let out_file = File::create(tmp.path())?;
        let mut group = out_file.create_group("features")?;
        to_h5(&fr_new, &mut group)?;

        // read it back in
        let fr_roundtrip = from_h5(&group)?;

        compare_feature_refs(&fr_roundtrip, &fr_new);
        Ok(())
    }

    #[test]
    fn gex_featuretest_features() -> Result<()> {
        test_feature_io(
            "test/feature/featuretest.csv",
            refdata_path("GRCh38-3.0.0/"),
            None,
        )
    }

    #[test]
    fn gex_crispr() -> Result<()> {
        if !refdata_available() {
            return Ok(());
        }

        // read from reference path
        let fr_orig = FeatureReference::from_paths(
            &refdata_path("GRCh38-3.0.0/"),
            Some(&FeatureReferenceFile::from(Path::new(
                "test/feature/crispr_features.csv",
            ))),
            None,
            None,
            None,
            None,
            None,
        )?;

        // write to h5
        let tmp = NamedTempFile::new_in(".")?;
        let out_file = File::create(tmp.path())?;
        let mut group = out_file.create_group("features")?;
        to_h5(&fr_orig, &mut group)?;

        // read it back in
        let fr_roundtrip = from_h5(&group)?;

        compare_feature_refs(&fr_roundtrip, &fr_orig);
        Ok(())
    }

    #[test]
    fn test_encode_ascii_xml() -> Result<()> {
        // mirrors test_(de|en)code_ascii_xml in lib/python/cellranger/test/test_hdf5.py
        fn decode(mut s: &str) -> Cow<'_, str> {
            let re = Regex::new(r"^&#([1-9]\d*);").unwrap();
            if re.is_match(s) {
                let mut t = String::new();
                while let Some(m) = re.captures(s) {
                    let m0 = m.get(0).unwrap();
                    let m1 = m.get(1).unwrap();
                    t.push_str(&s[0..m0.start()]);
                    let c = m1.as_str().parse::<u32>().unwrap();
                    t.push(char::try_from(c).unwrap());
                    s = &s[m0.end()..];
                }
                t.push_str(s);
                Cow::Owned(t)
            } else {
                Cow::Borrowed(s)
            }
        }

        assert_eq!("Hello, world!", encode_ascii_xml("Hello, world!"));
        assert_eq!("Hello, world!", decode(&encode_ascii_xml("Hello, world!")));
        assert_eq!("&#294;ello, world!", encode_ascii_xml("Ħello, world!"));
        assert_eq!("Ħello, world!", decode(&encode_ascii_xml("Ħello, world!")));

        Ok(())
    }
}

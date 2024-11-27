use anyhow::{anyhow, Result};
use barcode::Barcode;
use itertools::Itertools;
use martian_filetypes::json_file::JsonFile;
use martian_filetypes::LazyFileTypeIO;
use serde::{Deserialize, Serialize};
use std::cmp::Reverse;
use std::collections::{HashMap, HashSet};
use std::str::FromStr;
use vdj_ann::annotate::ContigAnnotation;
use vdj_types::{VdjChain, VdjRegion};

/// Represent an exact clonotype.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ExactClonotype {
    pub id: String,
    pub identity: Vec<ProductiveContig>,
    pub barcodes: Vec<String>,
}

#[derive(Hash, PartialEq, Eq, Debug)]
pub struct BarcodeContigUMI {
    pub barcode: Barcode,
    pub umi: usize,
}

/// Contig-level information for exact-clonotype grouping.
#[derive(Clone, Eq, PartialEq, Debug, Ord, PartialOrd, Serialize, Deserialize, Hash)]
pub struct ProductiveContig {
    pub chain: VdjChain,
    pub cdr3_nt: String,
    pub vdj_nt: String,
    pub c_ref_name: Option<String>,
    pub jc_delta: Option<i32>,
}

// TODO: lift these up to be inherent methods of ContigAnnotation.
pub trait XtraContigAnnotationMethods {
    fn get_jc_delta(&self) -> Option<i32>;
    fn get_vdj_seq(&self) -> Option<String>;
}

impl XtraContigAnnotationMethods for ContigAnnotation {
    fn get_vdj_seq(&self) -> Option<String> {
        if let (Some(v_gene), Some(j_gene)) =
            (self.get_region(VdjRegion::V), self.get_region(VdjRegion::J))
        {
            Some(self.sequence[v_gene.contig_match_start..j_gene.contig_match_end].to_string())
        } else {
            None
        }
    }

    /// Returns an i32 since J-region annotation end can overlap with C-region
    fn get_jc_delta(&self) -> Option<i32> {
        if let (Some(c_region), Some(j_region)) =
            (self.get_region(VdjRegion::C), self.get_region(VdjRegion::J))
        {
            Some(c_region.contig_match_start as i32 - j_region.contig_match_end as i32)
        } else {
            None
        }
    }
}

impl ProductiveContig {
    pub fn new(ann: &ContigAnnotation) -> Result<ProductiveContig> {
        Ok(ProductiveContig {
            chain: ann
                .chain_type()
                .ok_or_else(|| anyhow!("required chain_type not found"))?,
            cdr3_nt: ann
                .cdr3_seq
                .as_ref()
                .ok_or_else(|| anyhow!("required cdr3_nt not found"))?
                .to_string(),
            vdj_nt: ann
                .get_vdj_seq()
                .ok_or_else(|| anyhow!("required vdj_nt not found"))?,
            c_ref_name: ann.get_gene_name(VdjRegion::C).cloned(),
            jc_delta: ann.get_jc_delta(),
        })
    }
}

pub fn generate_exact_clonotypes(
    all_contigs_json: JsonFile<Vec<ContigAnnotation>>,
) -> Result<Vec<ExactClonotype>> {
    let mut prod_contigs_per_bc = HashMap::<String, Vec<ProductiveContig>>::new();
    for ann in all_contigs_json.lazy_reader()? {
        let ann: ContigAnnotation = ann?;
        if ann.is_cell && ann.is_productive() {
            prod_contigs_per_bc
                .entry(ann.barcode.clone())
                .or_default()
                .push(ProductiveContig::new(&ann)?);
        }
    }
    let mut exact_clonotypes_map: HashMap<Vec<ProductiveContig>, Vec<String>> = Default::default();
    for (bc, mut contigs) in prod_contigs_per_bc {
        contigs.sort(); // sort contigs within each barcode
        exact_clonotypes_map.entry(contigs).or_default().push(bc);
    }
    // sort the barcodes associated with each exact clonotype
    exact_clonotypes_map.values_mut().for_each(|bcs| bcs.sort());

    let exact_clonotypes_vec = exact_clonotypes_map
        .into_iter()
        .sorted_by(|(c0, _), (c1, _)| c0.len().cmp(&c1.len()).then_with(|| c0.cmp(c1)))
        .enumerate()
        .map(|(i, (contigs, bcs))| ExactClonotype {
            id: format!("exact_clonotype_{i}"),
            barcodes: bcs,
            identity: contigs,
        })
        .collect::<Vec<ExactClonotype>>();
    Ok(exact_clonotypes_vec)
}

pub fn generate_exact_contigs(
    all_contigs_json: &JsonFile<Vec<ContigAnnotation>>,
    cell_bcs: HashSet<String>,
) -> Result<HashMap<ProductiveContig, Vec<BarcodeContigUMI>>> {
    let mut bcs_per_prodcontig = HashMap::<ProductiveContig, Vec<BarcodeContigUMI>>::new();
    for ann in all_contigs_json.lazy_reader()? {
        let ann: ContigAnnotation = ann?;
        if ann.is_productive() && cell_bcs.contains(&ann.barcode) {
            bcs_per_prodcontig
                .entry(ProductiveContig::new(&ann)?)
                .or_default()
                .push(BarcodeContigUMI {
                    barcode: Barcode::from_str(&ann.barcode)?,
                    umi: ann.umi_count,
                });
        }
    }

    // sort the barcodes associated with each exact contig in descending order of umi counts
    bcs_per_prodcontig
        .values_mut()
        .for_each(|bcs| bcs.sort_by_key(|a| Reverse(a.umi)));

    Ok(bcs_per_prodcontig)
}

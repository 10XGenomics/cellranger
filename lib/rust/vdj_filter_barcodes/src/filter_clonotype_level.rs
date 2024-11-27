use crate::filter_log::{AsmCellFilter, FilterLogEntry, FilterLogger};
use serde::Serialize;
use std::collections::HashMap;
use strum::IntoEnumIterator;
use strum_macros::{Display, EnumIter};
use vdj_ann::annotate::ContigAnnotation;
use vdj_asm_utils::exact_clonotyping::ExactClonotype;

#[derive(Debug)]
pub struct WhitelistContamInfo {
    pub barcode: String,
    pub umis: usize,
}

#[derive(Serialize, EnumIter, Display)]
#[serde(rename_all = "snake_case")]
pub enum GelBeadBarcodeHalf {
    #[strum(to_string = "part_a")]
    PartA,
    #[strum(to_string = "part_b")]
    PartB,
}

impl WhitelistContamInfo {
    fn part(&self, which_part: &GelBeadBarcodeHalf) -> &str {
        // TODO (CELLRANGER-8643) exclude OH sequence from partA/B definition
        match which_part {
            GelBeadBarcodeHalf::PartA => &self.barcode[0..8],
            GelBeadBarcodeHalf::PartB => &self.barcode[8..16],
        }
    }
}

/// Construct a nested vector where each inner vector corresponds to an exact clonotype and contains
/// `WhitelistContamInfo` for each barcode. The `WhitelistContamInfo` includes the barcode and
/// total UMI counts associated with with productive contigs for that barcode.  
pub fn build_wlcontaminfo_per_exact_clonotype(
    contigs: &[ContigAnnotation],
    clonotypes: Vec<ExactClonotype>,
) -> Vec<Vec<WhitelistContamInfo>> {
    let contig_umis_per_bc = contigs
        .iter()
        .filter(|c| c.is_cell && c.is_productive())
        .fold(HashMap::new(), |mut acc, contig| {
            *acc.entry(&contig.barcode).or_insert(0) += contig.umi_count;
            acc
        });

    clonotypes
        .into_iter()
        .map(|c| {
            c.barcodes
                .iter()
                .map(|bc| WhitelistContamInfo {
                    barcode: bc.to_string(),
                    umis: *contig_umis_per_bc.get(bc).unwrap_or(&0),
                })
                .collect()
        })
        .collect()
}

pub fn whitelist_contamination_filter(
    barcodes: &Vec<WhitelistContamInfo>,
    mut filter_logger: Option<&mut FilterLogger>,
) -> Vec<String> {
    const GB_UMI_MULT: usize = 10;
    const GB_MIN_FRAC: f64 = 0.2;

    struct WLContaminants {
        barcode: String,
        max_group_umis: usize,
        umis: usize,
    }

    fn is_contam(group: &[&WhitelistContamInfo], exact_clonotype_size: f64) -> Vec<WLContaminants> {
        let max_umis = group.iter().map(|i| i.umis).max().unwrap();
        let contam_bcs: Vec<WLContaminants> = group
            .iter()
            .filter(|info| max_umis >= GB_UMI_MULT * info.umis)
            .map(|info| WLContaminants {
                barcode: info.barcode.clone(),
                max_group_umis: max_umis,
                umis: info.umis,
            })
            .collect();
        // TODO(CELLRANGER-8639): Re-evaluate the GB_MIN_FRAC threshold
        if contam_bcs.len() as f64 / exact_clonotype_size >= GB_MIN_FRAC {
            return contam_bcs;
        }
        Vec::new()
    }
    let exact_clonotype_size = barcodes.len() as f64;
    let mut filtered_bcs: Vec<String> = Vec::new();
    for bc_part in GelBeadBarcodeHalf::iter() {
        let mut same_part: HashMap<&str, Vec<&WhitelistContamInfo>> = HashMap::new();
        for bc in barcodes {
            same_part.entry(bc.part(&bc_part)).or_default().push(bc);
        }
        let contam: Vec<WLContaminants> = same_part
            .iter()
            .flat_map(|(_, v)| is_contam(v, exact_clonotype_size))
            .collect();
        for c in contam {
            filtered_bcs.push(c.barcode.clone());
            if let Some(ref mut logger) = filter_logger {
                logger.log(&FilterLogEntry::cell_calling(
                    c.barcode.clone(),
                    AsmCellFilter::GelBeadContamination {
                        shared_sequence: bc_part.to_string(),
                        max_group_umis: c.max_group_umis,
                        umis: c.umis,
                    },
                ));
            }
        }
    }
    filtered_bcs
}

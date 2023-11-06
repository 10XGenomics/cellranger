use crate::barcode_data::{BarcodeDataBrief, ContigChimeraData, ContigJunctionData};
use crate::filter_log::{AsmCellFilter, FilterLogEntry, FilterLogger, LowConfidenceReason};
use anyhow::Result;
use debruijn::dna_string::DnaString;
use std::cmp::{max, min};
use tenkit2::pack_dna::pack_bases_80;
use vdj_ann::annotate::ContigAnnotation;
use vdj_types::{VdjChain, VdjRegion};

pub fn confidence_filter(
    filter_params: &BarcodeFilteringParams,
    n50_n50_rpu: i32,
    low_confidence_reasons: &mut Vec<LowConfidenceReason>,
) -> bool {
    // Define confidence.
    // ◼ This is not carefully thought out.
    // ◼ Note that we could presumably improve by making calls after all
    // ◼ barcodes have been assembled.

    let mut ugly = false;
    if filter_params.num_productive_tra > 2
        || filter_params.num_productive_trb > 2
        || filter_params.num_good_contigs > 4
    {
        low_confidence_reasons.push(LowConfidenceReason::PutativeCellMultiplet {
            total_contigs: filter_params.num_good_contigs as usize,
            tra_trg_igh_contigs: filter_params.num_productive_tra as usize,
            trb_trd_igkl_contigs: filter_params.num_productive_trb as usize,
        });
        ugly = true;
    }
    if filter_params.nu3 < 3 && n50_n50_rpu > 2 {
        low_confidence_reasons.push(LowConfidenceReason::LowUmiSupport {
            n50_n50_rpu: n50_n50_rpu as usize,
            num_umis_min_3_reads: filter_params.nu3 as usize,
        });
        ugly = true;
    }

    if (filter_params.max_junct_supp <= 1
        && (filter_params.nu3 < 4 || filter_params.num_good_contigs > 2))
        || (filter_params.min_junct_supp <= 1 && filter_params.numn < 3)
    {
        low_confidence_reasons.push(LowConfidenceReason::LowJunctionSupport {
            min_junction_support_umis: filter_params.min_junct_supp as usize,
            max_junction_support_umis: filter_params.max_junct_supp as usize,
            n50_n50_rpu: n50_n50_rpu as usize,
            num_umis_min_3_reads: filter_params.nu3 as usize,
            num_umis_min_rpu_frac_reads: filter_params.numn as usize,
            total_contigs: filter_params.num_good_contigs as usize,
        });
        ugly = true;
    }
    !ugly
}

// Definition of cell.  There are separate conditions for TCR and BCR.
// The reason for having a second restriction for BCR is
// complicated.  On the dataset 90103_400_26_98, but presumably not restricted
// to that, we observed that a large set of barcodes were called cells, but had
// tiny UMI counts in the matched GEX run 90139, and therefore probably arose
// from GEMs having no cell.  Moreover, these barcodes seemed to correlate
// with ubiquitous sequences, possibly originating in plasma cells.  The second
// restriction on BCR is designed to address this.

#[allow(clippy::too_many_arguments)]
pub fn cell_filter(
    filter_params: &BarcodeFilteringParams,
    bc: &BarcodeCellInfo,
    denovo: bool,
    is_tcr: bool,
    is_bcr: bool,
    n50_n50_rpu: i32,
    mut filter_logger: Option<&mut FilterLogger>,
    low_confidence_reasons: Vec<LowConfidenceReason>,
) -> bool {
    const MIN_XUCOUNTS: usize = 3;
    const MIN_TOTAL_UCOUNTS: usize = 10;

    let mut is_cell = false;

    if is_tcr || denovo {
        is_cell = bc.xucounts.len() >= MIN_XUCOUNTS;
        if !is_cell {
            if let Some(ref mut logger) = filter_logger {
                logger.log(&FilterLogEntry::cell_calling(
                    bc.barcode.clone(),
                    AsmCellFilter::NotEnoughUmisTcrOrDenovo {
                        num_surviving_umis: bc.xucounts.len(),
                        param_min_num_surviving_umis: MIN_XUCOUNTS,
                    },
                ));
            }
        }
    }
    if is_bcr && !denovo {
        is_cell = bc.xucounts.len() >= MIN_XUCOUNTS && bc.total_ucounts >= MIN_TOTAL_UCOUNTS;
        if !is_cell {
            if let Some(ref mut logger) = filter_logger {
                logger.log(&FilterLogEntry::cell_calling(
                    bc.barcode.clone(),
                    AsmCellFilter::NotEnoughUmisBcr {
                        num_surviving_umis: bc.xucounts.len(),
                        param_min_num_surviving_umis: MIN_XUCOUNTS,
                        total_umis: bc.total_ucounts,
                        param_min_total_umis: MIN_TOTAL_UCOUNTS,
                    },
                ));
            }
        }
    }

    // Added condition.  To be a cell, there must be a contig having a V annotation.
    if filter_params.num_good_contigs + filter_params.num_reject_contigs == 0
        || (!denovo && !filter_params.have_v)
    {
        is_cell = false;
        if let Some(ref mut logger) = filter_logger {
            logger.log(&FilterLogEntry::cell_calling(
                bc.barcode.clone(),
                AsmCellFilter::NoContigWithVRegion {},
            ));
        }
    }

    // Added condition.  If a barcode has exactly one good contig,
    // and it has junction support one, then do not call it a cell.
    if filter_params.num_good_contigs == 1 && filter_params.max_junct_supp <= 1 {
        is_cell = false;
        if let Some(ref mut logger) = filter_logger {
            logger.log(&FilterLogEntry::cell_calling(
                bc.barcode.clone(),
                AsmCellFilter::NotEnoughJunctionSupport {},
            ));
        }
    }

    // Added condition.  Have to see at least one good confident contig to call a cell.
    if !bc.high_confidence || filter_params.num_good_contigs == 0 {
        is_cell = false;
        if let Some(ref mut logger) = filter_logger {
            logger.log(&FilterLogEntry::cell_calling(
                bc.barcode.clone(),
                AsmCellFilter::NoConfidentContig {
                    reasons: low_confidence_reasons,
                },
            ));
        }
    }

    // Added condition.
    let nun = bc.xucounts.len();
    if nun == 0 || (bc.xucounts[nun - 1] as f64) < 0.03 * n50_n50_rpu as f64 {
        is_cell = false;
        if let Some(ref mut logger) = filter_logger {
            logger.log(&FilterLogEntry::cell_calling(
                bc.barcode.clone(),
                AsmCellFilter::NotEnoughReadsPerUmi {
                    max_umi_reads: *bc.xucounts.last().unwrap_or(&0) as usize,
                    n50_n50_rpu: n50_n50_rpu as usize,
                },
            ));
        }
    }
    is_cell
}

#[derive(Default, Clone)]
pub struct Contigs {
    pub good_contigs: Vec<ContigAnnotation>,
    pub reject_contigs: Vec<ContigAnnotation>,
}

impl Contigs {
    pub fn build_chimdata(&self, is_cell: bool, denovo: bool) -> Result<Vec<ContigChimeraData>> {
        let mut chimdata = Vec::<ContigChimeraData>::new();
        let mut all_contigs = self.good_contigs.clone();
        all_contigs.extend(self.reject_contigs.clone());
        if !denovo {
            for contig in &all_contigs {
                if let Some(cdr3_seq) = &contig.cdr3_seq {
                    let v = if let Some(v_region) = contig.get_region(VdjRegion::V) {
                        v_region.feature.feature_id as i32
                    } else {
                        -1_i32
                    };
                    if v >= 0 {
                        chimdata.push(ContigChimeraData {
                            barcode: contig.barcode.to_string(),
                            cdr3: DnaString::from_acgt_bytes(cdr3_seq.as_bytes())
                                .clone()
                                .to_bytes(),
                            v_ref_id: v as usize,
                            umi_count: contig.umi_count,
                            productive: contig.productive.unwrap() && is_cell,
                        })
                    }
                }
            }
        }
        Ok(chimdata)
    }

    pub fn build_jundata(&self, high_confidence: bool) -> Result<Vec<ContigJunctionData>> {
        let mut jundata = Vec::<ContigJunctionData>::new();
        for contig in &self.good_contigs {
            const JREGION: usize = 80;
            let (j_stop, is_igh) = if let Some(j_region) = contig.get_region(VdjRegion::J) {
                let contig_match_len = j_region.contig_match_end - j_region.contig_match_start;
                let full_match = j_region.annotation_match_start + contig_match_len
                    == j_region.annotation_length;
                let min_length = j_region.contig_match_start + contig_match_len >= JREGION;
                match (full_match, min_length) {
                    (true, true) => (
                        Some(j_region.contig_match_end),
                        Some(j_region.feature.chain == VdjChain::IGH),
                    ),
                    (_, _) => (None, None),
                }
            } else {
                (None, None)
            };
            if let Some(j_stop) = j_stop {
                let tig = DnaString::from_acgt_bytes(contig.sequence.as_bytes()).to_ascii_vec();
                let p = j_stop - JREGION;
                let mut jxn_seq = [0_u8; 20];
                pack_bases_80(array_ref![tig, p, 80], &mut jxn_seq);
                let supp = min(65535, contig.junction_support.as_ref().unwrap().umis) as u16;
                jundata.push(ContigJunctionData {
                    jxn_seq,
                    umis: supp,
                    high_confidence,
                    is_igh: is_igh.unwrap(),
                });
            }
        }
        Ok(jundata)
    }
}

#[derive(PartialEq, Eq, PartialOrd, Ord, Clone, Serialize)]
pub struct BarcodeCellInfo {
    // the barcode
    pub barcode: String,

    // total number of umis
    pub total_ucounts: usize,

    // surviving nonsolo ucounts,
    pub xucounts: Vec<i32>,

    // do we have a pairing?
    pub paired: bool,

    // high confidence?
    pub high_confidence: bool,

    // is it a cell?
    pub now_a_cell: bool,

    pub jundata: Vec<ContigJunctionData>,

    pub chimdata: Vec<ContigChimeraData>,
}

impl From<BarcodeDataBrief> for BarcodeCellInfo {
    fn from(src: BarcodeDataBrief) -> Self {
        BarcodeCellInfo {
            barcode: src.barcode,
            total_ucounts: src.total_ucounts,
            xucounts: src.xucounts,
            paired: false,
            high_confidence: false,
            now_a_cell: false,
            jundata: Vec::<ContigJunctionData>::new(),
            chimdata: Vec::<ContigChimeraData>::new(),
        }
    }
}

#[derive(Default)]
pub struct BarcodeFilteringParams {
    // the barcode
    pub barcode: String,

    // num nonsolo surviving UMIs with >= 3 reads
    pub nu3: i32,

    // maximum reads supporting a nonsolo surviving UMI
    pub nmax: i32,

    // num nonsolo surviving UMIs with reads >= 0.5x n50_n50_rpu
    pub numn: i32,

    // lowest num of UMIs supporting the junction region of a good contig
    pub min_junct_supp: i32,

    // highest num of UMIs supporting the junction region of a good contig
    pub max_junct_supp: i32,

    // num good contigs
    pub num_good_contigs: i32,

    // num reject contigs
    pub num_reject_contigs: i32,

    // num productive TRA contigs
    // always 0 in denovo mode
    pub num_productive_tra: i32,

    // num productive TRB contigs
    // always 0 in denovo mode
    pub num_productive_trb: i32,

    // any contig with V annotation
    // always false in denovo mode
    pub have_v: bool,

    // is paired?
    // always false in denovo mode
    pub paired: bool,
}

impl BarcodeFilteringParams {
    pub fn build(
        contigs: &Contigs,
        bc: &BarcodeCellInfo,
        denovo: bool,
        n50_n50_rpu: i32,
        gd_mode: bool,
    ) -> Result<BarcodeFilteringParams> {
        // Reads per UMI related parameters
        let (mut nu3, mut nmax, mut numn) = (0, 0, 0);
        for x in &bc.xucounts {
            if *x >= 3 {
                nu3 += 1;
            }
            nmax = max(nmax, *x);
        }
        for x in &bc.xucounts {
            if *x as f64 >= 0.05_f64 * n50_n50_rpu as f64 {
                numn += 1;
            }
        }
        // Junction support related parameters
        let jsupp_umis = &contigs
            .good_contigs
            .iter()
            .map(|c| c.junction_support.as_ref().unwrap().umis)
            .collect::<Vec<i32>>();
        let (mut max_junct_supp, mut min_junct_supp) = (0, 1000000000);
        for umi in jsupp_umis {
            min_junct_supp = min(min_junct_supp, *umi);
            max_junct_supp = max(max_junct_supp, *umi);
        }
        // Num contigs
        let num_good_contigs = contigs.good_contigs.len();
        let num_reject_contigs = contigs.reject_contigs.len();

        // Num TRA/TRB and check for V-region
        let mut is_tra = Vec::new();
        let mut have_v = false;
        if !denovo {
            for contig in &contigs.good_contigs {
                let this_is_tra = match (
                    contig.get_region(VdjRegion::J).unwrap().feature.chain,
                    gd_mode,
                ) {
                    (VdjChain::IGH | VdjChain::TRA, _) => true,
                    (VdjChain::TRG, true) => true,
                    (_, _) => false,
                };
                is_tra.push(this_is_tra);
                if contig.get_region(VdjRegion::V).is_some() {
                    have_v = true;
                }
            }
            for contig in &contigs.reject_contigs {
                if contig.get_region(VdjRegion::V).is_some() {
                    have_v = true;
                }
            }
        }
        let num_productive_tra = is_tra.iter().filter(|&&b| b).count() as i32;
        let num_productive_trb = is_tra.iter().filter(|&&b| !b).count() as i32;

        // Pairing determination
        let mut paired = false;
        if num_productive_tra == 1 && num_productive_trb == 1 {
            paired = true;
            if n50_n50_rpu > 2 && max_junct_supp == 1 && nu3 < 4 {
                paired = false;
            }
            if min_junct_supp == 1 && numn < 3 {
                paired = false;
            }
            if n50_n50_rpu > 2 && nu3 < 3 {
                paired = false;
            }
        }
        Ok(BarcodeFilteringParams {
            barcode: bc.barcode.clone(),
            nu3,
            nmax,
            numn,
            min_junct_supp,
            max_junct_supp,
            num_good_contigs: num_good_contigs as i32,
            num_reject_contigs: num_reject_contigs as i32,
            num_productive_tra,
            num_productive_trb,
            have_v,
            paired,
        })
    }
}

use anyhow::Result;
use martian_derive::MartianStruct;
use martian_filetypes::json_file::JsonFile;
use martian_filetypes::lz4_file::Lz4;
use martian_filetypes::{LazyFileTypeIO, LazyWrite};
use serde::Deserialize;

#[derive(Clone, Copy, Serialize, Deserialize, Debug, MartianStruct)]
pub struct FilterSwitch {
    /// Turn on/off filters in the assembler that makes use of shared contigs to detect potential artifacts
    /// - Controls NonDominantJunction, WeakJunction, CommonCloneShadow and CommonCloneShadowSingleUmi
    pub asm_shared_contig: bool,
    /// Turn on/off enclone filters that makes use of shared contigs across clonotypes to detect
    /// potential artifacts
    /// - Sets NGRAPH_FILTER
    pub enclone_shared_contig: bool,
    /// Filters that remove putative doublets
    /// - Sets NWEAK_CHAINS, NFOURSIE_KILL, NDOUBLET, NSIG
    pub enclone_multiplet: bool,
    /// Turn on/off UMI thresholds in enclone
    /// - Sets NUMI, NUMI_RATIO
    pub enclone_umi: bool,
}

impl Default for FilterSwitch {
    fn default() -> Self {
        Self {
            asm_shared_contig: true,
            enclone_shared_contig: true,
            enclone_umi: true,
            enclone_multiplet: true,
        }
    }
}

pub struct FilterLogger {
    writer: <Lz4<JsonFile<Vec<FilterLogEntry>>> as LazyFileTypeIO<FilterLogEntry>>::LazyWriter,
}

impl FilterLogger {
    pub fn new(file: &Lz4<JsonFile<Vec<FilterLogEntry>>>) -> Result<Self> {
        Ok(FilterLogger {
            writer: file.lazy_writer()?,
        })
    }
    pub fn log(&mut self, entry: &FilterLogEntry) {
        self.writer.write_item(entry).unwrap();
    }
}

#[derive(Serialize, Deserialize)]
#[serde(tag = "category", content = "info")]
#[serde(rename_all = "snake_case")]
pub enum FilterLogEntry {
    CellCalling {
        barcode: String,
        filter: AsmCellFilter,
    },
}

impl FilterLogEntry {
    pub fn cell_calling(barcode: String, filter: AsmCellFilter) -> Self {
        FilterLogEntry::CellCalling { barcode, filter }
    }
}

#[derive(Serialize, Deserialize)]
#[serde(tag = "reason", content = "details")]
#[serde(rename_all = "snake_case")]
pub enum LowConfidenceReason {
    /// If we observe more than 2 conitgs of the same kind or observe more than 4 good contigs
    /// in a barcode, flag all the contigs as low confidence as it is likely to be a cell multiplet
    PutativeCellMultiplet {
        total_contigs: usize,
        tra_trg_igh_contigs: usize,
        trb_trd_igkl_contigs: usize,
    },
    /// If the N50 N50 reads per UMI is more than 2, we need to observe at least 3 surviving non-solo UMIs
    /// with a minimum of 3 reads in the barcode to be confident about the contigs.
    ///
    /// - `n50_n50_rpu`: For each barcode compute N50 reads per UMI. Compute the N50 of these values
    /// - `num_umis_min_3_reads`: No. of surviving non-solo UMIs with a minimum of 3 reads
    LowUmiSupport {
        n50_n50_rpu: usize,
        num_umis_min_3_reads: usize,
    },
    /// Junction segment is defined as 100 bases ending where the right end of a J region aligns to the contig
    ///
    /// - `min_junction_support_umis`: Minimum UMI support across the junction
    /// - `max_junction_support_umis`: Maximum UMI support across the junction
    /// - `n50_n50_rpu`: For each barcode compute N50 reads per UMI. Compute the N50 of these values
    /// - `num_umis_min_3_reads`: No. of surviving non-solo UMIs with a minimum of 3 reads
    /// - `num_umis_min_rpu_frac_reads`: No. of surviving non-solo UMIs with a minumum of `0.05*n50_n50_rpu` reads
    ///
    /// We flag all contigs in a barcode as low confidence if any of these conditions are met:
    /// - max_junction_support_umis <= 1 && (num_umis_min_3_reads < 4 || total_contigs > 2)
    /// - min_junction_support_umis <= 1 && num_umis_min_rpu_frac_reads < 3
    LowJunctionSupport {
        min_junction_support_umis: usize,
        max_junction_support_umis: usize,
        n50_n50_rpu: usize,
        num_umis_min_3_reads: usize,
        num_umis_min_rpu_frac_reads: usize,
        total_contigs: usize,
    },
}

#[derive(Serialize, Deserialize)]
#[serde(tag = "name", content = "details")]
#[serde(rename_all = "snake_case")]
pub enum AsmCellFilter {
    /// In TCR or denovo mode, each barcode needs at least `param_min_num_surviving_umis` surviving
    /// non-solo UMIs to be called as a cell.
    ///
    /// Definition of surviving umis.
    ///
    /// 1. Find "strong" paths in the graph that have an annotation, or in the denovo case,
    ///    match to a primer.
    /// 2. Find all the edges on a strong path.
    /// 3. Find all the reads on such a good edge.
    /// 4. Find the umis for these reads.
    /// 5. Restrict to those umis for which 50% of their kmers are contained in
    ///    good edges.
    /// 6. In the non-denovo case, if no strong path had a V annotation, kill all the
    ///    surviving umis.
    ///
    /// Non-solo => The UMI has more than 1 read or the N50 reads per UMI is 2
    ///
    NotEnoughUmisTcrOrDenovo {
        num_surviving_umis: usize,
        param_min_num_surviving_umis: usize,
    },

    /// In BCR each barcode needs at least `param_min_num_surviving_umis` surviving
    /// non-solo UMIs and have a total of at least `param_min_total_umis` UMIs to
    /// be called as a cell.
    NotEnoughUmisBcr {
        num_surviving_umis: usize,
        param_min_num_surviving_umis: usize,
        total_umis: usize,
        param_min_total_umis: usize,
    },

    /// To be a cell, there must be a contig having a V annotation. Not applicable in
    /// denovo.
    NoContigWithVRegion {},

    /// If a barcode has exactly one contig, and it has junction support one,
    /// then do not call it a cell.
    NotEnoughJunctionSupport {},

    /// Have to see a confident contig to call a cell.
    NoConfidentContig {
        /// Added after 7.0
        #[serde(default)]
        reasons: Vec<LowConfidenceReason>,
    },

    /// The maximum number of reads among the surviving UMIs needs to be at least
    /// 3% that of the n50_n50_rpu.
    /// - `n50_n50_rpu`: For each barcode compute N50 reads per UMI. Compute the N50 of these
    ///    values
    NotEnoughReadsPerUmi {
        max_umi_reads: usize,
        n50_n50_rpu: usize,
    },

    /// Consider the set of contigs sharing the same junction segment (100
    /// bases ending where the right end of a J region aligns to the contig). Let:
    ///
    /// - `cluster_size`: Number of contigs in the set
    /// - `cluster_median_junction_umis`: The median UMIs supporting the junction
    /// - `dominant_junction_umis`: Maximum UMIs supporting the junction
    /// - `dominant_contig`: Name of the contig with `dominant_junction_umis`
    ///
    /// If `cluster_size` is smaller than `param_min_cluster_size` and `cluster_median_junction_umis`
    /// is less than `param_max_median_junction_umis`, we deem barcodes with contigs whose `junction_umis`
    /// is such that `dominant_junction_umis >= param_min_umi_ratio * junction_umis` as non cell
    /// associated because it appears to arise from some sort of leakage from plasma cells.
    NonDominantJunction {
        contig: String,
        junction_umis: usize,
        param_min_umi_ratio: usize,
        dominant_contig: String,
        dominant_junction_umis: usize,
        cluster_size: usize,
        param_min_cluster_size: usize,
        cluster_median_junction_umis: u16,
        param_max_median_junction_umis: u16,
    },

    /// Consider pairs of barcodes sharing the same junction. If the following conditions are
    /// staisfied:
    /// - Junction UMIs for one of the barcodes ("dominant barcode") is greater than `param_min_dominant_umis`
    /// - Junction UMIs for the other barcode ("weak barcode") is 1
    /// - The dominant barcode has >=2 high confidence contigs
    /// - The weak barcode has >=3 high confidence contigs
    /// - There is only one common chain between the two cells allowing for some mutation
    ///
    /// then, we deem the "weak barcode" as non cell associated
    ///
    WeakJunction {
        contig: String,
        param_min_dominant_umis: usize,
        dominant_contig: String,
        dominant_junction_umis: usize,
    },

    /// Consider the the V segments that appear in the set of contigs sharing the same
    /// CDR3 nucleotide. If any V segment has collective UMI support (sum of UMIs across all
    /// barcodes) at least `param_chimera_ratio` times greater than another, then the barcodes
    /// containing weaker contigs are deemed non cell associated.
    ChimericContig {
        cdr3_nt: String,
        param_chimera_ratio: usize,
        contig_v_region_id: usize,
        dominant_v_region_id: usize,
        dominant_v_region_umis: usize,
    },

    /// Consider a barcode where we observe 2 or more high confidence contigs. Let
    /// - `J`: Set of junctions associated with these contigs
    /// - `multiplicity`: Number of barcodes where we observe the same set of junctions J.
    /// - `max_multiplicity`: The maximum number of barcodes that share a given set of (2 or more) high
    ///   confidence contigs where at least one contig has a junction in `J`.
    ///
    /// If `multiplicity <= param_max_kill` and `max_multiplicity >= param_min_ratio_big * multiplicity`, then
    /// we deem the barcode as non cell associated.
    ///
    /// Note: For barcodes with exactly 2 high confidence contigs, we account for SHMs
    /// in this check by making sure that the non-matching junctions differ in at least 10 bases.
    CommonCloneShadow {
        multiplicity: usize,
        max_multiplicity: usize,
        param_max_kill: usize,
        param_min_ratio_big: usize,
    },

    /// Same as `CommonCloneShadow`, but only applicable to barcodes with exactly 2 high confidence
    /// contigs, and have only a single umi supporting the junction of one of the contigs. In such a
    /// case we also do not account for SHMs
    CommonCloneShadowSingleUmi {
        multiplicity: usize,
        max_multiplicity: usize,
        param_max_kill: usize,
        param_min_ratio_big: usize,
    },

    /// The barcode is not called as a cell because its overhang does not match
    /// the overhang ids specified for this sample. Only applies to multiplexed VDJ.
    DifferentOverhang { overhang_id: String },
}

impl AsmCellFilter {
    pub fn label(&self) -> &'static str {
        match self {
            AsmCellFilter::NotEnoughUmisTcrOrDenovo { .. } => "LOW_UMI",
            AsmCellFilter::NotEnoughUmisBcr { .. } => "LOW_UMI",
            AsmCellFilter::NoContigWithVRegion {} => "NO_V_REGION",
            AsmCellFilter::NotEnoughJunctionSupport {} => "LOW_JUNCTION_SUPPORT",
            AsmCellFilter::NoConfidentContig { .. } => "NO_CONF_CONTIG",
            AsmCellFilter::NotEnoughReadsPerUmi { .. } => "LOW_RPU",
            AsmCellFilter::NonDominantJunction { .. } => "NON_DOMINANT_JUNCTION",
            AsmCellFilter::WeakJunction { .. } => "WEAK_JUNCTION",
            AsmCellFilter::ChimericContig { .. } => "CHIMERIC",
            AsmCellFilter::CommonCloneShadow { .. } => "COMMON_CLONE",
            AsmCellFilter::CommonCloneShadowSingleUmi { .. } => "COMMON_CLONE",
            AsmCellFilter::DifferentOverhang { .. } => "DIFFERENT_OVERHANG",
        }
    }
}

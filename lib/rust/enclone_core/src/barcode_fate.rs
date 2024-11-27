use serde::{Deserialize, Serialize};

/// Different reasons why a barcode which have productive contig(s)
/// are not called as cells by enclone
///
/// For more explanation, see <https://10xgenomics.github.io/enclone/pages/auto/default_filters.html>
/// and <https://10xgenomics.github.io/enclone/pages/auto/help.special.html>
#[derive(Serialize, Deserialize, Clone, Copy)]
pub enum BarcodeFate {
    Doublet,
    WeakChains,
    /// The barcode was not called as cell by the assembler
    NotAsmCell,
    FoursieKill,
    NotGexCell,
    Signature,
    /// Find and mark for deletion exact subclonotypes having a variant base in V..J that,
    /// accounting for all the cells in all the exact subclonotypes, never occurs as Q60
    /// doesn't occur as Q40 twice, and disagrees with the reference.
    Qual,
    Umi,
    UmiRatio,
    /// If a V..J segment appears in exactly one dataset, with frequency n, let x be the total
    /// number of productive pairs for that dataset, and let y be the total number of productive
    /// pairs for all datasets from the same origin.  If (x/y)^n <= 10^-6, i.e. the probability
    /// that assuming even distribution, all instances of that V..J ended up in that one dataset,
    /// delete all the productive pairs for that V..J segment that do not have at least 100
    /// supporting UMIs.  (Note no attempt to do Bonferroni correction.)
    ///
    /// For the case of two datasets for one origin, with equal numbers of productive pairs in
    /// each, this corresponds roughly to the case n = 20.
    ///
    /// Note that we could modify this to allow *some* occurrences in other datasets.
    ///
    /// There are only certain ways that these misdistribution events could happen:
    ///
    /// 1. A cell (and particularly a plasma cell or plasmablast) bursts after drawing cells to
    ///    make libraries, leaving behind cell fragments that seed separate GEMs
    ///    (probably most likely).
    /// 2. Multiple gel beads end up in one GEM.
    /// 3. Something involving like cells sticking together and subsequently separating.
    /// 4. Physical contamination of libraries.
    /// 5. Informatic mixup of libraries.
    /// 6. Nothing other than a low probability event (unlikely).
    ///
    /// Note that in case 1, we have evidence that a plasma cell or plasmablast existed in the
    /// original cells that were drawn (perhaps breaking up in the process of drawing), and was
    /// subsequently distintegrated.
    Cross,
    /// Filter out exact subclonotypes having more than one chain, but all of the same type.
    /// For example, the filter removes all exact subclonotypes having two TRA chains and
    /// no other chains
    Improper,
    GraphFilter,
    /// No productive contigs for this barcode. This will only happen
    /// when certain default filters are turned off
    NonProductive,
}

impl BarcodeFate {
    pub fn label(&self) -> &'static str {
        match self {
            BarcodeFate::Doublet => "DOUBLET",
            BarcodeFate::WeakChains => "WEAK_CHAINS",
            BarcodeFate::NotAsmCell => "CELL",
            BarcodeFate::FoursieKill => "FOURSIE_KILL",
            BarcodeFate::NotGexCell => "GEX",
            BarcodeFate::Signature => "SIGNATURE",
            BarcodeFate::Qual => "QUAL",
            BarcodeFate::Umi => "UMI",
            BarcodeFate::UmiRatio => "UMI_RATIO",
            BarcodeFate::Cross => "CROSS",
            BarcodeFate::Improper => "IMPROPER",
            BarcodeFate::GraphFilter => "GRAPH_FILTER",
            BarcodeFate::NonProductive => "PRODUCTIVE",
        }
    }
}

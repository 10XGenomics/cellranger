//! barcode_extensions::stubs
#![expect(missing_docs)]

use crate::select_barcode_corrector_common;
use anyhow::Result;
use barcode::{BarcodeConstruct, BarcodeCorrector, BcSegSeq, Whitelist, WhitelistSpec};
use cr_types::chemistry::BarcodeExtraction;
use metric::SimpleHistogram;
use std::ops::Range;

pub type CorrectionMap = ();
pub struct CorrectionMapBuilder;

impl CorrectionMapBuilder {
    pub fn build(&self, _: &SimpleHistogram<BcSegSeq>) -> Result<CorrectionMap> {
        unimplemented!();
    }

    pub fn new(_: &Whitelist) -> Self {
        unimplemented!();
    }
}

/// Select the appropriate barcode corrector for this configuration.
pub fn select_barcode_corrector(
    input: BarcodeConstruct<((&WhitelistSpec, Whitelist), SimpleHistogram<BcSegSeq>)>,
    barcode_extraction: Option<&BarcodeExtraction>,
    correction_map: Option<BarcodeConstruct<CorrectionMap>>,
) -> BarcodeConstruct<(BarcodeCorrector, Option<Range<usize>>)> {
    select_barcode_corrector_common(input, barcode_extraction, correction_map)
}

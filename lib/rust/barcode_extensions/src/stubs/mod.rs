use anyhow::Result;
use barcode::corrector::Posterior;
use barcode::{BarcodeConstruct, BarcodeCorrector, BcSegSeq, Whitelist};
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

pub fn select_barcode_corrector(
    input: BarcodeConstruct<(Whitelist, SimpleHistogram<BcSegSeq>)>,
    barcode_extraction: Option<&BarcodeExtraction>,
    correction_map: Option<BarcodeConstruct<CorrectionMap>>,
) -> BarcodeConstruct<(BarcodeCorrector, Option<Range<usize>>)> {
    assert!(correction_map.is_none());

    input.map(|(wl, bc_counts)| match barcode_extraction {
        Some(BarcodeExtraction::Independent) | None => (
            BarcodeCorrector::new(wl, bc_counts, Posterior::default()),
            None,
        ),
        Some(BarcodeExtraction::JointBc1Bc2 { .. }) => {
            unimplemented!();
        }
    })
}

//! barcode_extensions
#![deny(missing_docs)]

use barcode::corrector::Posterior;
use barcode::{
    BarcodeConstruct, BarcodeCorrector, BcSegSeq, CorrectSubNotIndel, Whitelist, WhitelistSpec,
};
use cr_types::chemistry::BarcodeExtraction;
use metric::SimpleHistogram;
use std::ops::Range;

#[cfg(feature = "tenx_internal")]
mod internal;
#[cfg(feature = "tenx_internal")]
pub use internal::*;
#[cfg(feature = "tenx_source_available")]
mod stubs;
#[cfg(feature = "tenx_source_available")]
pub use stubs::*;

/// Return a basic barcode corrector.
fn basic_corrector(
    wl: Whitelist,
    bc_counts: SimpleHistogram<BcSegSeq>,
) -> (BarcodeCorrector, Option<Range<usize>>) {
    (
        BarcodeCorrector::new(wl, bc_counts, Posterior::default()),
        None,
    )
}

/// Select the appropriate barcode corrector for this configuration.
fn select_barcode_corrector_common(
    input: BarcodeConstruct<((&WhitelistSpec, Whitelist), SimpleHistogram<BcSegSeq>)>,
    barcode_extraction: Option<&BarcodeExtraction>,
    correction_map: Option<BarcodeConstruct<CorrectionMap>>,
) -> BarcodeConstruct<(BarcodeCorrector, Option<Range<usize>>)> {
    assert!(correction_map.is_none());
    match input {
        BarcodeConstruct::GelBeadOnly(((spec, whitelist), counts)) => {
            assert!(barcode_extraction.is_none());
            match spec.whitelist_name() {
                Some("737K-flex-v2" | "737K-flex-v2-pd1") => BarcodeConstruct::GelBeadOnly((
                    BarcodeCorrector::new(whitelist, counts, CorrectSubNotIndel),
                    None,
                )),
                Some(_) | None => BarcodeConstruct::GelBeadOnly(basic_corrector(whitelist, counts)),
            }
        }
        BarcodeConstruct::GelBeadAndProbe(barcode_construct) => {
            BarcodeConstruct::GelBeadAndProbe(match barcode_extraction {
                None => barcode_construct
                    .map(|((_spec, whitelist), counts)| basic_corrector(whitelist, counts)),
                Some(BarcodeExtraction::VariableMultiplexingBarcode { .. }) => barcode_construct
                    .map(|((_spec, whitelist), counts)| {
                        (
                            BarcodeCorrector::new(whitelist, counts, CorrectSubNotIndel),
                            None,
                        )
                    }),
                Some(BarcodeExtraction::JointBc1Bc2 { .. }) => unreachable!(),
            })
        }
        BarcodeConstruct::Segmented(..) => unreachable!(),
    }
}

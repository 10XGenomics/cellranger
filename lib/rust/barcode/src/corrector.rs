//!
//! Corrects sequencing errors in barcodes, allowing up to one mismatch.
//!
use crate::whitelist::Whitelist;
use crate::{BarcodeSegment, BarcodeSegmentState, BcSegQual, BcSegSeq};
use metric::SimpleHistogram;

const BC_MAX_QV: u8 = 66; // This is the illumina quality value
pub(crate) const BASE_OPTS: [u8; 4] = [b'A', b'C', b'G', b'T'];
type Of64 = ordered_float::NotNan<f64>;

/// Implement the standard 10x barcode correction algorithm.
/// Requires the barcode whitelist (`whitelist`), and the observed counts
/// of each whitelist barcode, prior to correction (`bc_counts`).
pub struct BarcodeCorrector {
    whitelist: Whitelist,
    bc_counts: SimpleHistogram<BcSegSeq>,
    strategy: Box<dyn CorrectBarcode + Send + Sync>,
}

impl BarcodeCorrector {
    /// Load a barcode corrector from a whitelist path
    /// and a set of barcode count files.
    /// # Arguments
    /// * `whitelist` - Barcode whitelist
    /// * `bc_counts` - Prior counts of each barcode
    /// * `correction_strategy` - The barcode correction strategy to use
    pub fn new<S>(whitelist: Whitelist, bc_counts: SimpleHistogram<BcSegSeq>, strategy: S) -> Self
    where
        S: CorrectBarcode + Send + Sync + 'static,
    {
        Self {
            whitelist,
            bc_counts,
            strategy: Box::new(strategy),
        }
    }

    pub fn whitelist(&self) -> &Whitelist {
        &self.whitelist
    }

    /// Attempt to correct a non-whitelist barcode sequence onto the whitelist.
    ///
    /// For a description of the correction methods, see comments in BarcodeCorrectionMethod
    /// If there is no correction made, return None
    /// Otherwise return the distance to the corrected barcode
    pub fn correct_barcode(
        &self,
        observed_segment: &mut BarcodeSegment,
        qual: Option<BcSegQual>,
    ) -> Option<u16> {
        match self.correct_barcode_helper(*observed_segment, qual) {
            Some((segment, dist)) => {
                *observed_segment = segment;
                Some(dist)
            }
            None => None,
        }
    }

    // See comments in BarcodeCorrectionMethod
    fn correct_barcode_helper(
        &self,
        observed_segment: BarcodeSegment,
        qual: Option<BcSegQual>,
    ) -> Option<(BarcodeSegment, u16)> {
        self.strategy
            .correct_barcode(&self.whitelist, &self.bc_counts, observed_segment, qual)
    }
}

pub trait CorrectBarcode {
    fn correct_barcode(
        &self,
        whitelist: &Whitelist,
        bc_counts: &SimpleHistogram<BcSegSeq>,
        observed_segment: BarcodeSegment,
        qual: Option<BcSegQual>,
    ) -> Option<(BarcodeSegment, u16)>;
}

const BARCODE_CONFIDENCE_THRESHOLD: f64 = 0.975;

/// Compute a posterior distribution over whitelist barcodes given:
/// # A prior distribution over whitelist bacodes
/// # The likelihood of P(orig_barcode | wl_bc; qual) of the `observed_barcode` sequence
///   given a hidden whitelist barcode wl_bc, and the sequencer reported quality values qual.
/// Use these components to compute a posterior distribution of wl_bc candidates, and correct
/// the barcode if there's a whitelist barcode with posterior probability greater than
/// bc_confidence_threshold. If no quality scores are supplied then each base has the same
/// probability of error.
pub struct Posterior {
    /// threshold for sum of probability of error on barcode QVs. Barcodes exceeding
    /// this threshold will be marked as not valid.
    max_expected_barcode_errors: f64,
    /// if the posterior probability of a correction
    /// exceeds this threshold, the barcode will be corrected.
    bc_confidence_threshold: f64,
}

impl Default for Posterior {
    fn default() -> Self {
        Self {
            max_expected_barcode_errors: f64::MAX,
            bc_confidence_threshold: BARCODE_CONFIDENCE_THRESHOLD,
        }
    }
}

impl CorrectBarcode for Posterior {
    fn correct_barcode(
        &self,
        whitelist: &Whitelist,
        bc_counts: &SimpleHistogram<BcSegSeq>,
        observed_segment: BarcodeSegment,
        qual: Option<BcSegQual>,
    ) -> Option<(BarcodeSegment, u16)> {
        let mut a = observed_segment.sequence().seq().to_owned(); // Create a copy
        assert!(observed_segment.state == BarcodeSegmentState::Invalid);

        let mut best_option: Option<(Of64, BarcodeSegment)> = None;
        let mut total_likelihood = Of64::try_from(0.0).unwrap();

        for pos in 0..a.len() {
            let qv = qual.map_or(BC_MAX_QV, |q| q[pos].min(BC_MAX_QV));
            let existing = a[pos];
            for val in BASE_OPTS {
                if val == existing {
                    continue;
                }
                a[pos] = val;
                let mut trial_bc = BarcodeSegment::with_sequence(&a, observed_segment.state);

                if whitelist.check_and_update(&mut trial_bc) {
                    // Apply additive (Laplace) smoothing.
                    let raw_count = bc_counts.get(trial_bc.sequence());
                    let bc_count = 1 + raw_count;
                    let prob_edit = Of64::try_from(probability(qv)).unwrap();
                    let likelihood = prob_edit * Of64::try_from(bc_count as f64).unwrap();
                    match best_option {
                        None => best_option = Some((likelihood, trial_bc)),
                        Some(old_best) => best_option = Some(old_best.max((likelihood, trial_bc))),
                    }
                    total_likelihood += likelihood;
                }
            }
            a[pos] = existing;
        }
        drop(a);

        let thresh = Of64::try_from(self.bc_confidence_threshold).ok()?;

        let expected_errors: f64 = qual.map_or(0.0, |q| q.iter().copied().map(probability).sum());

        if let Some((best_like, best_bc)) = best_option {
            if expected_errors < self.max_expected_barcode_errors
                && best_like / total_likelihood >= thresh
            {
                return Some((best_bc, 1)); // 1-Hamming distance
            }
        }
        None
    }
}

pub fn probability(qual: u8) -> f64 {
    //33 is the illumina qual offset
    let q = f64::from(qual);
    (10_f64).powf(-(q - 33.0) / 10.0)
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::BcSegSeq;
    use metric::{SimpleHistogram, TxHashSet};
    use proptest::proptest;

    fn posterior(
        wl: Whitelist,
        bc_counts: SimpleHistogram<BcSegSeq>,
        max_expected_barcode_errors: f64,
        bc_confidence_threshold: f64,
    ) -> BarcodeCorrector {
        BarcodeCorrector::new(
            wl,
            bc_counts,
            Posterior {
                max_expected_barcode_errors,
                bc_confidence_threshold,
            },
        )
    }

    #[test]
    pub fn test_barcode_correction() {
        let mut wl = TxHashSet::default();

        let b1 = BcSegSeq::from_bytes(b"AAAAA");
        let b2 = BcSegSeq::from_bytes(b"AAGAC");
        let b3 = BcSegSeq::from_bytes(b"ACGAA");
        let b4 = BcSegSeq::from_bytes(b"ACGTT");

        wl.insert(b1);
        wl.insert(b2);
        wl.insert(b3);
        wl.insert(b4);

        let mut bc_counts = SimpleHistogram::default();
        bc_counts.insert(b1, 100);
        bc_counts.insert(b2, 11);
        bc_counts.insert(b3, 2);

        let corrector = posterior(Whitelist::Plain(wl), bc_counts, 1.0, 0.95);

        // Easy
        let mut t1 = BarcodeSegment::with_sequence(b"AAAAA", BarcodeSegmentState::Invalid);
        // Low quality
        assert_eq!(
            corrector
                .correct_barcode_helper(t1, Some(BcSegQual::from_bytes(&[34, 34, 34, 66, 66]))),
            None
        );
        assert!(corrector
            .correct_barcode(&mut t1, Some(BcSegQual::from_bytes(&[34, 34, 34, 66, 66])),)
            .is_none());
        assert_eq!(
            t1,
            BarcodeSegment::with_sequence(b"AAAAA", BarcodeSegmentState::Invalid)
        );

        // Trivial correction
        let mut t2 = BarcodeSegment::with_sequence(b"AAAAT", BarcodeSegmentState::Invalid);
        assert_eq!(
            corrector
                .correct_barcode_helper(t2, Some(BcSegQual::from_bytes(&[66, 66, 66, 66, 40])))
                .unwrap()
                .0,
            BarcodeSegment::new(b1.into(), BarcodeSegmentState::ValidAfterCorrection)
        );
        assert!(corrector
            .correct_barcode(&mut t2, Some(BcSegQual::from_bytes(&[66, 66, 66, 66, 40])),)
            .is_some());
        assert_eq!(
            t2,
            BarcodeSegment::new(b1.into(), BarcodeSegmentState::ValidAfterCorrection)
        );

        // Pseudo-count kills you
        let t3 = BarcodeSegment::with_sequence(b"ACGAT", BarcodeSegmentState::Invalid);
        assert_eq!(
            corrector
                .correct_barcode_helper(t3, Some(BcSegQual::from_bytes(&[66, 66, 66, 66, 66]))),
            None
        );

        // Quality help you
        let t4 = BarcodeSegment::with_sequence(b"ACGAT", BarcodeSegmentState::Invalid);
        assert_eq!(
            corrector
                .correct_barcode_helper(t4, Some(BcSegQual::from_bytes(&[66, 66, 66, 66, 40])))
                .unwrap()
                .0,
            BarcodeSegment::new(b3.into(), BarcodeSegmentState::ValidAfterCorrection)
        );

        // Counts help you
        let t5 = BarcodeSegment::with_sequence(b"ACAAA", BarcodeSegmentState::Invalid);
        assert_eq!(
            corrector
                .correct_barcode_helper(t5, Some(BcSegQual::from_bytes(&[66, 66, 66, 66, 40])))
                .unwrap()
                .0,
            BarcodeSegment::new(b1.into(), BarcodeSegmentState::ValidAfterCorrection)
        );
    }

    #[test]
    pub fn test_barcode_correction_no_valid_counts() {
        let mut wl = TxHashSet::default();

        let b1 = BcSegSeq::from_bytes(b"AAAAA");
        let b2 = BcSegSeq::from_bytes(b"AAGAC");
        let b3 = BcSegSeq::from_bytes(b"ACGAA");
        let b4 = BcSegSeq::from_bytes(b"ACGTT");

        wl.insert(b1);
        wl.insert(b2);
        wl.insert(b3);
        wl.insert(b4);

        let bc_counts = SimpleHistogram::default();

        let val = posterior(Whitelist::Plain(wl), bc_counts, 1.0, 0.95);

        // Easy
        let t1 = BarcodeSegment::with_sequence(b"AAAAA", BarcodeSegmentState::Invalid);
        // Low quality
        assert_eq!(
            val.correct_barcode_helper(t1, Some(BcSegQual::from_bytes(&[34, 34, 34, 66, 66]))),
            None
        );

        // Trivial correction
        let t2 = BarcodeSegment::with_sequence(b"AAAAT", BarcodeSegmentState::Invalid);
        assert_eq!(
            val.correct_barcode_helper(t2, Some(BcSegQual::from_bytes(&[66, 66, 66, 66, 40])))
                .unwrap()
                .0,
            BarcodeSegment::with_sequence(b1.seq(), BarcodeSegmentState::ValidAfterCorrection)
        );
    }

    proptest! {
        #[test]
        fn prop_test_n_in_barcode(
            n_pos in (0..16usize)
        ) {
            let mut wl = TxHashSet::default();
            let bc = BcSegSeq::from_bytes(b"GCGATTGACCCAAAGG");
            wl.insert(bc);

            let corrector = posterior(Whitelist::Plain(wl), SimpleHistogram::default(), 1.0, 0.975);

            let mut bc_seq_with_n = bc.seq().to_vec();
            bc_seq_with_n[n_pos] = b'N';
            let mut qual = vec![53; bc_seq_with_n.len()];
            qual[n_pos] = 35;
            let bc_with_n = BarcodeSegment::with_sequence(&bc_seq_with_n, BarcodeSegmentState::Invalid);

            assert_eq!(
                corrector.correct_barcode_helper(bc_with_n, Some(BcSegQual::from_bytes(&qual))),
                Some((BarcodeSegment::with_sequence(
                    bc.seq(),
                    BarcodeSegmentState::ValidAfterCorrection
                ), 1))
            );

        }
    }
}

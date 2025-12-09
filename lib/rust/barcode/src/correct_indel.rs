#![deny(missing_docs)]

use crate::corrector::{BARCODE_CONFIDENCE_THRESHOLD, CorrectBarcode};
use crate::{
    BarcodeSegment, BarcodeSegmentContent, BarcodeSegmentState, BcSegQual, BcSegSeq, Whitelist,
};
use fastq_set::sseq::{HammingIterOpt, InsertionIterOpt};
use itertools::{Itertools, chain};
use metric::{Histogram, SimpleHistogram};

/// Correct one edit in a barcode.
/// Set its state to `InvalidAmbiguous` if the correction is ambiguous.
struct CorrectOneEdit;

impl CorrectBarcode for CorrectOneEdit {
    fn correct_barcode(
        &self,
        whitelist: &Whitelist,
        bc_counts: &SimpleHistogram<BcSegSeq>,
        observed_segment: BarcodeSegment,
        qual: Option<BcSegQual>,
    ) -> Option<(BarcodeSegment, u16)> {
        Self::correct_barcode_with_edit_type(whitelist, bc_counts, observed_segment, qual)
            .map(|(segment, _edit_type)| (segment, 1))
    }
}

/// The possible single edits.
enum EditType {
    Substitution,
    Insertion,
    Deletion,
}

impl CorrectOneEdit {
    /// Attempt to correct a barcode.
    /// Return both the corrected sequence and the type of edit.
    fn correct_barcode_with_edit_type(
        whitelist: &Whitelist,
        bc_counts: &SimpleHistogram<BcSegSeq>,
        observed_segment: BarcodeSegment,
        _qual: Option<BcSegQual>,
    ) -> Option<(BarcodeSegment, EditType)> {
        assert!(!observed_segment.is_valid());
        let seq = observed_segment.sequence();
        assert!(!whitelist.contains(seq));

        let shorter_seq = BcSegSeq::from_bytes(&seq.as_bytes()[..seq.len() - 1]);
        let longer_seq_a: BcSegSeq = seq.iter().chain(b"A").collect();
        let longer_seq_c: BcSegSeq = seq.iter().chain(b"C").collect();
        let longer_seq_g: BcSegSeq = seq.iter().chain(b"G").collect();
        let longer_seq_t: BcSegSeq = seq.iter().chain(b"T").collect();

        let candidates = chain!(
            seq.one_hamming_iter(HammingIterOpt::MutateNBase)
                .map(|seq| (seq, EditType::Substitution)),
            shorter_seq
                .one_insertion_iter(InsertionIterOpt::ExcludeNBase)
                .map(|seq| (seq, EditType::Insertion)),
            longer_seq_a
                .one_deletion_iter()
                .map(|seq| (seq, EditType::Deletion)),
            longer_seq_c
                .one_deletion_iter()
                .map(|seq| (seq, EditType::Deletion)),
            longer_seq_g
                .one_deletion_iter()
                .map(|seq| (seq, EditType::Deletion)),
            longer_seq_t
                .one_deletion_iter()
                .map(|seq| (seq, EditType::Deletion)),
        )
        // NOTE: since unique_by keeps the first element, this implies that we
        // will always prefer describing a 1-hamming-distance edit as a subsitution,
        // rather than (say) a deletion of the last base. CorrectSubNotIndel
        // relies on this behavior.
        .unique_by(|(seq, _edit_type)| *seq);

        let mut sum_count = 0;
        let Some((corrected_count, corrected_seq, edit_type)) = candidates
            .filter_map(|(seq, edit_type)| {
                whitelist
                    .exact_match_with_translation(&seq)
                    .map(|corrected_seq| (corrected_seq, edit_type))
            })
            .map(|(seq, edit_type)| {
                // Apply additive (Laplace) smoothing.
                let count = 1 + bc_counts.get(&seq);
                sum_count += count;
                (count, seq, edit_type)
            })
            .max_by_key(|(count, seq, _edit_type)| (*count, *seq))
        else {
            // No match
            return None;
        };

        let segment = if corrected_count as f64 / sum_count as f64 >= BARCODE_CONFIDENCE_THRESHOLD {
            // Return the corrected barcode.
            BarcodeSegment {
                content: BarcodeSegmentContent::Sequence(corrected_seq),
                state: BarcodeSegmentState::ValidAfterCorrection,
            }
        } else {
            // The correction is ambiguous.
            BarcodeSegment {
                content: observed_segment.content,
                state: BarcodeSegmentState::InvalidAmbiguous,
            }
        };
        Some((segment, edit_type))
    }
}

/// Correct one substitution in a barcode, and detect but do not correct an indel.
/// Set its sequence to the corrected sequence if the correction is a unique substitution.
/// Set its state to `InvalidAmbiguous` if the correction is ambiguous.
/// Set its state to `InvalidIndel` if the correction is an indel.
pub struct CorrectSubNotIndel;

impl CorrectBarcode for CorrectSubNotIndel {
    fn correct_barcode(
        &self,
        whitelist: &Whitelist,
        bc_counts: &SimpleHistogram<BcSegSeq>,
        observed_segment: BarcodeSegment,
        qual: Option<BcSegQual>,
    ) -> Option<(BarcodeSegment, u16)> {
        assert!(!observed_segment.is_valid());

        let Some((corrected_seq, edit_type)) = CorrectOneEdit::correct_barcode_with_edit_type(
            whitelist,
            bc_counts,
            observed_segment,
            qual,
        ) else {
            // No match
            return None;
        };

        if corrected_seq.state == BarcodeSegmentState::InvalidAmbiguous {
            // The correction is ambiguous.
            return Some((corrected_seq, 1));
        }
        assert_eq!(
            corrected_seq.state,
            BarcodeSegmentState::ValidAfterCorrection
        );

        Some((
            BarcodeSegment {
                content: corrected_seq.content,
                state: match edit_type {
                    EditType::Substitution => BarcodeSegmentState::ValidAfterCorrection,
                    EditType::Deletion | EditType::Insertion => BarcodeSegmentState::InvalidIndel,
                },
            },
            1,
        ))
    }
}

use crate::{VdjHierarchy, VdjReference};
use fxhash::{FxBuildHasher, FxHashMap};
use std::collections::HashMap;

const DEFAULT_KMER_LENGTH: usize = 20;

/// Control the strategy used to classify reads in `KmerClassify`
/// There are two strategies:
/// * `Lazy` (Default): In this mode, the `choice` corrsponding to the first kmer match in the read
///     is returned. Because not all kmers are necessarily checked, this works out to be the
///     fastest strategy, but with a higher false positive rate. This is recommended if you are
///     using the classifier to pick between IG/TCR. There is barely any 20-mer common between
///     IG and TCR receptor reference sequences in human/mouse, so the classificiation here
///     is pretty accurate, unless there is some weird chimera that is formed.
/// * `Strict`: In this mode, a read is classified as a `choice` only if at least
///     one kmer in the read maps to `choice` and no kmers in the read maps to other choices.
#[derive(Debug, Copy, Clone)]
pub enum KmerClassifyStrategy {
    Lazy,
    Strict,
}

/// Kmer based classification of a read sequence into one of the `VdjHierarchy` `choices()`
/// based on a `VdjReference`.
///
/// This is useful when, for example, you want to check whether a read maps to TCR or IG.
//
/// Internally a hashmap from kmers to the `VdjHierarchy` is stored. These will only contain
/// kmers uniquely associated with a `choice`. i.e Any kmer which maps to multiple `choices`
/// is not stored. Due to this reason, it is expected that you use a fairly large kmer length.
/// The default kmer length is 20, but this can be adjusted via one of the contructors. The
/// strategy used for classification may be controlled by selecting the appropriate
/// `KmerClassifyStrategy`.
#[derive(Debug)]
pub struct KmerClassify<'a, T> {
    kmer_len: usize,
    kmer_hash: FxHashMap<&'a [u8], T>,
    strategy: KmerClassifyStrategy,
}

impl<'a, T> KmerClassify<'a, T>
where
    T: VdjHierarchy,
{
    /// Create a new `KmerClassify` object with a kmer length of 20 and `Lazy` strategy.
    pub fn new(reference: &'a VdjReference) -> Self {
        KmerClassify::with_kmer_len_and_strategy(
            reference,
            DEFAULT_KMER_LENGTH,
            KmerClassifyStrategy::Lazy,
        )
    }

    /// Create a new `KmerClassify` object with the given kmer length and `Lazy` strategy.
    pub fn with_kmer_len(reference: &'a VdjReference, kmer_len: usize) -> Self {
        KmerClassify::with_kmer_len_and_strategy(reference, kmer_len, KmerClassifyStrategy::Lazy)
    }

    /// Create a new `KmerClassify` object with a kmer length of 20 and the given strategy.
    pub fn with_strategy(reference: &'a VdjReference, strategy: KmerClassifyStrategy) -> Self {
        KmerClassify::with_kmer_len_and_strategy(reference, DEFAULT_KMER_LENGTH, strategy)
    }

    /// Create a new `KmerClassify` object with the given kmer length and strategy.
    pub fn with_kmer_len_and_strategy(
        reference: &'a VdjReference,
        kmer_len: usize,
        strategy: KmerClassifyStrategy,
    ) -> Self {
        // Hashmap from kmer to Option<T>. The value will be none if more than one value of `T` maps to that kmer.
        let mut kmer_hash_build = FxHashMap::default();
        for ref_entry in reference {
            let seq = ref_entry.seq();
            let choice = Some(T::from_entry(ref_entry));
            for kmer in seq.windows(kmer_len) {
                kmer_hash_build
                    .entry(kmer)
                    .and_modify(|val| {
                        // in-place mutable access to an occupied entry
                        if *val != choice {
                            *val = None;
                        }
                    })
                    .or_insert(choice);
            }
        }
        let mut kmer_hash =
            HashMap::with_capacity_and_hasher(kmer_hash_build.len(), FxBuildHasher::default());
        for (k, v) in kmer_hash_build {
            if let Some(val) = v {
                kmer_hash.insert(k, val);
            }
        }
        KmerClassify {
            kmer_hash,
            kmer_len,
            strategy,
        }
    }

    /// Classify the read into one of `T::choices()`. The choice is made based on the
    /// strategy `KmerClassifyStrategy`
    ///
    /// # Input
    /// * `read`: An iterator over ascii encoded byte sequence
    ///
    /// # Example
    /// ```rust
    /// use vdj_reference::{VdjReference, VdjReceptor, KmerClassify};
    /// use bio::io::fasta;
    /// let tiny_fasta = b">459|TRAV10 ENST00000390432|TRAV10|5'UTR|TR|TRA|None|00\nAGTCAACTTCTGGGAGCAGATCTCTGCAGAATAAAA\n>314|IGLJ2 ENST00000390322|IGLJ2|J-REGION|IG|IGL|None|00\nTGTGGTATTCGGCGGAGGGACCAAGCTGACCGTCCTAG";
    /// let vdj_ref = VdjReference::from_fasta_reader(fasta::Reader::new(&tiny_fasta[..])).unwrap();
    /// let classifier = KmerClassify::<VdjReceptor>::new(&vdj_ref);
    /// assert_eq!(classifier.classify(b"GGGAGCAGATCTCTGCAGAATAAA"), Some(VdjReceptor::TR));
    /// assert_eq!(classifier.classify(b"GTATTCGGCGGAGGGACCAAGCTGACCGTCCT"), Some(VdjReceptor::IG));
    /// assert_eq!(classifier.classify(b"AGACTGTCCTGGCTCTTGACACATCCCCAG"), None);
    /// ```
    pub fn classify(&self, read: &[u8]) -> Option<T> {
        use KmerClassifyStrategy::{Lazy, Strict};
        let mut choice = None;
        for kmer in read.windows(self.kmer_len) {
            if let Some(&val) = self.kmer_hash.get(&kmer) {
                match self.strategy {
                    Lazy => return Some(val), // Return the first matching choice
                    Strict => {
                        if let Some(c) = choice {
                            // We have found a match before.
                            if c != val {
                                // Found two distinct choices for this read. Return None.
                                return None;
                            }
                        } else {
                            choice = Some(val);
                        }
                    }
                }
            }
        }
        choice
    }

    /// Classify the read into one of `T::choices()`, reverse complementing the
    /// read prior to kmer search. Computing the reverse complement triggers an allocation
    /// of a new vector with the same length as the read
    pub fn classify_rc(&self, read: &[u8]) -> Option<T> {
        let read_rc = byteseq::revcomp(read);
        self.classify(&read_rc)
    }

    /// Classify the read into one of `T::choices()`, optionally reverse complementing the
    /// read prior to kmer search. Computing the reverse complement triggers an allocation
    /// of a new vector with the same length as the read
    pub fn classify_maybe_rc(&self, read: &[u8], rc: bool) -> Option<T> {
        if rc {
            self.classify_rc(read)
        } else {
            self.classify(read)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::VdjReceptor;
    use bio::io::fasta;
    use indoc::indoc;

    fn tiny_ref() -> VdjReference {
        let tiny_fasta = indoc!(
            "
            >102|IGHV2-70D ENST00000390634|IGHV2-70D|5'UTR|IG|IGH|None|00
            ATCTCCACCAGCTCCACCCTCCCCTGGGTTCAAAAGACGAGGACAGGGCCTCGCTCAGTGAATCCTGCTCTCCACC
            >433|TRAJ42 ENST00000390495|TRAJ42|J-REGION|TR|TRA|None|00
            TGAATTATGGAGGAAGCCAAGGAAATCTCATCTTTGGAAAAGGCACTAAACTCTCTGTTAAACCAA
        "
        )
        .as_bytes();
        VdjReference::from_fasta_reader(fasta::Reader::new(tiny_fasta)).unwrap()
    }
    #[test]
    fn test_kmer_classifier() {
        let vdj_reference = tiny_ref();
        let classifier: KmerClassify<'_, VdjReceptor> = KmerClassify::new(&vdj_reference);
        assert_eq!(classifier.kmer_hash.len(), 57 + 47);
        assert_eq!(
            classifier.classify(b"ACCCTCCCCTGGGTTCAAAAGACGAGGACAGGGCCTCGCTCAGTGAATCCT"),
            Some(VdjReceptor::IG)
        );
        assert_eq!(
            classifier.classify(b"AGCCAAGGAAATCTCATCTTTGGAAAAGGCACTAAACTCTCTGT"),
            Some(VdjReceptor::TR)
        );
        assert_eq!(
            classifier.classify(b"TCATATTCCCCCCGACGACATAGGCCTTTGCCTCGTCACTCACCCCTGCC"),
            None
        );
    }

    #[test]
    fn test_kmer_classifier_lazy() {
        let vdj_reference = tiny_ref();
        let classifier: KmerClassify<'_, VdjReceptor> = KmerClassify::new(&vdj_reference);

        assert_eq!(
            //                    <------------IG-----------><-----random------->
            classifier.classify(b"ACCCTCCCCTGGGTTCAAAAGACGAGGTTCGGGCTGCTAAGATCTCA"),
            Some(VdjReceptor::IG)
        );

        assert_eq!(
            classifier
                //          <-----random-------><------------IG-----------><-----random------->
                .classify(b"TAGTCCGACGAAGCACTATGACCCTCCCCTGGGTTCAAAAGACGAGGTTCGGGCTGCTAAGATCTCA"),
            Some(VdjReceptor::IG)
        );

        assert_eq!(
            classifier
                //          <-----random-------><------------TR-----------><-----random------->
                .classify(b"TAGTCCGACGAAGCACTATGAGCCAAGGAAATCTCATCTTTGGAAAATTCGGGCTGCTAAGATCTCA"),
            Some(VdjReceptor::TR)
        );

        assert_eq!(
            //                    <------------TR-----------><------------IG----------->
            classifier.classify(b"AGCCAAGGAAATCTCATCTTTGGAAAAACCCTCCCCTGGGTTCAAAAGACGAGG"),
            Some(VdjReceptor::TR)
        ); // False positive
    }

    #[test]
    fn test_kmer_classifier_strict() {
        let vdj_reference = tiny_ref();
        let classifier: KmerClassify<'_, VdjReceptor> =
            KmerClassify::with_strategy(&vdj_reference, KmerClassifyStrategy::Strict);

        assert_eq!(
            //                    <------------IG-----------><-----random------->
            classifier.classify(b"ACCCTCCCCTGGGTTCAAAAGACGAGGTTCGGGCTGCTAAGATCTCA"),
            Some(VdjReceptor::IG)
        );

        assert_eq!(
            classifier
                //          <-----random-------><------------IG-----------><-----random------->
                .classify(b"TAGTCCGACGAAGCACTATGACCCTCCCCTGGGTTCAAAAGACGAGGTTCGGGCTGCTAAGATCTCA"),
            Some(VdjReceptor::IG)
        );

        assert_eq!(
            classifier
                //          <-----random-------><------------TR-----------><-----random------->
                .classify(b"TAGTCCGACGAAGCACTATGAGCCAAGGAAATCTCATCTTTGGAAAATTCGGGCTGCTAAGATCTCA"),
            Some(VdjReceptor::TR)
        );

        assert_eq!(
            //                    <------------TR-----------><------------IG----------->
            classifier.classify(b"AGCCAAGGAAATCTCATCTTTGGAAAAACCCTCCCCTGGGTTCAAAAGACGAGG"),
            None
        ); // Ambiguous
    }

    #[test]
    fn test_kmer_classifier_rc() {
        let vdj_reference = tiny_ref();
        let classifier: KmerClassify<'_, VdjReceptor> = KmerClassify::new(&vdj_reference);

        assert_eq!(
            //                    <------------IG-----------><-----random------->
            classifier.classify(b"ACCCTCCCCTGGGTTCAAAAGACGAGGTTCGGGCTGCTAAGATCTCA"),
            Some(VdjReceptor::IG)
        );

        assert_eq!(
            //                       <-----random-------><------------IG-rc-------->
            classifier.classify_rc(b"TGAGATCTTAGCAGCCCGAACCTCGTCTTTTGAACCCAGGGGAGGGT"),
            Some(VdjReceptor::IG)
        );

        assert_eq!(
            //                             <------------IG-----------><-----random------->
            classifier.classify_maybe_rc(b"ACCCTCCCCTGGGTTCAAAAGACGAGGTTCGGGCTGCTAAGATCTCA", false),
            Some(VdjReceptor::IG)
        );

        assert_eq!(
            //                             <-----random-------><------------IG-rc-------->
            classifier.classify_maybe_rc(b"TGAGATCTTAGCAGCCCGAACCTCGTCTTTTGAACCCAGGGGAGGGT", true),
            Some(VdjReceptor::IG)
        );
    }
}

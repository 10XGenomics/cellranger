// Copyright (c) 2018 10x Genomics, Inc. All rights reserved.
#![expect(missing_docs)]

//! Sized, stack-allocated container for a short quality string.

use crate::array::{ArrayContent, ByteArray};
use std::iter::Iterator;
use std::str;

#[derive(Clone, Copy, PartialOrd, Ord, PartialEq, Eq)]
pub struct SQualityContents;

impl ArrayContent for SQualityContents {
    /// Ensure that the input byte slice contains only valid quality characters, and panic
    /// otherwise.
    fn validate_bytes(seq: &[u8]) {
        for (i, &c) in seq.iter().enumerate() {
            let q = c as i16 - 33;
            assert!(
                (0..94).contains(&q),
                "Invalid quality value {q} ASCII character {c} at position {i}"
            );
        }
    }
    fn expected_contents() -> &'static str {
        "A valid read quality value"
    }
}

/// Fixed-sized container for a short quality string, with capacity determined by the type `N`.
/// Used as a convenient container for a barcode or UMI quality string.
/// An `SQualityGen` is guaranteed to contain only valid quality characters.
pub type SQualityGen<const N: usize> = ByteArray<SQualityContents, N>;

/// Fixed-sized container for a short quality string, up to 23bp in length.
/// Used as a convenient container for a barcode or UMI quality string.
/// An `SQuality` is guaranteed to contain only valid quality characters.
pub type SQuality = SQualityGen<23>;

#[cfg(test)]
mod squality_test {
    use super::*;
    use bincode;
    use proptest::{prop_assert_eq, proptest};

    const VALID_CHARS: &[u8; 42] = br##"!"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJ"##;

    #[test]
    fn test_sseq_valid_quality() {
        assert_eq!(SQuality::from_bytes(&VALID_CHARS[0..23]).len(), 23);
        assert_eq!(SQuality::from_bytes(&VALID_CHARS[23..42]).len(), 19);
        assert_eq!(
            SQuality::from_bytes(&VALID_CHARS[0..23]).to_string(),
            str::from_utf8(&VALID_CHARS[0..23]).unwrap()
        );
    }

    #[test]
    #[should_panic]
    fn test_sseq_invalid_quality_1() {
        let _ = SQuality::from_bytes(b"GHIJ ");
    }

    #[test]
    #[should_panic]
    fn test_sseq_invalid_quality_2() {
        let _ = SQuality::from_bytes(b" HIJK");
    }

    #[test]
    fn test_serde() {
        let sseqs = vec![
            SQuality::from_bytes(&VALID_CHARS[0..23]),
            SQuality::from_bytes(&VALID_CHARS[23..41]),
        ];

        let mut buf = Vec::new();
        bincode::serialize_into(&mut buf, &sseqs).unwrap();
        let roundtrip: Vec<SQuality> = bincode::deserialize_from(&buf[..]).unwrap();
        assert_eq!(sseqs, roundtrip);
    }

    proptest! {
        #[test]
        fn prop_test_serde_squality(
            ref seq in "[!FGHIJ]{0, 23}",
        ) {
            let target = SQuality::from_bytes(seq.as_bytes());
            let encoded: Vec<u8> = bincode::serialize(&target).unwrap();
            let decoded: SQuality = bincode::deserialize(&encoded[..]).unwrap();
            prop_assert_eq!(target, decoded);
        }
    }
}

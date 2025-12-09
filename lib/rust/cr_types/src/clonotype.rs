#![expect(missing_docs)]
use anyhow::{Context, Result};

const CLONOTYPE_PREFIX: &str = "clonotype";
const CONSENSUS_PREFIX: &str = "consensus";
const CONCAT_REF: &str = "concat_ref";

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct ClonotypeId<'a> {
    pub id: usize,
    pub sample_id: Option<&'a str>,
}

impl<'a> ClonotypeId<'a> {
    /// Parse the provided string as a ClonotypeId.
    pub fn parse(s: &'a str) -> Result<Self> {
        Ok(
            if let Some((sample_id, clonotype_num)) = s.rsplit_once('_') {
                ClonotypeId {
                    id: parse_clonotype_number(clonotype_num)?,
                    sample_id: Some(sample_id),
                }
            } else {
                ClonotypeId {
                    id: parse_clonotype_number(s)?,
                    sample_id: None,
                }
            },
        )
    }
}

impl std::fmt::Display for ClonotypeId<'_> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if let Some(sample_id) = self.sample_id {
            write!(f, "{sample_id}_")?;
        }
        write!(f, "{CLONOTYPE_PREFIX}{}", self.id)
    }
}

/// Parse the number from a string that looks like clonotype123.
fn parse_clonotype_number(s: &str) -> Result<usize> {
    s[CLONOTYPE_PREFIX.len()..]
        .parse()
        .with_context(|| format!("invalid clonotype number: \"{s}\""))
}

impl ClonotypeId<'_> {
    pub fn consensus_name(&self, consensus_index: usize) -> String {
        format!("{self}_{CONSENSUS_PREFIX}_{consensus_index}")
    }
    pub fn concat_ref_name(&self, concat_ref_index: usize) -> String {
        format!("{self}_{CONCAT_REF}_{concat_ref_index}")
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_clonotype_id() {
        let clonotype_id = ClonotypeId {
            id: 1,
            sample_id: None,
        };
        assert_eq!(clonotype_id.to_string(), "clonotype1");
        assert_eq!(
            ClonotypeId::parse("clonotype1").unwrap(),
            ClonotypeId {
                id: 1,
                sample_id: None
            }
        );
        assert_eq!(clonotype_id.consensus_name(2), "clonotype1_consensus_2");
        assert_eq!(clonotype_id.concat_ref_name(3), "clonotype1_concat_ref_3");
    }

    #[test]
    fn test_clonotype_id_with_sample() {
        let clonotype_id = ClonotypeId {
            id: 1,
            sample_id: Some("test_sample"),
        };
        assert_eq!(clonotype_id.to_string(), "test_sample_clonotype1");
        assert_eq!(
            ClonotypeId::parse("test_sample_clonotype1").unwrap(),
            ClonotypeId {
                id: 1,
                sample_id: Some("test_sample")
            }
        );
        assert_eq!(
            clonotype_id.consensus_name(2),
            "test_sample_clonotype1_consensus_2"
        );
        assert_eq!(
            clonotype_id.concat_ref_name(3),
            "test_sample_clonotype1_concat_ref_3"
        );
    }

    proptest::proptest! {
        #[test]
        fn test_clonotype_id_roundtrip(id in 0..1_000_000usize) {
            let clonotype_id = ClonotypeId { id, sample_id: None };
            assert_eq!(ClonotypeId::parse(&clonotype_id.to_string()).unwrap(), clonotype_id);
        }

        #[test]
        fn test_clonotype_id_roundtrip_with_sample(id in 0..1_000_000usize) {
            let clonotype_id = ClonotypeId { id, sample_id: Some("test_sample") };
            assert_eq!(ClonotypeId::parse(&clonotype_id.to_string()).unwrap(), clonotype_id);
        }
    }
}

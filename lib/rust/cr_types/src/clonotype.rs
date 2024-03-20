pub const CLONOTYPE_PREFIX: &str = "clonotype";
pub const SAMPLE_PREFIX: &str = "sample";
pub const CONSENSUS_PREFIX: &str = "consensus";
pub const CONCAT_REF: &str = "concat_ref";

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub struct ClonotypeId {
    pub id: usize,
    pub sample_number: Option<usize>,
}

impl std::fmt::Display for ClonotypeId {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if let Some(sample_number) = self.sample_number {
            write!(f, "{SAMPLE_PREFIX}{sample_number}_")?;
        }
        write!(f, "{CLONOTYPE_PREFIX}{}", self.id)
    }
}

impl std::str::FromStr for ClonotypeId {
    type Err = anyhow::Error;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let (left, right) = if s.starts_with(CLONOTYPE_PREFIX) {
            (None, s)
        } else if s.starts_with(SAMPLE_PREFIX) {
            let (left, right) = s.split_once('_').unwrap();
            (Some(left), right)
        } else {
            anyhow::bail!("Invalid clonotype id: {}", s);
        };
        Ok(ClonotypeId {
            id: right[CLONOTYPE_PREFIX.len()..].parse()?,
            sample_number: left.map(|s| s[SAMPLE_PREFIX.len()..].parse()).transpose()?,
        })
    }
}

impl ClonotypeId {
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
            sample_number: None,
        };
        assert_eq!(clonotype_id.to_string(), "clonotype1");
        assert_eq!(
            "clonotype1".parse::<ClonotypeId>().unwrap(),
            ClonotypeId {
                id: 1,
                sample_number: None
            }
        );
        assert_eq!(clonotype_id.consensus_name(2), "clonotype1_consensus_2");
        assert_eq!(clonotype_id.concat_ref_name(3), "clonotype1_concat_ref_3");
    }

    #[test]
    fn test_clonotype_id_with_sample() {
        let clonotype_id = ClonotypeId {
            id: 1,
            sample_number: Some(3),
        };
        assert_eq!(clonotype_id.to_string(), "sample3_clonotype1");
        assert_eq!(
            "sample2_clonotype1".parse::<ClonotypeId>().unwrap(),
            ClonotypeId {
                id: 1,
                sample_number: Some(2)
            }
        );
        assert_eq!(
            clonotype_id.consensus_name(2),
            "sample3_clonotype1_consensus_2"
        );
        assert_eq!(
            clonotype_id.concat_ref_name(3),
            "sample3_clonotype1_concat_ref_3"
        );
    }

    proptest::proptest! {
        #[test]
        fn test_clonotype_id_roundtrip(id in 0..1_000_000usize) {
            let clonotype_id = ClonotypeId { id, sample_number: None };
            assert_eq!(clonotype_id.to_string().parse::<ClonotypeId>().unwrap(), clonotype_id);
        }

        #[test]
        fn test_clonotype_id_roundtrip_with_sample(id in 0..1_000_000usize, sample_number in 0..1_000_000usize) {
            let clonotype_id = ClonotypeId { id, sample_number: Some(sample_number) };
            assert_eq!(clonotype_id.to_string().parse::<ClonotypeId>().unwrap(), clonotype_id);
        }
    }
}

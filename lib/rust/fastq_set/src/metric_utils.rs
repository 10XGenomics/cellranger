#![expect(missing_docs)]
use bio::pattern_matching;

pub const ILLUMINA_QUAL_OFFSET: u8 = 33;

pub type Pattern = pattern_matching::bndm::BNDM;
pub struct PatternCheck {
    pattern: Pattern,
}

impl PatternCheck {
    pub fn new(pattern_seq: &[u8]) -> Self {
        PatternCheck {
            pattern: Pattern::new(pattern_seq),
        }
    }
    pub fn exists(&self, read: &[u8]) -> bool {
        self.pattern.find_all(read).next().is_some()
    }
}

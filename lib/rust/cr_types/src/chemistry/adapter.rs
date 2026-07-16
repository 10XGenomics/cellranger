//! Helpers for producting adapter trimmers.
use crate::LibraryType;
use crate::chemistry::{ChemistryDef, ChemistryName};
use crate::rna_read::AdapterTrimmingParameters;
use crate::types::AlignerParam;
use anyhow::{Result, bail};
use fastq_set::long_read_adapter_trimmer::{
    AdapterLayout, DEFAULT_MAX_R2_LENGTH, ReadSplittingStrategy, TrimmerConfig,
};

/// R2 length cap applied when aligning with STAR.
const STAR_MAX_R2_LENGTH: usize = 90;

/// Return a trimmer for the provided chemistry, library type, layout and params.
/// Return None if we can't uniquely determine one.
/// Return Err if we can't get necessary parameters from the chemistry.
pub fn trimmer_for_chemistry(
    chem: &ChemistryDef,
    library_type: LibraryType,
    layout: AdapterLayout,
    params: &AdapterTrimmingParameters,
) -> Result<Option<TrimmerConfig>> {
    if library_type.is_gex() && (chem.name.is_sc_3p_v3() || chem.name.is_sc_3p_v4()) {
        let Some(bc_umi_r1_offset) = chem.bc_umi_offset_into_r1() else {
            bail!("bc + umi not defined on R1 for {}", chem.name);
        };
        assert!(
            layout == layouts::FIVE_PRIME_GEX,
            "empirically determined adapter layout {layout:?} is incompatible with detected 3' GEX chemistry",
        );
        Ok(Some(TrimmerConfig {
            layout: layouts::THREE_PRIME_GEX,
            splitting_strategy: ReadSplittingStrategy::FixedOffset {
                offset: bc_umi_r1_offset,
            },
            allow_truncated_constructs: library_type.is_gex(),
            is_directional: params.is_directional,
            max_r2_length: max_r2_length_for_aligner(params.aligner),
        }))
    } else if chem.name.is_sc_5p_v3() {
        let Some(bc_umi_r1_offset) = chem.bc_umi_offset_into_r1() else {
            bail!("bc + umi not defined on R1 for {}", chem.name);
        };
        let splitting_strategy = if library_type.is_gex() {
            ReadSplittingStrategy::FixedOffset {
                offset: bc_umi_r1_offset,
            }
        } else {
            ReadSplittingStrategy::AdapterLocated {}
        };
        Ok(Some(TrimmerConfig {
            layout,
            splitting_strategy,
            allow_truncated_constructs: library_type.is_gex(),
            is_directional: params.is_directional,
            max_r2_length: max_r2_length_for_aligner(params.aligner),
        }))
    } else if chem.name.is_rtl() == Some(true) {
        Ok(Some(TrimmerConfig {
            layout,
            splitting_strategy: ReadSplittingStrategy::AdapterLocated {},
            allow_truncated_constructs: false,
            is_directional: params.is_directional,
            max_r2_length: max_r2_length_for_aligner(params.aligner),
        }))
    } else {
        Ok(None)
    }
}

/// Return a trimmer for 5' libraries.
pub fn trimmer_5p(
    params: &AdapterTrimmingParameters,
    library_type: LibraryType,
) -> Option<TrimmerConfig> {
    let bc_umi_r1_offset = ChemistryName::FivePrimeR2V3
        .def()
        .bc_umi_offset_into_r1()
        .expect("bc + umi offset should be defined for 5' v3 chemistry");
    let (layout, allow_truncated_constructs, splitting_strategy) = match library_type {
        LibraryType::Gex => (
            layouts::FIVE_PRIME_GEX,
            true,
            ReadSplittingStrategy::FixedOffset {
                offset: bc_umi_r1_offset,
            },
        ),
        LibraryType::Antibody => (
            layouts::FIVE_PRIME_ANTIBODY,
            false,
            ReadSplittingStrategy::AdapterLocated {},
        ),
        LibraryType::Crispr => (
            layouts::FIVE_PRIME_CRISPR,
            false,
            ReadSplittingStrategy::AdapterLocated {},
        ),
        _ => {
            return None;
        }
    };
    Some(TrimmerConfig {
        layout,
        splitting_strategy,
        allow_truncated_constructs,
        is_directional: params.is_directional,
        max_r2_length: max_r2_length_for_aligner(params.aligner),
    })
}

/// Return compatible trimmers for flex libraries.
pub fn trimmer_flex(
    params: &AdapterTrimmingParameters,
    library_type: LibraryType,
) -> Vec<TrimmerConfig> {
    let layouts: &[AdapterLayout] = match library_type {
        LibraryType::Gex => &[layouts::FLEX_GEX],
        LibraryType::Antibody => &[layouts::FLEX_ANTIBODY],
        LibraryType::Crispr => &[layouts::FLEX_CRISPR, layouts::FLEX_HYBRID_CRISPR],
        _ => return vec![],
    };
    let max_r2_length = max_r2_length_for_aligner(params.aligner);
    layouts
        .iter()
        .map(|layout| TrimmerConfig {
            layout: layout.clone(),
            splitting_strategy: ReadSplittingStrategy::AdapterLocated {},
            allow_truncated_constructs: false,
            is_directional: params.is_directional,
            max_r2_length,
        })
        .collect()
}

/// Returns the maximum R2 length to use for the given aligner.
fn max_r2_length_for_aligner(aligner: AlignerParam) -> usize {
    match aligner {
        AlignerParam::Star => STAR_MAX_R2_LENGTH,
        _ => DEFAULT_MAX_R2_LENGTH,
    }
}

/// Statically-defined adapter layouts.
pub mod layouts {

    use AdapterKind::{
        NexteraR2Bp12Rc, NonPolydtRc, R1Bp12, RtPrimerRc, SmallRnaR2Rc, Spacer, TruncatedR2Rc,
        X12Capture,
    };
    use fastq_set::adapters::AdapterKind;
    use fastq_set::long_read_adapter_trimmer::AdapterLayout;

    /// 3' gene expression libraries.
    pub const THREE_PRIME_GEX: AdapterLayout = AdapterLayout {
        terminal_5p: R1Bp12,
        probe_prefix: None,
        target_prefix: None,
        terminal_3p: RtPrimerRc,
    };

    /// 5' gene expression libraries.
    pub const FIVE_PRIME_GEX: AdapterLayout = AdapterLayout {
        terminal_5p: R1Bp12,
        probe_prefix: None,
        target_prefix: None,
        terminal_3p: NonPolydtRc,
    };

    /// All 5' antibody libraries.
    pub const FIVE_PRIME_ANTIBODY: AdapterLayout = AdapterLayout {
        terminal_5p: R1Bp12,
        probe_prefix: None,
        target_prefix: Some(Spacer),
        terminal_3p: NexteraR2Bp12Rc,
    };

    /// All 5' CRISPR libraries.
    pub const FIVE_PRIME_CRISPR: AdapterLayout = AdapterLayout {
        terminal_5p: R1Bp12,
        probe_prefix: None,
        target_prefix: Some(Spacer),
        terminal_3p: TruncatedR2Rc,
    };

    /// All flex gene expression libraries.
    pub const FLEX_GEX: AdapterLayout = AdapterLayout {
        terminal_5p: R1Bp12,
        probe_prefix: Some(X12Capture),
        target_prefix: Some(Spacer),
        terminal_3p: SmallRnaR2Rc,
    };

    /// All flex antibody libraries.
    pub const FLEX_ANTIBODY: AdapterLayout = AdapterLayout {
        terminal_5p: R1Bp12,
        probe_prefix: Some(X12Capture),
        target_prefix: Some(Spacer),
        terminal_3p: NexteraR2Bp12Rc,
    };

    /// Standard flex CRISPR libraries.
    pub const FLEX_CRISPR: AdapterLayout = AdapterLayout {
        terminal_5p: R1Bp12,
        probe_prefix: Some(X12Capture),
        target_prefix: Some(Spacer),
        terminal_3p: TruncatedR2Rc,
    };

    /// Hybrid flex CRISPR libraries.
    pub const FLEX_HYBRID_CRISPR: AdapterLayout = AdapterLayout {
        terminal_5p: R1Bp12,
        probe_prefix: None,
        target_prefix: Some(X12Capture),
        terminal_3p: TruncatedR2Rc,
    };

    #[cfg(test)]
    mod test {
        use super::*;

        #[test]
        fn test_short_3p_search() {
            assert!(!THREE_PRIME_GEX.terminal_3p_short_search());
            assert!(!FIVE_PRIME_GEX.terminal_3p_short_search());
            assert!(!FIVE_PRIME_ANTIBODY.terminal_3p_short_search());
            assert!(FIVE_PRIME_CRISPR.terminal_3p_short_search());
            assert!(!FLEX_GEX.terminal_3p_short_search());
            assert!(!FLEX_ANTIBODY.terminal_3p_short_search());
            assert!(FLEX_CRISPR.terminal_3p_short_search());
            assert!(FLEX_HYBRID_CRISPR.terminal_3p_short_search());
        }
    }
}

use anyhow::Result;
use martian::prelude::*;
use martian_derive::{make_mro, MartianStruct};
use serde::{Deserialize, Serialize};
use vdj_asm_utils::filter_log::FilterSwitch;

#[derive(Debug, Clone, Serialize, Deserialize, MartianStruct)]
pub struct StageInputs {
    pub disable_count: Option<bool>,
    pub is_antibody_only: Option<bool>,
    pub is_non_targeted_gex: Option<bool>,
    pub multiplet_filter: Option<bool>,
    pub shared_contig_filter: Option<bool>,
    pub umi_baseline_filter: Option<bool>,
}

#[derive(Debug, Clone, Serialize, Deserialize, MartianStruct)]
pub struct StageOutputs {
    pub filter_switch: FilterSwitch,
}

pub struct MakeFilterSwitch;

#[make_mro]
impl MartianMain for MakeFilterSwitch {
    type StageInputs = StageInputs;
    type StageOutputs = StageOutputs;
    fn main(&self, args: Self::StageInputs, _rover: MartianRover) -> Result<Self::StageOutputs> {
        let StageInputs {
            disable_count,
            is_antibody_only,
            is_non_targeted_gex,
            multiplet_filter,
            shared_contig_filter,
            umi_baseline_filter,
        } = args;

        let mut filter_switch = match (disable_count, is_antibody_only, is_non_targeted_gex) {
            (Some(false), Some(false), Some(true)) => FilterSwitch {
                asm_shared_contig: false,
                enclone_shared_contig: false,
                enclone_umi: false,
                enclone_multiplet: true,
            },
            _ => FilterSwitch::default(),
        };

        if let Some(filter) = multiplet_filter {
            filter_switch.enclone_multiplet = filter;
        }

        if let Some(filter) = shared_contig_filter {
            filter_switch.asm_shared_contig = filter;
            filter_switch.enclone_shared_contig = filter;
        }

        if let Some(filter) = umi_baseline_filter {
            filter_switch.enclone_umi = filter;
        }

        Ok(StageOutputs { filter_switch })
    }
}

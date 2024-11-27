//! Martian stage COPY_CHEMISTRY_SPEC
//! This is a shim stage to adapt older pipelines that only accept a single
//! chemistry spec from the user into the spec-per-library-type required by
//! DETECT_CHEMISTRY.
//!
//! This stage can also handle expanding a ChemistrySet into the individual
//! chemistries required for each library.

use anyhow::Result;
use cr_types::chemistry::{AutoOrRefinedChemistry, ChemistrySpecs};
use cr_types::sample_def::SampleDef;
use martian::prelude::*;
use martian_derive::{make_mro, MartianStruct};
use multi::config::ChemistryParam;
use serde::{Deserialize, Serialize};

#[derive(Clone, Deserialize, MartianStruct)]
pub struct CopyChemistrySpecStageInputs {
    pub sample_defs: Vec<SampleDef>,
    pub chemistry_spec: ChemistryParam,
}

#[derive(Clone, Serialize, Deserialize, MartianStruct)]
pub struct CopyChemistrySpecStageOutputs {
    pub chemistry_specs: ChemistrySpecs,
}

pub struct CopyChemistrySpec;

#[make_mro(volatile = strict)]
impl MartianMain for CopyChemistrySpec {
    type StageInputs = CopyChemistrySpecStageInputs;
    type StageOutputs = CopyChemistrySpecStageOutputs;

    fn main(&self, args: Self::StageInputs, _rover: MartianRover) -> Result<Self::StageOutputs> {
        Ok(Self::StageOutputs {
            chemistry_specs: args
                .sample_defs
                .iter()
                .map(|sample| {
                    let library_type = sample.library_type.unwrap_or_default();
                    let chem = match args.chemistry_spec {
                        ChemistryParam::AutoOrRefined(chem) => chem,
                        ChemistryParam::Set(set) => AutoOrRefinedChemistry::Refined(
                            set.chemistry_for_library_type(library_type)?,
                        ),
                    };
                    Ok((library_type, chem))
                })
                .collect::<Result<_>>()?,
        })
    }
}

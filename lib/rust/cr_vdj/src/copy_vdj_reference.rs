//! CopyVdjReference stage code
//! ```text
//! vdj_GRCh38_alts_ensembl-4.0.0/
//! ├── fasta
//! │   └── regions.fa
//! └── reference.json
//! ```

use anyhow::{Context, Result};
use martian::prelude::*;
use martian_derive::{make_mro, martian_filetype, MartianStruct};
use martian_filetypes::json_file::JsonFile;
use serde::{Deserialize, Serialize};
use std::path::PathBuf;

martian_filetype! {FaFile, "fa"}

#[derive(Debug, Clone, Serialize, Deserialize, MartianStruct)]
struct VdjRefFastaFolder {
    regions: FaFile,
}

#[derive(Debug, Clone, Serialize, Deserialize, MartianStruct)]
struct VdjRefFolder {
    fasta: VdjRefFastaFolder,
    reference: JsonFile<()>,
}

#[derive(Debug, Clone, Serialize, Deserialize, MartianStruct)]
pub struct CopyVdjReferenceStageInputs {
    vdj_reference_path: Option<PathBuf>,
}

#[derive(Debug, Clone, Serialize, Deserialize, MartianStruct)]
pub struct CopyVdjReferenceStageOutputs {
    vdj_reference: Option<VdjRefFolder>,
}

pub struct CopyVdjReference;

#[make_mro]
impl MartianMain for CopyVdjReference {
    type StageInputs = CopyVdjReferenceStageInputs;
    type StageOutputs = CopyVdjReferenceStageOutputs;
    fn main(&self, args: Self::StageInputs, rover: MartianRover) -> Result<Self::StageOutputs> {
        let vdj_reference = match args.vdj_reference_path {
            // NOTE: We don't need to worry about creating the folder structure. Martian will take
            // care of that when populating the outs folder by virtue of the struct nesting.
            Some(path) => Some(VdjRefFolder {
                fasta: {
                    VdjRefFastaFolder {
                        regions: {
                            let src = FaFile::new(path.join("fasta"), "regions");
                            let dest: FaFile = rover.make_path("regions");
                            std::fs::copy(&src, &dest).with_context(|| {
                                format!(
                                    "Error: unable to copy {} to {}",
                                    &src.as_ref().display(),
                                    &dest.as_ref().display()
                                )
                            })?;
                            dest
                        },
                    }
                },
                reference: {
                    let src: JsonFile<()> = JsonFile::new(&path, "reference");
                    let dest: JsonFile<()> = rover.make_path("reference");
                    std::fs::copy(&src, &dest).with_context(|| {
                        format!(
                            "Error: unable to copy {} to {}",
                            &src.as_ref().display(),
                            &dest.as_ref().display()
                        )
                    })?;
                    dest
                },
            }),
            None => None,
        };
        Ok(CopyVdjReferenceStageOutputs { vdj_reference })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use cr_lib::testing::tools::ensure_no_diff;
    use std::path::Path;

    #[test]
    fn test_no_vdj_reference() {
        assert!(CopyVdjReference
            .test_run_tmpdir(CopyVdjReferenceStageInputs {
                vdj_reference_path: None,
            })
            .unwrap()
            .vdj_reference
            .is_none());
    }

    #[test]
    fn test_with_vdj_reference() {
        let dir = tempfile::tempdir().unwrap();
        let vdj_ref = CopyVdjReference
            .test_run(
                &dir,
                CopyVdjReferenceStageInputs {
                    vdj_reference_path: Some(
                        "../dui_tests/test_resources/reference/vdj_GRCh38_alts_ensembl-4.0.0/"
                            .into(),
                    ),
                },
            )
            .unwrap()
            .vdj_reference
            .unwrap();
        println!("{vdj_ref:#?}");
        assert!(vdj_ref.reference.as_ref().exists());
        assert!(vdj_ref.fasta.regions.as_ref().exists());
        ensure_no_diff(
            &vdj_ref.reference,
            Path::new(
                "../dui_tests/test_resources/reference/vdj_GRCh38_alts_ensembl-4.0.0/reference.json",
            ),
        );
        ensure_no_diff(
            &vdj_ref.fasta.regions,
            Path::new(
                "../dui_tests/test_resources/reference/vdj_GRCh38_alts_ensembl-4.0.0/fasta/regions.fa",
            ),
        );
    }
}

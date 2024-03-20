//! WriteWsJson stage code

use crate::parse_aggr_csv::VdjAggrCsvLibrary;
use crate::setup_vdj_aggr::EncloneProtoMetaFormat;
use crate::websummary::annotation_card::VdjAggrAnnotationTable;
use crate::websummary::cdr3_table::{VdjAggrSharedCdr3, VdjAggrSharedCdr3Row};
use crate::websummary::cells_card::{VdjAggrCellsRow, VdjAggrCellsTable};
use crate::websummary::clonotype_hist::ClonotypeHist;
use crate::websummary::clonotype_table::{VdjAggrClonotypeRow, VdjAggrClonotypeTable};
use crate::websummary::hero_metrics::VdjAggrHeroMetrics;
use crate::websummary::{VdjAggrPipelineInfo, VdjAggrWsContent, VdjAggrWsSummaryTab};
use crate::write_contig_proto::ProtoFile;
use anyhow::Result;
use cr_types::clonotype::ClonotypeId;
use cr_websummary::{Percent, WsSample};
use enclone_proto::proto_io::ClonotypeIter;
use itertools::Itertools;
use martian::prelude::*;
use martian_derive::{make_mro, MartianStruct};
use martian_filetypes::json_file::JsonFile;
use martian_filetypes::tabular_file::CsvFile;
use martian_filetypes::{FileTypeRead, FileTypeWrite};
use metric::SimpleHistogram;
use serde::{Deserialize, Serialize};
use std::cmp::Reverse;
use std::collections::{BTreeMap, BinaryHeap, HashSet};
use std::convert::Into;
use std::path::PathBuf;
use vdj_reference::{VdjChain, VdjReceptor, VdjReferenceInfo};

const TOP_N: usize = 10;

// TODO: Sync with clonotyp_assigner::write_clonotype_outs
#[derive(Clone, Debug, Serialize, Deserialize, Default)]
struct ClonotypesCsvRow {
    clonotype_id: String,
    frequency: usize,
    proportion: f64,
    cdr3s_aa: String,
}

#[derive(Debug, Clone, Deserialize, MartianStruct)]
pub struct WriteWsJsonStageInputs {
    vdj_reference_path: PathBuf,
    libraries: Vec<VdjAggrCsvLibrary>,
    enclone_output: ProtoFile,
    // TODO: Read enclone_gem_well_meta from enclone_output
    enclone_gem_well_meta: EncloneProtoMetaFormat,
    sample_id: String,
    sample_desc: String,
    clonotypes_csv: CsvFile<ClonotypesCsvRow>,
    receptor: VdjReceptor,
}

#[derive(Debug, Clone, Serialize, Deserialize, MartianStruct)]
pub struct WriteWsJsonStageOutputs {
    web_summary_content: JsonFile<VdjAggrWsContent>,
    per_origin_hist: JsonFile<Vec<(String, SimpleHistogram<usize>)>>,
}

struct CellBarcodeSummary {
    gem_well: u32,
    clonotype_id: usize,
    productive_chains: HashSet<VdjChain>,
}

impl CellBarcodeSummary {
    fn is_paired(&self) -> bool {
        self.is_ighk_paired()
            || self.is_ighl_paired()
            || self.is_trab_paired()
            || self.is_trgd_paired()
    }
    fn is_ighk_paired(&self) -> bool {
        self.productive_chains.contains(&VdjChain::IGH)
            && self.productive_chains.contains(&VdjChain::IGK)
    }
    fn is_ighl_paired(&self) -> bool {
        self.productive_chains.contains(&VdjChain::IGH)
            && self.productive_chains.contains(&VdjChain::IGL)
    }
    fn is_trab_paired(&self) -> bool {
        self.productive_chains.contains(&VdjChain::TRA)
            && self.productive_chains.contains(&VdjChain::TRB)
    }
    fn is_trgd_paired(&self) -> bool {
        self.productive_chains.contains(&VdjChain::TRG)
            && self.productive_chains.contains(&VdjChain::TRD)
    }
}

fn top_shared_cdr3(
    per_cdr3_gw_hist: &BTreeMap<String, SimpleHistogram<u32>>,
    f: impl Fn(u32) -> String,
) -> Vec<VdjAggrSharedCdr3Row> {
    let mut rows = BinaryHeap::with_capacity(TOP_N + 1);
    for (cdr3_aa, hist) in per_cdr3_gw_hist {
        let mut per_cat_hist = SimpleHistogram::default();
        for (gw, n) in hist.distribution() {
            per_cat_hist.observe_by(&f(*gw), n.count());
        }
        rows.push(Reverse(VdjAggrSharedCdr3Row {
            num_categories: per_cat_hist.distribution().len(),
            num_cells: per_cat_hist.raw_counts().map(|c| *c as usize).sum(),
            cdr3_aa: cdr3_aa.to_string(),
        }));
        while rows.len() > TOP_N {
            rows.pop();
        }
    }
    rows.into_sorted_vec().into_iter().map(|r| r.0).collect()
}

// This is our stage struct
pub struct WriteWsJson;

#[make_mro(stage_name = WRITE_WEB_SUMMARY_JSON)]
impl MartianMain for WriteWsJson {
    type StageInputs = WriteWsJsonStageInputs;
    type StageOutputs = WriteWsJsonStageOutputs;
    fn main(&self, args: Self::StageInputs, rover: MartianRover) -> Result<Self::StageOutputs> {
        let mut cbc_summaries = Vec::new();
        let mut per_cdr3_gw_hist = BTreeMap::new(); // For each cdr3, the histogram of number of cells with that cdr3 per gem well
        for (i, cl) in ClonotypeIter::from_file(args.enclone_output)?.enumerate() {
            let chain_types: Vec<VdjChain> = cl
                .chains
                .into_iter()
                .map(|ch| ch.chain_type.parse().unwrap())
                .collect();
            for ex_cl in cl.exact_clonotypes {
                let ex_cl_chain_types: HashSet<_> = ex_cl
                    .chains
                    .iter()
                    .map(|c| chain_types[c.index as usize])
                    .collect();
                let cdr3_aas: Vec<_> = ex_cl
                    .chains
                    .iter()
                    .map(|c| {
                        format!(
                            "{}:{}",
                            chain_types[c.index as usize],
                            c.chain.cdr3_aa_string()
                        )
                    })
                    .collect();
                for bc in ex_cl.cell_barcodes {
                    let gem_well = bc.split('-').last().unwrap().parse()?;
                    cbc_summaries.push(CellBarcodeSummary {
                        gem_well,
                        clonotype_id: i + 1,
                        productive_chains: ex_cl_chain_types.clone(),
                    });
                    for cdr3_aa in &cdr3_aas {
                        per_cdr3_gw_hist
                            .entry(cdr3_aa.clone())
                            .or_insert_with(SimpleHistogram::default)
                            .observe(&gem_well);
                    }
                }
            }
        }

        let num_cells = cbc_summaries.len();
        let num_clonotypes = cbc_summaries.last().unwrap().clonotype_id;
        let num_paired = cbc_summaries.iter().filter(|s| s.is_paired()).count();
        let hero_metrics = VdjAggrHeroMetrics::new(num_cells, num_clonotypes, num_paired);

        let clonotypes = args.clonotypes_csv.read()?;
        let cl_table_rows: Vec<_> = clonotypes
            .into_iter()
            .take(TOP_N)
            .map(|r| VdjAggrClonotypeRow {
                id: r.clonotype_id.parse::<ClonotypeId>().unwrap().id,
                cdr3_aas: r.cdr3s_aa.split(';').map(Into::into).collect(),
                num_cells: r.frequency,
                proportion: Percent::Float(r.proportion),
            })
            .collect();

        let vdj_annotation = VdjAggrAnnotationTable {
            paired_cells_frac: (num_paired as f64) / (num_cells as f64),
            ighk_paired_cells_frac: match args.receptor {
                VdjReceptor::IG => Some(
                    (cbc_summaries.iter().filter(|s| s.is_ighk_paired()).count() as f64)
                        / (num_cells as f64),
                ),
                VdjReceptor::TR => None,
                VdjReceptor::TRGD => None,
            },
            ighl_paired_cells_frac: match args.receptor {
                VdjReceptor::IG => Some(
                    (cbc_summaries.iter().filter(|s| s.is_ighl_paired()).count() as f64)
                        / (num_cells as f64),
                ),
                VdjReceptor::TR => None,
                VdjReceptor::TRGD => None,
            },
            trab_paired_cells_frac: match args.receptor {
                VdjReceptor::IG => None,
                VdjReceptor::TR => Some(
                    (cbc_summaries.iter().filter(|s| s.is_trab_paired()).count() as f64)
                        / (num_cells as f64),
                ),
                VdjReceptor::TRGD => None,
            },
            trgd_paired_cells_frac: match args.receptor {
                VdjReceptor::IG => None,
                VdjReceptor::TR => None,
                VdjReceptor::TRGD => Some(
                    (cbc_summaries.iter().filter(|s| s.is_trgd_paired()).count() as f64)
                        / (num_cells as f64),
                ),
            },
        };

        // Compute the cells table
        // We need clonotype histogram per (donor, origin)
        let gw_info = args.enclone_gem_well_meta.read()?.per_gem_well_info;
        let mut per_origin_hist = BTreeMap::new();
        for summary in cbc_summaries {
            let info = &gw_info[&summary.gem_well];
            per_origin_hist
                .entry((
                    info.library_id.clone(),
                    info.donor.clone(),
                    info.origin.clone(),
                ))
                .or_insert_with(SimpleHistogram::default)
                .observe_owned(summary.clonotype_id);
        }

        let per_origin_hist_file: JsonFile<_> = rover.make_path("per_origin_hist");
        per_origin_hist_file.write(
            &per_origin_hist
                .iter()
                .map(|((_, _, o), h)| (o.clone(), h.clone()))
                .collect(),
        )?;

        let top_cl_ids: Vec<_> = cl_table_rows.iter().map(|c| c.id).collect();
        let proportions = per_origin_hist
            .iter()
            .map(|((_library_id, _donor, origin), hist)| {
                (
                    origin.clone(),
                    top_cl_ids
                        .iter()
                        .map(|id| (hist.get(id) as f64) / (num_cells as f64))
                        .collect(),
                )
            })
            .collect();
        let vdj_clonotype_hist = ClonotypeHist {
            ids: top_cl_ids,
            proportions,
        };

        let vdj_aggr_cells = VdjAggrCellsTable(
            per_origin_hist
                .into_iter()
                .map(|((library_id, donor, origin), hist)| VdjAggrCellsRow {
                    sample_id: library_id,
                    donor,
                    origin,
                    cells: hist.raw_counts().map(|c| *c as usize).sum(),
                    clonotypes: hist.distribution().len(),
                    diversity: hist.effective_diversity(),
                })
                .collect(),
        );

        let vdj_shared_cdr3 = VdjAggrSharedCdr3 {
            by_library: top_shared_cdr3(&per_cdr3_gw_hist, |gw| gw_info[&gw].library_id.clone()),
            by_donor: top_shared_cdr3(&per_cdr3_gw_hist, |gw| gw_info[&gw].donor.clone()),
            by_origin: top_shared_cdr3(&per_cdr3_gw_hist, |gw| gw_info[&gw].origin.clone()),
        };

        let ws_content = VdjAggrWsContent {
            sample: WsSample::aggr(args.sample_id.clone(), args.sample_desc.clone()),
            summary_tab: VdjAggrWsSummaryTab {
                hero_metrics,
                pipeline_info_table: VdjAggrPipelineInfo {
                    run_id: args.sample_id,
                    run_desc: args.sample_desc,
                    vdj_ref: {
                        let ref_info =
                            VdjReferenceInfo::from_reference_folder(&args.vdj_reference_path)?;
                        match ref_info.version {
                            Some(v) => format!("{}-{v}", ref_info.genomes),
                            None => ref_info.genomes,
                        }
                    },
                    pipeline_version: rover.pipelines_version(),
                    num_samples: args.libraries.len(),
                    num_donors: args.libraries.iter().map(|l| &l.donor).unique().count(),
                    num_origins: args.libraries.iter().map(|l| &l.origin).unique().count(),
                },
                vdj_clonotype: VdjAggrClonotypeTable(cl_table_rows),
                vdj_annotation,
                vdj_aggr_cells,
                vdj_shared_cdr3,
                vdj_clonotype_hist,
            },
        };

        let web_summary_content: JsonFile<_> = rover.make_path("ws_content");
        serde_json::to_writer_pretty(std::io::stdout().lock(), &ws_content)?;
        web_summary_content.write(&ws_content)?;

        Ok(WriteWsJsonStageOutputs {
            web_summary_content,
            per_origin_hist: per_origin_hist_file,
        })
    }
}

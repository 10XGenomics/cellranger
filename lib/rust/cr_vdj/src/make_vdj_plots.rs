#![expect(missing_docs)]
use crate::clonotype_hist::{ClonotypeHist, ProportionsType};
use crate::clonotype_table::{VdjAggrClonotypeRow, VdjAggrClonotypeTable};
use anyhow::Result;
use clonotype_assigner::write_clonotype_outs::ClonotypesCsvRow;
use cr_types::MetricsFile;
use cr_types::clonotype::ClonotypeId;
use cr_websummary::multi::websummary_vdj::ClonotypeInfo;
use cr_websummary::{CardWithTable, ChartWithHelp, GenericTable, Percent};
use martian::prelude::*;
use martian_derive::{MartianStruct, make_mro};
use martian_filetypes::json_file::JsonFile;
use martian_filetypes::tabular_file::CsvFile;
use martian_filetypes::{FileTypeRead, FileTypeWrite};
use metric::JsonReporter;
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Serialize, Deserialize, MartianStruct)]
pub struct MakeVdjPlotsStageInputs {
    per_sample: bool,
    sample_id: String,
    sample_desc: String,
    receptor: String,
    n50_n50_rpu: i64,
    // CMO or Hashtag multiplexing do not have total_read_pairs as a sample-level metric
    total_read_pairs: Option<i64>,
    clonotypes_csv: Option<CsvFile<ClonotypesCsvRow>>,
}
#[derive(Debug, Clone, Serialize, Deserialize, MartianStruct)]
pub struct MakeVdjPlotsStageOutputs {
    sample_info_summary: MetricsFile,
    clonotype_info_json: Option<JsonFile<ClonotypeInfo>>,
}
pub struct MakeVdjPlots;

fn make_clonotype_hist_plot(clonotypes: CsvFile<ClonotypesCsvRow>) -> ChartWithHelp {
    let (ids, proportions): (Vec<usize>, Vec<f64>) = clonotypes
        .read()
        .unwrap()
        .iter()
        .take(10)
        .map(|row| {
            (
                ClonotypeId::parse(&row.clonotype_id).unwrap().id,
                row.proportion,
            )
        })
        .unzip();
    ChartWithHelp::from(ClonotypeHist {
        ids,
        proportions: ProportionsType::SingleSample(proportions),
    })
}

fn make_clonotype_table(clonotypes: CsvFile<ClonotypesCsvRow>) -> GenericTable {
    let table_rows = clonotypes
        .read()
        .unwrap()
        .iter()
        .take(10)
        .map(|r| VdjAggrClonotypeRow {
            id: ClonotypeId::parse(&r.clonotype_id).unwrap().id,
            cdr3_aas: r.cdr3s_aa.split(';').map(Into::into).collect(),
            num_cells: r.frequency,
            proportion: Percent::Float(r.proportion),
        })
        .collect();
    let table: CardWithTable = VdjAggrClonotypeTable(table_rows).into();
    table.table
}

#[make_mro(mem_gb = 1)]
impl MartianMain for MakeVdjPlots {
    type StageInputs = MakeVdjPlotsStageInputs;
    type StageOutputs = MakeVdjPlotsStageOutputs;

    fn main(&self, args: Self::StageInputs, rover: MartianRover) -> Result<Self::StageOutputs> {
        let sample_info_summary: MetricsFile = rover.make_path("sample_info_summary");
        let mut sample_metrics = JsonReporter::default();
        sample_metrics.insert("sample_id".to_string(), args.sample_id);
        sample_metrics.insert("sample_desc".to_string(), args.sample_desc);
        sample_metrics.insert("chain_type".to_string(), args.receptor);
        // n50_n50_rpu is a library level metric and missing from the per sample
        if args.per_sample {
            sample_metrics.insert("n50_n50_rpu".to_string(), args.n50_n50_rpu);
            sample_metrics.insert("VDJ_total_read_pairs".to_string(), args.total_read_pairs);
        };
        sample_info_summary.write(&sample_metrics)?;
        let clonotype_info_json = if let Some(clonotype_csv) = args.clonotypes_csv {
            let clonotype_plot = make_clonotype_hist_plot(clonotype_csv.clone());

            let clonotype_info = ClonotypeInfo {
                table: make_clonotype_table(clonotype_csv),
                plot: clonotype_plot.plot,
                help: cr_websummary::TitleWithHelp {
                    help: clonotype_plot.help.help,
                    title: clonotype_plot.help.title,
                },
            };
            let json_path: JsonFile<ClonotypeInfo> = rover.make_path("clonotype_info");
            json_path.write(&clonotype_info)?;
            Some(json_path)
        } else {
            None
        };

        Ok(MakeVdjPlotsStageOutputs {
            sample_info_summary,
            clonotype_info_json,
        })
    }
}

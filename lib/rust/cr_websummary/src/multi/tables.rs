use crate::alert::Alert;
use crate::{
    card_with_table, CardWithTable, CountAndPercent, FloatAsInt, GenericTable, MakePretty, Percent,
    PercentF1, TableRow, TitleWithHelp,
};
use serde::{Deserialize, Serialize};
use websummary_derive::make_tables;

make_tables!("src/multi/tables.toml");

// TODO: Sync this with VDJ aggr table
card_with_table!(
    table_name: VdjClonotypeFrequencyTable,
    row_name: VdjClonotypeFrequencyRow,
    title: "Top 10 Clonotype CDR3 Sequences",
    help_text: "We list CDR3 sequences of the 10 most frequent clonotypes",
    columns: {
        id: usize : "Clonotype ID",
        cdr3_aas : Vec<String> : "CDR3s",
        num_cells: usize : "Frequency",
        proportion: Percent : "Proportion"
    }
);

impl Alert for VdjClonotypeFrequencyTable {}

// Logging options.

#[derive(Default)]
pub struct LogOpts {
    pub print_umi: bool,       // print umis
    pub show_supp: bool,       // show read support for each graph edge
    pub print_seq_edges: bool, // print sequence for each graph edge
    pub nreject: bool,         // don't show rejected contigs
    pub ngood: bool,           // don't show good contigs
    pub umi_seq: bool,         // print reads for each umi
    pub survive: bool,         // log details about surviving umis
    pub vis: bool,             // show visual align of V+J ref to contig
    pub json: bool,            // print json annotations for each contig
    pub strong_edges: bool,    // print edges in strong paths
    pub nucounts: bool,        // don't print ucounts
    pub trace_seq: String,     // trace the given DNA sequence through the graph
    pub trace_umis: bool,      // trace all UMIs through the graph
    pub print_strong: bool,    // print strong paths
    pub paths: bool,           // print paths for contigs
}

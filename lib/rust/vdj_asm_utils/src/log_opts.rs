use perf_stats::{elapsed, mem_usage_gb, peak_mem_usage_gb};
use std::io::Write;
use std::time::Instant;

// Logging options.

pub struct LogOpts {
    pub print_seq: bool,       // print dna sequence of contig
    pub print_qual: bool,      // print qual scores of contig
    pub print_umi: bool,       // print umis
    pub show_supp: bool,       // show read support for each graph edge
    pub print_seq_edges: bool, // print sequence for each graph edge
    pub nreject: bool,         // don't show rejected contigs
    pub ngood: bool,           // don't show good contigs
    pub ngraph: bool,          // don't show graph
    pub umi_seq: bool,         // print reads for each umi
    pub survive: bool,         // log details about surviving umis
    pub vis: bool,             // show visual align of V+J ref to contig
    pub clock: bool,           // print some timing info
    pub mem: bool,             // print some mem usage info
    pub npipeline: bool,       // don't show pipeline logging
    pub json: bool,            // print json annotations for each contig
    pub strong_edges: bool,    // print edges in strong paths
    pub keep_all: bool,        // keep all logs even if graph empty
    pub nucounts: bool,        // don't print ucounts
    pub trace_seq: String,     // trace the given DNA sequence through the graph
    pub trace_umis: bool,      // trace all UMIs through the graph
    pub print_strong: bool,    // print strong paths
    pub paths: bool,           // print paths for contigs
}

impl LogOpts {
    pub fn new() -> LogOpts {
        LogOpts {
            print_seq: false,
            print_qual: false,
            print_umi: false,
            show_supp: false,
            print_seq_edges: false,
            nreject: false,
            ngood: false,
            ngraph: true,
            umi_seq: false,
            survive: false,
            vis: false,
            clock: false,
            mem: false,
            npipeline: false,
            json: false,
            strong_edges: false,
            keep_all: false,
            nucounts: false,
            trace_seq: String::new(),
            trace_umis: false,
            print_strong: false,
            paths: false,
        }
    }

    pub fn report_perf_stats(&self, log: &mut Vec<u8>, t: &Instant, msg: &str) {
        if self.clock || self.mem {
            fwriteln!(
                log,
                "used {:.1} seconds {}, \
                 mem = {:.2} GB, peak = {:.2} GB",
                elapsed(t),
                msg,
                mem_usage_gb(),
                peak_mem_usage_gb()
            );
        }
    }

    pub fn report_perf_stats_now(&self, t: &Instant, msg: &str) {
        if self.clock {
            println!(
                "used {:.1} seconds {}, mem = {:.2} GB, peak = {:.2} GB",
                elapsed(t),
                msg,
                mem_usage_gb(),
                peak_mem_usage_gb()
            );
        }
    }
}

impl Default for LogOpts {
    fn default() -> Self {
        Self::new()
    }
}

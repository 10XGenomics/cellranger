use clap::builder::NonEmptyStringValueParser;
use clap::Parser;

#[derive(Parser, Debug, Clone)]
pub struct MrpArgs {
    /// Job manager to use. Valid options: local (default), sge,
    /// lsf, slurm or path to a .template file. Search for help on "Cluster Mode"
    /// at support.10xgenomics.com for more details on configuring
    /// the pipeline to use a compute cluster.
    #[clap(
        long,
        default_value = "local",
        value_name = "MODE",
        value_parser = NonEmptyStringValueParser::new(),
    )]
    jobmode: String,

    /// Set max cores the pipeline may request at one time. Only
    /// applies to local jobs.
    #[clap(long, value_name = "NUM")]
    localcores: Option<usize>,

    /// Set max GB the pipeline may request at one time. Only
    /// applies to local jobs.
    #[clap(long, value_name = "NUM")]
    localmem: Option<usize>,

    /// Set max virtual address space in GB for the pipeline.
    /// Only applies to local jobs.
    #[clap(long, value_name = "NUM")]
    localvmem: Option<usize>,

    /// Reserve enough threads for each job to ensure enough
    /// memory will be available, assuming each core on your
    /// cluster has at least this much memory available. Only
    /// applies to cluster jobmodes.
    #[clap(long, value_name = "NUM")]
    mempercore: Option<usize>,

    /// Set max jobs submitted to cluster at one time. Only
    /// applies to cluster jobmodes.
    #[clap(long, value_name = "NUM")]
    maxjobs: Option<usize>,

    /// Set delay between submitting jobs to cluster, in ms.
    /// Only applies to cluster jobmodes.
    #[clap(long, value_name = "NUM")]
    jobinterval: Option<usize>,

    /// The path to a JSON file that specifies stage-level
    /// overrides for cores and memory. Finer-grained
    /// than --localcores, --mempercore and --localmem.
    /// Consult https://support.10xgenomics.com/ for an example
    /// override file.
    #[clap(long, value_name = "PATH")]
    overrides: Option<String>,

    /// Enables Volatile Data Removal. Valid options:
    /// rolling (default), post, strict, or disable
    #[clap(
        long,
        hide = true,
        value_name = "MODE",
        value_parser = ["rolling", "post", "disable", "strict"],
    )]
    vdrmode: Option<String>,

    /// Output the results to this directory.
    #[clap(long, value_name = "PATH")]
    pub output_dir: Option<String>,

    /// Serve web UI at http://localhost:PORT
    #[clap(long, value_name = "PORT")]
    uiport: Option<usize>,

    /// Do not serve the web UI.
    #[clap(long)]
    disable_ui: bool,

    /// Keep web UI running after pipestance completes or fails.
    #[clap(long)]
    noexit: bool,

    /// Skip preflight checks.
    #[clap(long)]
    pub nopreflight: bool,
}

impl MrpArgs {
    /// Convert this struct into a vector of command line arguments.
    pub(crate) fn get_args(&self) -> Vec<String> {
        [
            Some(format!("--jobmode={}", self.jobmode)),
            self.jobinterval.map(|x| format!("--jobinterval={x}")),
            self.localcores.map(|x| format!("--localcores={x}")),
            self.localmem.map(|x| format!("--localmem={x}")),
            self.localvmem.map(|x| format!("--localvmem={x}")),
            self.maxjobs.map(|x| format!("--maxjobs={x}")),
            self.mempercore.map(|x| format!("--mempercore={x}")),
            self.overrides.as_ref().map(|x| format!("--overrides={x}")),
            self.uiport.map(|x| format!("--uiport={x}")),
            self.vdrmode.as_ref().map(|x| format!("--vdrmode={x}")),
            self.output_dir.as_ref().map(|x| format!("--psdir={x}")),
            self.disable_ui.then_some("--disable-ui".to_string()),
            self.noexit.then_some("--noexit".to_string()),
            self.nopreflight.then_some("--nopreflight".to_string()),
        ]
        .into_iter()
        .flatten()
        .collect()
    }
}

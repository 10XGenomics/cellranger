#![deny(missing_docs)]

struct PlannerOutput {
    /// An MRO source string
    mro: String,

    /// The version of the requested pipeline.
    version: MartianVersion,

    /// The name of the top-level call in the MRO, e.g. PHASER_SVCALLER_CS.
    call: String,

    /// Specifies the environment modifications required to run mrp.
    /// Note that this should not include MROFLAGS.  If it does, they
    /// will be ignored.
    env: EnvironmentSpec,

    /// Specifies the storage high-water-mark estimate for this pipeline,
    /// in bytes.
    sizealloc: i64,

    /// Specifies additional flags to pass to mrp.  Note that any of these may
    /// be overridden or ignored by the invocation requester.
    mroflags: Vec<String>,

    /// Specifies the location for a default mrp executable, if available, as
    /// well as its version.  This may be ignored by the invocation requester.
    /// If left unset, the latest version available to Enterprise will be used.
    mrp: MrpLocation,

    /// Specifies the list of pipelines available from this planner.
    pipeline: Vec<String>,
}

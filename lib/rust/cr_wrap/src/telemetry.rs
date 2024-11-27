use crate::env::PkgEnv;
use crate::IntoExitCode;
use std::ffi::OsString;
use std::os::unix::process::ExitStatusExt;
use std::path::PathBuf;
use std::process::{Child, Command, ExitStatus};
use std::sync::mpsc::sync_channel;
use std::thread;
use std::time::{Duration, Instant};

/// Path to the collector binary, relative to the tarball root.
pub const TELEMETRY_BIN_PATH: &str = "lib/bin/telemetry";

/// Command wrappers should implement this to determine if we should collect
/// telemetry for a particular executing command.
pub trait CollectTelemetry {
    /// Return true if we should collect telemetry for this command.
    fn should_collect_telemetry(&self) -> bool;
}

/// How long the telemetry collection should wait for all child processes to complete.
const JOIN_TIMEOUT: Duration = Duration::from_secs(5);

pub struct TelemetryCollector {
    /// The path to the telemetry collector binary.
    collector_binary: PathBuf,
    /// The CLI subcommand we're executing.
    subcommand: String,
    /// The CLI args, in a canonical --foo=bar form.
    cli_args: Vec<OsString>,
    /// True if the collector should run.
    should_run: bool,
    /// Handles to any collector child processes we've spawned, and the time
    /// at which we spawned them.
    procs: Vec<(Child, Instant)>,
}

impl TelemetryCollector {
    /// Initialize a telemetry collector for a particular command.
    pub fn new_for_command<T: CollectTelemetry>(
        cmd: &T,
        subcommand: String,
        cli_args: Vec<OsString>,
        pkg_env: &PkgEnv,
    ) -> Self {
        Self {
            collector_binary: pkg_env.subcmd_path(TELEMETRY_BIN_PATH),
            cli_args,
            subcommand,
            should_run: cmd.should_collect_telemetry(),
            procs: Default::default(),
        }
    }

    /// Start telemetry collection if configured to do so.
    pub fn collect(
        &mut self,
        mrp_args: Option<&str>,
        mrp_exit_status: Option<&ExitStatus>,
        wall_time: Option<Duration>,
    ) {
        if !self.should_run {
            return;
        }
        let mut command = Command::new(&self.collector_binary);
        command.arg("collect");
        if let Some(mrp_args) = mrp_args {
            command.arg("-mrp_flags").arg(mrp_args);
        }
        if let Some(mrp_exit_status) = mrp_exit_status {
            command
                .arg("-mrp_exit_code")
                .arg(mrp_exit_status.into_u8().to_string());
        }
        if let Some(wall_time) = wall_time {
            command
                .arg("-walltime")
                .arg(wall_time.as_secs_f64().to_string());
        }
        // Pass the canonicalized CLI args into the collector.
        command.arg("--").arg(&self.subcommand).args(&self.cli_args);

        match command.spawn() {
            Ok(child) => {
                self.procs.push((child, Instant::now()));
            }
            Err(err) => eprintln!("Failed to start telemetry collector: {err}"),
        }
    }

    /// Await completion of all telemetry collection processes, or time out.
    ///
    /// The timeout will be `JOIN_TIMEOUT` minus the time elapsed since we last
    /// spawned a telemetry collector process.
    ///
    /// This method is called when this type is dropped, so if the location of
    /// the wait isn't critical, this method does not need to be manually called.
    pub fn join(&mut self) {
        let Some((_, last_spawned_at)) = self.procs.iter().max_by_key(|(_, spawned_at)| spawned_at)
        else {
            return;
        };
        let timeout = JOIN_TIMEOUT.saturating_sub(last_spawned_at.elapsed());
        // Spawn a thread to wait on the child procs, then wait on that process
        // for at most the configured timeout.
        let children = std::mem::take(&mut self.procs);
        let (send, recv) = sync_channel::<()>(1);
        thread::spawn(move || {
            for (mut child, _) in children {
                let Ok(exit_status) = child.wait() else {
                    continue;
                };
                if let Some(sig) = exit_status.signal() {
                    eprintln!("telemetry collector terminated by signal: {sig}");
                }
            }
            let _ = send.send(());
        });
        let _ = recv.recv_timeout(timeout);
    }
}

impl Drop for TelemetryCollector {
    fn drop(&mut self) {
        self.join();
    }
}

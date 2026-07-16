#![deny(missing_docs)]
use anyhow::Result;
use regex::Regex;
use std::borrow::Cow;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

// Regexes to check for known incompatible OS versions.  If any of these match, the pipeline will fail with an error.
// The error can be suppressed for by setting TENX_IGNORE_DEPRECATED_OS=1 in your environment.
const REGEX_RELEASE_INCOMPATIBLE: &[(&str, &str)] = &[
    // Red Hat or derivatives 5, 6 or 7
    (r" [567]\.", "/etc/redhat-release"),
    (r" [567]\.", "/etc/rocks-release"),
    (r" [567]\.", "/etc/os-release"),
    (r" [567]\.", "/etc/system-release"),
    // Ubuntu 16 or earlier
    (r" 1[0-6]\.", "/etc/lsb-release"),
    (r" [1-9]\.", "/etc/lsb-release"),
    // Suse 10 or 11
    (r"SUSE.* 1[01]\b", "/etc/SuSE-release"),
    // Debian 7 or earlier
    (r#"PRETTY_NAME="[dD]ebian.*\b[1-7]\."#, "/etc/os-release"),
];

// Warn for OSes that are currently compatible but will be incompatible in the future, so users have time to upgrade before it becomes a hard error.
const REGEX_RELEASE_DEPRECATED: &[(&str, &str)] = &[
    // Red Hat or derivatives 5, 6 or 7
    (r" [567]\.", "/etc/redhat-release"),
    (r" [567]\.", "/etc/rocks-release"),
    (r" [567]\.", "/etc/os-release"),
    (r" [567]\.", "/etc/system-release"),
    // Ubuntu 18 or earlier
    (r" 1[0-8]\.", "/etc/lsb-release"),
    (r" [1-9]\.", "/etc/lsb-release"),
    // Suse 10 or 11
    (r"SUSE.* 1[01]\b", "/etc/SuSE-release"),
    // Debian 7 or earlier
    (r#"PRETTY_NAME="[dD]ebian.*\b[1-7]\."#, "/etc/os-release"),
];

fn search_file_for_re(path: &Path, regex: &str) -> Option<String> {
    if !path.exists() {
        return None;
    }

    let re = Regex::new(regex).unwrap();
    if let Ok(file) = File::open(path) {
        let mut reader = BufReader::new(file);
        let mut line = String::new();
        while reader.read_line(&mut line).unwrap_or(0) != 0 {
            if re.is_match(&line) {
                return Some(line);
            }
        }
    }

    None
}

enum PlatformResult {
    Ok,
    Warn(Cow<'static, str>),
    Err(Cow<'static, str>),
}

impl PlatformResult {
    fn into_result(self, mut handle_warning: impl FnMut(&str)) -> Result<()> {
        match self {
            PlatformResult::Ok => Ok(()),
            PlatformResult::Warn(message) => {
                handle_warning(&message);
                Ok(())
            }
            PlatformResult::Err(msg) => Err(anyhow::anyhow!(msg)),
        }
    }
}

fn check_known_os() -> PlatformResult {
    // check for a known lsb-release file, if we find a matching entry, error due to incompatibility.
    for (regex, release) in REGEX_RELEASE_INCOMPATIBLE {
        if let Some(line) = search_file_for_re(Path::new(release), regex) {
            return PlatformResult::Err(Cow::from(format!(
                "This operating system version is unsupported:
    {}
For more information, see support.10xgenomics.com/os-support.

To continue running this version for now, set TENX_IGNORE_DEPRECATED_OS=1
in your environment.",
                line.trim(),
            )));
        }
    }
    // check for a known lsb-release file, if we find a matching entry, warn about future incompatibility
    for (regex, release) in REGEX_RELEASE_DEPRECATED {
        if let Some(line) = search_file_for_re(Path::new(release), regex) {
            return PlatformResult::Warn(Cow::from(format!(
                "WARNING: The operating system version is deprecated and will no longer be supported in a future release:
    {}
Future releases of this pipeline will require a more current system.
For more information, see support.10xgenomics.com/os-support.",
                line.trim(),
            )));
        }
    }
    PlatformResult::Ok
}

fn check_kernel_version() -> PlatformResult {
    let uname = rustix::system::uname();
    let release = uname.release().to_str().unwrap();
    let re = Regex::new(r"\D*(\d+)\.(\d+)").unwrap();
    if let Some(caps) = re.captures(release) {
        let version: Result<Vec<i64>, _> = (1..=2).map(|i| caps[i].parse()).collect();
        if let Ok(&[major, minor, ..]) = version.as_deref() {
            if (major, minor) < (3, 10) {
                return PlatformResult::Err(Cow::from(format!(
                    "The kernel used by this operating system version is unsupported: {release}
This release requires kernel version 3.10.0 or higher. Future releases of this
pipeline will require kernel version 4.15.0 or higher.
For more information, see support.10xgenomics.com/os-support.

To continue running this version for now, set TENX_IGNORE_DEPRECATED_OS=1
in your environment."
                )));
            } else if (major, minor) < (4, 15) {
                return PlatformResult::Warn(Cow::from(format!(
                    "WARNING: The kernel used by this operating system version will no longer
be supported in a future release: {release}
Future releases of this pipeline will require kernel version 4.15.0 or
higher.  For more information, see support.10xgenomics.com/os-support."
                )));
            }
        }
    }
    PlatformResult::Ok
}

/// Check whether the GNU libc version is sufficient on Linux and do nothing on macOS.
#[cfg(not(target_os = "linux"))]
fn check_libc_version() -> PlatformResult {
    PlatformResult::Ok
}
#[cfg(target_os = "linux")]
fn check_libc_version() -> PlatformResult {
    let libc_version = unsafe { std::ffi::CStr::from_ptr(libc::gnu_get_libc_version()) }
        .to_str()
        .unwrap();
    let re = Regex::new(r"\D*(\d+)\.(\d+)").unwrap();
    if let Some(caps) = re.captures(libc_version) {
        let version: Result<Vec<i64>, _> = (1..=2).map(|i| caps[i].parse()).collect();
        if let Ok(&[major, minor, ..]) = version.as_deref() {
            if (major, minor) < (2, 17) {
                return PlatformResult::Err(Cow::from(format!(
                    "The glibc version of this operating system version is unsupported: {libc_version}
This release requires glibc version 2.17 or higher. Future releases of this
pipeline will require glibc version 2.28 or higher.
For more information, see support.10xgenomics.com/os-support.

To continue running this version for now, set TENX_IGNORE_DEPRECATED_OS=1
in your environment."
                )));
            } else if (major, minor) < (2, 28) {
                return PlatformResult::Warn(Cow::from(format!(
                    "WARNING: The glibc version of this operating system version will no longer
be supported in a future release: {libc_version}
Future releases of this pipeline will require glibc version 2.28 or
higher.  For more information, see support.10xgenomics.com/os-support."
                )));
            }
        }
    }
    PlatformResult::Ok
}

fn check_cpu_features() -> PlatformResult {
    #[cfg(target_arch = "x86_64")]
    if !is_x86_feature_detected!("avx") {
        return PlatformResult::Err(Cow::from(
            "This CPU does not support AVX, which is required. \
             Set TENX_IGNORE_DEPRECATED_OS=1 in your environment to suppress this error. \
             For more information, see \
             https://www.10xgenomics.com/support/software/cell-ranger/downloads/cr-system-requirements",
        ));
    }
    #[cfg(target_arch = "x86_64")]
    if !is_x86_feature_detected!("avx2") {
        return PlatformResult::Warn(Cow::from(
            "This CPU does not support AVX2, which will be required in the future. \
             For more information, see \
             https://www.10xgenomics.com/support/software/cell-ranger/downloads/cr-system-requirements",
        ));
    }
    PlatformResult::Ok
}

/// Checks for deprecated OS, kernel, libc and CPU features
pub fn oscheck(mut handle_warning: impl FnMut(&str)) -> Result<()> {
    check_known_os().into_result(&mut handle_warning)?;
    check_kernel_version().into_result(&mut handle_warning)?;
    check_libc_version().into_result(&mut handle_warning)?;
    check_cpu_features().into_result(&mut handle_warning)?;
    Ok(())
}

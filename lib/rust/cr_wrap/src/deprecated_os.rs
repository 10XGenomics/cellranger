use anyhow::Result;
use libc::{c_char, c_int};
use regex::Regex;
use std::borrow::Cow;
use std::ffi::CStr;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

const REGEX_RELEASE: &[(&str, &str)] = &[
    (r" [56]\.", "/etc/redhat-release"),
    (r" [56]\.", "/etc/rocks-release"),
    (r" [56]\.", "/etc/os-release"),
    (r" [56]\.", "/etc/system-release"),
    // Ubuntu 13 or earlier
    (r" 1[0-3]\.", "/etc/lsb-release"),
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
    fn into_result(self) -> Result<()> {
        match self {
            PlatformResult::Ok => Ok(()),
            PlatformResult::Warn(msg) => {
                eprintln!("{msg}");
                Ok(())
            }
            PlatformResult::Err(msg) => Err(anyhow::anyhow!(msg)),
        }
    }
}

fn check_known_os() -> PlatformResult {
    // check for a known lsb-release file, if we find a matching entry, we're good
    for (regex, release) in REGEX_RELEASE {
        if let Some(line) = search_file_for_re(Path::new(release), regex) {
            return PlatformResult::Err(Cow::from(format!(
                "This operating system version is unsupported or will soon be unsupported:
    {}
Future releases of this pipeline will require a more current system.
For more information, see support.10xgenomics.com/os-support.

To continue running this version for now, set TENX_IGNORE_DEPRECATED_OS=1
in your environment.",
                line.trim(),
            )));
        }
    }
    PlatformResult::Ok
}

fn to_str(s: &[c_char]) -> &str {
    unsafe {
        let bytes = CStr::from_ptr(s.as_ptr()).to_bytes();
        std::str::from_utf8_unchecked(bytes)
    }
}

fn check_kernel_version() -> PlatformResult {
    let mut utsname: libc::utsname = unsafe { std::mem::zeroed() };
    if unsafe { libc::uname(&mut utsname as *mut libc::utsname) } != 0 {
        return PlatformResult::Warn(Cow::from(
            "WARNING: Unable to determine kernel version using uname(2)",
        ));
    }
    let release = to_str(&utsname.release[..]);
    let re = Regex::new(r"\D*(\d+)\.(\d+)").unwrap();
    if let Some(caps) = re.captures(release) {
        let version: Result<Vec<i64>, _> = (1..=2).map(|i| caps[i].parse()).collect();
        if let Ok(&[major, minor, ..]) = version.as_deref() {
            if (major, minor) < (3, 10) {
                return PlatformResult::Err(Cow::from(format!(
                    "The kernel used by this operating system version is unsupported:
    {release}
This release requires kernel version 3.10.0 or higher.
For more information, see support.10xgenomics.com/os-support.

To continue running this version for now, set TENX_IGNORE_DEPRECATED_OS=1
in your environment."
                )));
            }
        }
    }
    PlatformResult::Ok
}

// See ${sysroot}/usr/include/bits/confname.h
const _CS_GNU_LIBC_VERSION: c_int = 2;

extern "C" {
    fn confstr(name: c_int, buf: *mut c_char, len: libc::size_t) -> libc::size_t;
}

fn check_libc_version() -> PlatformResult {
    let libc_version = unsafe {
        let len = confstr(_CS_GNU_LIBC_VERSION, std::ptr::null_mut(), 0);
        let mut buf = vec![0; len];
        let _ = confstr(_CS_GNU_LIBC_VERSION, buf.as_mut_ptr(), len);
        buf
    };
    let libc_version = to_str(&libc_version);
    let re = Regex::new(r"\D*(\d+)\.(\d+)").unwrap();
    if let Some(caps) = re.captures(libc_version) {
        let version: Result<Vec<i64>, _> = (1..=2).map(|i| caps[i].parse()).collect();
        if let Ok(&[major, minor, ..]) = version.as_deref() {
            if (major, minor) < (2, 17) {
                return PlatformResult::Err(Cow::from(format!(
                    "The glibc version of this operating system version is unsupported:
    {libc_version}
This release requires libc version 2.17 or higher.
For more information, see support.10xgenomics.com/os-support.

To continue running this version for now, set TENX_IGNORE_DEPRECATED_OS=1
in your environment."
                )));
            }
        }
    }
    PlatformResult::Ok
}

fn check_cpu_features() -> PlatformResult {
    if !is_x86_feature_detected!("sse4.2") {
        return PlatformResult::Err(Cow::from(
            "The current CPU does not support sse4.2 instructions, and is no longer supported.

For more information, see
https://support.10xgenomics.com/os-support.

To continue running this version for now, set TENX_IGNORE_DEPRECATED_OS=1
in your environment.",
        ));
    }
    if !is_x86_feature_detected!("popcnt") {
        return PlatformResult::Err(Cow::from(
            "The current CPU does not support popcnt instructions, and is no longer supported.

For more information, see
https://support.10xgenomics.com/os-support.

To continue running this version for now, set TENX_IGNORE_DEPRECATED_OS=1
in your environment.",
        ));
    }
    if !is_x86_feature_detected!("avx") {
        return PlatformResult::Warn(Cow::from(
            "The current CPU does not support avx instructions.

Future versions of 10X Genomics software will not support CPUs older than Intel Xeon E3 (Sandy Bridge) or AMD Opteron FX (circa 2011).

For more information, see
https://support.10xgenomics.com/os-support.",
        ));
    }
    PlatformResult::Ok
}

pub fn oscheck() -> Result<()> {
    check_known_os().into_result()?;
    check_kernel_version().into_result()?;
    check_libc_version().into_result()?;
    check_cpu_features().into_result()?;
    Ok(())
}

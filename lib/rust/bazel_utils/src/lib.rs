//!
//! Module for handling bazel specific quirks such as path resolution
//!

use anyhow::{bail, Result};
use std::env;
use std::path::{Path, PathBuf};

pub fn canonicalize(path: &Path) -> Result<PathBuf> {
    use std::path::Component::CurDir;
    let path = if path.is_absolute() {
        path.to_path_buf()
    } else {
        env::current_dir()?.join(path)
    };
    Ok(path.components().filter(|&x| x != CurDir).collect())
}

/// Determines whether a path is a file and has executable permissions.
pub fn is_executable_file<P: AsRef<Path>>(path: P) -> bool {
    use libc::{access, X_OK};
    use std::ffi::CString;
    use std::os::unix::ffi::OsStrExt;
    if path.as_ref().is_file() {
        if let Ok(path) = CString::new(path.as_ref().as_os_str().as_bytes()) {
            return unsafe { access(path.as_c_str().as_ptr(), X_OK) } == 0;
        }
    }
    false
}

pub fn current_exe() -> Result<PathBuf, std::io::Error> {
    let exe = PathBuf::from(env::args_os().next().unwrap());
    if let Some(p) = exe.parent() {
        // For paths containing a /, like ./exe or foo/bar/exe, and only for
        // such paths, there is no PATH lookup.
        // In order to work within a bazel build, canonicalize the parent
        // directory but not on the executable file itself, since that
        // will resolve symlinks, which breaks things in a runfiles tree.
        if !p.as_os_str().is_empty() {
            let base = exe.file_name().unwrap();
            return Ok(p.canonicalize()?.join(base));
        }
    }

    // This could fail if you use relative paths in PATH and change your working
    // directory before calling this. Don't do that.
    env::var_os("PATH")
        .and_then(|paths| env::split_paths(&paths).find(|path| is_executable_file(path.join(&exe))))
        .and_then(|path| canonicalize(&path).ok())
        .map_or_else(std::env::current_exe, |path| Ok(path.join(&exe)))
}

pub fn runfiles_dir() -> Result<PathBuf> {
    let exe = current_exe()?;
    let exe_name = exe
        .file_name()
        .expect("Could not get the exe name")
        .to_str()
        .unwrap();
    let runfiles = exe
        .parent()
        .expect("No parent directory for exe")
        .join(format!("{exe_name}.runfiles"));
    if runfiles.exists() {
        Ok(runfiles)
    } else {
        bail!(
            "Could not find runfiles directory: {}",
            runfiles.display().to_string()
        );
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_canonicalize() -> Result<()> {
        let cwd = std::env::current_dir()?;
        let r1 = canonicalize(Path::new("/an/absolute/path/to/./a/../thing"))?;
        assert_eq!(PathBuf::from("/an/absolute/path/to/a/../thing"), r1);
        let r2 = canonicalize(Path::new("a/relative/path/to/./a/../thing"))?;
        assert_eq!(cwd.join("a/relative/path/to/a/../thing"), r2);
        Ok(())
    }
}

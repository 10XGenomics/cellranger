//! Module for handling bazel specific quirks such as path resolution
#![deny(missing_docs)]

use anyhow::{Result, bail};
use std::env;
use std::path::{Path, PathBuf};

/// Canonicalize a path
pub fn canonicalize(path: &Path) -> Result<PathBuf> {
    use std::path::Component::CurDir;
    let path = if path.is_absolute() {
        path.to_path_buf()
    } else {
        env::current_dir()?.join(path)
    };
    Ok(path.components().filter(|&x| x != CurDir).collect())
}

/// Return whether a path is an executable file.
pub fn is_executable_file(path: &Path) -> bool {
    path.is_file() && rustix::fs::access(path, rustix::fs::Access::EXEC_OK).is_ok()
}

/// Return the path to the current executable
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
        .and_then(|paths| {
            env::split_paths(&paths).find(|path| is_executable_file(&path.join(&exe)))
        })
        .and_then(|path| canonicalize(&path).ok())
        .map_or_else(env::current_exe, |path| Ok(path.join(&exe)))
}

/// Return the path to the runfiles directory
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
        bail!("Could not find runfiles directory: {}", runfiles.display());
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_canonicalize() -> Result<()> {
        let cwd = env::current_dir()?;
        let r1 = canonicalize(Path::new("/an/absolute/path/to/./a/../thing"))?;
        assert_eq!(PathBuf::from("/an/absolute/path/to/a/../thing"), r1);
        let r2 = canonicalize(Path::new("a/relative/path/to/./a/../thing"))?;
        assert_eq!(cwd.join("a/relative/path/to/a/../thing"), r2);
        Ok(())
    }
}

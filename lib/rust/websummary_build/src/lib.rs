use anyhow::{Context, Result};
use std::path::Path;
use tenx_websummary::WebSummaryBuildFiles;

pub fn read_to_string(fname: &Path) -> Result<String> {
    std::fs::read_to_string(fname).with_context(|| format!("While opening {fname:?} file"))
}

pub fn build_files() -> Result<WebSummaryBuildFiles<'static>> {
    let dist_folder = {
        // Option 1: Load from runfiles dir
        if let Ok(runfiles_dir) = bazel_utils::runfiles_dir() {
            runfiles_dir.join("cellranger/lib/python/websummary/dist/")
        } else {
            // Executable in lib/bin/exe and dist folder in `lib/python/websummary/dist
            bazel_utils::current_exe()?
                .parent()
                .unwrap()
                .parent()
                .unwrap()
                .join("python/websummary/dist/")
        }
    };

    build_files_in(&dist_folder)
}

pub fn build_files_in(dist_folder: &Path) -> Result<WebSummaryBuildFiles<'static>> {
    let script_js = read_to_string(&dist_folder.join("tenx-websummary-script.min.js"))?;
    let styles_css = read_to_string(&dist_folder.join("tenx-websummary-styles.min.css"))?;
    let template_html = read_to_string(&dist_folder.join("template.html"))?;

    Ok(WebSummaryBuildFiles::new(
        script_js,
        styles_css,
        template_html,
    ))
}

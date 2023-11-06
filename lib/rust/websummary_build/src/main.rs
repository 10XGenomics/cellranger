use anyhow::Result;
use tenx_websummary::WebSummaryBuildFiles;
use websummary_build::build_files;

fn main() -> Result<()> {
    let WebSummaryBuildFiles {
        script_js,
        styles_css,
        template_html,
    } = build_files()?;

    assert!(!script_js.is_empty());
    assert!(!styles_css.is_empty());
    assert!(!template_html.is_empty());
    Ok(())
}

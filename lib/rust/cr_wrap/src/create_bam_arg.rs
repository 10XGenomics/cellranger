use anyhow::{bail, Result};
use clap::{ArgAction, Parser};

#[derive(Parser, Debug, Clone, PartialEq, Eq)]
pub struct CreateBam {
    /// Enable or disable BAM file generation.
    /// Setting --create-bam=false reduces the total computation time and the size of
    /// the output directory (BAM file not generated).
    /// We recommend setting --create-bam=true if unsure.
    /// See https://10xgen.com/create-bam for additional guidance.
    #[clap(long = "create-bam", value_name= "true|false", required = true, num_args = 1, action = ArgAction::Set)]
    create_bam: bool,

    // Maintain the previous argument so that we can cleanly guide the user to
    // the new version.
    #[clap(long = "no-bam", hide = true)]
    no_bam: bool,
}

impl CreateBam {
    /// Validate that only an accepted arg combo was provided.
    /// Return the collapsed value for whether or not we should produce a BAM file.
    pub fn validated(&self) -> Result<bool> {
        if self.no_bam {
            bail!(
                "--no-bam has been replaced with --create-bam=false to disable BAM file production"
            )
        }
        Ok(self.create_bam)
    }
}

#[cfg(test)]
mod test {
    use super::*;

    fn parse_args(args: &[&str]) -> Result<bool> {
        let mut full_args = vec!["program"];
        full_args.extend(args);
        CreateBam::try_parse_from(full_args)?.validated()
    }

    #[test]
    fn test_create_bam() {
        assert!(parse_args(&[]).is_err());
        assert!(parse_args(&["--create-bam", "--no-bam"]).is_err());
        assert!(parse_args(&["--no-bam"]).is_err());
        assert!(parse_args(&["--create-bam"]).is_err());
        assert!(parse_args(&["--create-bam=true"]).unwrap());
        assert!(!parse_args(&["--create-bam=false"]).unwrap());
        assert!(parse_args(&["--create-bam", "true"]).unwrap());
        assert!(!parse_args(&["--create-bam", "false"]).unwrap());
    }
}

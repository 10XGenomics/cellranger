use anyhow::{bail, Result};
use itertools::Itertools;
use serde::Serialize;
use std::fmt::{Debug, Formatter};
use std::str::FromStr;

// The help text when user inputs an invalid
// chemistry
enum UserHelp {
    Visible(&'static str),
    Hidden,
}

const ALLOWED_COUNT_CHEM_INPUTS: [(&str, UserHelp); 17] = [
    ("auto", UserHelp::Visible("auto detection (default)")),
    ("threeprime", UserHelp::Visible("Single Cell 3'")),
    ("fiveprime", UserHelp::Visible("Single Cell 5'")),
    ("SC3P_auto", UserHelp::Hidden),
    ("SC5P_auto", UserHelp::Hidden),
    ("SC3Pv1", UserHelp::Visible("Single Cell 3'v1")),
    ("SC3Pv2", UserHelp::Visible("Single Cell 3'v2")),
    ("SC3Pv3", UserHelp::Visible("Single Cell 3'v3")),
    ("SC3Pv3LT", UserHelp::Visible("Single Cell 3'v3 LT")),
    ("SC3Pv3HT", UserHelp::Visible("Single Cell 3'v3 HT")),
    ("SC5P-PE", UserHelp::Visible("Single Cell 5' paired end")),
    ("SC5P-R2", UserHelp::Visible("Single Cell 5' R2-only")),
    ("SC5PHT", UserHelp::Hidden),
    ("SC5P-R1", UserHelp::Hidden),
    (
        "SC-FB",
        UserHelp::Visible("Single Cell Antibody-only 3' v2 or 5'"),
    ),
    ("SFRP", UserHelp::Hidden),
    ("ARC-v1", UserHelp::Visible("GEX portion only of multiome")),
];

#[derive(Serialize, Clone, PartialEq, Eq)]
pub struct CountChemistryArg(pub String);

impl FromStr for CountChemistryArg {
    type Err = anyhow::Error;

    fn from_str(s: &str) -> Result<CountChemistryArg> {
        if ALLOWED_COUNT_CHEM_INPUTS.iter().any(|(chem, _)| *chem == s) {
            Ok(CountChemistryArg(s.to_string()))
        } else {
            bail!(
                "{s} is an invalid input to `--chemistry`. Supported options are:\n - {}",
                ALLOWED_COUNT_CHEM_INPUTS
                    .iter()
                    .filter_map(|(chem, help)| match help {
                        UserHelp::Visible(h) => Some(format!("{chem} for {h}")),
                        UserHelp::Hidden => None,
                    })
                    .join("\n - ")
            );
        }
    }
}

impl Debug for CountChemistryArg {
    fn fmt(&self, f: &mut Formatter<'_>) -> Result<(), std::fmt::Error> {
        Debug::fmt(&self.0, f)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use cr_types::chemistry::AutoOrRefinedChemistry;
    #[test]
    fn test_parse_allowed() {
        assert!(ALLOWED_COUNT_CHEM_INPUTS
            .iter()
            .all(|(chem, _)| chem.parse::<AutoOrRefinedChemistry>().is_ok()));
    }

    #[test]
    fn test_parse_count_chemistry_arg() {
        for chem in [
            "SC3Pv3", "SC3Pv3LT", "SC3Pv3HT", "auto", "SC5P-R2", "SC5PHT", "ARC-v1",
        ] {
            assert_eq!(
                chem.parse::<CountChemistryArg>().unwrap(),
                CountChemistryArg(chem.to_string())
            );
        }
        assert!("SCVDJ".parse::<CountChemistryArg>().is_err());
        assert_eq!(
            format!("{:?}", CountChemistryArg("fiveprime".into())),
            r#""fiveprime""#
        );
    }
}

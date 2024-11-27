//! Types shared between websummary derive macros and websummary construction.
use anyhow::{bail, ensure, Result};
use proc_macro2::TokenStream;
use quote::quote;
use serde::{Deserialize, Serialize};
use std::collections::HashSet;
use std::iter::zip;

#[derive(Serialize, Deserialize, Debug, Copy, Clone)]
#[serde(rename_all = "snake_case")]
pub enum AlertIfMetricIs {
    GreaterThanOrEqual,
    LessThanOrEqual,
}

impl quote::ToTokens for AlertIfMetricIs {
    fn to_tokens(&self, tokens: &mut TokenStream) {
        tokens.extend(match *self {
            Self::GreaterThanOrEqual => {
                quote![::cr_types::websummary::AlertIfMetricIs::GreaterThanOrEqual]
            }
            Self::LessThanOrEqual => {
                quote![::cr_types::websummary::AlertIfMetricIs::LessThanOrEqual]
            }
        });
    }
}

// Should be synced with the context in cr_websummary
#[derive(Serialize, Deserialize, Debug, Clone, Default)]
pub struct AlertConditions {
    pub is_rtl: Option<bool>,
    pub include_introns: Option<bool>,
}

impl AlertConditions {
    pub fn vec(&self) -> Vec<Option<bool>> {
        let AlertConditions {
            include_introns,
            is_rtl,
        } = self;
        vec![*include_introns, *is_rtl]
    }
}

impl quote::ToTokens for AlertConditions {
    fn to_tokens(&self, tokens: &mut TokenStream) {
        let Self {
            is_rtl,
            include_introns,
        } = self;
        let is_rtl = quote_option(is_rtl);
        let include_introns = quote_option(include_introns);
        tokens.extend(quote![
            ::cr_types::websummary::AlertConditions {
                is_rtl: #is_rtl,
                include_introns: #include_introns,
            }
        ]);
    }
}

#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct AlertConfig {
    pub error_threshold: Option<f64>,
    pub warn_threshold: Option<f64>,
    pub if_metric_is: Option<AlertIfMetricIs>,
    pub warn_title: Option<String>,
    // Use warn_title if missing
    pub error_title: Option<String>,
    pub detail: String,
    #[serde(default)]
    pub conditions: AlertConditions,
}

impl quote::ToTokens for AlertConfig {
    fn to_tokens(&self, tokens: &mut TokenStream) {
        let Self {
            error_threshold,
            warn_threshold,
            if_metric_is,
            warn_title,
            error_title,
            detail,
            conditions,
        } = self;
        let error_threshold = quote_option(error_threshold);
        let warn_threshold = quote_option(warn_threshold);
        let if_metric_is = quote_option(if_metric_is);
        let warn_title = quote_string_option(warn_title);
        let error_title = quote_string_option(error_title);
        tokens.extend(quote![
            ::cr_types::websummary::AlertConfig {
                error_threshold: #error_threshold,
                warn_threshold: #warn_threshold,
                if_metric_is: #if_metric_is,
                warn_title: #warn_title,
                error_title: #error_title,
                detail: #detail.to_string(),
                conditions: #conditions,
            }
        ]);
    }
}

impl AlertConfig {
    /// Return the comparison operation to determine if this alert should fire.
    pub fn alert_if(&self) -> Result<AlertIfMetricIs> {
        match (self.error_threshold, self.warn_threshold, self.if_metric_is) {
            (Some(e), Some(w), None) => {
                if e < w {
                    Ok(AlertIfMetricIs::LessThanOrEqual)
                } else if e > w {
                    Ok(AlertIfMetricIs::GreaterThanOrEqual)
                } else {
                    bail!(
                        "error threshold ({e}) and warning threshold({w}) do not \
                        have a strict < or > relation"
                    );
                }
            }
            (Some(_), None, Some(if_metric_is)) | (None, Some(_), Some(if_metric_is)) => {
                Ok(if_metric_is)
            }
            (Some(_), Some(_), Some(_)) => {
                bail!(
                    "do not specify `if_metric_is` when both error and warn thresholds are \
                    specified as it is automatically inferred"
                );
            }
            (Some(_), None, None) | (None, Some(_), None) => {
                bail!(
                    "please specify `if_metric_is` in the alert as one of \
                    \"greater_than_or_equal\" or \"less_than_or_equal\"; with only a single threshold \
                    specified, it cannot be automatically inferred"
                );
            }
            (None, None, _) => bail!(
                "at least one of error_threshold or warn_threshold needs to be specified \
                in the alert"
            ),
        }
    }

    pub fn warn_title(&self, name: &str) -> &str {
        if let Some(ref w) = self.warn_title {
            return w;
        }
        if let Some(ref e) = self.error_title {
            return e;
        }
        panic!("At least one of warn_title or error_title is required for the alerts under {name}")
    }
    pub fn error_title(&self, name: &str) -> &str {
        if let Some(ref e) = self.error_title {
            return e;
        }
        if let Some(ref w) = self.warn_title {
            return w;
        }
        panic!("At least one of warn_title or error_title is required for the alerts under {name}")
    }
}

#[derive(Serialize, Deserialize, Debug, Clone)]
#[serde(deny_unknown_fields)]
pub struct MetricConfig {
    pub header: String,
    pub json_key: Option<String>,
    #[serde(rename = "type")]
    pub ty: String,
    #[serde(default)] // false
    pub optional: bool,
    pub help: Option<String>,
    #[serde(default)] // Empty vec
    pub alerts: Vec<AlertConfig>,
}

impl MetricConfig {
    pub fn validate(&self) -> Result<()> {
        // Ensure that expected transformed type is compatible with alerts.
        // FIXME extract transformer keys as constants.
        ensure!(
            self.alerts.is_empty() || self.ty != "String",
            "metric \"{}\" cannot be used in alerts",
            self.header,
        );

        let conditions: Vec<_> = self.alerts.iter().map(|a| &a.conditions).collect();
        ensure!(
            check_exclusive_conditions(&conditions),
            "The conditions for for the metric {} are not exclusive:\n{conditions:#?}.\
            \n Please specify mutually exclusive conditions for each alert.",
            &self.header
        );
        Ok(())
    }
}

impl quote::ToTokens for MetricConfig {
    fn to_tokens(&self, tokens: &mut TokenStream) {
        let Self {
            header,
            json_key,
            ty,
            optional,
            help,
            alerts,
        } = self;
        let json_key = quote_string_option(json_key);
        let help = quote_string_option(help);
        tokens.extend(quote![
            ::cr_types::websummary::MetricConfig {
                header: #header.to_string(),
                json_key: #json_key,
                ty: #ty.to_string(),
                optional: #optional,
                help: #help,
                alerts: vec![
                    #(#alerts),*
                ]
            }
        ]);
    }
}

fn check_exclusive_conditions(conditions: &[&AlertConditions]) -> bool {
    if conditions.len() < 2 {
        return true;
    }
    let first = conditions[0].vec();
    let mut seen_ids = HashSet::with_capacity(conditions.len());
    for cond in conditions {
        let cond_vec = cond.vec();
        let mut this_id = Vec::new();
        for (left, right) in zip(&first, cond_vec) {
            // Both should be None or Some
            match (left, right) {
                (Some(_), Some(r)) => this_id.push(r),
                (None, None) => {}
                _ => return false,
            }
        }
        if seen_ids.contains(&this_id) {
            return false;
        }
        seen_ids.insert(this_id);
    }
    true
}

fn quote_option<T: quote::ToTokens>(arg: &Option<T>) -> proc_macro2::TokenStream {
    match arg {
        None => quote![None],
        Some(x) => quote![Some(#x)],
    }
}

fn quote_string_option(arg: &Option<String>) -> proc_macro2::TokenStream {
    match arg {
        None => quote![None],
        Some(x) => quote![Some(#x.to_string())],
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_check_exclusive_contexts() {
        assert!(check_exclusive_conditions(&[]));
        assert!(check_exclusive_conditions(&[&AlertConditions {
            is_rtl: None,
            include_introns: None
        }]));
        assert!(check_exclusive_conditions(&[&AlertConditions {
            is_rtl: None,
            include_introns: None
        }]));

        assert!(!check_exclusive_conditions(&[
            &AlertConditions {
                is_rtl: None,
                include_introns: None
            },
            &AlertConditions {
                include_introns: None,
                is_rtl: None
            }
        ]));

        assert!(check_exclusive_conditions(&[
            &AlertConditions {
                is_rtl: Some(true),
                include_introns: None
            },
            &AlertConditions {
                is_rtl: Some(false),
                include_introns: None
            }
        ]));

        assert!(!check_exclusive_conditions(&[
            &AlertConditions {
                is_rtl: Some(true),
                include_introns: None
            },
            &AlertConditions {
                is_rtl: Some(true),
                include_introns: None
            }
        ]));
    }
}

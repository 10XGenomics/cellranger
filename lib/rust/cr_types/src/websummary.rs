//! Types shared between websummary derive macros and websummary construction.
use proc_macro2::TokenStream;
use quote::quote;
use serde::{Deserialize, Serialize};

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

impl AlertIfMetricIs {
    fn symbol(&self) -> proc_macro2::TokenStream {
        match *self {
            AlertIfMetricIs::GreaterThanOrEqual => quote![>=],
            AlertIfMetricIs::LessThanOrEqual => quote![<=],
        }
    }
}

// Should be synced with the context in cr_websummary
#[derive(Serialize, Deserialize, Debug, Clone, Default)]
pub struct AlertConditions {
    pub is_hybrid_capture: Option<bool>,
    pub is_rtl: Option<bool>,
    pub is_lt_chemistry: Option<bool>,
    pub is_arc_chemistry: Option<bool>,
    pub include_introns: Option<bool>,
}

impl AlertConditions {
    pub fn vec(&self) -> Vec<Option<bool>> {
        let AlertConditions {
            is_hybrid_capture,
            is_lt_chemistry,
            is_arc_chemistry,
            include_introns,
            is_rtl,
        } = self;
        vec![
            *is_hybrid_capture,
            *is_lt_chemistry,
            *is_arc_chemistry,
            *include_introns,
            *is_rtl,
        ]
    }
}

impl quote::ToTokens for AlertConditions {
    fn to_tokens(&self, tokens: &mut TokenStream) {
        let Self {
            is_hybrid_capture,
            is_rtl,
            is_lt_chemistry,
            is_arc_chemistry,
            include_introns,
        } = self;
        let is_hybrid_capture = quote_option(is_hybrid_capture);
        let is_rtl = quote_option(is_rtl);
        let is_lt_chemistry = quote_option(is_lt_chemistry);
        let is_arc_chemistry = quote_option(is_arc_chemistry);
        let include_introns = quote_option(include_introns);
        tokens.extend(quote![
            ::cr_types::websummary::AlertConditions {
                is_hybrid_capture: #is_hybrid_capture,
                is_rtl: #is_rtl,
                is_lt_chemistry: #is_lt_chemistry,
                is_arc_chemistry: #is_arc_chemistry,
                include_introns: #include_introns,
            }
        ]);
    }
}

#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct AlertConfig {
    #[serde(default)]
    pub rank: u8,
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
            rank,
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
                rank: #rank,
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
    pub fn symbol(&self, name: &str) -> proc_macro2::TokenStream {
        match (self.error_threshold, self.warn_threshold) {
            (Some(e), Some(w)) => {
                assert!(
                    self.if_metric_is.is_none(),
                    "Do not specify `if_metric_is` in the alert for {name}. When both error and warn \
                    thresholds are specified, it is automatically inferred."
                );
                if e < w {
                    quote![<=]
                } else if e > w {
                    quote![>=]
                } else {
                    panic!(
                        "ERROR: Error threshold ({e}) and warning threshold({w}) do not have \
                        strict < or > relation for {name}"
                    );
                }
            }
            (Some(_), None) => {
                assert!(
                    self.if_metric_is.is_some(),
                    "Please specify `if_metric_is` in the alert for {name} as one of \
                    \"greater_than_or_equal\" or \"less_than_or_equal\". With only an error threshold \
                    specified, it cannot be automatically inferred."
                );
                self.if_metric_is.unwrap().symbol()
            }
            (None, Some(_)) => {
                assert!(
                    self.if_metric_is.is_some(),
                    "Please specify `if_metric_is` in the alert for {name} as one of \
                    \"greater_than_or_equal\" or \"less_than_or_equal\". With only a warn threshold \
                    specified, it cannot be automatically inferred."
                );
                self.if_metric_is.unwrap().symbol()
            }
            (None, None) => panic!(
                "At least one of error_threshold or warn_threshold needs to be specified \
                in the alert for {name}"
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

use cr_types::websummary::{AlertConditions, MetricConfig};
use heck::ToUpperCamelCase;
use itertools::Itertools;
use proc_macro::TokenStream;
use quote::{format_ident, quote};
use serde::Deserialize;
use std::collections::{BTreeMap, HashSet};
use std::iter::zip;
use std::path::PathBuf;
use syn::{parse_macro_input, Fields, ItemStruct, LitStr};

#[derive(Deserialize, Debug)]
struct CardWithTableToml {
    title: String,
    help: Option<String>,
    entries: Vec<String>,
    #[serde(flatten)]
    entry_info: BTreeMap<String, MetricConfig>,
}

impl CardWithTableToml {
    fn validate(&self) {
        for e in &self.entries {
            assert!(
                self.entry_info.contains_key(e),
                "Information for {} does not exist in table '{}'",
                e,
                self.title
            );
        }
        assert!(
            self.entries.len() == self.entry_info.len(),
            "Duplicate entries listed under table '{}'",
            &self.title
        );
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

#[proc_macro]
pub fn make_tables(item: TokenStream) -> TokenStream {
    let input = parse_macro_input!(item as LitStr);

    let cwd = PathBuf::from(std::env::var("CARGO_MANIFEST_DIR").unwrap());

    let file_path = cwd.join(input.value());
    assert!(
        file_path.exists(),
        "File {} does not exist. Relative path: {}, cwd: {}.",
        file_path.display(),
        input.value(),
        cwd.display(),
    );
    let file_content = std::fs::read_to_string(&file_path).unwrap();

    let parsed: BTreeMap<String, CardWithTableToml> = toml::from_str(&file_content).unwrap();

    // eprintln!("{:#?}", parsed);

    let mut q = quote![];
    for (k, v) in parsed {
        v.validate();
        let table_struct_name = format_ident!("{}Table", k.to_upper_camel_case());
        let row_struct_name = format_ident!("{}Row", k.to_upper_camel_case());
        let title = &v.title;

        let mut row_struct_fields = quote![];
        let mut table_rows = quote![
            let mut rows = Vec::new();
        ];
        let mut metric_headers_and_help = quote![
            let mut headers = Vec::new();
            let mut help_data = Vec::new();
        ];
        if let Some(ref h) = v.help {
            metric_headers_and_help = quote![
                #metric_headers_and_help
                help_data.push(::cr_websummary::TermDesc::with_one_desc("", #h.to_string()));
            ];
        }
        let mut row_struct_field_from_metrics = quote![];
        let mut row_struct_to_json_summary_items = quote![];
        for e in &v.entries {
            let name = format_ident!("{}", e);
            let info = &v.entry_info[e];
            let ty = format_ident!("{}", info.ty);
            let json_key = info.json_key.as_ref().unwrap_or(e);
            row_struct_fields = quote![
                #row_struct_fields
                #[serde(alias = #json_key)]
                pub #name: ::std::option::Option<#ty>,
            ];
            row_struct_field_from_metrics = quote![
                #row_struct_field_from_metrics
                #name: match val.get(#json_key) {
                    Some(v) => #ty::try_from_value(v)?,
                    None => None,
                },
            ];
            row_struct_to_json_summary_items = quote![
                #row_struct_to_json_summary_items
                ::cr_websummary::multi::websummary::JsonMetricSummary {
                    key: stringify!(#name).to_string(),
                    value: serde_json::json!(&self.#name),
                    category: String::new(),
                    library_type: String::new(),
                    config: #info,
                },
            ];
            table_rows = if info.optional {
                quote![
                    #table_rows
                    if let Some(v) = src.#name {
                        rows.push(v.make_pretty());
                    }
                ]
            } else {
                quote![
                    #table_rows
                    if let Some(v) = src.#name {
                        rows.push(v.make_pretty());
                    } else {
                        rows.push("---".to_string());
                    }
                ]
            };

            let header = &info.header;
            let this_help = match info.help {
                Some(ref h) => quote![
                    help_data.push(::cr_websummary::TermDesc::with_one_desc(#header, #h));
                ],
                None => quote![],
            };

            metric_headers_and_help = if info.optional {
                quote![
                    #metric_headers_and_help
                    let visible = src.0.iter().any(|row| row.#name.is_some());
                    if visible {
                        headers.push(#header.to_string());
                        #this_help
                    }
                ]
            } else {
                quote![
                    #metric_headers_and_help
                    headers.push(#header.to_string());
                    #this_help
                ]
            }
        }

        let metric_alerts: Vec<(_, _, Vec<_>)> = v
            .entries
            .iter()
            .flat_map(|e| -> Vec<_> {
                v.entry_info[e]
                    .alerts
                    .iter()
                    .sorted_by_key(|a| a.rank)
                    .group_by(|a| a.rank)
                    .into_iter()
                    .map(|(_rank, group)| (e, v.entry_info[e].optional, group.collect()))
                    .collect()
            })
            .collect();

        let mut alert_quote = quote![
            use ::cr_websummary::MakePretty;
            let mut alert_specs = Vec::new();
        ];
        for (name, _, alerts) in metric_alerts {
            let name_ident = format_ident!("{}", name);

            let conditions: Vec<_> = alerts.iter().map(|a| &a.conditions).collect();
            assert!(
                check_exclusive_conditions(&conditions),
                "The conditions for alerts with rank {} for the metric {} are not exclusive:\n {:#?}.\
                \n If these alerts are independent, specify a different rank for the alerts. \
                Otherwise specify mutually exclusive conditions for each alert.",
                alerts[0].rank,
                name,
                conditions
            );

            let alert_val = quote![
                let #name_ident = self.#name_ident.map(|m| (m.as_f64(), m.make_pretty()));
            ];

            for alert in alerts {
                let error_title = alert.error_title(name);
                let warn_title = alert.warn_title(name);
                let detail = &alert.detail;
                let symbol = alert.symbol(name);
                let mut alert_specs = quote![
                    let mut has_error = false;
                ];
                if let Some(error_threshold) = alert.error_threshold {
                    alert_specs = quote![
                        #alert_specs
                        if val #symbol #error_threshold {
                            alert_specs.push(::cr_websummary::alert::AlertSpec {
                                level: ::cr_websummary::alert::AlertLevel::Error,
                                title: #error_title.to_string(),
                                formatted_value: formatted_value.clone(),
                                message: #detail.to_string()
                            });
                            has_error = true;
                        }
                    ];
                }
                if let Some(warn_threshold) = alert.warn_threshold {
                    alert_specs = quote![
                        #alert_specs
                        if !has_error && (val #symbol #warn_threshold) {
                            alert_specs.push(::cr_websummary::alert::AlertSpec {
                                level: ::cr_websummary::alert::AlertLevel::Warn,
                                title: #warn_title.to_string(),
                                formatted_value,
                                message: #detail.to_string()
                            });
                        }
                    ];
                }
                let AlertConditions {
                    is_hybrid_capture,
                    is_lt_chemistry,
                    is_arc_chemistry,
                    include_introns,
                    is_rtl,
                } = alert.conditions;
                let mut conditions_quote = quote![
                    let mut conditions_are_met = true;
                ];
                if let Some(is_hybrid_capture) = is_hybrid_capture {
                    conditions_quote = quote![
                        #conditions_quote
                        conditions_are_met = conditions_are_met && ctx.is_hybrid_capture == #is_hybrid_capture;
                    ];
                }
                if let Some(include_introns) = include_introns {
                    conditions_quote = quote![
                        #conditions_quote
                        conditions_are_met = conditions_are_met && ctx.include_introns == #include_introns;
                    ];
                }
                if let Some(lt_chemistry) = is_lt_chemistry {
                    conditions_quote = quote![
                        #conditions_quote
                        conditions_are_met = conditions_are_met && ctx.is_lt_chemistry == #lt_chemistry;
                    ];
                }
                if let Some(arc_chemistry) = is_arc_chemistry {
                    conditions_quote = quote![
                        #conditions_quote
                        conditions_are_met = conditions_are_met && ctx.is_arc_chemistry == #arc_chemistry;
                    ];
                }
                if let Some(rtl) = is_rtl {
                    conditions_quote = quote![
                        #conditions_quote
                        conditions_are_met = conditions_are_met && ctx.is_rtl == #rtl;
                    ];
                }
                alert_quote = quote![
                    #alert_quote
                    #alert_val
                    if let Some((val, formatted_value)) = #name_ident {
                        #conditions_quote
                        if conditions_are_met {
                            #alert_specs
                        }
                    }
                ];
            }
        }

        q = quote![
            #q
            #[automatically_derived]
            #[derive(Debug, Default, Clone, Serialize, Deserialize, PartialEq)]
            pub struct #table_struct_name(pub Vec<#row_struct_name>);

            impl #table_struct_name {
                /// Create a table with a single row from the json
                ///
                /// NOTE: Clones the val
                pub fn from_json_value(val: &::serde_json::value::Value) -> Self {
                    #table_struct_name(vec![::serde_json::from_value(val.clone()).unwrap()])
                }

                /// Create a table with a single row from a metric hashmap
                pub fn from_metrics<S: ::std::hash::BuildHasher>(val: &::std::collections::HashMap<::std::string::String, ::serde_json::value::Value, S>) -> Result<Self, ::anyhow::Error> {
                    Ok(#table_struct_name(vec![#row_struct_name::from_metrics(val)?]))
                }

            }

            #[automatically_derived]
            impl ::cr_websummary::alert::Alert for #table_struct_name {
                fn alerts(&self, ctx: &::cr_websummary::alert::AlertContext) -> Vec<::cr_websummary::alert::AlertSpec> {
                    self.0.iter().map(|r| r.alerts(ctx)).flatten().collect()
                }
            }

            #[automatically_derived]
            impl ::cr_websummary::multi::websummary::ToJsonSummary for #table_struct_name {
                fn to_json_summary(&self) -> Vec<::cr_websummary::multi::websummary::JsonMetricSummary> {
                    self.0.iter().map(::cr_websummary::multi::websummary::ToJsonSummary::to_json_summary).flatten().collect()
                }
            }

            #[automatically_derived]
            impl From<#table_struct_name> for ::cr_websummary::CardWithMetric {
                fn from(src: #table_struct_name) -> ::cr_websummary::CardWithMetric {
                    #metric_headers_and_help
                    let table = ::cr_websummary::GenericTable {
                        header: Some(headers),
                        rows: src.0.into_iter().map(|row| row.into()).collect(),
                    };
                    let help = ::cr_websummary::TitleWithTermDesc {
                        title: #title.to_string(),
                        data: help_data,
                    };

                    ::cr_websummary::CardWithMetric { table, help }
                }
            }

            // One row within the table
            #[automatically_derived]
            #[allow(clippy::derive_partial_eq_without_eq)]
            #[derive(Debug, Default, Clone, Serialize, Deserialize, PartialEq)]
            pub struct #row_struct_name {
                #row_struct_fields
            }

            impl #row_struct_name {
                /// Create a row from the json
                pub fn from_metrics<S: ::std::hash::BuildHasher>(val: &::std::collections::HashMap<::std::string::String, ::serde_json::value::Value, S>) -> Result<Self, ::anyhow::Error> {
                    use ::cr_websummary::value::TryFromValue;
                    Ok(#row_struct_name {
                        #row_struct_field_from_metrics
                    })
                }
            }

            #[automatically_derived]
            impl ::cr_websummary::multi::websummary::ToJsonSummary for #row_struct_name {
                fn to_json_summary(&self) -> Vec<::cr_websummary::multi::websummary::JsonMetricSummary> {
                    vec![
                        #row_struct_to_json_summary_items
                    ]
                }
            }

            #[automatically_derived]
            impl ::cr_websummary::alert::Alert for #row_struct_name {
                fn alerts(&self, ctx: &::cr_websummary::alert::AlertContext) -> Vec<::cr_websummary::alert::AlertSpec> {
                    #alert_quote
                    alert_specs
                }
            }

            #[automatically_derived]
            impl From<#row_struct_name> for ::cr_websummary::TableRow {
                fn from(src: #row_struct_name) -> ::cr_websummary::TableRow {
                    use ::cr_websummary::MakePretty;
                    #table_rows
                    ::cr_websummary::TableRow(rows)
                }
            }
        ];
    }

    let file_str = file_path.display().to_string();
    let ts: TokenStream = quote![
        pub const TABLES_TOML: &str = include_str!(#file_str);
        #q
    ]
    .into();
    ts
}

const ALERT_NOT_ON_NAMED_STRUCT_ERROR: &str =
    r#"#[derive(Alert)] can only be used on structs with named fields."#;

#[proc_macro_derive(Alert)]
pub fn websummary_alert(item: proc_macro::TokenStream) -> proc_macro::TokenStream {
    // ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    // STEP 1
    // ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    // Make sure that #[derive(Alert)] is used on struct
    // The way we achieve it is to try parsing the input TokenStrean as `ItemStruct`
    // and checking the parse result
    let Ok(item_struct) = syn::parse::<ItemStruct>(item.clone()) else {
        let span = proc_macro2::TokenStream::from(item);
        return syn::Error::new_spanned(span, ALERT_NOT_ON_NAMED_STRUCT_ERROR)
            .to_compile_error()
            .into();
    };

    // ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    // STEP 2
    // ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    // Collect the fields of the struct. If we get a unit struct or tuple struct
    // produce a compiler error.
    let fields = match item_struct.fields.clone() {
        Fields::Named(f) => f.named,
        _ => {
            let span = proc_macro2::TokenStream::from(item);
            return syn::Error::new_spanned(span, ALERT_NOT_ON_NAMED_STRUCT_ERROR)
                .to_compile_error()
                .into();
        }
    };

    // ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    // STEP 3
    // ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    // Generate tokenstream for `alerts` calls for each field
    let mut vec_inner = Vec::new();
    for field in fields {
        let name = field.ident.clone().unwrap();
        vec_inner.push(quote![
            self.#name.alerts(ctx)
        ]);
    }

    // ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    // STEP 4
    // ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    // Generate the `impl Alert` token stream
    // Handle generics in the struct
    let (impl_generics, ty_generics, where_clause) = item_struct.generics.split_for_impl();
    let item_ident = item_struct.ident.clone();

    quote![
        #[automatically_derived]
        impl #impl_generics ::cr_websummary::alert::Alert for #item_ident #ty_generics #where_clause {
            fn alerts(&self, ctx: &::cr_websummary::alert::AlertContext) -> Vec<::cr_websummary::alert::AlertSpec> {
                let mut alert_specs = Vec::new();
                #(alert_specs.extend(#vec_inner);)*
                alert_specs
            }
        }
    ].into()
}

/// ToCsvRows derive macro for generating metrics CSV from Websummary JSON Rust data structure
/// The websummary structure is a bunch of nested structs holding metrics in tables
/// We use this proc macro to derive for structs representing the Websummary Tabs the ability to convert it into a vector of CSV rows holding the metrics
/// The implementation provided marches over the fields of the websummary tab, and converts into CSV those fields that are MetricCards
/// The first two columns of the final metrics CSV ("Sample or Library" and "Library Type" are not filled out by this implementation),
/// And rather are implemented in the LibraryWebSummary and SampleWebSummary structs
/// and those columns are filled out there for each websummary tab.
#[proc_macro_derive(ToCsvRows)]
pub fn derive_to_csv_rows(item: TokenStream) -> TokenStream {
    let item_struct = parse_macro_input!(item as ItemStruct);
    let ident = item_struct.ident;

    let mut append_rows = quote![let mut csv_rows: Vec<Vec<String>> = vec![];]; // this chunk of code appends the CSV rows for each struct field

    for field_name in item_struct
        .fields
        .iter()
        .map(|f| f.ident.as_ref().expect("No identifier for struct field."))
    {
        append_rows = quote![
            #append_rows
            csv_rows.append(&mut self
                .#field_name
                .to_csv_rows());
        ];
    }
    quote! [
        #[automatically_derived]
        impl ToCsvRows for #ident {
            fn to_csv_rows(self) -> Vec<Vec<String>> {
                #append_rows;
                csv_rows
            }
        }
    ]
    .into()
}

/// Derive an implementation of ToJsonSummary by reflecting all struct fields.
/// Every struct member must implement ToJsonSummary.
#[proc_macro_derive(ToJsonSummary)]
pub fn derive_to_json_summary(item: TokenStream) -> TokenStream {
    let item_struct = parse_macro_input!(item as ItemStruct);
    let ident = item_struct.ident;

    let mut append_rows = quote![let mut vals: Vec<JsonMetricSummary> = vec![];];

    for field_name in item_struct
        .fields
        .iter()
        .map(|f| f.ident.as_ref().expect("No identifier for struct field."))
    {
        append_rows = quote![
            #append_rows
            vals.append(&mut self
                .#field_name
                .to_json_summary());
        ];
    }
    quote! [
        #[automatically_derived]
        impl ToJsonSummary for #ident {
            fn to_json_summary(&self) -> Vec<JsonMetricSummary> {
                #append_rows;
                vals
            }
        }
    ]
    .into()
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_check_exclusive_contexts() {
        assert!(check_exclusive_conditions(&[]));
        assert!(check_exclusive_conditions(&[&AlertConditions {
            is_hybrid_capture: None,
            is_lt_chemistry: None,
            is_arc_chemistry: None,
            is_rtl: None,
            include_introns: None
        }]));
        assert!(check_exclusive_conditions(&[&AlertConditions {
            is_hybrid_capture: Some(true),
            is_lt_chemistry: None,
            is_arc_chemistry: None,
            is_rtl: None,
            include_introns: None
        }]));

        assert!(!check_exclusive_conditions(&[
            &AlertConditions {
                is_hybrid_capture: None,
                is_lt_chemistry: None,
                is_arc_chemistry: None,
                is_rtl: None,
                include_introns: None
            },
            &AlertConditions {
                is_hybrid_capture: None,
                is_lt_chemistry: None,
                is_arc_chemistry: None,
                include_introns: None,
                is_rtl: None
            }
        ]));

        assert!(check_exclusive_conditions(&[
            &AlertConditions {
                is_hybrid_capture: Some(true),
                is_lt_chemistry: None,
                is_arc_chemistry: None,
                is_rtl: None,
                include_introns: None
            },
            &AlertConditions {
                is_hybrid_capture: Some(false),
                is_lt_chemistry: None,
                is_arc_chemistry: None,
                is_rtl: None,
                include_introns: None
            }
        ]));

        assert!(!check_exclusive_conditions(&[
            &AlertConditions {
                is_hybrid_capture: None,
                is_lt_chemistry: None,
                is_arc_chemistry: None,
                is_rtl: None,
                include_introns: None
            },
            &AlertConditions {
                is_hybrid_capture: Some(false),
                is_lt_chemistry: None,
                is_arc_chemistry: None,
                is_rtl: None,
                include_introns: None
            }
        ]));

        assert!(check_exclusive_conditions(&[
            &AlertConditions {
                is_hybrid_capture: None,
                is_lt_chemistry: Some(true),
                is_arc_chemistry: None,
                is_rtl: None,
                include_introns: None
            },
            &AlertConditions {
                is_hybrid_capture: None,
                is_lt_chemistry: Some(false),
                is_arc_chemistry: None,
                is_rtl: None,
                include_introns: None
            }
        ]));

        assert!(!check_exclusive_conditions(&[
            &AlertConditions {
                is_hybrid_capture: None,
                is_lt_chemistry: Some(true),
                is_arc_chemistry: None,
                is_rtl: None,
                include_introns: None
            },
            &AlertConditions {
                is_hybrid_capture: None,
                is_lt_chemistry: Some(true),
                is_arc_chemistry: None,
                is_rtl: None,
                include_introns: None
            }
        ]));
    }
}

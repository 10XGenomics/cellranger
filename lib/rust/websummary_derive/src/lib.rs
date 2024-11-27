use cr_types::websummary::{AlertConfig, MetricConfig};
use cr_types::LibraryType;
use heck::ToUpperCamelCase;
use proc_macro::TokenStream;
use quote::{format_ident, quote, ToTokens};
use serde::Deserialize;
use std::collections::BTreeMap;
use std::path::PathBuf;
use syn::{parse_macro_input, Fields, ItemStruct, LitStr};

#[allow(unused)]
#[derive(Default, Debug, Deserialize)]
struct Conditions {
    library_types: Vec<LibraryType>,
    is_cmo_multiplexed: Option<bool>,
    /// Covers both RTL and OCM multiplexing.
    is_read_multiplexed: Option<bool>,
    is_rtl: Option<bool>,
    has_gdna: Option<bool>,
}

#[allow(unused)]
#[derive(Deserialize, Debug)]
struct CardWithTableToml {
    title: String,
    tier: String,
    #[serde(default)]
    conditions: Conditions,
    entries: Vec<String>,
    /// Tables may specify a column in the
    /// table that can be used to group the metrics from the individual rows
    /// together. This is a stopgap solution until we refactor tables to not be
    /// table-structured.
    ///
    /// The value associated with the key column will be extracted and added to
    /// every metric in the row when we output flattened JSON metrics.
    ///
    /// This mechanism is also used to populate the first two columns in the
    /// metrics CSV output, so most tables should specify this key even if
    /// they are only going to have one row in the websummary.
    #[serde(default)]
    group_by_key: Option<String>,
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
        if let Some(group_by_key) = &self.group_by_key {
            assert!(
                self.entries.contains(group_by_key),
                "group_by_key {group_by_key} is not present among the entries for table '{}'",
                self.title
            );
        }

        for config in self.entry_info.values() {
            config.validate().unwrap();
        }
    }
}

struct AlertConfigSlice<'a>(&'a [AlertConfig]);

impl<'a> ToTokens for AlertConfigSlice<'a> {
    fn to_tokens(&self, tokens: &mut proc_macro2::TokenStream) {
        let mut quoted = quote![];
        for cfg in self.0 {
            quoted = quote! [
                #quoted
                #cfg,
            ];
        }
        tokens.extend(quote![
            &[#quoted]
        ]);
    }
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
        let mut row_struct_field_from_metrics = quote![];
        let mut row_struct_to_json_summary_items = quote![];
        let (grouping_key, metric_csv_grouping_header) = if let Some(group_by_key) = v.group_by_key
        {
            let grouping_entry = &v.entry_info[&group_by_key];
            let metric_csv_grouping_header = grouping_entry.header.as_str();
            let grouping_header_optional = grouping_entry.optional;
            let group_by_key = format_ident!("{group_by_key}");
            let grouping_key = quote![Some(serde_json::json!(&self.#group_by_key))];
            (
                grouping_key,
                quote![Some(
                    ::cr_websummary::GroupingHeader{
                        header: #metric_csv_grouping_header.to_string(),
                        optional: #grouping_header_optional,
                    })
                ],
            )
        } else {
            (quote![None], quote![None])
        };
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
                ::cr_websummary::multi::websummary::JsonMetricSummary::new(
                    stringify!(#name).to_string(),
                    self.#name.as_ref(),
                    String::new(),
                    String::new(),
                    #info,
                    #grouping_key,
                    ctx,
                ),
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

        let mut alert_quote = quote![
            use ::cr_websummary::MakePretty;
            let mut alert_specs = Vec::new();
        ];
        for name in &v.entries {
            let alerts = AlertConfigSlice(v.entry_info[name].alerts.as_slice());
            let name_ident = format_ident!("{}", name);

            alert_quote = quote![
                #alert_quote
                alert_specs.extend(::cr_websummary::multi::websummary::JsonMetricSummary::construct_alerts(
                    #name,
                    self.#name_ident.as_ref(),
                    #alerts,
                    ctx,
                ));
            ];
        }

        q = quote![
            #q
            #[automatically_derived]
            #[derive(Debug, Default, Clone, Serialize, Deserialize, PartialEq)]
            pub struct #table_struct_name(pub Vec<#row_struct_name>);

            impl #table_struct_name {
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
                fn to_json_summary(&self, ctx: &::cr_websummary::alert::AlertContext) -> Vec<::cr_websummary::multi::websummary::JsonMetricSummary> {
                    self.0.iter().map(|item| item.to_json_summary(ctx)).flatten().collect()
                }
            }

            #[automatically_derived]
            impl From<#table_struct_name> for ::cr_websummary::CardWithMetric {
                fn from(src: #table_struct_name) -> ::cr_websummary::CardWithMetric {
                    #metric_headers_and_help
                    let table = ::cr_websummary::GenericTable {
                        header: Some(headers),
                        rows: src.0.into_iter().map(|row| row.into()).collect(),
                        grouping_header: #metric_csv_grouping_header,
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
                fn to_json_summary(&self, ctx: &::cr_websummary::alert::AlertContext) -> Vec<::cr_websummary::multi::websummary::JsonMetricSummary> {
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
            fn to_csv_rows(&self) -> Vec<Vec<String>> {
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

    let mut append_rows =
        quote![let mut vals: Vec<::cr_websummary::multi::websummary::JsonMetricSummary> = vec![];];

    for field_name in item_struct
        .fields
        .iter()
        .map(|f| f.ident.as_ref().expect("No identifier for struct field."))
    {
        append_rows = quote![
            #append_rows
            vals.append(&mut self
                .#field_name
                .to_json_summary(ctx));
        ];
    }
    quote! [
        #[automatically_derived]
        impl ToJsonSummary for #ident {
            fn to_json_summary(&self, ctx: &AlertContext) -> Vec<::cr_websummary::multi::websummary::JsonMetricSummary> {
                #append_rows;
                vals
            }
        }
    ]
    .into()
}

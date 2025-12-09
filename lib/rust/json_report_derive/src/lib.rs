//! This crate defines a procedural macro for deriving the
//! `JsonReport` trait. `JsonReport` trait can be auto derived for a
//! struct if all the fields of the structs implement
//! the `JsonReport` trait
//!
//! # How to use
//! Add `json_report_derive` as a dependency to you `Cargo.toml`.
//! Just apply `#[derive(JsonReport)]` to a struct.
//!
//! # What you write
//!
//! ``` ignore
//! #[derive(JsonReport)]
//! struct Foo {
//!     bar1: CountMetric,
//!     bar2: PercentMetric,
//! }
//! ```
//!
//! # What you get (just for illustration)
//!
//! ``` ignore
//! struct Foo {
//!     bar1: CountMetric,
//!     bar2: PercentMetric,
//! }
//!
//! #[automatically_derived]
//! impl JsonReport for Foo {
//!     fn to_json_reporter(&self) -> JsonReporter {
//!         let mut full_reporter = JsonReporter::default();
//!         full_reporter.merge(self.bar1.to_json_reporter().add_prefix("bar1"));
//!         full_reporter.merge(self.bar2.to_json_reporter().add_prefix("bar2"));
//!         full_reporter
//!     }
//! }
//! ```
#![deny(missing_docs)]

use itertools::multizip;
use proc_macro::TokenStream;
use quote::quote;
use syn::{Attribute, Data, DeriveInput, Fields, Ident};

/// The function that gets called when we `#[derive(JsonReport)]`
///
/// # Arguments
///
/// * `input` - Rust code correspondint to the struct
///
/// # Returns
/// * a `TokenStream`, which is the Rust code corresponding to the
///   implementation of the trait
///
#[proc_macro_derive(JsonReport, attributes(json_report))]
pub fn derive_json_report_trait(input: TokenStream) -> TokenStream {
    // Parse the input tokens into a syntax tree.
    let input: DeriveInput = syn::parse(input).unwrap();

    // Used in the quasi-quotation below as `#name`.
    let name = input.ident;

    // Handle generics in the struct
    let (impl_generics, ty_generics, where_clause) = input.generics.split_for_impl();

    // Collect the struct fields
    let fieldinfo = collect_struct_fields(&input.data);

    let idents = fieldinfo
        .iter()
        .filter(|f| !f.skip_report && !f.block_report)
        .map(|f| &f.ident);

    let names = fieldinfo
        .iter()
        .filter(|f| !f.skip_report && !f.block_report)
        .map(|f| &f.name);

    let inline = fieldinfo
        .iter()
        .filter(|f| !f.skip_report && !f.block_report)
        .map(|f| &f.inline_report);

    let block_idents = fieldinfo
        .iter()
        .filter(|f| f.block_report)
        .map(|f| &f.ident);

    let block_names = fieldinfo.iter().filter(|f| f.block_report).map(|f| &f.name);

    let extend = extend_tokens(&input.attrs);

    let expanded = quote! {
        // The generated impl
        #[automatically_derived]
        impl #impl_generics ::metric::JsonReport for #name #ty_generics #where_clause {
            fn to_json_reporter(&self) -> ::metric::JsonReporter {
                let mut full_reporter = ::metric::JsonReporter::default();
                #(
                    let reporter = self.#idents.to_json_reporter();
                    ::metric::Metric::merge(&mut full_reporter,
                        if #inline {
                            reporter
                        } else {
                            reporter.add_prefix(#names)
                        }
                    );
                )*

                // For fields using #[json_report(block)]
                #(
                    let value = ::serde_json::to_value(&self.#block_idents).unwrap();
                    if !value.is_null() {
                        full_reporter.insert(#block_names, value);
                    }
                )*

                // Call the custom function
                #extend

                full_reporter
            }
        }
    };

    // Hand the output tokens back to the compiler.
    TokenStream::from(expanded)
}

struct FieldInfo {
    ident: Ident,
    name: String,
    skip_report: bool,
    inline_report: bool,
    block_report: bool,
}

// Collect the names of struct fields from the input data
// Panics if called on a tuple struct or unit struct or enum
fn collect_struct_fields(data: &Data) -> Vec<FieldInfo> {
    let idents: Vec<Ident>;
    let names: Vec<String>;
    let skip_reports: Vec<bool>;
    let inline_reports: Vec<bool>;
    let block_reports: Vec<bool>;
    match *data {
        Data::Struct(ref data) => match data.fields {
            Fields::Named(ref fields) => {
                idents = fields
                    .named
                    .iter()
                    .filter_map(|f| f.ident.clone())
                    .collect();
                names = fields
                    .named
                    .iter()
                    .filter_map(|x| x.ident.as_ref())
                    .map(ToString::to_string)
                    .collect();

                let field_attrs: Vec<(bool, bool, bool)> = fields
                    .named
                    .iter()
                    .map(|f| parse_field_attributes(&f.attrs))
                    .collect();

                skip_reports = field_attrs.iter().map(|x| x.0).collect();
                inline_reports = field_attrs.iter().map(|x| x.1).collect();
                block_reports = field_attrs.iter().map(|x| x.2).collect();
            }
            Fields::Unnamed(_) | Fields::Unit => {
                panic!("You can only derive JsonReport for normal structs!");
            }
        },
        Data::Enum(_) | Data::Union(_) => {
            panic!("You can only derive JsonReport for structs as of now!");
        }
    }
    let mut result = Vec::new();
    for (ident, name, skip_report, inline_report, block_report) in
        multizip((idents, names, skip_reports, inline_reports, block_reports))
    {
        result.push(FieldInfo {
            ident,
            name,
            skip_report,
            inline_report,
            block_report,
        });
    }
    result
}

fn parse_field_attributes(attrs: &[Attribute]) -> (bool, bool, bool) {
    let mut ident = None;
    for attr in attrs {
        if !attr.path().is_ident("json_report") {
            continue;
        };
        let syn::Meta::List(list) = &attr.meta else {
            panic!("Expected #[json_report(...)]");
        };
        list.parse_nested_meta(|meta| {
            assert!(ident.is_none(), "Expected at most one #[json_report(...)]");
            ident = Some(meta.path.get_ident().unwrap().to_string());
            Ok(())
        })
        .unwrap();
    }

    match ident.as_deref() {
        None => (false, false, false),
        Some("skip") => (true, false, false),
        Some("inline") => (false, true, false),
        Some("block") => (false, false, true),
        Some(s) => panic!("One of block, inline, or skip expected in #[json_report({s})]"),
    }
}

// Parse the extend function name argument from the container attributes.
// #[json_report(extend="function_name")]
fn parse_extend_argument(outer_attrs: &[Attribute]) -> Option<String> {
    let mut extend = None;
    for attr in outer_attrs {
        if !attr.path().is_ident("json_report") {
            continue;
        }
        let syn::Meta::List(list) = &attr.meta else {
            panic!("Expected #[json_report(extend=\"...\")]");
        };
        list.parse_nested_meta(|meta| {
            assert!(
                meta.path.is_ident("extend"),
                "Expected #[json_report(extend=\"...\")]",
            );
            assert!(
                extend.is_none(),
                "Expected at most one #[json_report(extend=\"...\")]"
            );
            let litstr: syn::LitStr = meta.value().unwrap().parse().unwrap();
            extend = Some(litstr.value());
            Ok(())
        })
        .unwrap();
    }
    extend
}

fn extend_tokens(attrs: &[Attribute]) -> proc_macro2::TokenStream {
    let Some(fname) = parse_extend_argument(attrs) else {
        return quote! {};
    };
    let fname_ident = Ident::new(&fname, proc_macro2::Span::call_site());
    quote! {
        let extend_reporter = self.#fname_ident();
        full_reporter.merge(extend_reporter);
    }
}

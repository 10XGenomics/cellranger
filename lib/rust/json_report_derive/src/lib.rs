#![deny(
    missing_copy_implementations,
    non_upper_case_globals,
    trivial_casts,
    trivial_numeric_casts,
    unsafe_code,
    unstable_features,
    unused_import_braces,
    unused_qualifications
)]

//!
//! This is documentation for the `json-metric_derive` crate.
//!
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
//! #[macro_use]
//! extern crate json_report_derive;
//!
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
//! #[macro_use]
//! extern crate json_report_derive;
//!
//! struct Foo {
//!     bar1: CountMetric,
//!     bar2: PercentMetric,
//! }
//! #[automatically_derived]
//! impl JsonReport for Foo {
//!     fn to_json_reporter(&self) -> JsonReporter {
//!         let mut full_reporter = JsonReporter::new();
//!         let mut reporter = bar1.to_json_reporter();
//!         reporter.add_prefix("bar1");
//!         full_reporter.merge(reporter);
//!         let mut reporter = bar2.to_json_reporter();
//!         reporter.add_prefix("bar2");
//!         full_reporter.merge(reporter);
//!         full_reporter
//!     }
//! }
//! ```

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
/// implementation of the trait
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
        .filter(|f| !f.skip_report)
        .map(|f| &f.ident);
    let names = fieldinfo.iter().filter(|f| !f.skip_report).map(|f| &f.name);
    let names_copy = names.clone();
    let inline = fieldinfo
        .iter()
        .filter(|f| !f.skip_report)
        .map(|f| &f.inline_report);

    let block = fieldinfo
        .iter()
        .filter(|f| !f.skip_report)
        .map(|f| &f.block_report);

    let extend = extend_tokens(&input.attrs);

    let expanded = quote! {
        // The generated impl
        #[automatically_derived]
        #[allow(unused_attributes)]
        impl #impl_generics ::metric::JsonReport for #name #ty_generics #where_clause {
            fn to_json_reporter(&self) -> ::metric::JsonReporter {
                let mut full_reporter = ::metric::JsonReporter::new();
                #(
                    let mut reporter = self.#idents.to_json_reporter();
                    if (#block) {
                        let block_val = reporter.block_value();
                        full_reporter.insert(#names_copy, block_val);
                    } else {
                        if !(#inline) {
                            reporter.add_prefix(#names);
                        }
                        ::metric::Metric::merge(&mut full_reporter, reporter)
                    }

                )*
                // Call the custom function
                #extend

                full_reporter
            }
        }
    };

    // Hand the output tokens back to the compiler.
    proc_macro::TokenStream::from(expanded)
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
                    .map(std::string::ToString::to_string)
                    .collect();

                let inner_attrs: Vec<(bool, bool, bool)> = fields
                    .named
                    .iter()
                    .map(|f| check_inner_attrs(&f.attrs))
                    .collect();

                skip_reports = inner_attrs.iter().map(|x| x.0).collect();
                inline_reports = inner_attrs.iter().map(|x| x.1).collect();
                block_reports = inner_attrs.iter().map(|x| x.2).collect();
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

fn check_inner_attrs(inner_attrs: &[Attribute]) -> (bool, bool, bool) {
    let mut skip = false;
    let mut inline = false;
    let mut found = false;
    let mut block = false;

    for attr in inner_attrs {
        // Cast the attribute as syn::Meta
        // match attr.interpret_meta().unwrap() {
        // Things like #[json_report(extend="function_name")] is parsed as a Meta::List
        if let syn::Meta::List(list) = attr.parse_meta().unwrap() {
            if !list.path.is_ident("json_report") {
                continue;
            }
            // Iterate through stuff in the paranthesis
            for nestedmeta in list.nested {
                match nestedmeta {
                    syn::NestedMeta::Meta(syn::Meta::Path(ref path)) => {
                        // skip is parsed as a Path
                        assert!(
                            path.is_ident("skip") || path.is_ident("inline") || path.is_ident("block"),
                            "Expecting only #[json_report(skip)] or #[json_report(inline)] or #[json_report(block)] as inner attribute"
                        );
                        if found {
                            panic!("Found multiple #[json_report()] for one field")
                        }
                        if path.is_ident("skip") {
                            skip = true;
                        }
                        if path.is_ident("inline") {
                            inline = true;
                        }
                        if path.is_ident("block") {
                            block = true;
                        }
                        found = true;
                    }
                    _ => panic!("Error parsing the inner #[json_report()] attribute"),
                }
            }
        }
    }
    assert!(
        !(inline && block),
        "You cannot #[json_report(inline)] and #[json_report(block)] on the same field"
    );
    (skip, inline, block)
}

// Parse function_name from the outer attribute
// #[json_report(extend="function_name")]
fn function_name_from_outer_attrs(outer_attrs: &[Attribute]) -> Option<String> {
    let mut result = None;
    let mut found = false;
    // Iterate through attributes
    for attr in outer_attrs {
        // Cast the attribute as syn::Meta
        // Things like #[json_report(extend="function_name")] is parsed as a Meta::List
        if let syn::Meta::List(list) = attr.parse_meta().unwrap() {
            if !list.path.is_ident("json_report") {
                continue;
            }
            // Iterate through stuff in the paranthesis
            for nestedmeta in list.nested {
                match nestedmeta {
                    syn::NestedMeta::Meta(syn::Meta::NameValue(ref name_value)) => {
                        // extend="function_name" is parsed as a Name value pair
                        assert!(
                            name_value.path.is_ident("extend"),
                            "Expecting only #[json_report(extend=\"name\")] as outer attribute"
                        );
                        match name_value.lit {
                            // Phew! Finally we reached the function name
                            syn::Lit::Str(ref s) => {
                                if found {
                                    panic!("Multiple #[json_report(extend = \"*\")] found");
                                }
                                result = Some(s.value());
                                found = true;
                            }
                            _ => panic!("Error parsing the outer #[json_report] attribute"),
                        };
                    }
                    _ => panic!("Error parsing the outer #[json_report] attribute"),
                }
            }
        }
    }
    result
}

fn extend_tokens(outer_attrs: &[Attribute]) -> proc_macro2::TokenStream {
    if let Some(fname) = function_name_from_outer_attrs(outer_attrs) {
        let fname_ident = Ident::new(&fname, proc_macro2::Span::call_site());
        quote! {
            let extend_reporter = self.#fname_ident();
            full_reporter.merge(extend_reporter);
        }
    } else {
        quote! {}
    }
}

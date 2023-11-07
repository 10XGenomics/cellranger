#![deny(
    missing_copy_implementations,
    non_upper_case_globals,
    trivial_casts,
    trivial_numeric_casts,
    unsafe_code,
    unstable_features,
    unused_extern_crates,
    unused_import_braces,
    unused_qualifications
)]
#![recursion_limit = "128"]

//!
//! This is documentation for the `metric_derive` crate.
//!
//! This crate defines a procedural macro for deriving the
//! Metric Trait. Metric trait can be auto derived for a
//! struct if all the fields of the structs implement
//! the Metric trait
//!
//! # How to use
//! Add `metric_derive` as a dependency to you `Cargo.toml`.
//! Just apply `#[derive(Metric)]` to a struct.
//!
//! # What you write
//!
//! ``` ignore
//! #[macro_use]
//! extern crate metric_derive;
//!
//! #[derive(Metric)]
//! struct Foo {
//!     bar1: CountMetric,
//!     bar2: PercentMetric,
//! }
//! ```
//!
//! # What you get
//!
//! ``` ignore
//! #[macro_use]
//! extern crate metric_derive;
//!
//! struct Foo {
//!     bar1: CountMetric,
//!     bar2: PercentMetric,
//! }
//! #[automatically_derived]
//! impl Metric for Foo {
//!     fn new() -> Self {
//!         Foo {
//!             bar1: Metric::new(),
//!             bar2: Metric::new()
//!         }
//!     }
//!     fn merge(&mut self, other: &Self) {
//!         self.bar1.merge(&other.bar1);
//!         self.bar2.merge(&other.bar2);
//!     }
//! }
//! ```

#[macro_use]
extern crate quote;

use proc_macro::TokenStream;
use syn::spanned::Spanned;
use syn::{Data, DeriveInput, Fields, Index};

/// The function that gets called when we #[derive(Metric)]
///
/// # Arguments
///
/// * `input` - Rust code correspondint to the struct
///
/// # Returns
/// * a `TokenStream`, which is the Rust code corresponding to the
/// implementation of the metric trait
///
#[proc_macro_derive(Metric)]
pub fn derive_metric_trait(input: TokenStream) -> TokenStream {
    // Parse the input tokens into a syntax tree.
    let input: DeriveInput = syn::parse(input).unwrap();

    // Used in the quasi-quotation below as `#name`.
    let name = input.ident;

    // Handle generics in the struct
    let (impl_generics, ty_generics, where_clause) = input.generics.split_for_impl();

    let expanded = match input.data {
        Data::Struct(ref data) => match data.fields {
            Fields::Named(ref fields) => {
                let fnames = fields
                    .named
                    .iter()
                    .map(|f| f.ident.as_ref().unwrap())
                    .collect::<Vec<_>>();
                let fnames_copy1 = fnames.clone();
                let fnames_copy2 = fnames.clone();

                quote! {
                    // The generated impl
                    #[automatically_derived]
                    #[allow(unused_attributes)]
                    impl #impl_generics Metric for #name #ty_generics #where_clause {
                        fn new() -> Self {
                            #name {
                                #(
                                    #fnames: Metric::new()
                                ),*
                            }
                        }
                        fn merge(&mut self, other: Self) {
                            #(
                                self.#fnames_copy1.merge(other.#fnames_copy2);
                            )*
                        }
                    }

                    #[automatically_derived]
                    #[allow(unused_attributes)]
                    impl #impl_generics ::std::ops::Add for #name #ty_generics #where_clause {
                        type Output = #name #ty_generics;
                        fn add(mut self, other: Self::Output) -> Self::Output {
                            self.merge(other);
                            self
                        }
                    }

                    #[automatically_derived]
                    #[allow(unused_attributes)]
                    impl #impl_generics ::std::ops::AddAssign for #name #ty_generics #where_clause {
                        fn add_assign(&mut self, other: #name #ty_generics) {
                            self.merge(other);
                        }
                    }

                    #[automatically_derived]
                    #[allow(unused_attributes)]
                    impl #impl_generics ::std::iter::Sum for #name #ty_generics #where_clause {
                        fn sum<IteratorType: Iterator<Item=#name #ty_generics>>(iter: IteratorType) -> #name #ty_generics {
                            iter.fold(Metric::new(), |a, b| a + b)
                        }
                    }
                }
            }
            Fields::Unnamed(ref fields) => {
                let indices = fields
                    .unnamed
                    .iter()
                    .enumerate()
                    .map(|(i, f)| Index {
                        index: i as u32,
                        span: f.span(),
                    })
                    .collect::<Vec<_>>();
                let indices_copy = indices.clone();
                let new_tokens = (0..indices.len())
                    .map(|_| quote! { Metric::new() })
                    .collect::<Vec<_>>();
                quote! {
                    // The generated impl
                    #[automatically_derived]
                    #[allow(unused_attributes)]
                    impl #impl_generics Metric for #name #ty_generics #where_clause {
                        fn new() -> Self {
                            #name(#(#new_tokens),*)
                        }
                        fn merge(&mut self, other: Self) {
                            #(
                                self.#indices.merge(other.#indices_copy);
                            )*
                        }
                    }

                    #[automatically_derived]
                    #[allow(unused_attributes)]
                    impl #impl_generics ::std::ops::Add for #name #ty_generics #where_clause {
                        type Output = #name #ty_generics;
                        fn add(mut self, other: Self::Output) -> Self::Output {
                            self.merge(other);
                            self
                        }
                    }

                    #[automatically_derived]
                    #[allow(unused_attributes)]
                    impl #impl_generics ::std::ops::AddAssign for #name #ty_generics #where_clause {
                        fn add_assign(&mut self, other: #name #ty_generics) {
                            self.merge(other);
                        }
                    }

                    #[automatically_derived]
                    #[allow(unused_attributes)]
                    impl #impl_generics ::std::iter::Sum for #name #ty_generics #where_clause {
                        fn sum<IteratorType: Iterator<Item=#name #ty_generics>>(iter: IteratorType) -> #name #ty_generics {
                            iter.fold(Metric::new(), |a, b| a + b)
                        }
                    }
                }
            }
            Fields::Unit => {
                panic!("You can only derive Metric for normal structs!");
            }
        },
        Data::Enum(_) | Data::Union(_) => {
            panic!("You can only derive Metric for structs as of now!");
        }
    };

    // Hand the output tokens back to the compiler.
    expanded.into()
}

#![deny(missing_docs)]
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
//! struct Foo {
//!     bar1: CountMetric,
//!     bar2: PercentMetric,
//! }
//! #[automatically_derived]
//! impl Metric for Foo {
//!     fn merge(&mut self, other: &Self) {
//!         self.bar1.merge(&other.bar1);
//!         self.bar2.merge(&other.bar2);
//!     }
//! }
//! ```

use proc_macro::TokenStream;
use quote::quote;
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
///   implementation of the metric trait
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
                    impl #impl_generics Metric for #name #ty_generics #where_clause {
                        fn merge(&mut self, other: Self) {
                            #(
                                self.#fnames_copy1.merge(other.#fnames_copy2);
                            )*
                        }
                    }

                    #[automatically_derived]
                    impl #impl_generics ::std::ops::Add for #name #ty_generics #where_clause {
                        type Output = #name #ty_generics;
                        fn add(mut self, other: Self::Output) -> Self::Output {
                            self.merge(other);
                            self
                        }
                    }

                    #[automatically_derived]
                    impl #impl_generics ::std::ops::AddAssign for #name #ty_generics #where_clause {
                        fn add_assign(&mut self, other: #name #ty_generics) {
                            self.merge(other);
                        }
                    }

                    #[automatically_derived]
                    impl #impl_generics ::std::iter::Sum for #name #ty_generics #where_clause {
                        fn sum<IteratorType: Iterator<Item=#name #ty_generics>>(iter: IteratorType) -> #name #ty_generics {
                            iter.fold(Default::default(), |a, b| a + b)
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
                quote! {
                    // The generated impl
                    #[automatically_derived]
                    impl #impl_generics Metric for #name #ty_generics #where_clause {
                        fn merge(&mut self, other: Self) {
                            #(
                                self.#indices.merge(other.#indices_copy);
                            )*
                        }
                    }

                    #[automatically_derived]
                    impl #impl_generics ::std::ops::Add for #name #ty_generics #where_clause {
                        type Output = #name #ty_generics;
                        fn add(mut self, other: Self::Output) -> Self::Output {
                            self.merge(other);
                            self
                        }
                    }

                    #[automatically_derived]
                    impl #impl_generics ::std::ops::AddAssign for #name #ty_generics #where_clause {
                        fn add_assign(&mut self, other: #name #ty_generics) {
                            self.merge(other);
                        }
                    }

                    #[automatically_derived]
                    impl #impl_generics ::std::iter::Sum for #name #ty_generics #where_clause {
                        fn sum<IteratorType: Iterator<Item=#name #ty_generics>>(iter: IteratorType) -> #name #ty_generics {
                            iter.fold(Default::default(), |a, b| a + b)
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

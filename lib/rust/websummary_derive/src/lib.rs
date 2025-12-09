//! websummary_derive
#![deny(missing_docs)]

use proc_macro::TokenStream;
use quote::quote;
use syn::{Fields, ItemStruct};

const ALERT_NOT_ON_NAMED_STRUCT_ERROR: &str =
    r"#[derive(Alert)] can only be used on structs with named fields.";

/// Implement #[derive(Alert)] macro
#[proc_macro_derive(Alert)]
pub fn websummary_alert(item: TokenStream) -> TokenStream {
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

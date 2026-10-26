//! Derive macros for geometric traits. Re-exported by `monstertruck-traits`
//! (feature `"derive"`).
//!
//! # Examples
//!
//! The example is `ignore`d rather than run: it needs `monstertruck-traits`,
//! which DEPENDS on this crate, so this crate can never dev-depend on it to
//! link the doctest. These derives are exercised for real in
//! `monstertruck-traits/tests/derives.rs` (behind that crate's `derive` and
//! `polynomial` features), where the dependency direction works.
//!
//! ```rust,ignore
//! use monstertruck_traits::prelude::*;
//!
//! /// An enum of curve types -- derive macros delegate trait methods
//! /// to the inner type via match arms.
//! #[derive(Clone, ParametricCurve, BoundedCurve)]
//! pub enum MyCurve {
//!     Line(Line<Point3>),
//!     Nurbs(NurbsCurve<Vector4>),
//! }
//!
//! let curve: MyCurve = MyCurve::Line(/* ... */);
//! let pt = curve.evaluate(0.5); // dispatches to Line::evaluate
//! ```
//!
//! Users do not need to depend on this crate directly:
//!
//! ```toml
//! monstertruck-traits = { version = "0.3", features = ["derive"] }
//! ```
#![cfg_attr(not(debug_assertions), deny(warnings))]
#![deny(clippy::all, rust_2018_idioms)]
#![warn(
    missing_docs,
    missing_debug_implementations,
    trivial_casts,
    trivial_numeric_casts,
    unsafe_code,
    unstable_features,
    unused_import_braces,
    unused_qualifications
)]

use proc_macro::TokenStream;
use proc_macro_crate::{FoundCrate, crate_name};
use proc_macro_error::proc_macro_error;
use proc_macro2::{Span, TokenStream as TokenStream2};
use quote::*;
use syn::*;

macro_rules! derive_try {
    ($result:expr) => {
        match $result {
            Ok(value) => value,
            Err(error) => return error.into_compile_error().into(),
        }
    };
}

fn single_variant_field(variant: &Variant) -> Result<&syn::Field> {
    let mut fields = variant.fields.iter();
    let field = fields.next().ok_or_else(|| {
        Error::new_spanned(
            variant,
            "derive variants must contain exactly one delegated field",
        )
    })?;
    if fields.next().is_some() {
        Err(Error::new_spanned(
            variant,
            "derive variants must contain exactly one delegated field",
        ))
    } else {
        Ok(field)
    }
}

fn enum_field_type<'a>(
    variants: &'a punctuated::Punctuated<Variant, token::Comma>,
    ty: &Ident,
) -> Result<&'a Type> {
    let mut variants = variants.iter();
    let first_variant = variants
        .next()
        .ok_or_else(|| Error::new_spanned(ty, "derive requires at least one enum variant"))?;
    let first = single_variant_field(first_variant)?;
    for variant in variants {
        single_variant_field(variant)?;
    }
    Ok(&first.ty)
}

fn tuple_struct_field(fields: &Fields) -> Result<&Type> {
    let Fields::Unnamed(fields) = fields else {
        return Err(Error::new_spanned(
            fields,
            "derive supports only single-field tuple structs",
        ));
    };
    let mut fields_iter = fields.unnamed.iter();
    let field = fields_iter.next().ok_or_else(|| {
        Error::new_spanned(fields, "derive supports only single-field tuple structs")
    })?;
    if fields_iter.next().is_some() {
        Err(Error::new_spanned(
            fields,
            "derive supports only single-field tuple structs",
        ))
    } else {
        Ok(&field.ty)
    }
}

fn unsupported_input(ty: &Ident) -> TokenStream2 {
    Error::new_spanned(
        ty,
        "derive supports only enums and single-field tuple structs",
    )
    .into_compile_error()
}

fn dependency_path(package: &str) -> Result<TokenStream2> {
    match crate_name(package) {
        Ok(FoundCrate::Itself) => Ok(quote! { crate }),
        Ok(FoundCrate::Name(name)) => {
            let ident = Ident::new(&name, Span::call_site());
            Ok(quote! { ::#ident })
        }
        Err(error) => Err(Error::new(
            Span::call_site(),
            format!("cannot resolve `{package}` for generated derive code: {error}"),
        )),
    }
}

fn dependency_scope(
    package: &str,
    alias: &str,
    implementation: TokenStream2,
) -> Result<TokenStream2> {
    let path = dependency_path(package)?;
    let alias = Ident::new(alias, Span::call_site());
    Ok(quote! {
        const _: () = {
            use #path as #alias;
            #implementation
        };
    })
}

fn traits_scope(implementation: TokenStream2) -> Result<TokenStream2> {
    dependency_scope("monstertruck-traits", "monstertruck_traits", implementation)
}

fn io_scope(implementation: TokenStream2) -> Result<TokenStream2> {
    dependency_scope("monstertruck-io", "monstertruck_io", implementation)
}

fn enumerate_impl_return_something<'a>(
    variants: impl IntoIterator<Item = &'a Variant> + 'a,
    method: TokenStream2,
    method_variants: TokenStream2,
) -> TokenStream2 {
    let impls = variants
        .into_iter()
        .map(|variant| {
            let field = single_variant_field(variant)?;
            let variant_name = &variant.ident;
            Ok(match &field.ident {
                Some(ident) => {
                    quote! { Self::#variant_name { #ident } => #method(#ident, #method_variants) }
                }
                None => quote! { Self::#variant_name(got) => #method(got, #method_variants) },
            })
        })
        .collect::<Result<Vec<_>>>();
    match impls {
        Ok(impls) => quote! { match self { #(#impls),* } },
        Err(error) => error.into_compile_error(),
    }
}

fn enumerate_impl_return_self<'a>(
    variants: impl IntoIterator<Item = &'a Variant> + 'a,
    method: TokenStream2,
    method_variants: TokenStream2,
) -> TokenStream2 {
    let impls = variants
        .into_iter()
        .map(|variant| {
            let field = single_variant_field(variant)?;
            let variant_name = &variant.ident;
            Ok(match &field.ident {
                Some(ident) => quote! {
                    Self::#variant_name { #ident } => Self::#variant_name { #ident: #method(#ident, #method_variants) }
                },
                None => quote! { Self::#variant_name(got) => Self::#variant_name(#method(got, #method_variants)) },
            })
        })
        .collect::<Result<Vec<_>>>();
    match impls {
        Ok(impls) => quote! { match self { #(#impls),* } },
        Err(error) => error.into_compile_error(),
    }
}

#[derive(Clone, Debug)]
struct Field {
    var: TokenStream2,
    ty: TokenStream2,
}

macro_rules! fields {
    ($($var: tt : $ty: tt),*) => {
        vec![$(Field {
            var: quote! { $var },
            ty: quote! { $ty },
        }),*]
    };
}

#[derive(Clone, Debug)]
struct Method<I> {
    name: TokenStream2,
    generics: Option<TokenStream2>,
    self_field: TokenStream2,
    fields: Vec<Field>,
    return_type: TokenStream2,
    variants: I,
    trait_name: TokenStream2,
}

macro_rules! methods {
    (
        $variants: ident, $trait_name: ident,
        $(fn $name: ident (
            $self_field: expr,
            $($var: ident: $ty: ty),*$(,)?
        ) -> $return_type: ty),*$(,)?
    ) => {
        vec![$(Method {
            name: quote! { $name },
            generics: None,
            self_field: quote! { $self_field, },
            fields: fields!($($var: $ty),*),
            return_type: quote! { $return_type },
            variants: $variants,
            trait_name: $trait_name.clone(),
        }
        .to_token_stream()),*]
    };
    (
        $variants: ident, $trait_name: ident,
        $(fn $name: ident <$($(($const: tt))? $gen: ident: $path: path),*> (
            $self_field: expr,
            $($var: ident: $ty: ty),*$(,)?
        ) -> $return_type: ty),*$(,)?
    ) => {
        vec![$(Method {
            name: quote! { $name },
            generics: Some(quote! { <$($($const)? $gen: $path),*> }),
            self_field: quote! { $self_field, },
            fields: fields!($($var: $ty),*),
            return_type: quote! { $return_type },
            variants: $variants,
            trait_name: $trait_name.clone(),
        }
        .to_token_stream()),*]
    };
}

impl<'a, I> Method<I>
where I: IntoIterator<Item = &'a Variant> + 'a + Copy
{
    fn to_token_stream(&'a self) -> TokenStream2 {
        let method_name = &self.name;
        let generics = &self.generics;
        let trait_name = &self.trait_name;
        let self_field = &self.self_field;
        let fields = self
            .fields
            .iter()
            .map(|f| {
                let var = &f.var;
                let ty = &f.ty;
                quote! { #var: #ty }
            })
            .collect::<Vec<_>>();
        let vals = self.fields.iter().map(|f| &f.var).collect::<Vec<_>>();
        let return_type = &self.return_type;
        let implement = if return_type.to_string() == "Self" {
            enumerate_impl_return_self::<'_>(
                self.variants,
                quote! { #trait_name::#method_name },
                quote! { #(#vals),* },
            )
        } else {
            enumerate_impl_return_something::<'a>(
                self.variants,
                quote! { #trait_name::#method_name },
                quote! { #(#vals),* },
            )
        };
        quote! {
            fn #method_name #generics (#self_field #(#fields),*) -> #return_type { #implement }
        }
    }
}

/// Derive macro generating an impl of the trait `BoundedCurve` for Enums or single field tuple structs.
#[proc_macro_error]
#[proc_macro_derive(BoundedCurve)]
pub fn derive_bounded_curve(input: TokenStream) -> TokenStream {
    let input = parse_macro_input!(input as DeriveInput);
    let trait_name = quote! { monstertruck_traits::BoundedCurve };
    let ty = input.ident;
    let generics_params = input.generics;
    let where_predicates = generics_params
        .where_clause
        .iter()
        .flat_map(|x| &x.predicates);
    let implementation = match input.data {
        Data::Enum(DataEnum { ref variants, .. }) => {
            let top_ty = derive_try!(enum_field_type(variants, &ty));
            let methods = methods! {
                variants, trait_name,
                fn range_tuple(&self,) -> (f64, f64),
                fn front(&self,) -> Self::Point,
                fn back(&self,) -> Self::Point,
            };
            quote! {
                #[automatically_derived]
                impl #generics_params #trait_name for #ty #generics_params
                where
                    #(#where_predicates,)*
                    #top_ty: #trait_name,

                    Self: monstertruck_traits::ParametricCurve<Point = <#top_ty as monstertruck_traits::ParametricCurve>::Point>, {
                    #(#methods)*
                }
            }
        }
        Data::Struct(DataStruct { ref fields, .. }) => {
            let field_type = derive_try!(tuple_struct_field(fields));
            quote! {
                #[automatically_derived]
                impl #generics_params #trait_name for #ty #generics_params
                where
                    #(#where_predicates,)*
                    #field_type: #trait_name,
                    Self: monstertruck_traits::ParametricCurve<Point
                        = <#field_type as monstertruck_traits::ParametricCurve>::Point>, {
                    fn range_tuple(&self) -> (f64, f64) { self.0.range_tuple() }
                    fn front(&self) -> Self::Point { self.0.front() }
                    fn back(&self) -> Self::Point { self.0.back() }
                }
            }
        }
        _ => unsupported_input(&ty),
    };
    derive_try!(traits_scope(implementation)).into()
}

/// Derive macro generating an impl of the trait `BoundedSurface` for Enums or single field tuple structs.
#[proc_macro_error]
#[proc_macro_derive(BoundedSurface)]
pub fn derive_bounded_surface(input: TokenStream) -> TokenStream {
    let input = parse_macro_input!(input as DeriveInput);
    let trait_name = quote! { monstertruck_traits::BoundedSurface };
    let ty = input.ident;
    let generics_params = input.generics;
    let where_predicates = generics_params
        .where_clause
        .iter()
        .flat_map(|x| &x.predicates);
    let implementation = match input.data {
        Data::Enum(DataEnum { ref variants, .. }) => {
            let top_ty = derive_try!(enum_field_type(variants, &ty));
            let methods = methods! {
                variants, trait_name,
                fn range_tuple(&self,) -> ((f64, f64), (f64, f64)),
            };
            quote! {
                #[automatically_derived]
                impl #generics_params #trait_name for #ty #generics_params
                where
                    #(#where_predicates,)*
                    #top_ty: #trait_name,
                    {
                    #(#methods)*
                }
            }
        }
        Data::Struct(DataStruct { ref fields, .. }) => {
            let field_type = derive_try!(tuple_struct_field(fields));
            quote! {
                #[automatically_derived]
                impl #generics_params #trait_name for #ty #generics_params
                where
                    #(#where_predicates,)*
                    #field_type: #trait_name {
                    fn range_tuple(&self) -> ((f64, f64), (f64, f64)) {
                        self.0.range_tuple()
                    }
                }
            }
        }
        _ => unsupported_input(&ty),
    };
    derive_try!(traits_scope(implementation)).into()
}

/// Derive macro generating an impl of the trait `Cut` for Enums or single field tuple structs.
#[proc_macro_error]
#[proc_macro_derive(Cut)]
pub fn derive_cut(input: TokenStream) -> TokenStream {
    let input = parse_macro_input!(input as DeriveInput);
    let trait_name = quote! { monstertruck_traits::Cut };
    let ty = input.ident;
    let generics_params = input.generics;
    let where_predicates = generics_params
        .where_clause
        .iter()
        .flat_map(|x| &x.predicates);
    let implementation = match input.data {
        Data::Enum(DataEnum { ref variants, .. }) => {
            let top_ty = derive_try!(enum_field_type(variants, &ty));
            let methods = methods! {
                variants, trait_name,
                fn cut(&mut self, t: f64) -> Self,
            };
            quote! {
                #[automatically_derived]
                impl #generics_params #trait_name for #ty #generics_params
                where
                    #(#where_predicates,)*
                    #top_ty: #trait_name,
                    {
                    #(#methods)*
                }
            }
        }
        Data::Struct(DataStruct { ref fields, .. }) => {
            let field_type = derive_try!(tuple_struct_field(fields));
            quote! {
                #[automatically_derived]
                impl #generics_params #trait_name for #ty #generics_params
                where
                    #(#where_predicates,)*
                    #field_type: #trait_name {
                    fn cut(&mut self, t: f64) -> Self { Self(self.0.cut(t)) }
                }
            }
        }
        _ => unsupported_input(&ty),
    };
    derive_try!(traits_scope(implementation)).into()
}

/// Derive macro generating an impl of the trait `Invertible` for Enums or single field tuple structs.
#[proc_macro_error]
#[proc_macro_derive(Invertible)]
pub fn derive_invertible(input: TokenStream) -> TokenStream {
    let input = parse_macro_input!(input as DeriveInput);
    let trait_name = quote! { monstertruck_traits::Invertible };
    let ty = input.ident;
    let generics_params = input.generics;
    let where_predicates = generics_params
        .where_clause
        .iter()
        .flat_map(|x| &x.predicates);
    let implementation = match input.data {
        Data::Enum(DataEnum { ref variants, .. }) => {
            let top_ty = derive_try!(enum_field_type(variants, &ty));
            let methods = methods! {
                variants, trait_name,
                fn invert(&mut self,) -> (),
                fn inverse(&self,) -> Self,
            };
            quote! {
                #[automatically_derived]
                impl #generics_params #trait_name for #ty #generics_params
                where
                    #(#where_predicates,)*
                    #top_ty: #trait_name,
                    #top_ty: #trait_name,

                    Self: Clone {
                    #(#methods)*
                }
            }
        }
        Data::Struct(DataStruct { ref fields, .. }) => {
            let field_type = derive_try!(tuple_struct_field(fields));
            quote! {
                #[automatically_derived]
                impl #generics_params #trait_name for #ty #generics_params
                where
                    #(#where_predicates,)*
                    #field_type: #trait_name,
                    Self: Clone {
                    fn invert(&mut self) { self.0.invert() }
                    fn inverse(&self) -> Self { Self(self.0.inverse()) }
                }
            }
        }
        _ => unsupported_input(&ty),
    };
    derive_try!(traits_scope(implementation)).into()
}

/// Derive macro generating an impl of the trait `ParameterDivision1D` for Enums or single field tuple structs.
#[proc_macro_error]
#[proc_macro_derive(ParameterDivision1D)]
pub fn derive_parameter_division_1d(input: TokenStream) -> TokenStream {
    let input = parse_macro_input!(input as DeriveInput);
    let trait_name = quote! { monstertruck_traits::ParameterDivision1D };
    let ty = input.ident;
    let generics_params = input.generics;
    let where_predicates = generics_params
        .where_clause
        .iter()
        .flat_map(|x| &x.predicates);
    let implementation = match input.data {
        Data::Enum(DataEnum { ref variants, .. }) => {
            let top_ty = derive_try!(enum_field_type(variants, &ty));
            let methods = methods! {
                variants, trait_name,
                fn try_parameter_division(&self, range: (f64, f64), tol: f64) -> Option<(Vec<f64>, Vec<Self::Point>)>,
                fn parameter_division(&self, range: (f64, f64), tol: f64) -> (Vec<f64>, Vec<Self::Point>),
            };
            quote! {
                #[automatically_derived]
                impl #generics_params #trait_name for #ty #generics_params
                where
                    #(#where_predicates,)*
                    #top_ty: #trait_name,
                    {
                    type Point = <#top_ty as #trait_name>::Point;
                    #(#methods)*
                }
            }
        }
        Data::Struct(DataStruct { ref fields, .. }) => {
            let field_type = derive_try!(tuple_struct_field(fields));
            quote! {
                #[automatically_derived]
                impl #generics_params #trait_name for #ty #generics_params
                where
                    #(#where_predicates,)*
                    #field_type: #trait_name {
                    type Point = <#field_type as #trait_name>::Point;
                    fn try_parameter_division(&self, range: (f64, f64), tol: f64) -> Option<(Vec<f64>, Vec<Self::Point>)> {
                        self.0.try_parameter_division(range, tol)
                    }
                    fn parameter_division(&self, range: (f64, f64), tol: f64) -> (Vec<f64>, Vec<Self::Point>) {
                        self.0.parameter_division(range, tol)
                    }
                }
            }
        }
        _ => unsupported_input(&ty),
    };
    derive_try!(traits_scope(implementation)).into()
}

/// Derive macro generating an impl of the trait `ParameterDivision2D` for Enums or single field tuple structs.
#[proc_macro_error]
#[proc_macro_derive(ParameterDivision2D)]
pub fn derive_parameter_division_2d(input: TokenStream) -> TokenStream {
    let input = parse_macro_input!(input as DeriveInput);
    let trait_name = quote! { monstertruck_traits::ParameterDivision2D };
    let ty = input.ident;
    let generics_params = input.generics;
    let where_predicates = generics_params
        .where_clause
        .iter()
        .flat_map(|x| &x.predicates);
    let implementation = match input.data {
        Data::Enum(DataEnum { ref variants, .. }) => {
            let top_ty = derive_try!(enum_field_type(variants, &ty));
            let methods = methods! {
                variants, trait_name,
                fn parameter_division(&self, range: ((f64, f64), (f64, f64)), tol: f64) -> (Vec<f64>, Vec<f64>),
            };
            quote! {
                #[automatically_derived]
                impl #generics_params #trait_name for #ty #generics_params
                where
                    #(#where_predicates,)*
                    #top_ty: #trait_name,
                    {
                    #(#methods)*
                }
            }
        }
        Data::Struct(DataStruct { ref fields, .. }) => {
            let field_type = derive_try!(tuple_struct_field(fields));
            quote! {
                #[automatically_derived]
                impl #generics_params #trait_name for #ty #generics_params
                where
                    #(#where_predicates,)*
                    #field_type: #trait_name {
                    fn parameter_division(&self, range: ((f64, f64), (f64, f64)), tol: f64) -> (Vec<f64>, Vec<f64>) {
                        self.0.parameter_division(range, tol)
                    }
                }
            }
        }
        _ => unsupported_input(&ty),
    };
    derive_try!(traits_scope(implementation)).into()
}

/// Derive macro generating an impl of the trait `ParametricCurve` for Enums or single field tuple structs.
#[proc_macro_error]
#[proc_macro_derive(ParametricCurve)]
pub fn derive_parametric_curve(input: TokenStream) -> TokenStream {
    let input = parse_macro_input!(input as DeriveInput);
    let trait_name = quote! { monstertruck_traits::ParametricCurve };
    let ty = input.ident;
    let generics_params = input.generics;
    let where_predicates = generics_params
        .where_clause
        .iter()
        .flat_map(|x| &x.predicates);
    let implementation = match input.data {
        Data::Enum(DataEnum { ref variants, .. }) => {
            let top_ty = derive_try!(enum_field_type(variants, &ty));
            let methods = methods!(
                variants,
                trait_name,
                fn evaluate(&self, t: f64) -> Self::Point,
                fn derivative(&self, t: f64) -> Self::Vector,
                fn derivative_2(&self, t: f64) -> Self::Vector,
                fn derivative_n(&self, n: usize, t: f64) -> Self::Vector,
                fn derivatives(&self, n: usize, t: f64) -> monstertruck_traits::__private::CurveDerivatives<Self::Vector>,
                fn parameter_range(&self,) -> monstertruck_traits::ParameterRange,
                fn period(&self,) -> Option<f64>,
            );
            quote! {
                #[automatically_derived]
                impl #generics_params #trait_name for #ty #generics_params
                where #(#where_predicates,)*
                      #top_ty: #trait_name,
                      #top_ty: #trait_name,

                      Self: Clone {
                    type Point = <#top_ty as #trait_name>::Point;
                    type Vector = <#top_ty as #trait_name>::Vector;
                    #(#methods)*
                }
            }
        }
        Data::Struct(DataStruct { ref fields, .. }) => {
            let field_type = derive_try!(tuple_struct_field(fields));
            quote! {
                #[automatically_derived]
                impl #generics_params #trait_name for #ty #generics_params
                where #(#where_predicates,)*
                      #field_type: #trait_name,
                      Self: Clone, {
                    type Point = <#field_type as #trait_name>::Point;
                    type Vector = <#field_type as #trait_name>::Vector;
                    fn evaluate(&self, t: f64) -> Self::Point { self.0.evaluate(t) }
                    fn derivative(&self, t: f64) -> Self::Vector { self.0.derivative(t) }
                    fn derivative_2(&self, t: f64) -> Self::Vector { self.0.derivative_2(t) }
                    fn derivative_n(&self, n: usize, t: f64) -> Self::Vector { self.0.derivative_n(n, t) }
                    fn derivatives(&self, n: usize, t: f64) -> monstertruck_traits::__private::CurveDerivatives<Self::Vector> { self.0.derivatives(n, t) }
                    fn parameter_range(&self) -> monstertruck_traits::ParameterRange { self.0.parameter_range() }
                    fn period(&self) -> Option<f64> { self.0.period() }
                }
            }
        }
        _ => unsupported_input(&ty),
    };
    derive_try!(traits_scope(implementation)).into()
}

/// Derive macro generating an impl of the trait `ParametricSurface` for Enums or single field tuple structs.
#[proc_macro_error]
#[proc_macro_derive(ParametricSurface)]
pub fn derive_parametric_surface(input: TokenStream) -> TokenStream {
    let input = parse_macro_input!(input as DeriveInput);
    let trait_name = quote! { monstertruck_traits::ParametricSurface };
    let ty = input.ident;
    let generics_params = input.generics;
    let where_predicates = generics_params
        .where_clause
        .iter()
        .flat_map(|x| &x.predicates);
    let implementation = match input.data {
        Data::Enum(DataEnum { ref variants, .. }) => {
            let top_ty = derive_try!(enum_field_type(variants, &ty));
            let methods = methods!(
                variants,
                trait_name,
                fn evaluate(&self, s: f64, t: f64) -> Self::Point,
                fn derivative_u(&self, s: f64, t: f64) -> Self::Vector,
                fn derivative_v(&self, s: f64, t: f64) -> Self::Vector,
                fn derivative_uu(&self, s: f64, t: f64) -> Self::Vector,
                fn derivative_uv(&self, s: f64, t: f64) -> Self::Vector,
                fn derivative_vv(&self, s: f64, t: f64) -> Self::Vector,
                fn derivative_mn(&self, m: usize, n: usize, s: f64, t: f64) -> Self::Vector,
                fn derivatives(&self, max_order: usize, u: f64, v: f64) -> monstertruck_traits::__private::SurfaceDerivatives<Self::Vector>,
                fn parameter_range(&self,) -> (monstertruck_traits::ParameterRange, monstertruck_traits::ParameterRange),
                fn period_u(&self,) -> Option<f64>,
                fn period_v(&self,) -> Option<f64>,
            );
            quote! {
                #[automatically_derived]
                impl #generics_params #trait_name for #ty #generics_params
                where
                    #(#where_predicates,)*
                    #top_ty: #trait_name,
                    {
                    type Point = <#top_ty as #trait_name>::Point;
                    type Vector = <#top_ty as #trait_name>::Vector;
                    #(#methods)*
                }
            }
        }
        Data::Struct(DataStruct { ref fields, .. }) => {
            let field_type = derive_try!(tuple_struct_field(fields));
            quote! {
                #[automatically_derived]
                impl #generics_params #trait_name for #ty #generics_params
                where
                    #(#where_predicates,)*
                    #field_type: #trait_name, {
                    type Point = <#field_type as #trait_name>::Point;
                    type Vector = <#field_type as #trait_name>::Vector;
                    #[inline(always)]
                    fn evaluate(&self, s: f64, t: f64) -> Self::Point { self.0.evaluate(s, t) }
                    #[inline(always)]
                    fn derivative_u(&self, s: f64, t: f64) -> Self::Vector { self.0.derivative_u(s, t) }
                    #[inline(always)]
                    fn derivative_v(&self, s: f64, t: f64) -> Self::Vector { self.0.derivative_v(s, t) }
                    #[inline(always)]
                    fn derivative_uu(&self, s: f64, t: f64) -> Self::Vector { self.0.derivative_uu(s, t) }
                    #[inline(always)]
                    fn derivative_uv(&self, s: f64, t: f64) -> Self::Vector { self.0.derivative_uv(s, t) }
                    #[inline(always)]
                    fn derivative_vv(&self, s: f64, t: f64) -> Self::Vector { self.0.derivative_vv(s, t) }
                    #[inline(always)]
                    fn derivative_mn(&self, m: usize, n: usize, s: f64, t: f64) -> Self::Vector { self.0.derivative_mn(m, n, s, t) }
                    #[inline(always)]
                    fn derivatives(&self, max_order: usize, s: f64, t: f64) -> monstertruck_traits::__private::SurfaceDerivatives<Self::Vector> {
                        self.0.derivatives(max_order, s, t)
                    }
                    #[inline(always)]
                    fn parameter_range(&self,) -> ((std::ops::Bound<f64>, std::ops::Bound<f64>), (std::ops::Bound<f64>, std::ops::Bound<f64>)) {
                        self.0.parameter_range()
                    }
                    #[inline(always)]
                    fn period_u(&self) -> Option<f64> { self.0.period_u() }
                    #[inline(always)]
                    fn period_v(&self) -> Option<f64> { self.0.period_v() }
                }
            }
        }
        _ => unsupported_input(&ty),
    };
    derive_try!(traits_scope(implementation)).into()
}

/// Derive macro generating an impl of the trait `ParametricSurface3D` for Enums or single field tuple structs.
#[proc_macro_error]
#[proc_macro_derive(ParametricSurface3D)]
pub fn derive_parametric_surface3d(input: TokenStream) -> TokenStream {
    let input = parse_macro_input!(input as DeriveInput);
    let trait_name0 = quote! { monstertruck_traits::ParametricSurface };
    let trait_name1 = quote! { monstertruck_traits::ParametricSurface3D };
    let ty = input.ident;
    let generics_params = input.generics;
    let where_predicates = generics_params
        .where_clause
        .iter()
        .flat_map(|x| &x.predicates)
        .collect::<Vec<_>>();
    let implementation = match input.data {
        Data::Enum(DataEnum { ref variants, .. }) => {
            let top_ty = derive_try!(enum_field_type(variants, &ty));
            let methods0 = methods!(
                variants,
                trait_name0,
                fn evaluate(&self, s: f64, t: f64) -> Self::Point,
                fn derivative_u(&self, s: f64, t: f64) -> Self::Vector,
                fn derivative_v(&self, s: f64, t: f64) -> Self::Vector,
                fn derivative_uu(&self, s: f64, t: f64) -> Self::Vector,
                fn derivative_uv(&self, s: f64, t: f64) -> Self::Vector,
                fn derivative_vv(&self, s: f64, t: f64) -> Self::Vector,
                fn derivative_mn(&self, m: usize, n: usize, s: f64, t: f64) -> Self::Vector,
                fn derivatives(&self, max_order: usize, u: f64, v: f64) -> monstertruck_traits::__private::SurfaceDerivatives<Self::Vector>,
                fn parameter_range(&self,) -> (monstertruck_traits::ParameterRange, monstertruck_traits::ParameterRange),
                fn period_u(&self,) -> Option<f64>,
                fn period_v(&self,) -> Option<f64>,
            );
            let methods1 = methods!(
                variants,
                trait_name1,
                fn normal(&self, u: f64, v: f64) -> monstertruck_traits::__private::Vector3,
            );
            quote! {
                #[automatically_derived]
                impl #generics_params #trait_name0 for #ty #generics_params
                where
                    #(#where_predicates,)*
                    #top_ty: #trait_name0,
                    {
                    type Point = monstertruck_traits::__private::Point3;
                    type Vector = monstertruck_traits::__private::Vector3;
                    #(#methods0)*
                }

                #[automatically_derived]
                impl #generics_params #trait_name1 for #ty #generics_params
                where
                    #(#where_predicates,)*
                    #top_ty: #trait_name0,
                    {
                    #(#methods1)*
                }
            }
        }
        Data::Struct(DataStruct { ref fields, .. }) => {
            let field_type = derive_try!(tuple_struct_field(fields));
            quote! {
                #[automatically_derived]
                impl #generics_params #trait_name0 for #ty #generics_params
                where
                    #(#where_predicates,)*
                    #field_type: #trait_name0, {
                    type Point = monstertruck_traits::__private::Point3;
                    type Vector = monstertruck_traits::__private::Vector3;
                    fn evaluate(&self, s: f64, t: f64) -> Self::Point { self.0.evaluate(s, t) }
                    fn derivative_u(&self, s: f64, t: f64) -> Self::Vector { self.0.derivative_u(s, t) }
                    fn derivative_v(&self, s: f64, t: f64) -> Self::Vector { self.0.derivative_v(s, t) }
                    fn derivative_uu(&self, s: f64, t: f64) -> Self::Vector { self.0.derivative_uu(s, t) }
                    fn derivative_uv(&self, s: f64, t: f64) -> Self::Vector { self.0.derivative_uv(s, t) }
                    fn derivative_vv(&self, s: f64, t: f64) -> Self::Vector { self.0.derivative_vv(s, t) }
                    fn derivative_mn(&self, m: usize, n: usize, s: f64, t: f64) -> Self::Vector { self.0.derivative_mn(m, n, s, t) }
                    fn derivatives(&self, max_order: usize, u: f64, v: f64) -> monstertruck_traits::__private::SurfaceDerivatives<Self::Vector> {
                        self.0.derivatives(max_order, u, v)
                    }
                    fn parameter_range(&self,) -> (monstertruck_traits::ParameterRange, monstertruck_traits::ParameterRange) {
                        self.0.parameter_range()
                    }
                    fn period_u(&self) -> Option<f64> { self.0.period_u() }
                    fn period_v(&self) -> Option<f64> { self.0.period_v() }
                }
                #[automatically_derived]
                impl #generics_params #trait_name1 for #ty #generics_params
                where
                    #(#where_predicates,)*
                    #field_type: #trait_name0, {
                    fn normal(&self, u: f64, v: f64) -> monstertruck_traits::__private::Vector3 { self.0.normal(u, v) }
                }
            }
        }
        _ => unsupported_input(&ty),
    };
    derive_try!(traits_scope(implementation)).into()
}

/// Derive macro generating an impl of the trait `SearchNearestParameter<CurveParameter>` for Enums or single field tuple structs.
#[proc_macro_error]
#[proc_macro_derive(SearchNearestParameterD1)]
pub fn derive_snp_d1(input: TokenStream) -> TokenStream {
    let input = parse_macro_input!(input as DeriveInput);
    let trait_name = quote! { monstertruck_traits::SearchNearestParameter::<monstertruck_traits::CurveParameter> };
    let ty = input.ident;
    let generics_params = input.generics;
    let where_predicates = generics_params
        .where_clause
        .iter()
        .flat_map(|x| &x.predicates);
    let implementation = match input.data {
        Data::Enum(DataEnum { ref variants, .. }) => {
            let top_ty = derive_try!(enum_field_type(variants, &ty));
            let methods = methods!(
                variants,
                trait_name,
                fn search_nearest_parameter<H: Into<monstertruck_traits::SearchParameterHint1D>>(
                    &self,
                    pt: Self::Point,
                    hint: H,
                    trials: usize,
                ) -> Option<f64>,
            );
            quote! {
                #[automatically_derived]
                impl #generics_params #trait_name for #ty #generics_params
                where
                    #(#where_predicates,)*
                    #top_ty: #trait_name,
                    {
                    type Point = <#top_ty as #trait_name>::Point;
                    #(#methods)*
                }
            }
        }
        Data::Struct(DataStruct { ref fields, .. }) => {
            let field_type = derive_try!(tuple_struct_field(fields));
            quote! {
                #[automatically_derived]
                impl #generics_params #trait_name for #ty #generics_params
                where
                    #(#where_predicates,)*
                    #field_type: #trait_name, {
                    type Point = <#field_type as #trait_name>::Point;
                    fn search_nearest_parameter<H: Into<monstertruck_traits::SearchParameterHint1D>>(
                        &self,
                        pt: Self::Point,
                        hint: H,
                        trials: usize,
                    ) -> Option<f64> {
                        self.0.search_nearest_parameter(pt, hint, trials)
                    }
                }
            }
        }
        _ => unsupported_input(&ty),
    };
    derive_try!(traits_scope(implementation)).into()
}

/// Derive macro generating an impl of the trait `SearchNearestParameter<SurfaceParameter>` for Enums or single field tuple structs.
#[proc_macro_error]
#[proc_macro_derive(SearchNearestParameterD2)]
pub fn derive_snp_d2(input: TokenStream) -> TokenStream {
    let input = parse_macro_input!(input as DeriveInput);
    let trait_name = quote! { monstertruck_traits::SearchNearestParameter::<monstertruck_traits::SurfaceParameter> };
    let ty = input.ident;
    let generics_params = input.generics;
    let where_predicates = generics_params
        .where_clause
        .iter()
        .flat_map(|x| &x.predicates);
    let implementation = match input.data {
        Data::Enum(DataEnum { ref variants, .. }) => {
            let top_ty = derive_try!(enum_field_type(variants, &ty));
            let methods = methods!(
                variants,
                trait_name,
                fn search_nearest_parameter<H: Into<monstertruck_traits::SearchParameterHint2D>>(
                    &self,
                    pt: Self::Point,
                    hint: H,
                    trials: usize,
                ) -> Option<(f64, f64)>,
            );
            quote! {
                #[automatically_derived]
                impl #generics_params #trait_name for #ty #generics_params
                where
                    #(#where_predicates,)*
                    #top_ty: #trait_name,
                    {
                    type Point = <#top_ty as #trait_name>::Point;
                    #(#methods)*
                }
            }
        }
        Data::Struct(DataStruct { ref fields, .. }) => {
            let field_type = derive_try!(tuple_struct_field(fields));
            quote! {
                #[automatically_derived]
                impl #generics_params #trait_name for #ty #generics_params
                where
                    #(#where_predicates,)*
                    #field_type: #trait_name, {
                    type Point = <#field_type as #trait_name>::Point;
                    fn search_nearest_parameter<H: Into<monstertruck_traits::SearchParameterHint2D>>(
                        &self,
                        pt: Self::Point,
                        hint: H,
                        trials: usize,
                    ) -> Option<(f64, f64)> {
                        self.0.search_nearest_parameter(pt, hint, trials)
                    }
                }
            }
        }
        _ => unsupported_input(&ty),
    };
    derive_try!(traits_scope(implementation)).into()
}

/// Derive macro generating an impl of the trait `SearchParameter<CurveParameter>` for Enums or single field tuple structs.
#[proc_macro_error]
#[proc_macro_derive(SearchParameterD1)]
pub fn derive_sp_d1(input: TokenStream) -> TokenStream {
    let input = parse_macro_input!(input as DeriveInput);
    let trait_name =
        quote! { monstertruck_traits::SearchParameter::<monstertruck_traits::CurveParameter> };
    let ty = input.ident;
    let generics_params = input.generics;
    let where_predicates = generics_params
        .where_clause
        .iter()
        .flat_map(|x| &x.predicates);
    let implementation = match input.data {
        Data::Enum(DataEnum { ref variants, .. }) => {
            let top_ty = derive_try!(enum_field_type(variants, &ty));
            let methods = methods!(
                variants,
                trait_name,
                fn search_parameter<H: Into<monstertruck_traits::SearchParameterHint1D>>(
                    &self,
                    pt: Self::Point,
                    hint: H,
                    trials: usize,
                ) -> Option<f64>,
            );
            quote! {
                #[automatically_derived]
                impl #generics_params #trait_name for #ty #generics_params
                where
                    #(#where_predicates,)*
                    #top_ty: #trait_name,
                    {
                    type Point = <#top_ty as #trait_name>::Point;
                    #(#methods)*
                }
            }
        }
        Data::Struct(DataStruct { ref fields, .. }) => {
            let field_type = derive_try!(tuple_struct_field(fields));
            quote! {
                #[automatically_derived]
                impl #generics_params #trait_name for #ty #generics_params
                where
                    #(#where_predicates,)*
                    #field_type: #trait_name, {
                    type Point = <#field_type as #trait_name>::Point;
                    fn search_parameter<H: Into<monstertruck_traits::SearchParameterHint1D>>(
                        &self,
                        pt: Self::Point,
                        hint: H,
                        trials: usize,
                    ) -> Option<f64> {
                        self.0.search_parameter(pt, hint, trials)
                    }
                }
            }
        }
        _ => unsupported_input(&ty),
    };
    derive_try!(traits_scope(implementation)).into()
}

/// Derive macro generating an impl of the trait `SearchParameter<SurfaceParameter>` for Enums or single field tuple structs.
#[proc_macro_error]
#[proc_macro_derive(SearchParameterD2)]
pub fn derive_sp_d2(input: TokenStream) -> TokenStream {
    let input = parse_macro_input!(input as DeriveInput);
    let trait_name =
        quote! { monstertruck_traits::SearchParameter::<monstertruck_traits::SurfaceParameter> };
    let ty = input.ident;
    let generics_params = input.generics;
    let where_predicates = generics_params
        .where_clause
        .iter()
        .flat_map(|x| &x.predicates);
    let implementation = match input.data {
        Data::Enum(DataEnum { ref variants, .. }) => {
            let top_ty = derive_try!(enum_field_type(variants, &ty));
            let methods = methods!(
                variants,
                trait_name,
                fn search_parameter<H: Into<monstertruck_traits::SearchParameterHint2D>>(
                    &self,
                    pt: Self::Point,
                    hint: H,
                    trials: usize,
                ) -> Option<(f64, f64)>,
            );
            quote! {
                #[automatically_derived]
                impl #generics_params #trait_name for #ty #generics_params
                where
                    #(#where_predicates,)*
                    #top_ty: #trait_name,
                    {
                    type Point = <#top_ty as #trait_name>::Point;
                    #(#methods)*
                }
            }
        }
        Data::Struct(DataStruct { ref fields, .. }) => {
            let field_type = derive_try!(tuple_struct_field(fields));
            quote! {
                #[automatically_derived]
                impl #generics_params #trait_name for #ty #generics_params
                where
                    #(#where_predicates,)*
                    #field_type: #trait_name, {
                    type Point = <#field_type as #trait_name>::Point;
                    fn search_parameter<H: Into<monstertruck_traits::SearchParameterHint2D>>(
                        &self,
                        pt: Self::Point,
                        hint: H,
                        trials: usize,
                    ) -> Option<(f64, f64)> {
                        self.0.search_parameter(pt, hint, trials)
                    }
                }
            }
        }
        _ => unsupported_input(&ty),
    };
    derive_try!(traits_scope(implementation)).into()
}

/// Derive macro generating an impl of the trait `Transformed<Matrix3>` for enums or single field tuple structs.
#[proc_macro_error]
#[proc_macro_derive(TransformedM3)]
pub fn derive_transform_m3(input: TokenStream) -> TokenStream {
    let input = parse_macro_input!(input as DeriveInput);
    let trait_name =
        quote! { monstertruck_traits::Transformed::<monstertruck_traits::__private::Matrix3> };
    let ty = input.ident;
    let generics_params = input.generics;
    let where_predicates = generics_params
        .where_clause
        .iter()
        .flat_map(|x| &x.predicates);
    let implementation = match input.data {
        Data::Enum(DataEnum { ref variants, .. }) => {
            let top_ty = derive_try!(enum_field_type(variants, &ty));
            let methods = methods!(
                variants,
                trait_name,
                fn transform_by(&mut self, mat: monstertruck_traits::__private::Matrix3) -> (),
                fn transformed(&self, mat: monstertruck_traits::__private::Matrix3) -> Self,
            );
            quote! {
                #[automatically_derived]
                impl #generics_params #trait_name for #ty #generics_params
                where
                    #(#where_predicates,)*
                    #top_ty: #trait_name,
                    {
                    #(#methods)*
                }
            }
        }
        Data::Struct(DataStruct { ref fields, .. }) => {
            let field_type = derive_try!(tuple_struct_field(fields));
            quote! {
                #[automatically_derived]
                impl #generics_params #trait_name for #ty #generics_params
                where
                    #(#where_predicates,)*
                    #field_type: #trait_name, {
                    fn transform_by(&mut self, mat: monstertruck_traits::__private::Matrix3) {
                        self.0.transform_by(mat)
                    }
                    fn transformed(&self, mat: monstertruck_traits::__private::Matrix3) -> Self {
                        Self(self.0.transformed(mat))
                    }
                }
            }
        }
        _ => unsupported_input(&ty),
    };
    derive_try!(traits_scope(implementation)).into()
}

/// Derive macro generating an impl of the trait `Transformed<Matrix4>` for enums or single field tuple structs.
#[proc_macro_error]
#[proc_macro_derive(TransformedM4)]
pub fn derive_transform_m4(input: TokenStream) -> TokenStream {
    let input = parse_macro_input!(input as DeriveInput);
    let trait_name =
        quote! { monstertruck_traits::Transformed::<monstertruck_traits::__private::Matrix4> };
    let ty = input.ident;
    let generics_params = input.generics;
    let where_predicates = generics_params
        .where_clause
        .iter()
        .flat_map(|x| &x.predicates);
    let implementation = match input.data {
        Data::Enum(DataEnum { ref variants, .. }) => {
            let top_ty = derive_try!(enum_field_type(variants, &ty));
            let methods = methods!(
                variants,
                trait_name,
                fn transform_by(&mut self, mat: monstertruck_traits::__private::Matrix4) -> (),
                fn transformed(&self, mat: monstertruck_traits::__private::Matrix4) -> Self,
            );
            quote! {
                #[automatically_derived]
                impl #generics_params #trait_name for #ty #generics_params
                where
                    #(#where_predicates,)*
                    #top_ty: #trait_name,
                    {
                    #(#methods)*
                }
            }
        }
        Data::Struct(DataStruct { ref fields, .. }) => {
            let field_type = derive_try!(tuple_struct_field(fields));
            quote! {
                #[automatically_derived]
                impl #generics_params #trait_name for #ty #generics_params
                where
                    #(#where_predicates,)*
                    #field_type: #trait_name, {
                    fn transform_by(&mut self, mat: monstertruck_traits::__private::Matrix4) {
                        self.0.transform_by(mat)
                    }
                    fn transformed(&self, mat: monstertruck_traits::__private::Matrix4) -> Self {
                        Self(self.0.transformed(mat))
                    }
                }
            }
        }
        _ => unsupported_input(&ty),
    };
    derive_try!(traits_scope(implementation)).into()
}
/// Derive macro generating an impl of the trait `StepLength` for enums or single field tuple structs.
#[proc_macro_error]
#[proc_macro_derive(StepLength)]
pub fn derive_step_length(input: TokenStream) -> TokenStream {
    let input = parse_macro_input!(input as DeriveInput);
    let trait_name = quote! { monstertruck_io::step::save::StepLength };
    let ty = input.ident;
    let generics_params = input.generics;
    let where_predicates = generics_params
        .where_clause
        .iter()
        .flat_map(|x| &x.predicates);
    let implementation = match input.data {
        Data::Enum(DataEnum { ref variants, .. }) => {
            let top_ty = derive_try!(enum_field_type(variants, &ty));
            let methods = methods!(
                variants,
                trait_name,
                fn step_length(&self,) -> usize,
            );
            quote! {
                #[automatically_derived]
                impl #generics_params #trait_name for #ty #generics_params
                where
                    #(#where_predicates,)*
                    #top_ty: #trait_name,
                    {
                    #(#methods)*
                }
            }
        }
        Data::Struct(DataStruct { ref fields, .. }) => {
            let field_type = derive_try!(tuple_struct_field(fields));
            quote! {
                #[automatically_derived]
                impl #generics_params #trait_name for #ty #generics_params
                where
                    #(#where_predicates,)*
                    #field_type: #trait_name, {
                    fn step_length(&self) -> usize { self.0.step_length() }
                }
            }
        }
        _ => unsupported_input(&ty),
    };
    derive_try!(io_scope(implementation)).into()
}

// Internal generator shared by the `StepFormat` derive and its
// upstream-compatible `DisplayByStep` alias.
fn derive_step_format_impl(input: TokenStream) -> TokenStream {
    let input = parse_macro_input!(input as DeriveInput);
    let trait_name = quote! { monstertruck_io::step::save::StepFormat };
    let ty = input.ident;
    let generics_params = input.generics;
    let where_predicates = generics_params
        .where_clause
        .iter()
        .flat_map(|x| &x.predicates);
    let implementation = match input.data {
        Data::Enum(DataEnum { ref variants, .. }) => {
            let top_ty = derive_try!(enum_field_type(variants, &ty));
            let methods = methods!(
                variants,
                trait_name,
                fn fmt(&self, idx: usize, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result,
            );
            quote! {
                #[automatically_derived]
                impl #generics_params #trait_name for #ty #generics_params
                where
                    #(#where_predicates,)*
                    #top_ty: #trait_name,
                    {
                    #(#methods)*
                }
            }
        }
        Data::Struct(DataStruct { ref fields, .. }) => {
            let field_type = derive_try!(tuple_struct_field(fields));
            quote! {
                #[automatically_derived]
                impl #generics_params #trait_name for #ty #generics_params
                where
                    #(#where_predicates,)*
                    #field_type: #trait_name, {
                    fn fmt(&self, idx: usize, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
                        #trait_name::fmt(&self.0, idx, f)
                    }
                }
            }
        }
        _ => unsupported_input(&ty),
    };
    derive_try!(io_scope(implementation)).into()
}

/// Derive macro generating an impl of the trait `StepFormat` for enums or single field tuple structs.
#[proc_macro_error]
#[proc_macro_derive(StepFormat)]
pub fn derive_step_format(input: TokenStream) -> TokenStream { derive_step_format_impl(input) }

// `DisplayByStep` is the upstream `truck-stepio` derive name. We renamed
// the trait to `StepFormat` because "display by step" reads as
// "display step-by-step" rather than "format as STEP". The alias is
// kept registered so external code ported from `truck` still compiles
// against `monstertruck`; new code should use `StepFormat`. Slated for
// removal once downstream callers are off the old name.
/// Deprecated alias for [`StepFormat`].
#[proc_macro_error]
#[proc_macro_derive(DisplayByStep)]
pub fn derive_display_by_step(input: TokenStream) -> TokenStream { derive_step_format_impl(input) }

/// Derive macro generating an impl of the trait `StepCurve` for enums or single field tuple structs.
#[proc_macro_error]
#[proc_macro_derive(StepCurve)]
pub fn derive_step_curve(input: TokenStream) -> TokenStream {
    let input = parse_macro_input!(input as DeriveInput);
    let trait_name = quote! { monstertruck_io::step::save::StepCurve };
    let ty = input.ident;
    let generics_params = input.generics;
    let where_predicates = generics_params
        .where_clause
        .iter()
        .flat_map(|x| &x.predicates);
    let implementation = match input.data {
        Data::Enum(DataEnum { ref variants, .. }) => {
            let top_ty = derive_try!(enum_field_type(variants, &ty));
            let methods = methods!(
                variants,
                trait_name,
                fn same_sense(&self,) -> bool,
            );
            quote! {
                #[automatically_derived]
                impl #generics_params #trait_name for #ty #generics_params
                where
                    #(#where_predicates,)*
                    #top_ty: #trait_name,
                    {
                    #(#methods)*
                }
            }
        }
        Data::Struct(DataStruct { ref fields, .. }) => {
            let field_type = derive_try!(tuple_struct_field(fields));
            quote! {
                #[automatically_derived]
                impl #generics_params #trait_name for #ty #generics_params
                where
                    #(#where_predicates,)*
                    #field_type: #trait_name, {
                    fn same_sense(&self) -> bool { #trait_name::same_sense(&self.0) }
                }
            }
        }
        _ => unsupported_input(&ty),
    };
    derive_try!(io_scope(implementation)).into()
}

/// Derive macro generating an impl of the trait `StepSurface` for enums or single field tuple structs.
#[proc_macro_error]
#[proc_macro_derive(StepSurface)]
pub fn derive_step_surface(input: TokenStream) -> TokenStream {
    let input = parse_macro_input!(input as DeriveInput);
    let trait_name = quote! { monstertruck_io::step::save::StepSurface };
    let ty = input.ident;
    let generics_params = input.generics;
    let where_predicates = generics_params
        .where_clause
        .iter()
        .flat_map(|x| &x.predicates);
    let implementation = match input.data {
        Data::Enum(DataEnum { ref variants, .. }) => {
            let top_ty = derive_try!(enum_field_type(variants, &ty));
            let methods = methods!(
                variants,
                trait_name,
                fn same_sense(&self,) -> bool,
            );
            quote! {
                #[automatically_derived]
                impl #generics_params #trait_name for #ty #generics_params
                where
                    #(#where_predicates,)*
                    #top_ty: #trait_name,
                    {
                    #(#methods)*
                }
            }
        }
        Data::Struct(DataStruct { ref fields, .. }) => {
            let field_type = derive_try!(tuple_struct_field(fields));
            quote! {
                #[automatically_derived]
                impl #generics_params #trait_name for #ty #generics_params
                where
                    #(#where_predicates,)*
                    #field_type: #trait_name, {
                    fn same_sense(&self) -> bool { #trait_name::same_sense(&self.0) }
                }
            }
        }
        _ => unsupported_input(&ty),
    };
    derive_try!(io_scope(implementation)).into()
}

/// Derive macro generating an impl of the trait `ToSameGeometry<Self>` for structs.
#[proc_macro_error]
#[proc_macro_derive(SelfSameGeometry)]
pub fn derive_self_same_geometry(input: TokenStream) -> TokenStream {
    let input = parse_macro_input!(input as DeriveInput);
    let ident = input.ident;
    let generics_params = input.generics;
    let implementation = quote! {
        #[automatically_derived]
        impl #generics_params monstertruck_traits::ToSameGeometry<Self> for #ident #generics_params where Self: Clone {
            #[inline(always)]
            fn to_same_geometry(&self) -> Self { self.clone() }
        }
    };
    derive_try!(traits_scope(implementation)).into()
}

#[cfg(test)]
mod validation_tests {
    use super::*;

    #[test]
    fn enum_validation_checks_every_variant_shape() {
        let valid: DeriveInput = parse_quote! { enum Example { First(u8), Second(u16) } };
        let Data::Enum(data) = &valid.data else {
            unreachable!()
        };
        let Ok(first) = enum_field_type(&data.variants, &valid.ident) else {
            unreachable!()
        };
        assert_eq!(quote! { #first }.to_string(), "u8");

        let invalid: DeriveInput = parse_quote! { enum Invalid { First(u8), Second(u16, u32) } };
        let Data::Enum(data) = &invalid.data else {
            unreachable!()
        };
        assert!(enum_field_type(&data.variants, &invalid.ident).is_err());
    }

    #[test]
    fn invalid_delegation_shapes_return_syn_errors() {
        let empty: DeriveInput = parse_quote! { enum Empty {} };
        let Data::Enum(data) = &empty.data else {
            unreachable!()
        };
        assert!(enum_field_type(&data.variants, &empty.ident).is_err());

        let named: DeriveInput = parse_quote! { struct Named { value: u8 } };
        let Data::Struct(data) = &named.data else {
            unreachable!()
        };
        assert!(tuple_struct_field(&data.fields).is_err());
    }
}

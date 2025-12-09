//! io_utils
// Copyright (c) 2021 10x Genomics, Inc. All rights reserved.
#![deny(missing_docs)]

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓
// WRITE STUFF
// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

/// fwriteln! is just like writeln! except that it has an expect call tacked on.
/// fwrite! similar.
#[macro_export]
macro_rules! fwriteln {
    ($f:expr, $u:expr) => {
        writeln!( $f, $u ).expect("writeln! failed")
    };
    ($f:expr, $u:expr, $($x:tt)*) => {
        writeln!( $f, $u, $($x)* )
            .unwrap_or_else(|_| panic!( "writeln! failed while writing \"{}\"", $u ) )
    };
}

/// fwrite! is just like write! except that it has an expect call tacked on.
#[macro_export]
macro_rules! fwrite {
    ($f:expr, $u:expr) => {
        write!( $f, $u ).expect( "write! failed" )
    };
    ($f:expr, $u:expr, $($x:tt)*) => {
        write!( $f, $u, $($x)* )
            .unwrap_or_else(|_| panic!( "write! failed while writing \"{}\"", $u ) )
    };
}

#[test]
fn test_fwrite() {
    use std::io::Write;
    let mut buf = vec![];
    fwrite!(&mut buf, "please delete this entire crate");
}

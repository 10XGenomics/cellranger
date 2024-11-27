// Copyright (c) 2021 10x Genomics, Inc. All rights reserved.

// This file contains miscellaneous utilities for input and output.

use bincode::{deserialize_from, serialize_into};
use serde::de::DeserializeOwned;
use serde::Serialize;
use std::ffi::OsStr;
use std::fmt::Debug;
use std::fs;
use std::fs::File;
use std::io::prelude::*;
use std::path::Path;

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓
// GET CONTENTS OF DIRECTORY
// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

pub fn dir_list(d: &str) -> Vec<String> {
    let x = fs::read_dir(d).unwrap_or_else(|_| panic!("failed to read directory {d}"));
    let mut y = Vec::<String>::new();
    for f in x {
        let s: String = f.unwrap().file_name().into_string().unwrap();
        y.push(s);
    }
    y.sort();
    y
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓
// TEST FOR EXISTENCE OF FILE
// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

pub fn path_exists(p: impl AsRef<Path>) -> bool {
    p.as_ref().exists()
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓
// WRITE STUFF
// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// fwriteln! is just like writeln! except that it has an expect call tacked on.
// fwrite! similar.

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

// fwrite! is just like write! except that it has an expect call tacked on.

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

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓
// OPEN FILES
// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

#[macro_export]
macro_rules! open_for_read {
    ($filename:expr) => {
        ::std::io::BufReader::new(
            ::std::fs::File::open(::core::convert::AsRef::<::std::path::Path>::as_ref(
                $filename,
            ))
            .unwrap_or_else(|_| {
                panic!(
                    "Could not open file \"{}\"",
                    ::core::convert::AsRef::<::std::path::Path>::as_ref($filename)
                        .to_string_lossy(),
                )
            }),
        )
    };
}

#[macro_export]
macro_rules! open_for_write_new {
    ($filename:expr) => {
        ::std::io::BufWriter::new(
            ::std::fs::File::create(::core::convert::AsRef::<::std::path::Path>::as_ref(
                $filename,
            ))
            .unwrap_or_else(|_| {
                panic!(
                    "Could not create file \"{}\"",
                    ::core::convert::AsRef::<::std::path::Path>::as_ref($filename)
                        .to_string_lossy()
                )
            }),
        )
    };
}

pub fn open_lz4<P: AsRef<Path>>(filename: P) -> lz4::Decoder<File> {
    let f = File::open(filename).expect("Failed to open file for reading");
    lz4::Decoder::new(f).expect("Failed to create lz4 decoder")
}

// If you accidentally pass a gzipped file to this it will succeed in opening the file,
// but then when you try to run read_line, the read will return !is_ok().  This seems horrible.

pub fn open_maybe_compressed<P: AsRef<Path>>(filename: P) -> Box<dyn Read> {
    match filename.as_ref().extension().and_then(OsStr::to_str) {
        Some("lz4") => Box::new(open_lz4(filename)) as Box<dyn Read>,
        _ => Box::new(File::open(filename).expect("Failed to open file for reading"))
            as Box<dyn Read>,
    }
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓
// READ A FILE
// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓
// CODE TO DO READS AND WRITES USING SERDE
// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

pub fn write_obj<T: Serialize, P: AsRef<Path> + Debug>(g: &T, filename: P) {
    let f = match std::fs::File::create(&filename) {
        Err(err) => panic!("couldn't create file {filename:?}: {err}"),
        Ok(f) => f,
    };
    let mut writer = std::io::BufWriter::new(f);
    serialize_into(&mut writer, &g)
        .unwrap_or_else(|_| panic!("write_obj of file {filename:?} failed"));
}

pub fn read_obj<T: DeserializeOwned, P: AsRef<Path> + Debug>(filename: P) -> T {
    let f = match std::fs::File::open(&filename) {
        Err(err) => panic!("couldn't open file {filename:?}: {err}"),
        Ok(f) => f,
    };
    let mut reader = std::io::BufReader::new(f);
    deserialize_from(&mut reader).unwrap_or_else(|_| panic!("read_obj of file {filename:?} failed"))
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓
// PRINT MACRO
// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// Print a list of things, useful for debugging.
// Example: if x is 3 and y is y, then printme!(x,y) yields
// x = 3, y = 7,

#[allow(unused_macros)]
#[macro_export]
macro_rules! printme {
        ( $( $x:expr ),* ) => {
            println!(concat!( $( stringify!($x), " = {}, ", )* ), $($x,)*)
        }
    }

#[allow(unused_macros)]
#[macro_export]
macro_rules! eprintme {
        ( $( $x:expr ),* ) => {
            eprintln!(concat!( $( stringify!($x), " = {}, ", )* ), $($x,)*)
        }
    }

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓
// READ FILE TO STRING AND PRINT FILE NAME IF IT DOESN'T EXIST
// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

pub fn read_to_string_safe<P: AsRef<Path>>(path: P) -> String {
    fs::read_to_string(&path).unwrap_or_else(|_| {
        panic!(
            "Could not open file \"{}\".",
            path.as_ref().to_str().unwrap()
        )
    })
}

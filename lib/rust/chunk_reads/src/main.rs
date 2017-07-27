// `error_chain!` can recurse deeply
#![recursion_limit = "1024"]

#[macro_use]
extern crate error_chain;

#[macro_use]
extern crate serde_derive;

#[macro_use]
extern crate serde_json;

extern crate docopt;
extern crate fastq;
extern crate rustc_serialize;
extern crate rayon;
extern crate tempdir;

#[cfg(test)]
extern crate fastq_10x;

use rayon::prelude::*;

use std::io::Read;
use std::path::{Path, PathBuf};

use std::fs::File;
use std::io;
use docopt::Docopt;
use fastq::{parse_path, Record};
use serde_json::Value;
use serde_json::map::Map;

// We'll put our errors in an `errors` module, and other modules in
// this crate will `use errors::*;` to get access to everything
// `error_chain!` creates.
mod errors {
    // Create the Error, ErrorKind, ResultExt, and Result types
    error_chain! {
        foreign_links {
            Io(::std::io::Error) #[doc = "Link to a `std::error::Error` type."];
            Json(::serde_json::Error);
        }
    }
}

use errors::*;


const USAGE: &'static str = "
Usage:
  chunk_reads [options] <output-path> <prefix> [--r1=<fn>] [--r2=<fn>] [--i1=<fn>] [--i2=<fn>] [--martian-args=<args-file>]
  chunk_reads (-h | --help)

Options:
  -h --help            Show this screen.
  --reads-per-fastq=N  Number of reads per FASTQ chunk [default: 50000000]
";

#[derive(Serialize, Deserialize)]
pub struct MartianArgs {
    read_chunks: ReadChunk,
    reads_interleaved: bool,
}

#[derive(Serialize, Deserialize)]
pub struct ReadChunk {
    R1: Option<String>,
    R2: Option<String>,
    I1: Option<String>,
    I2: Option<String>,
}

#[derive(Debug, RustcDecodable, Clone)]
pub struct Args {
    arg_output_path: String,
    arg_prefix: String,
    flag_reads_per_fastq: usize,
    flag_martian_args: Option<String>,
    flag_r1: Option<String>,
    flag_r2: Option<String>,
    flag_i1: Option<String>,
    flag_i2: Option<String>,
}


/// Read JSON from a chunk file
fn read_file<P: AsRef<Path>>(filename: P) -> Result<String> {
    let mut f = try!(File::open(filename.as_ref()));
    let mut buf = String::new();
    let s = try!(f.read_to_string(&mut buf));
    Ok(buf)
}

fn write_json<P: AsRef<Path>>(filename: P, val: &serde_json::Value) -> Result<()>
{
    let mut f = File::create(&filename)?;
    serde_json::to_writer_pretty(&mut f, val)?;
    Ok(())
}

fn main() {
    if let Err(ref e) = run() {
        println!("error: {}", e);

        for e in e.iter().skip(1) {
            println!("caused by: {}", e);
        }

        // The backtrace is not always generated. Try to run this example
        // with `RUST_BACKTRACE=1`.

        if let Some(backtrace) = e.backtrace() {
            println!("------------");
            println!("If you believe this is a bug in bamtofastq, please report a bug to software@10xgenomics.com.");
            println!("{:?}", backtrace);
        }

        ::std::process::exit(1);
    }
}


fn run() -> Result<Vec<(Option<String>, Option<String>, Option<String>, Option<String>)>> {

    println!("chunk_reads v{}", "VERSION");
    let args: Args = Docopt::new(USAGE)
                         .and_then(|d| d.decode())
                         .unwrap_or_else(|e| e.exit());

    run_args(args)
}


fn run_args(args: Args) -> Result<Vec<(Option<String>, Option<String>, Option<String>, Option<String>)>> {

    // Load arguments
    let (read_chunks, interleaved) =
        if args.flag_martian_args.is_some() {
            let path = args.flag_martian_args.as_ref().unwrap();
            let martian_args_string = read_file(&Path::new(path))?;
            let martian_args: MartianArgs = serde_json::from_str(&martian_args_string)?;
            (martian_args.read_chunks, martian_args.reads_interleaved)
        } else {
            (ReadChunk {
                R1: args.flag_r1.clone(),
                R2: args.flag_r2.clone(),
                I1: args.flag_i1.clone(),
                I2: args.flag_i2.clone(),
            }, false)
        };

    let mut chunks = Vec::new();

    if read_chunks.R1.is_some() {
        chunks.push(("R1", read_chunks.R1.unwrap().clone(), interleaved));
    }

    if read_chunks.R2.is_some() {
        chunks.push(("R2", read_chunks.R2.unwrap().clone(), false));
    }

    if read_chunks.I1.is_some() {
        chunks.push(("I1", read_chunks.I1.unwrap().clone(), false));
    }

    if read_chunks.I2.is_some() {
        chunks.push(("I2", read_chunks.I2.unwrap().clone(), false));
    }

    println!("chunks: {:?}", chunks);
    let out_path = Path::new(&args.arg_output_path);

    let mut _out_chunks = Vec::new();
    rayon::initialize(rayon::Configuration::new().num_threads(2));
    chunks.into_par_iter().map(|(name, in_file, interleaved)| {

        let read_per_chunk =
            if interleaved {
                args.flag_reads_per_fastq * 2
            } else {
                args.flag_reads_per_fastq
            };

        chunk_fastq(Path::new(&in_file), read_per_chunk, out_path, name)
    }).collect_into(&mut _out_chunks);

    let mut out_chunks = vec![];
    for oc in _out_chunks {
        let o = oc?;
        out_chunks.push(o);
    }

    // Write out the new read_chunk results
    let mut out_read_chunks = vec![];
    let mut out_read_sets = vec![];

    for i in 0 .. out_chunks[0].1.len() {

        let mut bt = Map::new();
        for &(ref read_name, ref files) in out_chunks.iter() {
            bt.insert(read_name.clone(), Value::String(files[i].to_str().unwrap().to_string()));
        }

        {
            let get_key = |name| {
                bt.get(name).map(|v| v.as_str().unwrap().to_string())
            };

            // Tuple of output filenames -- used for testing
            let read_set = (get_key("R1"), get_key("R2"), get_key("I1"), get_key("I2"));
            out_read_sets.push(read_set);
        }

        out_read_chunks.push(Value::Object(bt));
    }

    let mut out_chunks_path = out_path.to_path_buf();
    out_chunks_path.push(Path::new("read_chunks.json"));
    write_json(out_chunks_path, &Value::Array(out_read_chunks));

    println!("{:?}", &out_chunks);
    Ok(out_read_sets)
}


fn make_fastq_path(out_path: &Path, out_prefix: &str, chunk_number: usize) -> PathBuf {
    let mut p = out_path.to_path_buf();
    let filename = format!("{}-{:04}.fastq", out_prefix, chunk_number);
    p.push(Path::new(&filename));
    p
}

fn chunk_fastq(p: &Path, nrecord: usize, out_path: &Path, out_prefix: &str) -> Result<(String, Vec<PathBuf>)> {

    println!("opening: {:?}", p);
    let rr = parse_path(Some(p), |mut parser| {
        let mut paths: Vec<PathBuf> = vec![];
        let mut total = 0;
        let mut this_chunk = 0;
        let mut total_chunks = 0;


        let path = make_fastq_path(out_path, &out_prefix, total_chunks);
        let output = File::create(&path).unwrap();
        paths.push(path);
        let mut bufwrite = std::io::BufWriter::new(output);

        let ret = parser.each(|rec| {
            total += 1;
            this_chunk += 1;
            rec.write(&mut bufwrite);

            if this_chunk == nrecord {
                total_chunks += 1;
                let path = make_fastq_path(out_path, &out_prefix, total_chunks);
                let output = File::create(&path).unwrap();
                paths.push(path);
                bufwrite = std::io::BufWriter::new(output);
                this_chunk = 0;
            }

            true
        }).map(|r| paths);

        println!("got recs: {}", total);
        ret
    });


    match rr {
        Ok(Ok(paths)) => Ok((out_prefix.to_string(), paths)),
        Ok(Err(v)) => Err(v).chain_err(|| "fastq parsing error"),
        Err(v) => Err(v).chain_err(|| "fastq file opening error"),
    }
}




#[cfg(test)]
mod tests {
    use tempdir;
    use super::*;
    use fastq_10x::*;
    use std::collections::HashMap;

    type ReadSet = HashMap<Vec<u8>, RawReadSet>;

    pub fn load_fastq_set<I: Iterator<Item=RawReadSet>>(reads: &mut ReadSet, iter: I) {
        for r in iter {
            reads.insert((r.0).0.clone(), r);
        }
    }


    pub fn strict_compare_read_sets(orig_set: ReadSet, new_set: ReadSet) {

        assert_eq!(orig_set.len(), new_set.len());

        let mut keys1: Vec<Vec<u8>> = orig_set.keys().cloned().collect();
        keys1.sort();

        let mut keys2: Vec<Vec<u8>> = new_set.keys().cloned().collect();
        keys2.sort();

        for (k1, k2) in keys1.iter().zip(keys2.iter()) {
            assert_eq!(k1, k2);
            assert_eq!(orig_set.get(k1), new_set.get(k2));
        }
    }


    // Ensure reads are equivalent through mkfastq code path
    #[test]
    fn test_mkfastq() {
        let tempdir = tempdir::TempDir::new("chunk_reads_test").expect("create temp dir");
        let tmp_path = tempdir.path();

        let args = Args {
            arg_output_path: tmp_path.to_str().unwrap().to_string(),
            arg_prefix: "p1".to_string(),
            flag_reads_per_fastq: 1000,
            flag_martian_args: None,
            flag_r1: Some("test/mkfastq/pbmc8k_S1_L007_R1_001.fastq.gz".to_string()),
            flag_r2: Some("test/mkfastq/pbmc8k_S1_L007_R2_001.fastq.gz".to_string()),
            flag_i1: Some("test/mkfastq/pbmc8k_S1_L007_I1_001.fastq.gz".to_string()),
            flag_i2: None,
        };

        let out_path_sets = super::run_args(args).unwrap();

        println!("opening orig");
        let original_read_iter = open_fastq_pair_iter(
            "test/mkfastq/pbmc8k_S1_L007_R1_001.fastq.gz",
            "test/mkfastq/pbmc8k_S1_L007_R2_001.fastq.gz",
            Some("test/mkfastq/pbmc8k_S1_L007_I1_001.fastq.gz"));

        let mut orig_reads = ReadSet::new();
        load_fastq_set(&mut orig_reads, original_read_iter);

        println!("opening chunks");
        let mut output_reads = ReadSet::new();
        for (r1, r2, i1, i2) in out_path_sets {
            load_fastq_set(&mut output_reads, open_fastq_pair_iter(r1.unwrap(), r2.unwrap(), i1));
        }

        println!("comparing {} reads", orig_reads.len());
        strict_compare_read_sets(orig_reads, output_reads);
    }

    // Ensure reads are equivalent through bcl_processor code path
    #[test]
    fn test_bclprocessor() {
        let tempdir = tempdir::TempDir::new("outs1").expect("create temp dir");
        let tmp_path = tempdir.path();

        let args = Args {
            arg_output_path: tmp_path.to_str().unwrap().to_string(),
            arg_prefix: "p1".to_string(),
            flag_reads_per_fastq: 1000,
            flag_martian_args: Some("test/bcl_processor/chunk.json".to_string()),
            flag_r1: None,
            flag_r2: None,
            flag_i1: None,
            flag_i2: None,
        };

        let out_path_sets = super::run_args(args).unwrap();

        println!("opening orig");
        let original_read_iter = open_interleaved_fastq_pair_iter(
            "test/bcl_processor/read-RA_si-GTAATTGC_lane-007-chunk-005.fastq.gz",
            Some("test/bcl_processor/read-I1_si-GTAATTGC_lane-007-chunk-005.fastq.gz"));

        let mut orig_reads = ReadSet::new();
        load_fastq_set(&mut orig_reads, original_read_iter);

        println!("opening chunks");
        let mut output_reads = ReadSet::new();
        for (r1, r2, i1, i2) in out_path_sets {
            load_fastq_set(&mut output_reads, open_interleaved_fastq_pair_iter(r1.unwrap(), i1));
        }

        println!("comparing {} reads", orig_reads.len());
        strict_compare_read_sets(orig_reads, output_reads);
    }

    // Ensure reads are equivalent through bcl_processor code path
    #[test]
    fn test_invalid() {
        let tempdir = tempdir::TempDir::new("outs1").expect("create temp dir");
        let tmp_path = tempdir.path();

        let args = Args {
            arg_output_path: tmp_path.to_str().unwrap().to_string(),
            arg_prefix: "p1".to_string(),
            flag_reads_per_fastq: 1000,
            flag_martian_args: Some("test/invalid_fastq/chunk.json".to_string()),
            flag_r1: None,
            flag_r2: None,
            flag_i1: None,
            flag_i2: None,
        };

        let out_path_sets = super::run_args(args);

        match out_path_sets {
            Ok(_) => assert!(false, "expected an error"),
            _ => ()
        }
    }

}

//! This tool is used the python code in lib/python/cellranger/reference.py to
//! quickly load properties of the Genes and Transcripts from the GTF/ FASTA
//! files.

use std::path::Path;
use transcriptome::python_gene_index;

pub fn main() {
    let args: Vec<String> = std::env::args().collect();

    if args.len() != 3 {
        println!(
            "Compute Gene Index data from a 10x Genomics reference directory and write to JSON"
        );
        println!("If an error occurs, the process will return 1 and write a message to stderr");
        println!("usage: gtf_to_gene_index <reference_path> <output_json_filename>");
        std::process::exit(1);
    }

    let reference_path = Path::new(&args[1]);
    let json_file = Path::new(&args[2]);
    let res = python_gene_index::write_gene_index(reference_path, json_file);

    if let Err(e) = res {
        // write message and stack trace, exit code = 1;
        eprintln!("{e}");
        for c in e.chain().skip(1) {
            eprintln!("\tCaused by: {c}");
        }
        std::process::exit(1);
    }
}

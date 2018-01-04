//
// Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
//

use std::fs::File;
use std::io;
use std::io::BufWriter;
use std::path::Path;
use std::str;
use serde;
use serde_json;
use rust_htslib::bam;

pub const PROC_BC_SEQ_TAG: &'static [u8]     = b"CB";
pub const PROC_UMI_SEQ_TAG: &'static [u8]    = b"UB";
pub const GENE_ID_TAG: &'static [u8]         = b"GX";
pub const EXTRA_FLAGS_TAG: &'static [u8]     = b"xf";

//pub const HIGH_CONF_MAPQ: u8 = 255;

bitflags! {
    #[derive(Default)]
    pub struct ExtraFlags: u32 {
        // Confidently mapped to transcriptome
        const CONF_MAPPED = 1u32;
        const LOW_SUPPORT_UMI = 2u32;
        const GENE_DISCORDANT = 4u32;
        // This read is representative for a molecule and can be treated as a UMI count.
        const UMI_COUNT = 8u32;
    }
}

/// Write an object to a json file
pub fn write_json_file<P: AsRef<Path>, T: serde::Serialize>(path: P, data: &T)
                                                        -> Result<(), io::Error> {
    serde_json::to_writer_pretty(&mut BufWriter::new(File::create(path)?), data)?;
    Ok(())
}

/// Read an object from a json file
pub fn read_json_file<P: AsRef<Path>, T: serde::de::DeserializeOwned>(path: P)
                                                                      -> Result<T, io::Error> {
    Ok(serde_json::from_reader(File::open(path)?)?)
}


/// Get an alignment record's cell barcode
pub fn get_read_barcode(record: &bam::Record) -> Option<Vec<u8>> {
    record.aux(PROC_BC_SEQ_TAG).map(|x| x.string()).map(|x| x.to_owned())
}

/// Get an alignment record's processed UMI sequence
pub fn get_read_umi(record: &bam::Record) -> Option<Vec<u8>> {
    record.aux(PROC_UMI_SEQ_TAG).map(|x| x.string()).map(|x| x.to_owned())
}

/// Get an alignment record's list of aligned gene ID strings
pub fn get_read_gene_ids(record: &bam::Record) -> Vec<String> {
    record.aux(GENE_ID_TAG).map_or(vec![], |x| str::from_utf8(x.string()).unwrap()
                             .split(';')
                             .map(|y| y.to_owned())
                             .collect()
    )
}

/// Get an alignment record's extra flags
pub fn get_read_extra_flags(record: &bam::Record) -> ExtraFlags {
    record.aux(EXTRA_FLAGS_TAG).map_or(Default::default(), |x| ExtraFlags::from_bits_truncate(x.integer() as u32))
}

/// Return true if a read should be considered for PCR-duplicate marking
pub fn is_read_dup_candidate(record: &bam::Record) -> bool {
    let is_primary_alignment = !record.is_secondary();
    let is_conf_mapped = get_read_extra_flags(&record).intersects(ExtraFlags::CONF_MAPPED);
    is_primary_alignment && is_conf_mapped
}

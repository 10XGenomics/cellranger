use std::fs;
use std::io;
use std::fs::File;
use std::cmp::{min, max};
use std::path::Path;
use std::mem::drop;
use itertools::Itertools;

use bincode::deserialize_from;

use rust_htslib::bam;
use rust_htslib::bam::Read;

#[derive(Deserialize, Debug)]
pub struct GzipHeader {
    pub id1: u8, pub id2: u8, pub cm: u8, pub flg: u8, 
    pub mtime: u32, pub xfl: u8, pub os: u8, pub xlen: u16,
}
#[derive(Deserialize, Debug)]
pub struct GzipExtra {
    pub id1: u8, pub id2: u8, pub slen: u16, pub bsize: u16,
}
#[derive(Deserialize, Debug)]
pub struct GzipISIZE {
    pub isize: u32,
}

fn parse_bgzf_header(mut bgzf: &mut File, total_bytes: u64) -> Option<u64> {
    use self::io::Seek;

    let cur_pos = bgzf.seek(io::SeekFrom::Current(0))
        .expect("Failed to get current position in BAM file");
    if cur_pos == total_bytes {
        return None;
    }

    let header: GzipHeader = deserialize_from(&mut bgzf)
        .expect("Failed to deserialize gzip header - invalid BAM file");

    if header.id1 != 31 || header.id2 != 139 {
        panic!("Invalid gzip header in BAM file");
    }

    // Determine the BGZF block size
    let extra: GzipExtra = deserialize_from(&mut bgzf)
        .expect("Failed to deserialize gzip extra field - invalid BAM file");
    if extra.id1 != 66 || extra.id2 != 67 || extra.slen != 2 {
        panic!("BSIZE field not found in BGZF header - invalid BAM file");
    }

    let next_pos = cur_pos + (extra.bsize as u64) + 1;
    bgzf.seek(io::SeekFrom::Start(next_pos))
        .expect("Failed to seek to next block position in BAM file");
    Some(next_pos as u64)
}


/// Iterate over BAM records while simultaneously tracking position of the current record's start
/// This is to get around a borrow-violation when using tell() inside a records() loop
pub struct BamTellIter<'a> {
    pub bam: &'a mut bam::Reader,
}
impl<'a> Iterator for BamTellIter<'a> {
    type Item = Result<(i64, bam::record::Record), bam::ReadError>;

    fn next(&mut self) -> Option<Result<(i64, bam::record::Record), bam::ReadError>> {
        let mut record = bam::record::Record::new();
        let pos = self.bam.tell();
        match self.bam.read(&mut record) {
            Err(bam::ReadError::NoMoreRecord) => None,
            Ok(())   => Some(Ok((pos, record))),
            Err(err) => Some(Err(err)),
        }
    }
}

/// Find a valid virtual offset for iteration by scanning
/// through reads from a given real byte offset.
fn find_valid_virtual_offset<K: Eq>(bam: &mut bam::Reader,
                                    pos: u64,
                                    chunk_bound_key: &Fn(bam::Record) -> Option<K>)
                                    -> Option<i64> {
    use self::bam::Read;
    bam.seek((pos << 16) as i64)
        .expect("Failed to seek BAM file");
    let mut n = 0;

    let mut prev_key = None;
    let tell_iter = BamTellIter { bam: bam };
    for (cur_voff, read) in tell_iter.map(|x| x.unwrap()) {
        n += 1;
        if n == 1000 {
            warn!("Taking a long time to find a chunk boundary. Check your chunk_bound_key function?");
        }

        let new_key = chunk_bound_key(read);

        // - If boundary func returned None, this is a good boundary
        // - If the key has changed, this is a good boundary
        if new_key.is_none() || prev_key.is_some() && prev_key != new_key {
            return Some(cur_voff);
        } else {
            prev_key = new_key;
        }
    }
    // Chunk has no valid start position or has no reads.
    None
}

/// Return block offsets of a bam file.
pub fn bam_block_offsets<P: AsRef<Path>>(bam_path : &P) -> Vec<u64> {
    let bam_bytes = fs::metadata(bam_path)
        .expect("Failed to get filesize of input BAM").len();
    // Get first block w/ alignment info.
    let bam = bam::Reader::from_path(bam_path).expect("Failed to open input BAM");
    // Real offset of first block
    let first_block = (bam.tell() >> 16) as u64;
    info!("Found first block at real offset {}", first_block);

    // Get BGZF block positions. Open as a regular file.
    let mut bgzf = File::open(bam_path).expect("Failed to open input BAM");
    let mut block_offsets: Vec<u64> = Vec::new();
    while let Some(offset) = parse_bgzf_header(&mut bgzf, bam_bytes) {
        if offset >= first_block {
            block_offsets.push(offset);
        }
    }
    drop(bgzf);
    info!("Found {} blocks", block_offsets.len());

    block_offsets
}


/// Split a BAM along chunk boundaries
/// Returns a list of tuples of (start, end) virtual offsets,
/// where end=None means EOF
pub fn chunk_bam_records<P: AsRef<Path>, K: Eq>(bam_path: &P,
                                                block_offsets: &Vec<u64>,
                                                chunk_bound_key: &Fn(bam::Record) -> Option<K>,
                                                chunk_size_gb: f64,
                                                max_chunks: u64) -> Vec<(i64, Option<i64>)> {
    use self::bam::Read;
    
    // Compute chunk size
    let bam_bytes = fs::metadata(bam_path)
        .expect("Failed to get filesize of input BAM").len();
    let size_gb = bam_bytes as f64 / 1e9_f64;
    let num_chunks = min(max_chunks, max(1, (size_gb as f64 / chunk_size_gb).ceil() as u64));
    info!("size_gb={}, chunk_size_gb={}, num_chunks={}", size_gb, chunk_size_gb, num_chunks);

    // Get first block w/ alignment info.
    let mut bam = bam::Reader::from_path(bam_path).expect("Failed to open input BAM");
    // Real offset of first block
    let first_block = (bam.tell() >> 16) as u64;
    info!("Found first block at real offset {}", first_block);

    let block_step = max(1, block_offsets.len() / num_chunks as usize);

    // Real block offsets to start boundary-searches at
    let start_offsets = block_offsets.into_iter().step(block_step).skip(1);

    // Virtual offsets to start chunks at
    let mut start_voffs: Vec<i64> = start_offsets
        .filter_map(|x| find_valid_virtual_offset(&mut bam, *x, chunk_bound_key))
        .collect();
    start_voffs.insert(0, (first_block << 16) as i64);

    // Remove duplicate chunk starts
    start_voffs.dedup();
    // Determine chunk ends
    let mut end_voffs: Vec<Option<i64>> = start_voffs.clone().into_iter().skip(1)
        .map(|x| Some(x)).collect();
    end_voffs.push(None);

    start_voffs.into_iter().zip(end_voffs.into_iter()).collect()
}

/// Iterate over BAM records starting at virtual offset A and ending at virtual offset B
/// (excluding the record starting at offset B)
pub struct BamChunkIter<'a> {
    tell_iter: BamTellIter<'a>,
    end: Option<i64>,
}
impl<'a> BamChunkIter<'a> {
    /// Takes a bam::Reader and a tuple of virtual offsets (start, end)
    /// where end==None signifies EOF
    pub fn new(bam: &'a mut bam::Reader, range: (i64, Option<i64>)) -> BamChunkIter<'a> {
        use self::bam::Read;

        bam.seek(range.0).expect("Failed to seek to start of BAM chunk");

        BamChunkIter {
            tell_iter: BamTellIter { bam: bam },
            end: range.1,
        }
    }
}
impl<'a> Iterator for BamChunkIter<'a> {
    type Item = Result<bam::record::Record, bam::ReadError>;

    fn next(&mut self) -> Option<Result<bam::record::Record, bam::ReadError>> {
        match self.tell_iter.next() {
            // Reached EOF
            None => None,
            // Read error
            Some(Err(err)) => Some(Err(err)),
            Some(Ok((pos, record))) => match self.end {
                // No chunk end specified, so always return the record
                None => Some(Ok(record)),
                Some(end) => match pos < end {
                    // Still inside chunk
                    true => Some(Ok(record)),
                    // Reached end of chunk
                    false => None,
                }
            }
        }
    }
}



#[cfg(test)]
mod tests {
    extern crate env_logger;
    extern crate rand;

    use std::path::{Path, PathBuf};
    use std::fs;
    use self::rand::Rng;
    use rust_htslib::bam;
    use std::collections::{HashMap};

    use chunk_bam;

    /// Generate a random nucleotide sequence of length len
    fn random_sequence<R: Rng>(len: usize, rng: &mut R) -> Vec<u8> {
        let nucs = b"ACGT";
        (0..len).map(|_| nucs[rng.gen_range(0, 4)]).collect()
    }

    fn create_sorted_bam<P: AsRef<Path>>(filename: P, num_tids: u64, reads_per_tid: u64) {
        use rust_htslib::bam::header::{Header, HeaderRecord};

        let mut rng = rand::thread_rng();

        let read_len = 100;
        let qual = vec![b'I'; read_len];
        let cigar = bam::record::CigarString(vec![bam::record::Cigar::Match(read_len as u32)]);

        let mut header = Header::new();
        for i in 0..num_tids {
            header.push_record(HeaderRecord::new(b"SQ").push_tag(b"SN", &format!("chr{}", 1+i)).push_tag(b"LN", &10000000));
        }
        let mut writer = bam::Writer::from_path(filename, &header)
            .expect("Failed to create bam writer");

        for tid in 0..num_tids {
            for pos in 0..reads_per_tid {
                let mut rec = bam::Record::new();
                let seq = random_sequence(read_len, &mut rng);
                rec.set(b"1", &cigar, &seq, &qual);
                rec.set_tid(tid as i32);
                rec.set_pos(pos as i32);
                writer.write(&rec).expect("Failed to write BAM record");
            }
        }
    }

    #[test]
    fn test_empty_bam() {
        let _ = env_logger::init();

        let mut test_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        test_dir.push("test_output/chunk_bam");
        fs::create_dir_all(&test_dir)
            .expect("Failed to create test output dir");

        let test_bam = test_dir.join("empty.bam");
        create_sorted_bam(&test_bam, 6, 0);

        // Chunk it
        let block_offsets = chunk_bam::bam_block_offsets(&test_bam);
        let chunks = chunk_bam::chunk_bam_records(&test_bam, &block_offsets, &|r| Some(r.tid()), 0.001, 256);
        println!("{:?}", chunks);

        // Read the chunks
        let mut bam = bam::Reader::from_path(test_bam).unwrap();

        let mut tidpos_count: HashMap<(i32, i32), u64> = HashMap::new();
        for range in chunks {
            for rec in chunk_bam::BamChunkIter::new(&mut bam, range).map(|x| x.unwrap()) {
                // Record the (tid, pos)
                *tidpos_count.entry((rec.tid(), rec.pos())).or_insert(0) += 1;
            }
        }
        assert!(tidpos_count.len() == 0);
    }

    #[test]
    fn test_chunk_bam() {
        //let _ = env_logger::init();

        // Create test BAM
        let mut test_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        test_dir.push("test_output/chunk_bam");
        fs::create_dir_all(&test_dir)
            .expect("Failed to create test output dir");

        let test_bam = test_dir.join("pos_sorted.bam");
        let num_tids = 6;
        let reads_per_tid = 10000;
        create_sorted_bam(&test_bam, num_tids, reads_per_tid);

        // Chunk it
        let block_offsets = chunk_bam::bam_block_offsets(&test_bam);
        let chunks = chunk_bam::chunk_bam_records(&test_bam, &block_offsets, &|r| Some(r.tid()), 0.001, 256);
        println!("{:?}", chunks);

        // Read the chunks
        let mut bam = bam::Reader::from_path(test_bam).unwrap();

        let mut tidpos_count: HashMap<(i32, i32), u64> = HashMap::new();

        for range in chunks {
            for rec in chunk_bam::BamChunkIter::new(&mut bam, range).map(|x| x.unwrap()) {
                // Record the (tid, pos)
                *tidpos_count.entry((rec.tid(), rec.pos())).or_insert(0) += 1;
            }
        }

        // Saw all the reads
        assert!(tidpos_count.len() == (num_tids*reads_per_tid) as usize);

        // Saw each read once
        for count in tidpos_count.values() {
            assert!(count == &1);
        }
    }
}

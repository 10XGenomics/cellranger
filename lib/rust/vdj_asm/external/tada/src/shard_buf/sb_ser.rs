//! Efficiently write Rust structs to shard files from multiple threads.
/*
use std::fs::File;
use std::io::{Write, BufWriter, BufReader, Read};
use std::collections::HashMap;
use std::sync::mpsc::sync_channel;
use std::sync::mpsc::{SyncSender, Receiver};
use std::path::{Path, PathBuf};

use std::thread;
use std::thread::JoinHandle;
use std::slice;
use std::mem;

use std::fs::create_dir;
use std::fs::{self, DirEntry};

use std::io;
use std::io::Seek;

use std::marker::PhantomData;

use bincode;
use rustc_serialize::{Encodable, Decodable};
use bincode::rustc_serialize::{EncodingResult, DecodingResult, encode_into, decode_from};


pub struct BincodeChunkWriter<T: Encodable, R: Fn(&T) -> u64> {
    writer: BufWriter<File>,
    key: R,
    buffers: HashMap<u64, Vec<T>>,
    entries: Vec<(u64, u64)>,
    max_len: usize,
}

// Write chunks w/ bincode
impl<T: Encodable, R: Fn(&T) -> u64> BincodeChunkWriter<T,R> {
    pub fn new<P: AsRef<Path>>(path: P, key: R, max_len: usize) -> BincodeChunkWriter<T, R> {
        let f = File::create(path).unwrap();
        let w = BufWriter::new(f);
        BincodeChunkWriter {
            writer: w,
            key: key,
            buffers: HashMap::new(),
            entries: Vec::new(),
            max_len: max_len
         }
    }

    pub fn send(&mut self, v: T) {
        let k = (self.key)(&v);
        let write = {
            let b = self.buffers.entry(k).or_insert_with(|| Vec::new());
            b.push(v);
            b.len() == self.max_len
        };

        // if a buffer fills up, write it out
        if write {
            let buf = self.buffers.remove(&k).unwrap();
            let current_pos = self.writer.seek(io::SeekFrom::Current(0)).unwrap();
            self.write(&buf);
            self.entries.push((k, current_pos));
        }
    }

    pub fn write(&mut self, t: &Vec<T>) {
        encode_into(t, &mut self.writer, bincode::SizeLimit::Infinite);
    }
}

impl<T: Encodable, R: Fn(&T) -> u64> Drop for BincodeChunkWriter<T,R> {
    fn drop(&mut self) {
        let _r = self.writer.flush();
    }
}
*/

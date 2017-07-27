//! Efficiently write Rust structs to shard files from multiple threads.

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

use std::marker::PhantomData;

mod sb_ser;


pub struct PodWriter {
    writer: BufWriter<File>,
}

/// Bit-blit arrays of Copy types to a file
impl PodWriter {
    pub fn new<P: AsRef<Path>>(path: P) -> PodWriter {
        let f = File::create(path).unwrap();
        let w = BufWriter::new(f);
        PodWriter { writer: w }
    }

    pub fn write<T: Copy>(&mut self, t: &Vec<T>) {

        unsafe {
            let data_ptr = t.as_ptr() as *mut u8;
            let len = t.len() * mem::size_of::<T>();
            let write_buf = slice::from_raw_parts(data_ptr, len);
            self.writer.write(write_buf).unwrap();
        }
        // self.writer.flush();
        // println!("flush res: {:?}", fr)
    }
}

impl Drop for PodWriter {
    fn drop(&mut self) {
        let _r = self.writer.flush();
    }
}


pub struct PodReader<T: Copy> {
    reader: BufReader<File>,
    phantom: PhantomData<T>,
}

/// Bit-blit arrays of Copy types to a file
impl<T: Copy> PodReader<T> {
    pub fn new<P: AsRef<Path>>(path: P) -> PodReader<T> {
        let f = File::open(path).expect("couldn't open");
        let w = BufReader::new(f);
        PodReader {
            reader: w,
            phantom: PhantomData,
        }
    }

    #[inline(never)]
    pub fn read(&mut self) -> Vec<T> {
        unsafe {
            let sz = mem::size_of::<T>();

            let mut byte_buf = Vec::new();
            let _ = self.reader.read_to_end(&mut byte_buf).unwrap();

            let items = byte_buf.len() / sz;
            let data_ptr = byte_buf.as_ptr() as *mut T;
            let vec = slice::from_raw_parts(data_ptr, items).to_vec();
            vec
        }
    }
}


// one possible implementation of fs::walk_dir only visiting files
fn visit_dirs(dir: &Path, cb: &Fn(&DirEntry)) -> io::Result<()> {
    if try!(fs::metadata(dir)).is_dir() {
        for entry in try!(fs::read_dir(dir)) {
            let entry = try!(entry);
            if try!(fs::metadata(entry.path())).is_dir() {
                try!(visit_dirs(&entry.path(), cb));
            } else {
                cb(&entry);
            }
        }
    }
    Ok(())
}


struct ShardReader<'a, T: Copy> {
    path: &'a Path,
    phantom: PhantomData<T>,
}

impl<'a, T: Copy> ShardReader<'a, T> {
    fn new(path: &'a Path) -> ShardReader<'a, T> {
        ShardReader {
            path: path,
            phantom: PhantomData,
        }
    }

    fn num_shards(&self) -> usize {
        let mut shard_count = 0;
        for entry in fs::read_dir(self.path).unwrap() {
            let e = entry.unwrap();
            let file_name_str = e.file_name();
            let file_name = file_name_str.to_str().expect("not a string");

            if file_name.starts_with("shard_file") {
                shard_count += 1;
            }
        }
        shard_count
    }

    fn get_reader(&self, idx: usize) -> PodReader<T> {
        let p = Path::join(self.path, format!("shard_file_{}.bin", idx));
        PodReader::new(p)
    }
}


struct ShardWriterThread<T: Copy> {
    num_shards: usize,
    rx: Receiver<(usize, Vec<T>)>,
    writers: HashMap<usize, PodWriter>,
}

impl<'a, T: 'static + Send + Copy> ShardWriterThread<T> {
    fn spawn(rx: Receiver<Option<(usize, Vec<T>)>>,
             mut writers: HashMap<usize, PodWriter>,
             num_shards: usize)
             -> thread::JoinHandle<HashMap<usize, usize>> {
        let thread_closure = move || {
            let mut counts: HashMap<usize, usize> = HashMap::new();

            loop {
                match rx.recv() {
                    Ok(Some((idx, v))) => {
                        let shard = idx % num_shards;
                        let cc = counts.entry(shard).or_insert(0);
                        *cc += v.len();

                        // Write to the correct output file
                        let writer = writers.get_mut(&shard)
                                            .expect("got shard index higher, but don't have \
                                                     writer");
                        writer.write(&v);

                    }

                    Ok(None) => {
                        // The writers are not having drop called
                        // so we manually flush them
                        for (_, w) in writers.iter_mut() {
                            let _r = w.writer.flush();
                        }
                        break;
                    }

                    Err(_) => break,
                }
            }
            counts
        };

        thread::spawn(thread_closure)
    }
}

pub struct Shard<'a, T: Copy + Send> {
    path: &'a Path,
    num_shards: usize,
    threads: Vec<JoinHandle<HashMap<usize, usize>>>,
    txs: Vec<SyncSender<Option<(usize, Vec<T>)>>>,
    filenames: Vec<PathBuf>
}

impl<'a, T: 'static + Copy + Send> Shard<'a, T> {
    pub fn new(path: &'a Path, num_shards: usize, num_threads: usize) -> Shard<'a, T> {
        let mut txs = Vec::new();
        let mut threads = Vec::new();
        let mut filenames = Vec::new();

        create_dir(path).unwrap();

        // Create num_threads to service num_shards
        // data for shard i is handled by thread i % num_threads
        for p in (0..num_threads).into_iter() {
            // Create communication channel
            let (tx, rx) = sync_channel::<Option<(usize, Vec<T>)>>(10);
            txs.push(tx);

            // Create output
            let writer_idxs: Vec<usize> = (0..num_shards)
                                              .filter(|i| i % num_threads == p)
                                              .collect();

            let mut writers = HashMap::new();

            for idx in writer_idxs {
                let filename = format!("shard_file_{}.bin", idx);
                let wpath = path.join(Path::new(&filename));
                filenames.push(wpath.clone());
                let w = PodWriter::new(wpath);
                writers.insert(idx, w);
            }

            let thread = ShardWriterThread::spawn(rx, writers, num_shards);
            threads.push(thread);
        }

        Shard {
            path: path,
            num_shards: num_shards,
            threads: threads,
            txs: txs,
            filenames: filenames,
        }
    }
    /// Shutdown
    pub fn finish(self) -> Vec<PathBuf> {
        for tx in self.txs.iter() {
            tx.send(None).unwrap();
        }

        for t in self.threads {
            t.join().unwrap();
        }

        self.filenames
    }
}

pub struct ShardProducer<T: Copy + Send> {
    tx_channels: Vec<SyncSender<Option<(usize, Vec<T>)>>>,
    buffer: HashMap<usize, Vec<T>>,
    buf_size: usize,
    num_shards: usize,
}

impl<'a, T: Copy + Send> ShardProducer<T> {
    pub fn new(shard: &Shard<'a, T>) -> ShardProducer<T> {
        let mut new_txs = Vec::new();
        for t in shard.txs.iter() {
            new_txs.push(t.clone())
        }

        ShardProducer {
            tx_channels: new_txs,
            buffer: HashMap::new(),
            buf_size: 4096,
            num_shards: shard.num_shards,
        }
    }

    pub fn send(&mut self, raw_idx: usize, item: T) {
        let shard_idx = raw_idx % self.num_shards;
        let send = {
            let entry = self.buffer.entry(shard_idx).or_insert_with(|| Vec::new());
            entry.push(item);
            entry.len() == self.buf_size
        };

        if send {
            let send_entry = self.buffer.remove(&shard_idx).unwrap();
            let out_ch = self.tx_channels.get(shard_idx % self.tx_channels.len()).unwrap();
            out_ch.send(Some((shard_idx, send_entry))).unwrap();
        }
    }

    pub fn send_many(&mut self, items: Vec<(usize, T)>) {
        for (idx, item) in items {
            self.send(idx, item)
        }
    }

    pub fn finished(&mut self) {
        let idxs: Vec<usize> = self.buffer.keys().cloned().collect();

        for idx in idxs {
            let vec = self.buffer.remove(&idx).unwrap();
            let out_ch = self.tx_channels.get(idx % self.tx_channels.len()).unwrap();
            out_ch.send(Some((idx, vec))).unwrap();
        }
    }
}

impl<'a, T: Copy + Send> Drop for ShardProducer<T> {
    fn drop(&mut self) {
        self.finished();
    }
}


#[cfg(test)]
mod pod_tests {
    use std::path::Path;
    use std::fs::remove_dir_all;

    #[derive(Copy, Clone, Eq, PartialEq)]
    struct T1 {
        a: u64,
        b: u32,
        c: u16,
        d: u8,
    }

    fn write_data(filename: &str, d: &Vec<T1>) {
        // let filename = "test_pod.bin";
        let mut w = super::PodWriter::new(filename);
        w.write(d);
    }

    #[test]
    fn pod_round_trip() {
        let filename = "test_pod.bin";
        let mut write_vec = Vec::new();
        for i in 0..100 {
            let tt = T1 {
                a: i as u64,
                b: i as u32,
                c: i as u16,
                d: i as u8,
            };
            write_vec.push(tt);
        }
        write_data(filename, &write_vec);


        let mut reader = super::PodReader::new(filename);
        let read_vec = reader.read();

        assert!(write_vec == read_vec)
    }

    #[test]
    fn sharded_writing() {
        let n = 100000;
        let shards = 16;

        let test_path = Path::new("shard_test_path");
        if test_path.exists() {
            remove_dir_all(test_path).unwrap();
        }
        let s = super::Shard::new(test_path, shards, 2);

        let mut w = super::ShardProducer::new(&s);

        for i in 0..n {
            let tt = T1 {
                a: i as u64,
                b: i as u32,
                c: i as u16,
                d: i as u8,
            };
            w.send(i % shards, tt);
        }

        w.finished();
        s.finish();

        let shard_reader = super::ShardReader::new(test_path);
        let mut s1 = shard_reader.get_reader(0);
        let outs: Vec<T1> = s1.read();

        assert!(shard_reader.num_shards() == shards);
        assert!(outs.len() == n / shards);

        for i in outs {
            assert!(i.a % shards as u64 == 0)
        }
    }
}

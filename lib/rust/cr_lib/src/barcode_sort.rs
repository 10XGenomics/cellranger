use anyhow::Result;
use barcode::{Barcode, BarcodeConstruct, BarcodeConstructMetric, HasBarcode, Whitelist};
use cr_bam::constants::{ALN_BC_DISK_CHUNK_SZ, ALN_BC_ITEM_BUFFER_SZ, ALN_BC_SEND_BUFFER_SZ};
use cr_types::rna_read::{RnaChunk, RnaProcessor, RnaRead};
use fastq_set::read_pair::ReadPair;
use fastq_set::{FastqProcessor, ProcessResult};
use metric::{CountMetric, Metric};
use serde::Serialize;
use shardio::{ShardSender, ShardWriter, SortKey};
use std;
use std::borrow::Cow;
use std::marker::PhantomData;
use std::path::Path;

// 1GiB if above holds
const BG_FASTA_READAHEAD: usize = 1024;

pub trait ReadVisitor {
    type ReadType;
    fn visit_processed_read(&mut self, read: &mut Self::ReadType) -> Result<()>;
    fn visit_unprocessed_read(&mut self, read: ReadPair, reason: String) -> Result<()>;
}

struct DefaultVisitor<T> {
    phantom: PhantomData<T>,
}
impl<T> ReadVisitor for DefaultVisitor<T> {
    type ReadType = T;
    fn visit_processed_read(&mut self, _: &mut Self::ReadType) -> Result<()> {
        Ok(())
    }
    fn visit_unprocessed_read(&mut self, _: ReadPair, _: String) -> Result<()> {
        Ok(())
    }
}

pub struct BarcodeSortWorkflow {
    sorter: BarcodeAwareSorter<RnaRead, BarcodeOrder>,
    processor: RnaProcessor,
}

impl BarcodeSortWorkflow {
    pub fn new(
        chunk: RnaChunk,
        valid_path: &Path,
        invalid_path: &Path,
        whitelist: BarcodeConstruct<Whitelist>,
    ) -> Result<Self> {
        let sorter = BarcodeAwareSorter::new(valid_path, invalid_path)?;
        Ok(BarcodeSortWorkflow {
            processor: RnaProcessor::new(chunk, whitelist),
            sorter,
        })
    }

    pub fn execute_workflow_with_visitor<V>(
        &mut self,
        max_iters: Option<usize>,
        visitor: &mut V,
    ) -> Result<()>
    where
        V: ReadVisitor<ReadType = RnaRead>,
    {
        let mut nreads = 0;
        for read_result in self
            .processor
            .iter_background(BG_FASTA_READAHEAD)?
            .take(max_iters.unwrap_or(std::usize::MAX))
        {
            nreads += 1;
            match read_result? {
                ProcessResult::Processed(mut read) => {
                    visitor.visit_processed_read(&mut read)?;
                    self.sorter.process(read)?;
                }
                ProcessResult::Unprocessed { read, reason } => {
                    visitor.visit_unprocessed_read(read, reason)?;
                }
            }
        }
        println!("Processed {nreads} reads");
        self.sorter.finish()
    }

    pub fn num_valid_items(&self) -> i64 {
        self.sorter.valid_items
    }

    pub fn num_valid_segments(&self) -> BarcodeConstructMetric<CountMetric> {
        self.sorter.valid_items_in
    }

    pub fn num_invalid_items(&self) -> i64 {
        self.sorter.invalid_items
    }
}

pub struct BarcodeOrder;
impl<T> SortKey<T> for BarcodeOrder
where
    T: HasBarcode,
{
    type Key = Barcode;
    fn sort_key(v: &T) -> Cow<'_, Barcode> {
        Cow::Owned(v.barcode())
    }
}

struct BarcodeAwareSorter<T, Order = BarcodeOrder>
where
    T: 'static + HasBarcode + Send + Serialize,
    Order: SortKey<T>,
    Order::Key: 'static + Send + Serialize,
{
    // note: sender must appear before writers, because senders
    // must be dropped before writers.
    valid_sender: ShardSender<T, Order>,
    valid_writer: ShardWriter<T, Order>,
    invalid_sender: ShardSender<T, Order>,
    invalid_writer: ShardWriter<T, Order>,
    valid_items: i64,
    valid_items_in: BarcodeConstructMetric<CountMetric>,
    invalid_items: i64,
}

impl<T, Order> BarcodeAwareSorter<T, Order>
where
    T: 'static + HasBarcode + Send + Serialize,
    Order: SortKey<T>,
    <Order as SortKey<T>>::Key: 'static + Send + Ord + Serialize + Clone,
{
    fn new(valid_path: &Path, invalid_path: &Path) -> Result<Self> {
        let valid_writer: ShardWriter<T, Order> = ShardWriter::new(
            valid_path,
            ALN_BC_SEND_BUFFER_SZ,
            ALN_BC_DISK_CHUNK_SZ,
            ALN_BC_ITEM_BUFFER_SZ,
        )?;
        let valid_sender = valid_writer.get_sender();

        let invalid_writer: ShardWriter<T, Order> = ShardWriter::new(
            invalid_path,
            ALN_BC_SEND_BUFFER_SZ,
            ALN_BC_DISK_CHUNK_SZ,
            ALN_BC_ITEM_BUFFER_SZ,
        )?;
        let invalid_sender = invalid_writer.get_sender();

        Ok(BarcodeAwareSorter {
            valid_writer,
            valid_sender,
            invalid_writer,
            invalid_sender,
            valid_items: 0,
            valid_items_in: BarcodeConstructMetric::<CountMetric>::new(),
            invalid_items: 0,
        })
    }

    fn process(&mut self, read: T) -> Result<()> {
        self.valid_items_in.merge(
            read.segmented_barcode()
                .segments_valid()
                .map(CountMetric::from)
                .into(),
        );

        if read.barcode().is_valid() {
            self.valid_items += 1;
            self.valid_sender.send(read)?;
        } else {
            self.invalid_items += 1;
            self.invalid_sender.send(read)?;
        }
        Ok(())
    }

    fn finish(&mut self) -> Result<()> {
        self.valid_sender.finished()?;
        self.invalid_sender.finished()?;
        self.valid_writer.finish()?;
        self.invalid_writer.finish()?;
        Ok(())
    }
}

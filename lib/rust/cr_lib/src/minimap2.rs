#![deny(missing_docs)]

use barcode::whitelist::ReqStrand;
use minimap2::ffi::{MM_F_OUT_MD, MM_F_SPLICE_FOR, MM_F_SPLICE_REV};
use minimap2::{Aligner, Built, Preset};
use rust_htslib::bam::{Header, HeaderView, Record};
use std::path::Path;
use std::sync::Arc;

pub struct MinimapReference {
    pub aligner: MinimapAligner,
    pub header: Header,
}

impl MinimapReference {
    pub fn load(preset: Preset, strandedness: ReqStrand, index: &Path, tx_bed: &Path) -> Self {
        let mut aligner = Aligner::builder().preset(preset);

        match strandedness {
            ReqStrand::Forward => aligner.mapopt.flag &= !MM_F_SPLICE_REV as i64,
            ReqStrand::Reverse => aligner.mapopt.flag &= !MM_F_SPLICE_FOR as i64,
        }

        // Return the MD FLAG
        aligner.mapopt.flag |= MM_F_OUT_MD as i64;

        let aligner = aligner
            .with_cigar()
            .with_sam_out()
            .with_index(index, None)
            .unwrap();

        aligner.read_junction_lr(tx_bed.to_str().unwrap()).unwrap();

        let mut header = Header::new();
        aligner.populate_header(&mut header);

        Self {
            aligner: MinimapAligner {
                inner: Arc::new(aligner),
                header_view: HeaderView::from_header(&header),
            },
            header,
        }
    }

    pub fn get_aligner(&self) -> MinimapAligner {
        self.aligner.clone()
    }
}

#[derive(Clone)]
pub struct MinimapAligner {
    inner: Arc<Aligner<Built>>,
    header_view: HeaderView,
}

// Manual implementation of Send marker trait is safe because header_view is only
// read and not modified and Aligner within the inner field already implements Send.
unsafe impl Send for MinimapAligner {}

impl MinimapAligner {
    pub fn align_read(&mut self, name: &[u8], read: &[u8], qual: &[u8]) -> Vec<Record> {
        // Minimap2 will throw an error on empty reads - so just construct an empty record.
        if read.is_empty() {
            // Make an unmapped record and return it
            let rec = new_record(name, read, qual);
            return vec![rec];
        }

        self.inner
            .map_to_sam(read, Some(qual), Some(name), &self.header_view, None, None)
            .unwrap()
    }
}

pub fn new_record(name: &[u8], seq: &[u8], qual: &[u8]) -> Record {
    let mut record = Record::new();
    record.set(name, None, seq, qual);
    record.set_tid(-1);
    record.set_pos(-1);
    record.set_unmapped();
    record.set_mtid(-1);
    record.set_mpos(-1);
    record
}

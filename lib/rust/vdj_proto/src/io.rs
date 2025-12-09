//! Reading/Writing the vdj contig proto file which is used for VDJ aggr. This file contains
//! metadata about the VDJ analysis, the VDJ reference and a list of annotated contigs.
//!
//! The messages are defined in `types.proto`. The file format is as follows
//! ```text
//!
//! +-------------------------+
//! | Number of messages(k)   |
//! |       [4 bytes]         |
//! +-------------------------+
//!
//! +------------+----------------------------+
//! | Length     |      Message 0             |
//! | [4 bytes]  |      [Length Bytes]        |
//! +------------+----------------------------+
//! ...
//! ...
//! +------------+----------------------------+
//! | Length     |      Message k-1           |
//! | [4 bytes]  |      [Length Bytes]        |
//! +------------+----------------------------+
//! ```
//! where each `Message` is of type `VdjProtoMessage`, which can be one of metadata, reference
//! or an annotated contig. For convenience and ease of access, the metadata and the reference are
//! written as the first two messages followed by a list of contigs.
//!
//! This follows the recommendation here: https://developers.google.com/protocol-buffers/docs/techniques#streaming
//!
//! The number of messages and the length are written in Big endian.
#![expect(missing_docs)]

use crate::types::vdj_proto_message::MessageContent;
use crate::types::{BarcodeData, MetricsSummary, VdjMetadata, VdjProtoMessage, VdjReferenceRaw};
use anyhow::{Context, Result};
use byteorder::{BigEndian, ReadBytesExt, WriteBytesExt};
use enclone_proto::proto_io::{ProtoReader, ProtoWriter};
use serde_json::{Map, Value};
use std::convert::TryInto;
use std::fs::File;
use std::io::{BufReader, BufWriter, Seek};
use std::path::Path;

pub const PROTOBUF_VERSION: &str = "1.1.0";

/// Helper struct for writing vdj contig proto file.
///
/// Note: We use the `File` interface instead of `Write` because we need to seek to the
/// beginning of the file and write the number of messages on finish/drop.
pub struct VdjProtoWriter {
    // It's an option so that we can take ownership of inner
    // value on finish() with a mutable reference
    writer: Option<ProtoWriter<BufWriter<File>>>,
    num_messages: u32,
}

impl VdjProtoWriter {
    pub fn new(
        file_name: &Path,
        metadata: VdjMetadata,
        reference: VdjReferenceRaw,
        metrics: MetricsSummary,
    ) -> Result<Self> {
        let mut buf_writer = BufWriter::new(
            File::create(file_name)
                .with_context(|| format!("While creating file: {}", file_name.display()))?,
        );

        // We will update this at the end
        buf_writer.write_u32::<BigEndian>(0u32)?;

        let mut writer = ProtoWriter::with_writer(buf_writer);

        // IMPORTANT: Any changes made here need to be evaluated in sync with VdjProtoReader.
        // Some APIs in the reader explicitly assumes that these two messages exist always
        writer.encode_and_write(VdjProtoMessage::from(metadata))?;
        writer.encode_and_write(VdjProtoMessage::from(reference))?;
        writer.encode_and_write(VdjProtoMessage::from(metrics))?;

        Ok(VdjProtoWriter {
            writer: Some(writer),
            num_messages: 3, // Metadata, Reference and Metrics
        })
    }

    // DO NOT CALL THIS DIRECTLY. USE finish()
    fn _finish(&mut self) -> Result<()> {
        if let Some(w) = self.writer.take() {
            let mut writer = w.finish();
            writer.rewind()?;
            writer.write_u32::<BigEndian>(self.num_messages)?;
        }
        Ok(())
    }

    /// Finish writing the proto file.
    ///
    /// We need to seek to the beginning of the file and write the number of messages.
    /// Call this function to explicitly catch any errors as a result. On dropping we attempt
    /// to do it ignoring any errors
    pub fn finish(mut self) -> Result<()> {
        self._finish()
    }

    /// Write a contig annotation to the proto file
    pub fn write_annotation(
        &mut self,
        annotation: vdj_ann::annotate::ContigAnnotation,
    ) -> Result<()> {
        self.writer
            .as_mut()
            .unwrap() // Guaranteed to exist by construction
            .encode_and_write(VdjProtoMessage::from(annotation))?;
        self.num_messages += 1;
        Ok(())
    }

    /// Write a barcodedata to the proto file
    pub fn write_barcode_data(
        &mut self,
        barcode_brief: vdj_asm_utils::barcode_data::BarcodeDataBrief,
    ) -> Result<()> {
        self.writer
            .as_mut()
            .unwrap() // Guaranteed to exist by construction
            .encode_and_write(VdjProtoMessage::from(barcode_brief))?;
        self.num_messages += 1;
        Ok(())
    }
}

impl Drop for VdjProtoWriter {
    fn drop(&mut self) {
        let _ = self._finish();
    }
}

pub struct VdjProtoReader {
    reader: ProtoReader<BufReader<File>>,
    num_messages: u32,
    messages_read: u32,
}

impl VdjProtoReader {
    pub fn new(file_name: &Path) -> Result<VdjProtoReader> {
        let mut buf_reader =
            BufReader::new(File::open(file_name).with_context(|| {
                format!("While opening file for reading: {}", file_name.display())
            })?);
        let num_messages = buf_reader.read_u32::<BigEndian>()?;
        Ok(VdjProtoReader {
            reader: ProtoReader::from_reader(buf_reader),
            num_messages,
            messages_read: 0,
        })
    }

    pub fn len(&self) -> usize {
        self.num_messages as usize
    }
    pub fn is_empty(&self) -> bool {
        self.num_messages == 0
    }

    pub fn read_metadata(file_name: &Path) -> Result<VdjMetadata> {
        let reader = VdjProtoReader::new(file_name)?;
        for message in reader {
            if let MessageContent::Metadata(m) = message?.content() {
                return Ok(m);
            }
        }
        unreachable!("Metadata should always be present in the proto file")
    }

    pub fn read_reference(file_name: &Path) -> Result<VdjReferenceRaw> {
        let reader = VdjProtoReader::new(file_name)?;
        for message in reader {
            if let MessageContent::Reference(r) = message?.content() {
                return Ok(r);
            }
        }
        unreachable!("Reference should always be present in the proto file")
    }

    fn _read_metrics(file_name: &Path) -> Result<MetricsSummary> {
        let reader = VdjProtoReader::new(file_name)?;
        for message in reader {
            if let MessageContent::Metrics(r) = message?.content() {
                return Ok(r);
            }
        }
        unreachable!("Metrics should always be present in the proto file")
    }

    pub fn read_metrics(file_name: &Path) -> Result<Map<String, Value>> {
        match serde_json::from_str(&Self::_read_metrics(file_name)?.raw_json)? {
            Value::Object(o) => Ok(o),
            _ => unreachable!("Metrics summary JSON should be a map"),
        }
    }

    /// Returns an iterator over contig annotations. Everything is wrapped in a result to bubble up
    /// any IO errors
    pub fn read_annotations(
        file_name: &Path,
    ) -> Result<impl Iterator<Item = Result<vdj_ann::annotate::ContigAnnotation>> + use<>> {
        Ok(VdjProtoReader::new(file_name)?.filter_map(|message| {
            // We have a Result<VdjProtoMessage>
            match message {
                // We want to filter all except the annotation
                Ok(m) => match m.content() {
                    MessageContent::Annotation(ann) => Some(Ok(ann.into())),
                    _ => None,
                },
                // We want to bubble up the Error
                Err(e) => Some(Err(e)),
            }
        }))
    }

    pub fn read_barcode_data(
        file_name: &Path,
    ) -> Result<Option<impl Iterator<Item = Result<BarcodeData>> + use<>>> {
        if Self::read_metadata(file_name)?.protobuf_version.is_empty() {
            Ok(None)
        } else {
            let iterator = VdjProtoReader::new(file_name)?.filter_map(|message| match message {
                Ok(m) => match m.content() {
                    MessageContent::BarcodeData(bc) => Some(Ok(bc)),
                    _ => None,
                },
                Err(e) => Some(Err(e)),
            });
            Ok(Some(iterator))
        }
    }
}

impl Iterator for VdjProtoReader {
    type Item = Result<VdjProtoMessage>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.messages_read < self.num_messages {
            self.messages_read += 1;
            Some(self.reader.read_and_decode().map_err(anyhow::Error::from))
        } else {
            None
        }
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        let remaining = if self.messages_read < self.num_messages {
            let remaining = (self.num_messages - self.messages_read).try_into();
            if let Ok(remaining) = remaining {
                remaining
            } else {
                return (0, None);
            }
        } else {
            0
        };
        (remaining, Some(remaining))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::io::PROTOBUF_VERSION;
    use itertools::Itertools;
    use proptest::prelude::*;

    // Create an arbitrary meaningless VdjMetadata for testing purposes
    fn arbitrary_metadata() -> impl Strategy<Value = VdjMetadata> {
        (
            "[a-z0-9]*",
            "[a-z0-9-.]*",
            0..2i32,
            0..8u32,
            1..100_000u32,
            any::<String>(),
            any::<String>(),
            any::<String>(),
            PROTOBUF_VERSION,
            0..2i32,
        )
            .prop_map(
                |(
                    reference_fasta_hash,
                    pipeline_version,
                    receptor,
                    num_gem_wells,
                    number_of_cells,
                    sample_id,
                    sample_desc,
                    multi_config_sha,
                    protobuf_version,
                    multiplexing_method,
                )| VdjMetadata {
                    reference_fasta_hash,
                    pipeline_version,
                    receptor,
                    gem_wells: (1..=num_gem_wells).collect(),
                    number_of_cells,
                    sample_id,
                    sample_desc,
                    multi_config_sha,
                    protobuf_version,
                    multiplexing_method,
                },
            )
            .boxed()
    }

    // Create an arbitrary meaningless VdjReference for testing purposes
    fn arbitrary_reference() -> impl Strategy<Value = VdjReferenceRaw> {
        (any::<String>(), any::<String>())
            .prop_map(|(regions, ref_json)| VdjReferenceRaw { regions, ref_json })
            .boxed()
    }

    proptest! {
        #[test]
        fn test_simple_roundtrip(
            metadata in arbitrary_metadata(),
            reference in arbitrary_reference(),
            finish in any::<bool>(),
        ) {
            let file = tempfile::NamedTempFile::new().unwrap().into_temp_path();
            let metrics = MetricsSummary {raw_json: "{}".into()};
            let writer = VdjProtoWriter::new(
                &file,
                metadata.clone(),
                reference.clone(),
                metrics.clone(),
            ).unwrap();
            if finish {
                writer.finish().unwrap();
            } else {
                drop(writer);
            }

            let reader = VdjProtoReader::new(&file).unwrap();
            assert_eq!(reader.len(), 3);
            let messages: Vec<_> = reader.try_collect().unwrap();
            assert_eq!(messages, vec![
                metadata.clone().into(),
                reference.clone().into(),
                metrics.clone().into(),
            ]);

            assert_eq!(VdjProtoReader::read_metadata(&file).unwrap(), metadata);
            assert_eq!(VdjProtoReader::read_reference(&file).unwrap(), reference);
            assert_eq!(VdjProtoReader::_read_metrics(&file).unwrap(), metrics);
        }
    }

    #[test]
    fn test_write_annotations() -> Result<()> {
        let metadata = VdjMetadata {
            reference_fasta_hash: "abcd".into(),
            pipeline_version: "4.1".into(),
            receptor: 0,
            gem_wells: vec![1],
            number_of_cells: 100,
            sample_id: "100_100".into(),
            sample_desc: "Human PBMC".into(),
            multi_config_sha: "e96907faf862b9269bdd4649597778d2c44954dfa877095ceb934b8ce043bbca"
                .into(),
            protobuf_version: String::new(),
            multiplexing_method: 0,
        };

        let reference = VdjReferenceRaw {
            regions: ">fake\nACGT".into(),
            ref_json: "{}".into(),
        };

        let file = tempfile::NamedTempFile::new()?.into_temp_path();
        let metrics = MetricsSummary {
            raw_json: r#"{"vdj_filtered_bcs": 100}"#.into(),
        };
        let mut writer =
            VdjProtoWriter::new(&file, metadata.clone(), reference.clone(), metrics.clone())?;
        let anns: Vec<vdj_ann::annotate::ContigAnnotation> =
            serde_json::from_reader(File::open("../vdj_asm_asm/test.json")?)?;

        let num_messages = 3 + anns.len();
        for ann in &anns {
            writer.write_annotation(ann.clone())?;
        }
        writer.finish()?;

        {
            let reader = VdjProtoReader::new(&file)?;
            assert_eq!(reader.len(), num_messages);
            let anns_read: Vec<_> = VdjProtoReader::read_annotations(&file)?.try_collect()?;
            assert_eq!(anns_read, anns);
        }

        assert_eq!(VdjProtoReader::read_metadata(&file)?, metadata);
        assert_eq!(VdjProtoReader::read_reference(&file)?, reference);
        assert_eq!(VdjProtoReader::_read_metrics(&file)?, metrics);

        Ok(())
    }

    #[test]
    fn test_cellranger7_compatibility() -> Result<()> {
        let file = Path::new("test/cellranger7-1_contig_info.pb");
        let metadata = VdjProtoReader::read_metadata(file)?;
        assert!(metadata.protobuf_version.is_empty());
        assert!(metadata.multiplexing_method == 0);

        let barcode_data = VdjProtoReader::read_barcode_data(file)?;
        assert!(barcode_data.is_none());

        Ok(())
    }
}

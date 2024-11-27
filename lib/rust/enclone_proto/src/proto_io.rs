// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

//!
//! Code for reading and writing custom proto files
//!
//! Protobuf format by default is not suitable for streaming. We need to write our own logic to
//! track where a message begins and ends. The documentation recommends using a length delimited
//! format (https://developers.google.com/protocol-buffers/docs/techniques#streaming).
//!
//! ## Message Format
//! We use the following format to write a single protobuf message:
//! ```text
//! +------------+----------------------------+
//! | Length     |      Message               |
//! | [4 bytes]  |      [Length Bytes]        |
//! +------------+----------------------------+
//! ```
//! In the above diagram,
//! - `Length` is an unsigned 32 bit integer stored in **Big endian** order.
//! - If there are multiple messages, they are stored consecutively following the same format.

use crate::types::{Clonotype, EncloneOutputs};
use byteorder::{BigEndian, ReadBytesExt, WriteBytesExt};
use prost::Message;
use std::fs::File;
use std::io::{BufReader, BufWriter, Read, Write};
use std::path::Path;

const BUFFER_CAPACITY: usize = 1_000_000;

#[derive(Debug, thiserror::Error)]
pub enum Error {
    #[error("Expected to get {expected} bytes from the reader. Got {got} bytes!")]
    Truncated { expected: usize, got: usize },

    #[error(transparent)]
    Io(#[from] std::io::Error),

    #[error(transparent)]
    ProtoDecode(#[from] prost::DecodeError),

    #[error(transparent)]
    ProtoEncode(#[from] prost::EncodeError),
}

/// A helper struct to write a length delimited protobuf encoded message into the inner `writer`.
pub struct ProtoWriter<W: Write> {
    // Buffer space for storing the encoded message
    encode_buffer: Vec<u8>,
    writer: W,
}

impl<W: Write> ProtoWriter<W> {
    pub fn with_writer(writer: W) -> Self {
        ProtoWriter {
            encode_buffer: Vec::with_capacity(BUFFER_CAPACITY),
            writer,
        }
    }
    /// Writes a single message in length delimited format. Returns the total number
    /// of bytes written
    pub fn encode_and_write<M>(&mut self, message: M) -> Result<usize, Error>
    where
        M: Message,
    {
        // Write the message length as a big endian unsigned 32 bit integer
        let encoded_len = message.encoded_len();
        self.writer.write_u32::<BigEndian>(encoded_len as u32)?;

        // Encode the message in protobuf format and write it to the underlying writer
        message.encode(&mut self.encode_buffer)?;
        self.writer.write_all(&self.encode_buffer)?;
        self.encode_buffer.clear();
        Ok(encoded_len + 4) // +4 because of the u32
    }

    /// Consume self and return the inner writer
    pub fn finish(self) -> W {
        self.writer
    }
}

/// A helper struct to read a length delimited protobuf encoded message from the inner `reader`.
pub struct ProtoReader<R: Read> {
    decode_buffer: Vec<u8>,
    reader: R,
}

impl<R: Read> ProtoReader<R> {
    pub fn from_reader(reader: R) -> Self {
        ProtoReader {
            decode_buffer: Vec::with_capacity(BUFFER_CAPACITY),
            reader,
        }
    }
    // Clear the decode_buffer and fill it with `num_bytes` bytes from the reader
    fn read_exact(&mut self, num_bytes: usize) -> Result<(), Error> {
        self.decode_buffer.clear();
        self.reader
            .by_ref()
            .take(num_bytes as u64)
            .read_to_end(&mut self.decode_buffer)?;
        // If we did not get num_bytes bytes, return an error
        if self.decode_buffer.len() != num_bytes {
            return Err(Error::Truncated {
                expected: num_bytes,
                got: self.decode_buffer.len(),
            });
        }
        Ok(())
    }
    // Skip a message from the underlying reader
    pub fn skip(&mut self) -> Result<(), Error> {
        self.read_exact(4)?;
        let decoded_len = self.decode_buffer.as_slice().read_u32::<BigEndian>()?;
        self.read_exact(decoded_len as usize)
    }
    pub fn read_and_decode<M>(&mut self) -> Result<M, Error>
    where
        M: Message + Default,
    {
        // Attempt to take 4 bytes from the buffer for length
        self.read_exact(4)?;
        // decode the 4 bytes as a big endian u32
        let decoded_len = self.decode_buffer.as_slice().read_u32::<BigEndian>()?;

        // Decode the message
        self.read_exact(decoded_len as usize)?;
        let decoded_message = M::decode(self.decode_buffer.as_slice())?;
        Ok(decoded_message)
    }
}

/// The enclone outputs are stored in the protobuf file as follows:
/// ```text
/// +------------+----------------------------+
/// | Length     |      Version String        |....
/// | [4 bytes]  |      [Length Bytes]        |....
/// +------------+----------------------------+
///
/// +------------+----------------------------+
/// | Length     |      Metadata (aggr)       |....
/// | [4 bytes]  |      [Length Bytes]        |....
/// +------------+----------------------------+
///
/// +------------+----------------------------+------------+------------------------+
/// | Length     |      Universal Reference   | Length     |      Donor Reference   |....
/// | [4 bytes]  |      [Length Bytes]        | [4 bytes]  |      [Length Bytes]    |....
/// +------------+----------------------------+------------+------------------------+
///
/// +------------+-------------------------------+------------+-------------------------------+
/// | Length     |      Number of clonotypes (N) | Length     |      Clonotype 0              |....
/// | [4 bytes]  |      [Length Bytes]           | [4 bytes]  |      [Length Bytes]           |....
/// +------------+-------------------------------+------------+-------------------------------+
///
///     +------------+----------------------------------------+
/// ... | Length     |      Clonotype N-1                     |
/// ... | [4 bytes]  |      [Length Bytes]                    |
///     +------------+----------------------------------------+
/// ```
/// The newlines are only showed for illustration
pub fn write_proto(enclone_outputs: EncloneOutputs, path: impl AsRef<Path>) -> Result<(), Error> {
    let writer = BufWriter::new(File::create(path)?);
    let mut proto_writer = ProtoWriter::with_writer(writer);

    // Write the version
    proto_writer.encode_and_write(enclone_outputs.version)?;
    // Write the metadata
    proto_writer.encode_and_write(enclone_outputs.metadata)?;
    // Write the universal reference
    proto_writer.encode_and_write(enclone_outputs.universal_reference)?;
    // Write the donor reference
    proto_writer.encode_and_write(enclone_outputs.donor_reference)?;
    // Write the number of clonotypes. Not bothering to write this raw
    proto_writer.encode_and_write(enclone_outputs.clonotypes.len() as u32)?;
    for cl in enclone_outputs.clonotypes {
        proto_writer.encode_and_write(cl)?;
    }
    Ok(())
}

/// A read that mirrors the write above. The fields until the list of clonotypes are read here.
/// The clonotypes are assigned an empty vector.
pub fn read_proto_until_clonotypes(
    path: impl AsRef<Path>,
) -> Result<(EncloneOutputs, ProtoReader<impl Read>), Error> {
    let reader = BufReader::new(File::open(path)?);
    let mut proto_reader = ProtoReader::from_reader(reader);

    // Read the version
    let version = proto_reader.read_and_decode()?;
    // Read the metadata
    let metadata = proto_reader.read_and_decode()?;
    // Read the universal reference
    let universal_reference = proto_reader.read_and_decode()?;
    // Read the donor reference
    let donor_reference = proto_reader.read_and_decode()?;
    // Number of clonotypes
    let num_clonotypes: u32 = proto_reader.read_and_decode()?;

    Ok((
        EncloneOutputs {
            version,
            metadata,
            universal_reference,
            donor_reference,
            num_clonotypes,
            clonotypes: Vec::new(),
        },
        proto_reader,
    ))
}

/// A read that mirrors the write above. It is possible to stream through the
/// clonotypes instead of loading everything into memory.
pub fn read_proto(path: impl AsRef<Path>) -> Result<EncloneOutputs, Error> {
    let (mut output, mut proto_reader) = read_proto_until_clonotypes(path)?;
    let mut clonotypes = Vec::new();
    for _ in 0..output.num_clonotypes {
        clonotypes.push(proto_reader.read_and_decode()?);
    }
    output.clonotypes = clonotypes;
    Ok(output)
}

/// Iterator over clonotypes
pub struct ClonotypeIter<R: Read> {
    index: u32,
    num_clonotypes: u32,
    proto_reader: ProtoReader<R>,
}

impl ClonotypeIter<BufReader<File>> {
    pub fn from_file(file: impl AsRef<Path>) -> Result<Self, Error> {
        ClonotypeIter::from_reader(BufReader::new(File::open(file)?))
    }
}

impl<R: Read> ClonotypeIter<R> {
    pub fn from_reader(reader: R) -> Result<Self, Error> {
        let mut proto_reader = ProtoReader::from_reader(reader);
        // Skip version, metadata, universal reference, donor reference
        for _ in 0..4 {
            proto_reader.skip()?;
        }
        let num_clonotypes: u32 = proto_reader.read_and_decode()?;
        Ok(ClonotypeIter {
            index: 0,
            num_clonotypes,
            proto_reader,
        })
    }
}

impl<R: Read> Iterator for ClonotypeIter<R> {
    type Item = Clonotype;

    fn next(&mut self) -> Option<Self::Item> {
        if self.index < self.num_clonotypes {
            let cl = match self.proto_reader.read_and_decode() {
                Ok(c) => c,
                Err(e) => panic!("Failed to decode clonotype due to {e}"),
            };
            self.index += 1;
            Some(cl)
        } else {
            None
        }
    }
}

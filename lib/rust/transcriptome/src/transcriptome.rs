use crate::parse_gtf::{parse_gtf_line, validate_gtf_line};
use crate::Gene;
use anyhow::{anyhow, bail, Context, Result};
use bio::io::fasta::IndexedReader;
use bio_types::strand::ReqStrand;
use flate2::read::MultiGzDecoder;
use std::collections::hash_map::Entry;
use std::collections::{BTreeMap, HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader, ErrorKind, Read, Seek};
use std::path::Path;

#[derive(Hash, Eq, PartialEq, Debug, Clone, Ord, PartialOrd, Copy)]
pub struct GeneIdx(pub u32);

#[derive(Hash, Eq, PartialEq, Debug, Clone, Ord, PartialOrd, Copy)]
pub struct TranscriptIdx(pub u32);

#[derive(Eq, PartialEq, Debug, Clone, Hash)]
pub struct TranscriptomeGene {
    pub idx: GeneIdx,
    pub id: String,
    pub name: String,
    pub properties: Vec<(String, String)>,
}

impl TranscriptomeGene {
    pub fn to_gene(&self) -> Gene {
        Gene {
            id: self.id.clone(),
            name: self.name.clone(),
        }
    }
}

#[derive(Debug)]
pub struct Transcript {
    pub idx: TranscriptIdx,
    pub gene_idx: GeneIdx,
    pub id: String,
    pub chrom: String,
    pub strand: ReqStrand,
    pub exons: Vec<Exon>,
    pub properties: Vec<(String, String)>,
}

impl Transcript {
    pub fn get_sequence<R: Read + Seek>(
        &self,
        fasta_reader: &mut IndexedReader<R>,
    ) -> Result<Vec<u8>> {
        crate::transcript_sequence::get_transcript_sequence(
            fasta_reader,
            &self.chrom,
            &self.exons,
            self.strand,
        )
    }

    pub fn len(&self) -> u64 {
        self.exons.iter().map(Exon::len).sum()
    }

    pub fn gc_content<R: Read + Seek>(&self, fasta_reader: &mut IndexedReader<R>) -> Result<f64> {
        let buf = self.get_sequence(fasta_reader)?;
        if buf.is_empty() {
            Ok(0.0)
        } else {
            let n = buf.len();
            let gc = buf
                .into_iter()
                .filter(|&c| matches!(c, b'c' | b'C' | b'g' | b'G'))
                .count();
            Ok((gc as f64) / (n as f64))
        }
    }

    pub fn is_empty(&self) -> bool {
        !self.exons.iter().any(|e| !e.is_empty())
    }

    pub fn start(&self) -> u64 {
        if self.exons.is_empty() {
            return 0;
        }
        self.exons[0].start
    }

    pub fn end(&self) -> u64 {
        if self.exons.is_empty() {
            return 0;
        }
        self.exons[self.exons.len() - 1].end
    }
}

pub struct Transcriptome {
    pub genes: Vec<TranscriptomeGene>,
    pub gene_id_to_idx: HashMap<String, GeneIdx>,
    pub transcripts: Vec<Transcript>,
    pub transcript_id_to_idx: HashMap<String, TranscriptIdx>,
    pub gene_to_transcripts: BTreeMap<GeneIdx, Vec<TranscriptIdx>>, // map from gene_id to list of transcript_ids
}

impl Transcriptome {
    /// Read a possibly-compressed GTF file from a reader.
    pub fn from_reader(mut reader: impl BufRead) -> Result<Transcriptome> {
        load_from_gtf_reader(&mut reader)
    }

    /// Read a possibly-compressed GTF file.
    pub fn from_gtf_path(path: &Path) -> Result<Transcriptome> {
        match File::open(path) {
            Ok(file) => Transcriptome::from_reader(BufReader::new(file)),
            Err(err) if err.kind() == ErrorKind::NotFound => {
                Transcriptome::from_reader(BufReader::new(MultiGzDecoder::new(
                    File::open(path.with_extension("gtf.gz"))
                        .with_context(|| path.display().to_string())?,
                )))
            }
            Err(err) => Err(err).with_context(|| path.display().to_string()),
        }
    }

    /// Read a possibly-compressed GTF file from a path to a reference transcriptome.
    pub fn from_reference_path(path: &Path) -> Result<Transcriptome> {
        Transcriptome::from_gtf_path(&path.join("genes/genes.gtf"))
    }

    pub fn dummy() -> Transcriptome {
        Transcriptome {
            genes: Vec::new(),
            gene_id_to_idx: HashMap::new(),
            transcripts: Vec::new(),
            transcript_id_to_idx: HashMap::new(),
            gene_to_transcripts: BTreeMap::new(),
        }
    }

    /// Check if there are any transcripts with no exons and return Err(<informative message>) if we find any.
    pub fn check_for_transcripts_with_no_exons(&self) -> Result<()> {
        let empty_transcripts: String = self
            .transcripts
            .iter()
            .filter_map(|tx| {
                if tx.exons.is_empty() {
                    Some(format!(
                        "Transcript ID: {}, Gene Name: {}\n",
                        tx.id, &self.genes[tx.gene_idx.0 as usize].name
                    ))
                } else {
                    None
                }
            })
            .collect();
        if empty_transcripts.is_empty() {
            Ok(())
        } else {
            anyhow::bail!(
                "Transcripts with no exons were detected in the GTF.\n\
                 Please ensure there is at least one `exon` entry for each `transcript` in the GTF.\n\
                 The following transcripts are lacking an exon entry:\n{empty_transcripts}",
            );
        }
    }
}

#[derive(Hash, Eq, PartialEq, Debug, Clone, Ord, PartialOrd)]
pub struct Exon {
    pub start: u64,
    pub end: u64,
}

impl Exon {
    pub fn len(&self) -> u64 {
        self.end - self.start
    }
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
}

/// Read a `Transcriptome` from a GTF file. Handles `gene`, `transcript` and `exon` GTF entries.
fn load_from_gtf_reader(in_gtf: &mut dyn BufRead) -> Result<Transcriptome> {
    let mut transcripts: Vec<Transcript> = Vec::new();
    let mut transcript_id_to_idx: HashMap<String, TranscriptIdx> = HashMap::new();

    let mut genes = Vec::new();
    let mut gene_id_to_idx: HashMap<String, GeneIdx> = HashMap::new();

    // To keep track of genes that are added by parsing exons/transcripts but
    // later updated by parsing a gene
    let mut genes_not_from_file = HashSet::new();

    for (line_num, line) in in_gtf.lines().enumerate() {
        let line = line?;
        if line.starts_with('#') {
            continue;
        }
        let make_err =
            |msg: &str| anyhow!("Parsing GTF on line {line_num}: {msg}\nLine = '{line}'");

        let Ok((_, rec)) = parse_gtf_line(line.as_bytes()) else {
            // Since parsing failed, validation should fail.
            // Parsing errors are unreadable, so run full (slow) validation.
            // If for some reason our human-readable validator doesn't fail,
            // return a generic error message.
            bail!(make_err(
                &match validate_gtf_line(line.as_bytes()) {
                    Ok(()) =>
                        anyhow!("please check this line of your GTF file for formatting errors"),
                    Err(err) => err,
                }
                .to_string()
            ));
        };

        // start is 1-based, convert to 0-based
        let start = rec.start - 1;
        // end is 1-based, inclusive. same value as 0-based exclusive.
        let end = rec.end;

        let strand = match rec.strand {
            b"+" => ReqStrand::Forward,
            b"-" => ReqStrand::Reverse,
            _ => bail!(make_err("unknown strand")),
        };

        let chrom = std::str::from_utf8(rec.seqname)?;

        if rec.feature_type == b"gene" {
            let id = rec.get_attr("gene_id")?;
            let name = rec.get_attr("gene_name").unwrap_or_else(|_| id.clone());
            let idx = GeneIdx(genes.len() as u32);

            let gene = TranscriptomeGene {
                idx,
                id: id.clone(),
                name,
                properties: rec.all_attributes(),
            };
            // Make sure we don't let duplicates through
            match gene_id_to_idx.entry(id) {
                Entry::Occupied(entry) => {
                    if !genes_not_from_file.remove(entry.key()) {
                        bail!("Duplicate Gene ID found in GTF: {}", entry.key());
                    }
                    // If a `transcript` appears before a `gene`, we generate an entry for that gene
                    // but we should override that with the actual gene once we observe it.
                    let idx_int = entry.get().0 as usize;
                    genes[idx_int] = gene;
                }
                Entry::Vacant(entry) => {
                    entry.insert(idx);
                    genes.push(gene);
                }
            }
        } else if rec.feature_type == b"transcript" {
            let id = rec.get_attr("transcript_id")?;
            let gene_id = rec.get_attr("gene_id")?;

            match transcript_id_to_idx.get(&id) {
                Some(idx) => {
                    let tx = &mut transcripts[idx.0 as usize];

                    // we already made a place-holder transcript. If so, make sure things
                    // match, and update the attributes

                    // verify that the chrom and strand of the exon entries match the transcript entries.
                    anyhow::ensure!(
                        tx.chrom == chrom,
                        "Previous 'exon' entry for transcript_id {} had chrom == {}, but 'transcript' entry has chrom == {}",
                        id,
                        tx.chrom,
                        chrom
                    );
                    anyhow::ensure!(
                        tx.strand == strand,
                        "Previous 'exon' entry for transcript_id {} had strand == {}, but 'transcript' entry has strand == {}",
                        id,
                        tx.strand,
                        strand
                    );

                    // update the attribute of transcript entry.
                    tx.properties = rec.all_attributes();
                }
                None => {
                    // make a fresh transcript
                    let idx = TranscriptIdx(transcripts.len() as u32);

                    // handle a missing gene for this tx
                    let gene_idx = *gene_id_to_idx.entry(gene_id.clone()).or_insert_with(|| {
                        let gene_name = rec
                            .get_attr("gene_name")
                            .unwrap_or_else(|_| gene_id.clone());
                        let new_gene_idx = GeneIdx(genes.len() as u32);

                        let gene = TranscriptomeGene {
                            idx: new_gene_idx,
                            id: gene_id.clone(),
                            name: gene_name,
                            properties: vec![],
                        };
                        genes_not_from_file.insert(gene_id.clone());
                        genes.push(gene);
                        new_gene_idx
                    });

                    let transcript = Transcript {
                        idx,
                        id: id.clone(),
                        gene_idx,
                        chrom: chrom.to_string(),
                        strand,
                        exons: vec![],
                        properties: rec.all_attributes(),
                    };

                    transcript_id_to_idx.insert(id, idx);
                    transcripts.push(transcript);
                }
            }
        } else if rec.feature_type == b"exon" {
            let exon = Exon { start, end };
            let transcript_id = rec.get_attr("transcript_id")?;

            // the transcript hasn't been declared -- make a  dummy
            if transcript_id_to_idx.get(&transcript_id).is_none() {
                let gene_id = rec.get_attr("gene_id")?;
                let idx = TranscriptIdx(transcripts.len() as u32);

                // handle a missing gene for this tx
                if gene_id_to_idx.get(&gene_id).is_none() {
                    let gene_name = rec
                        .get_attr("gene_name")
                        .unwrap_or_else(|_| gene_id.clone());
                    let new_gene_idx = GeneIdx(genes.len() as u32);

                    let gene = TranscriptomeGene {
                        idx: new_gene_idx,
                        id: gene_id.clone(),
                        name: gene_name,
                        properties: vec![],
                    };
                    genes_not_from_file.insert(gene_id.clone());
                    gene_id_to_idx.insert(gene_id.clone(), new_gene_idx);
                    genes.push(gene);
                }

                let gene_idx = *gene_id_to_idx.get(&gene_id).expect("missing gene");

                let transcript = Transcript {
                    idx,
                    id: transcript_id.clone(),
                    gene_idx,
                    chrom: chrom.to_string(),
                    strand,
                    exons: vec![],
                    properties: vec![],
                };

                transcript_id_to_idx.insert(transcript_id.clone(), idx);
                transcripts.push(transcript);
            }

            let tx_idx = transcript_id_to_idx.get(&transcript_id).ok_or_else(|| {
                make_err(&format!(
                    "this row references transcript_id={transcript_id} but this \
                     transcript has no preceding 'transcript' row in the GTF"
                ))
            })?;

            transcripts[tx_idx.0 as usize].exons.push(exon);
        }
    }

    let mut gene_to_transcripts = BTreeMap::new();
    for tx in &mut transcripts {
        // Sort the exons so they're in coordinate order on the genome
        tx.exons.sort();

        // Tabulate the set of transcripts for each gene
        gene_to_transcripts
            .entry(tx.gene_idx)
            .or_insert_with(Vec::new)
            .push(tx.idx);
    }

    Ok(Transcriptome {
        genes,
        gene_id_to_idx,
        transcripts,
        transcript_id_to_idx,
        gene_to_transcripts,
    })
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::transcript_index::TranscriptIndex;
    use anyhow::Result;
    use martian_filetypes::tabular_file::TsvFile;

    // manual test that the GTF path matches the gene_index_tab created in the
    // legacy python code.
    #[test]
    fn test_tx_index() -> Result<()> {
        let f =
            File::open("test/gtf/GRCh38_small.gtf.gz").expect("Failed to open file for reading");
        let rdr = flate2::read::GzDecoder::new(f);
        let txome = Transcriptome::from_reader(BufReader::new(rdr))?;
        let old = TranscriptIndex::new(&TsvFile::from("test/gtf/GRCh38_small_gene_index.tsv"))?;
        let new = TranscriptIndex::from_transcriptome(&txome);
        assert_eq!(old, new);
        Ok(())
    }

    #[test]
    fn test_no_gene_name() -> Result<()> {
        let f = File::open("test/gtf/no_gene_names.gtf").expect("Failed to open file for reading");
        let _txome = Transcriptome::from_reader(BufReader::new(f))?;
        Ok(())
    }

    // test duplicates can't be parsed.
    #[test]
    fn test_no_dups() -> Result<()> {
        let gtf = "\
NC_000086.8	BestRefSeq	gene	168758039	168761913	.	-	.	gene_id \"G530011O06Rik\"
NC_000087.8	BestRefSeq	gene	90762409	90766319	.	-	.	gene_id \"G530011O06Rik\"; ";
        let txome = Transcriptome::from_reader(BufReader::new(gtf.as_bytes()));
        let expected = "Duplicate Gene ID found in GTF: G530011O06Rik";
        assert!(txome.is_err());
        if let Err(e) = txome {
            assert_eq!(format!("{e}"), expected);
        }
        // This should succeed
        let alt_gtf = "\
NC_000087.8	BestRefSeq	exon	90765690	90766319	.	-	.	gene_id \"G530011O06Rik\"; transcript_id \"NR_137283.1\";
NC_000086.8	BestRefSeq	gene	168758039	168761913	.	-	.	gene_id \"G530011O06Rik\";";
        let txome2 = Transcriptome::from_reader(BufReader::new(alt_gtf.as_bytes()));
        assert!(txome2.is_ok());
        Ok(())
    }

    #[test]
    fn exon_then_transcript() -> Result<()> {
        // Make sure an exon that appears before it's corresponding transcript is actually retained

        let gtf = r#"\MD043-011_TRB1	Genbank	gene	1	644	.	+	.	gbkey "Gene"; gene_name "gn_MD043-011_TRB1"; gene_biotype "protein_coding"; gene_id "gi_MD043-011_TRB1";
MD043-011_TRB1	Genbank	CDS	1	644	.	+	0	exon_number "1"; gbkey "CDS"; gene_name "gn_MD043-011_TRB1"; gene_id "gi_MD043-011_TRB1"; product "E6"; protein_id "pi_MD043-011_TRB1"; transcript_id "ti_MD043-011_TRB1";
MD043-011_TRB1	Genbank	exon	1	644	.	+	.	exon_number "1"; gbkey "CDS"; gene_name "gn_MD043-011_TRB1"; gene_id "gi_MD043-011_TRB1"; product "E6"; protein_id "pi_MD043-011_TRB1"; transcript_id "ti_MD043-011_TRB1";
MD043-011_TRB1	Genbank	transcript	1	644	.	+	.	exon_number "1"; gbkey "CDS"; gene_name "gn_MD043-011_TRB1"; gene_id "gi_MD043-011_TRB1"; product "E6"; protein_id "pi_MD043-011_TRB1"; transcript_id "ti_MD043-011_TRB1";"#;

        let txome = Transcriptome::from_reader(BufReader::new(gtf.as_bytes()))?;
        assert_eq!(txome.transcripts[0].exons.len(), 1);
        assert!(txome.transcripts[0]
            .properties
            .iter()
            .any(|(name, _)| name == "protein_id"));
        Ok(())
    }
}

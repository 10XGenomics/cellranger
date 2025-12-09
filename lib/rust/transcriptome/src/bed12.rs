use crate::transcriptome::{Exon, Transcript};
use itertools::Itertools;
use martian_derive::martian_filetype;
use martian_filetypes::table_config;
use martian_filetypes::tabular_file::{DelimitedFormat, TableConfig};
use serde::{Deserialize, Serialize, Serializer};

martian_filetype!(Bed12File, "bed");
table_config! { BedFormatNoHeader, b'\t', "bed", false, None }
/// A type alias for the BED12 file format, representing a tabular file with no header.
pub type Bed12Format = DelimitedFormat<Bed12Transcript, Bed12File, BedFormatNoHeader>;

fn comma_separated<S, T>(vec: &[T], serializer: S) -> Result<S::Ok, S::Error>
where
    S: Serializer,
    T: std::fmt::Display,
{
    let mut s = vec.iter().join(",");
    s.push(',');
    serializer.serialize_str(&s)
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct Bed12Transcript {
    pub chrom: String,
    pub start: u64,
    pub end: u64,
    pub name: String,
    pub score: u16,
    pub strand: String,
    pub thick_start: u64,
    pub thick_end: u64,
    pub item_rgb: String,
    pub block_count: usize,
    #[serde(serialize_with = "comma_separated")]
    pub block_sizes: Vec<u64>,
    #[serde(serialize_with = "comma_separated")]
    pub block_starts: Vec<u64>,
}

impl Bed12Transcript {
    pub fn from_tx(tx: &Transcript) -> Self {
        assert!(
            !tx.exons.is_empty(),
            "Transcript must have at least one exon"
        );

        Bed12Transcript {
            chrom: tx.chrom.clone(),
            start: tx.start(),
            end: tx.end(),
            name: format!("{}|{}|{}", tx.id, tx.transcript_type(), tx.gene_name()),
            score: 1000,
            strand: tx.strand.strand_symbol().to_owned(),
            thick_start: tx.thick_start(),
            thick_end: tx.thick_end(),
            item_rgb: match tx.transcript_type().as_str() {
                "protein_coding" => "0,128,255",
                _ => "196,196,196",
            }
            .to_string(),
            block_count: tx.exons.len(),
            block_sizes: tx.exons.iter().map(Exon::len).collect(),
            block_starts: tx.exons.iter().map(|e| e.start - tx.start()).collect(),
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::Transcriptome;
    use anyhow::Result;
    use std::io::BufReader;
    #[test]
    fn tx_noncoding_fwd() -> Result<()> {
        let gtf = r#"chr1	HAVANA	transcript	11869	14409	.	+	.	gene_id "ENSG00000290825"; gene_version "1"; transcript_id "ENST00000456328"; transcript_version "2"; gene_type "lncRNA"; gene_name "DDX11L2"; transcript_type "lncRNA"; transcript_name "DDX11L2-202"; level 2; transcript_support_level "1"; tag "basic"; tag "Ensembl_canonical"; havana_transcript "OTTHUMT00000362751.1";
chr1	HAVANA	exon	11869	12227	.	+	.	gene_id "ENSG00000290825"; gene_version "1"; transcript_id "ENST00000456328"; transcript_version "2"; gene_type "lncRNA"; gene_name "DDX11L2"; transcript_type "lncRNA"; transcript_name "DDX11L2-202"; exon_number 1; exon_id "ENSE00002234944"; exon_version "1"; level 2; transcript_support_level "1"; tag "basic"; tag "Ensembl_canonical"; havana_transcript "OTTHUMT00000362751.1";
chr1	HAVANA	exon	12613	12721	.	+	.	gene_id "ENSG00000290825"; gene_version "1"; transcript_id "ENST00000456328"; transcript_version "2"; gene_type "lncRNA"; gene_name "DDX11L2"; transcript_type "lncRNA"; transcript_name "DDX11L2-202"; exon_number 2; exon_id "ENSE00003582793"; exon_version "1"; level 2; transcript_support_level "1"; tag "basic"; tag "Ensembl_canonical"; havana_transcript "OTTHUMT00000362751.1";
chr1	HAVANA	exon	13221	14409	.	+	.	gene_id "ENSG00000290825"; gene_version "1"; transcript_id "ENST00000456328"; transcript_version "2"; gene_type "lncRNA"; gene_name "DDX11L2"; transcript_type "lncRNA"; transcript_name "DDX11L2-202"; exon_number 3; exon_id "ENSE00002312635"; exon_version "1"; level 2; transcript_support_level "1"; tag "basic"; tag "Ensembl_canonical"; havana_transcript "OTTHUMT00000362751.1";"#;

        let txome = Transcriptome::from_reader(BufReader::new(gtf.as_bytes()))?;
        let bed12: Vec<Bed12Transcript> = txome.convert_to_bed12().collect();
        assert_eq!(bed12.len(), 1);

        let expected = Bed12Transcript {
            chrom: "chr1".to_string(),
            start: 11868,
            end: 14409,
            name: "ENST00000456328|lncRNA|DDX11L2".to_string(),
            score: 1000,
            strand: "+".to_string(),
            thick_start: 11868,
            thick_end: 14409,
            item_rgb: "196,196,196".to_string(),
            block_count: 3,
            block_sizes: vec![359, 109, 1189],
            block_starts: vec![0, 744, 1352],
        };
        assert_eq!(expected, bed12[0]);
        Ok(())
    }

    #[test]
    fn tx_noncoding_rev() -> Result<()> {
        let gtf = r#"chr1	HAVANA	transcript	129081	133566	.	-	.	gene_id "ENSG00000238009"; gene_version "6"; transcript_id "ENST00000453576"; transcript_version "2"; gene_type "lncRNA"; gene_name "ENSG00000238009"; transcript_type "lncRNA"; transcript_name "ENST00000453576"; level 2; transcript_support_level "2"; havana_gene "OTTHUMG00000001096.2"; havana_transcript "OTTHUMT00000003689.1";
chr1	HAVANA	exon	133374	133566	.	-	.	gene_id "ENSG00000238009"; gene_version "6"; transcript_id "ENST00000453576"; transcript_version "2"; gene_type "lncRNA"; gene_name "ENSG00000238009"; transcript_type "lncRNA"; transcript_name "ENST00000453576"; exon_number 1; exon_id "ENSE00001737600"; exon_version "2"; level 2; transcript_support_level "2"; havana_gene "OTTHUMG00000001096.2"; havana_transcript "OTTHUMT00000003689.1";
chr1	HAVANA	exon	129081	129223	.	-	.	gene_id "ENSG00000238009"; gene_version "6"; transcript_id "ENST00000453576"; transcript_version "2"; gene_type "lncRNA"; gene_name "ENSG00000238009"; transcript_type "lncRNA"; transcript_name "ENST00000453576"; exon_number 2; exon_id "ENSE00001827073"; exon_version "1"; level 2; transcript_support_level "2"; havana_gene "OTTHUMG00000001096.2"; havana_transcript "OTTHUMT00000003689.1";"#;

        let txome = Transcriptome::from_reader(BufReader::new(gtf.as_bytes()))?;
        let bed12: Vec<Bed12Transcript> = txome.convert_to_bed12().collect();
        assert_eq!(bed12.len(), 1);

        let expected = Bed12Transcript {
            chrom: "chr1".to_string(),
            start: 129080,
            end: 133566,
            name: "ENST00000453576|lncRNA|ENSG00000238009".to_string(),
            score: 1000,
            strand: "-".to_string(),
            thick_start: 129080,
            thick_end: 133566,
            item_rgb: "196,196,196".to_string(),
            block_count: 2,
            block_sizes: vec![143, 193],
            block_starts: vec![0, 4293],
        };
        assert_eq!(expected, bed12[0]);
        Ok(())
    }

    #[test]
    fn tx_coding_rev() -> Result<()> {
        let gtf = r#"chr1	HAVANA	transcript	450740	451678	.	-	.	gene_id "ENSG00000284733"; gene_version "2"; transcript_id "ENST00000426406"; transcript_version "4"; gene_type "protein_coding"; gene_name "OR4F29"; transcript_type "protein_coding"; transcript_name "OR4F29-201"; level 2; protein_id "ENSP00000409316.1"; transcript_support_level "NA"; hgnc_id "HGNC:31275"; tag "basic"; tag "Ensembl_canonical"; tag "MANE_Select"; tag "appris_principal_1"; tag "CCDS"; ccdsid "CCDS72675.1"; havana_gene "OTTHUMG00000002860.3"; havana_transcript "OTTHUMT00000007999.3";
chr1	HAVANA	exon	450740	451678	.	-	.	gene_id "ENSG00000284733"; gene_version "2"; transcript_id "ENST00000426406"; transcript_version "4"; gene_type "protein_coding"; gene_name "OR4F29"; transcript_type "protein_coding"; transcript_name "OR4F29-201"; exon_number 1; exon_id "ENSE00003989331"; exon_version "1"; level 2; protein_id "ENSP00000409316.1"; transcript_support_level "NA"; hgnc_id "HGNC:31275"; tag "basic"; tag "Ensembl_canonical"; tag "MANE_Select"; tag "appris_principal_1"; tag "CCDS"; ccdsid "CCDS72675.1"; havana_gene "OTTHUMG00000002860.3"; havana_transcript "OTTHUMT00000007999.3";
chr1	HAVANA	CDS	450743	451678	.	-	0	gene_id "ENSG00000284733"; gene_version "2"; transcript_id "ENST00000426406"; transcript_version "4"; gene_type "protein_coding"; gene_name "OR4F29"; transcript_type "protein_coding"; transcript_name "OR4F29-201"; exon_number 1; exon_id "ENSE00003989331"; exon_version "1"; level 2; protein_id "ENSP00000409316.1"; transcript_support_level "NA"; hgnc_id "HGNC:31275"; tag "basic"; tag "Ensembl_canonical"; tag "MANE_Select"; tag "appris_principal_1"; tag "CCDS"; ccdsid "CCDS72675.1"; havana_gene "OTTHUMG00000002860.3"; havana_transcript "OTTHUMT00000007999.3";
chr1	HAVANA	start_codon	451676	451678	.	-	0	gene_id "ENSG00000284733"; gene_version "2"; transcript_id "ENST00000426406"; transcript_version "4"; gene_type "protein_coding"; gene_name "OR4F29"; transcript_type "protein_coding"; transcript_name "OR4F29-201"; exon_number 1; exon_id "ENSE00003989331"; exon_version "1"; level 2; protein_id "ENSP00000409316.1"; transcript_support_level "NA"; hgnc_id "HGNC:31275"; tag "basic"; tag "Ensembl_canonical"; tag "MANE_Select"; tag "appris_principal_1"; tag "CCDS"; ccdsid "CCDS72675.1"; havana_gene "OTTHUMG00000002860.3"; havana_transcript "OTTHUMT00000007999.3";
chr1	HAVANA	stop_codon	450740	450742	.	-	0	gene_id "ENSG00000284733"; gene_version "2"; transcript_id "ENST00000426406"; transcript_version "4"; gene_type "protein_coding"; gene_name "OR4F29"; transcript_type "protein_coding"; transcript_name "OR4F29-201"; exon_number 1; exon_id "ENSE00003989331"; exon_version "1"; level 2; protein_id "ENSP00000409316.1"; transcript_support_level "NA"; hgnc_id "HGNC:31275"; tag "basic"; tag "Ensembl_canonical"; tag "MANE_Select"; tag "appris_principal_1"; tag "CCDS"; ccdsid "CCDS72675.1"; havana_gene "OTTHUMG00000002860.3"; havana_transcript "OTTHUMT00000007999.3";
chr1	HAVANA	UTR	450740	450742	.	-	.	gene_id "ENSG00000284733"; gene_version "2"; transcript_id "ENST00000426406"; transcript_version "4"; gene_type "protein_coding"; gene_name "OR4F29"; transcript_type "protein_coding"; transcript_name "OR4F29-201"; exon_number 1; exon_id "ENSE00003989331"; exon_version "1"; level 2; protein_id "ENSP00000409316.1"; transcript_support_level "NA"; hgnc_id "HGNC:31275"; tag "basic"; tag "Ensembl_canonical"; tag "MANE_Select"; tag "appris_principal_1"; tag "CCDS"; ccdsid "CCDS72675.1"; havana_gene "OTTHUMG00000002860.3"; havana_transcript "OTTHUMT00000007999.3";"#;

        let txome = Transcriptome::from_reader(BufReader::new(gtf.as_bytes()))?;
        let bed12: Vec<Bed12Transcript> = txome.convert_to_bed12().collect();
        assert_eq!(bed12.len(), 1);

        let expected = Bed12Transcript {
            chrom: "chr1".to_string(),
            start: 450739,
            end: 451678,
            name: "ENST00000426406|protein_coding|OR4F29".to_string(),
            score: 1000,
            strand: "-".to_string(),
            thick_start: 450742,
            thick_end: 451678,
            item_rgb: "0,128,255".to_string(),
            block_count: 1,
            block_sizes: vec![939],
            block_starts: vec![0],
        };
        assert_eq!(expected, bed12[0]);
        Ok(())
    }
}

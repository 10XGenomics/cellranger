use anyhow::{anyhow, bail, ensure, Result};
use nom::branch::alt;
use nom::bytes::complete::{is_not, tag, take_until, take_while, take_while1};
use nom::character::complete::{char, digit1};
use nom::combinator::{all_consuming, map_res, opt};
use nom::error::{ErrorKind, ParseError};
use nom::sequence::{delimited, pair, separated_pair, terminated, tuple};
use nom::{self, IResult};
use smallvec::SmallVec;
use std::cmp::Ordering;
use std::num::{ParseFloatError, ParseIntError};
use std::str;
use std::str::FromStr;

/// A parsed GTF record. The string
/// fields are borrowed from a buffer
/// containing the GTF line.
#[derive(Debug, Clone, PartialEq)]
pub struct Record<'a> {
    pub seqname: &'a [u8],
    pub source: &'a [u8],
    pub feature_type: &'a [u8],
    pub start: u64,
    pub end: u64,
    pub score: Option<f64>,
    pub strand: &'a [u8],
    pub frame: &'a [u8],
    pub attributes: AttrVec<'a>,
}

impl<'a> Record<'a> {
    pub fn all_attributes(&self) -> Vec<(String, String)> {
        let mut vec: Vec<(String, String)> = Vec::new();

        for (k, v) in &self.attributes {
            let k = std::str::from_utf8(k).unwrap().to_string();
            let v = std::str::from_utf8(v).unwrap().to_string();
            vec.push((k, v));
        }

        vec
    }

    pub fn get_attr(&self, attribute: &str) -> Result<String> {
        for (k, v) in &self.attributes {
            if k == &attribute.as_bytes() {
                return Ok(std::str::from_utf8(v)?.to_string());
            }
        }

        bail!("attribute not found: {}", attribute)
    }
}

/// Parse one line of a GTF file into a `Record<'a>`. The
/// record will borrow slices from the input line.
pub fn parse_gtf_line(line: &[u8]) -> IResult<&[u8], Record<'_>> {
    let gtf_line = tuple((
        is_not("\t\r\n "), // seqname
        char('\t'),
        is_not("\t\r\n "), // source
        char('\t'),
        is_not("\t\r\n "), // feature
        char('\t'),
        parse_u64, // start
        char('\t'),
        parse_u64, // end
        char('\t'),
        score, // score
        char('\t'),
        is_not("\t\r\n "), // strand -- FIXME!!
        char('\t'),
        is_not("\t\r\n "), // frame -- FIXME 0,1,2,or '.'
        char('\t'),
        gtf_attributes, // attributes
    ));

    let v = map_res(gtf_line, convert_to_record);
    all_consuming(v)(line)
}

impl<'a> PartialOrd for Record<'a> {
    fn partial_cmp(&self, other: &Record<'_>) -> Option<Ordering> {
        let r = self.seqname.cmp(other.seqname);
        if r != Ordering::Equal {
            return Some(r);
        }

        let r = self.start.cmp(&other.start);
        if r != Ordering::Equal {
            return Some(r);
        }

        let r = self.end.cmp(&other.end);
        if r != Ordering::Equal {
            return Some(r);
        }

        let r = self.feature_type.cmp(other.feature_type);
        if r != Ordering::Equal {
            return Some(r);
        }

        Some(Ordering::Equal)
    }
}

/// convert ascii byte slice contain a decimal integer to u64
fn u64_from_str(input: &[u8]) -> Result<u64, ParseIntError> {
    u64::from_str(str::from_utf8(input).unwrap())
}

/// parse an integer from the input
fn parse_u64(input: &[u8]) -> IResult<&[u8], u64> {
    map_res(digit1, u64_from_str)(input)
}

/// return None unconditionally
fn empty_f64_option(_: &[u8]) -> Result<Option<f64>, ParseFloatError> {
    Ok(None)
}

/// convert ascii byte from
fn f64_option(input: &[u8]) -> Result<Option<f64>, ParseFloatError> {
    f64::from_str(str::from_utf8(input).unwrap()).map(Some)
}

/// parse a GTF score field to Option<64>. An empty field '.',
/// will return None, otherwise Some(f64).
fn score(input: &[u8]) -> IResult<&[u8], Option<f64>> {
    alt((
        map_res(tag("."), empty_f64_option),
        map_res(nom::number::complete::recognize_float, f64_option),
    ))(input)
}

/// Is character a valid token for a GTF
/// attribute
#[inline]
fn is_token(c: u8) -> bool {
    !matches!(c,
        128..=255 |
        0..=31 |
        b' ' |
        b'"' |
        b'(' | b')' |
        b',' |
        b'/' |
        b':' |
        b';' |
        b'<' |
        b'=' |
        b'>' |
        b'?' |
        b'@' |
        b'[' |
        b'\\' |
        b']' |
        b'{' | b'}'
    )
}

type AttrVec<'a> = SmallVec<[(&'a [u8], &'a [u8]); 16]>;

fn gtf_attributes(input: &[u8]) -> IResult<&[u8], AttrVec<'_>> {
    terminated(
        separated_list_smallvec(
            pair(tag(";"), take_while1(|c| c == b' ')),
            separated_pair(
                take_while1(is_token),
                take_while1(|c| c == b' '),
                alt((
                    delimited(char('"'), take_until("\""), char('"')),
                    take_while1(is_token),
                )),
            ),
        ),
        opt(pair(opt(tag(";")), take_while(|c| c == b' '))),
    )(input)
}

/// raw fields of a GTF line, separated by tab characters
/// this is used transiently and will be converted to a
/// `Record<'a>`.
type RecInnerSep<'a> = (
    &'a [u8],
    char,
    &'a [u8],
    char,
    &'a [u8],
    char,
    u64,
    char,
    u64,
    char,
    Option<f64>,
    char,
    &'a [u8],
    char,
    &'a [u8],
    char,
    AttrVec<'a>,
);

fn convert_to_record(inp: RecInnerSep<'_>) -> Result<Record<'_>, ParseFloatError> {
    Ok(Record {
        seqname: inp.0,
        source: inp.2,
        feature_type: inp.4,
        start: inp.6,
        end: inp.8,
        score: inp.10,
        strand: inp.12,
        frame: inp.14,
        attributes: inp.16,
    })
}

/// Replacement for `separated_list` in nom, that returns items in a `SamllVec`
/// to avoid allocations in the tight inner loop of GTF attribute parsing.
fn separated_list_smallvec<I, O, O2, E, F, G>(
    mut sep: G,
    mut f: F,
) -> impl FnMut(I) -> IResult<I, SmallVec<[O; 16]>, E>
where
    I: Clone + PartialEq,
    F: FnMut(I) -> IResult<I, O, E>,
    G: FnMut(I) -> IResult<I, O2, E>,
    E: ParseError<I>,
{
    use nom::Err;

    move |mut i: I| {
        let mut res = SmallVec::new();

        match f(i.clone()) {
            Err(Err::Error(_)) => return Ok((i, res)),
            Err(e) => return Err(e),
            Ok((i1, o)) => {
                if i1 == i {
                    return Err(Err::Error(E::from_error_kind(i1, ErrorKind::SeparatedList)));
                }

                res.push(o);
                i = i1;
            }
        }

        loop {
            match sep(i.clone()) {
                Err(Err::Error(_)) => return Ok((i, res)),
                Err(e) => return Err(e),
                Ok((i1, _)) => {
                    if i1 == i {
                        return Err(Err::Error(E::from_error_kind(i1, ErrorKind::SeparatedList)));
                    }

                    match f(i1.clone()) {
                        Err(Err::Error(_)) => return Ok((i, res)),
                        Err(e) => return Err(e),
                        Ok((i2, o)) => {
                            if i2 == i {
                                return Err(Err::Error(E::from_error_kind(
                                    i2,
                                    ErrorKind::SeparatedList,
                                )));
                            }

                            res.push(o);
                            i = i2;
                        }
                    }
                }
            }
        }
    }
}

/// Check that a line of a GTF file is valid.
/// This function is intended to be used on a line that failed parsing in order
/// to provide a user-facing error message.
pub fn validate_gtf_line(line: &[u8]) -> Result<()> {
    let pieces: Vec<_> = line.split(|c| *c == b'\t').collect();
    ensure!(
        pieces.len() == 9,
        "expected 9 tab-separated elements but found {}",
        pieces.len()
    );
    #[allow(clippy::type_complexity)]
    let validators: [(&str, fn(&[u8]) -> Result<()>); 9] = [
        ("seqname", validate_no_space),
        ("source", validate_no_space),
        ("feature", validate_no_space),
        ("start", validate_u64),
        ("end", validate_u64),
        ("score", validate_score),
        ("strand", validate_no_space),
        ("frame", validate_no_space),
        ("attributes", validate_gtf_attributes),
    ];
    for (piece, (item_name, validator)) in std::iter::zip(pieces, validators) {
        if let Err(err) = validator(piece) {
            bail!("{item_name}: {err}");
        }
    }
    Ok(())
}

fn validate_no_space(input: &[u8]) -> Result<()> {
    ensure!(!input.iter().any(|c| *c == b' '), "cannot contain spaces");
    Ok(())
}

fn validate_u64(input: &[u8]) -> Result<()> {
    let input_str = str::from_utf8(input).map_err(|_| anyhow!("invalid UTF-8 character(s)"))?;
    u64::from_str(input_str).map_err(|_| anyhow!("expected an integer, not \"{input_str}\""))?;
    Ok(())
}

fn validate_score(input: &[u8]) -> Result<()> {
    score(input).map_err(|_| anyhow!("expected \".\" or a number"))?;
    Ok(())
}

fn validate_gtf_attributes(input: &[u8]) -> Result<()> {
    println!("{}", str::from_utf8(input).unwrap());
    all_consuming(gtf_attributes)(input).map_err(|_| anyhow!("invalid attributes format"))?;
    Ok(())
}

#[cfg(test)]
mod test {
    use super::*;

    const ORIG: &[u8] = br#"1	havana	exon	29554	30039	.	+	.	gene_id "ENSG00000243485"; gene_version "5"; transcript_id "ENST00000473358"; transcript_version "1"; exon_number "1"; gene_name "MIR1302-2HG"; gene_source "havana"; gene_biotype "lincRNA"; transcript_name "MIR1302-2HG-202"; transcript_source "havana"; transcript_biotype "lincRNA"; exon_id "ENSE00001947070"; exon_version "1"; tag "basic"; transcript_support_level "5""#;
    const NO_QUOTES: &[u8] = br#"1	havana	exon	29554	30039	.	+	.	gene_id ENSG00000243485; gene_version "5"; transcript_id "ENST00000473358"; transcript_version "1"; exon_number "1"; gene_name "MIR1302-2HG"; gene_source "havana"; gene_biotype "lincRNA"; transcript_name "MIR1302-2HG-202"; transcript_source "havana"; transcript_biotype "lincRNA"; exon_id "ENSE00001947070"; exon_version "1"; tag "basic"; transcript_support_level "5""#;
    const BAD1: &[u8] = br#"1	havana	exon	29554	30039	.	+	.	gene    _id ENSG00000243485; gene_version "5"; transcript_id "ENST00000473358"; transcript_version "1"; exon_number "1"; gene_name "MIR1302-2HG"; gene_source "havana"; gene_biotype "lincRNA"; transcript_name "MIR1302-2HG-202"; transcript_source "havana"; transcript_biotype "lincRNA"; exon_id "ENSE00001947070"; exon_version "1"; tag "basic"; transcript_support_level "5""#;
    const TWO_SPACE: &[u8] = br#"1	havana	exon	29554	30039	.	+	.	gene_id  ENSG00000243485; gene_version "5"; transcript_id "ENST00000473358"; transcript_version "1"; exon_number "1"; gene_name "MIR1302-2HG"; gene_source "havana"; gene_biotype "lincRNA"; transcript_name "MIR1302-2HG-202"; transcript_source "havana"; transcript_biotype "lincRNA"; exon_id "ENSE00001947070"; exon_version "1"; tag "basic"; transcript_support_level "5""#;
    const GRCH_120: &[u8] = br#"1	havana	gene	29554	31109	.	+	.	gene_id "ENSG00000243485"; gene_version "3"; gene_name "RP11-34P13.3"; gene_source "havana"; gene_biotype "lincRNA"; havana_gene "OTTHUMG00000000959"; havana_gene_version "2";"#;
    const TRAILING_SEMI_SPACE: &[u8] = br#"1	havana	gene	29554	31109	.	+	.	gene_id "ENSG00000243485"; gene_version "3"; gene_name "RP11-34P13.3"; gene_source "havana"; gene_biotype "lincRNA"; havana_gene "OTTHUMG00000000959"; havana_gene_version "2"; "#;
    const TRAILING_SPACE: &[u8] = br#"1	havana	gene	29554	31109	.	+	.	gene_id "ENSG00000243485"; gene_version "3"; gene_name "RP11-34P13.3"; gene_source "havana"; gene_biotype "lincRNA"; havana_gene "OTTHUMG00000000959"; havana_gene_version "2" "#;
    const EMPTY_ATTR_VAL: &[u8] = br#"chr21	HAVANA	gene	34073578	34106260	.	+	.	gene_id ""; gene_type "protein_coding"; gene_name ""; hgnc_id "HGNC:11038"; tag "overlapping_locus"; tag "retrogene"; havana_gene "OTTHUMG00000065821.4";"#;
    const SPACE_IN_FIELD: &[u8] = br#" 1	havana	exon	29554	30039	.	+	.	gene_id "ENSG00000243485"; gene_version "5"; transcript_id "ENST00000473358"; transcript_version "1"; exon_number "1"; gene_name "MIR1302-2HG"; gene_source "havana"; gene_biotype "lincRNA"; transcript_name "MIR1302-2HG-202"; transcript_source "havana"; transcript_biotype "lincRNA"; exon_id "ENSE00001947070"; exon_version "1"; tag "basic"; transcript_support_level "5""#;

    #[test]
    fn test_gtf_spaces() {
        let original = parse_gtf_line(ORIG);
        let no_quotes = parse_gtf_line(NO_QUOTES);
        assert_eq!(original, no_quotes);

        let two_space = parse_gtf_line(TWO_SPACE);
        assert_eq!(original, two_space);

        assert!(parse_gtf_line(SPACE_IN_FIELD).is_err());
        assert_eq!(
            validate_gtf_line(SPACE_IN_FIELD).unwrap_err().to_string(),
            "seqname: cannot contain spaces"
        );
    }

    #[test]
    fn bad_in_attr_name() {
        assert!(parse_gtf_line(BAD1).is_err());
        assert_eq!(
            validate_gtf_line(BAD1).unwrap_err().to_string(),
            "attributes: invalid attributes format"
        );
    }

    #[test]
    fn trailing_semi() {
        let b = parse_gtf_line(GRCH_120);
        println!("{b:?}");
        assert!(b.is_ok());
    }

    #[test]
    fn trailing_space() {
        let b = parse_gtf_line(TRAILING_SPACE);
        println!("{b:?}");
        assert!(b.is_ok());
    }

    #[test]
    fn trailing_semi_space() {
        let b = parse_gtf_line(TRAILING_SEMI_SPACE);
        println!("{b:?}");
        assert!(b.is_ok());
    }

    #[test]
    fn empty_attr_value() {
        let b = parse_gtf_line(EMPTY_ATTR_VAL);
        println!("{b:?}");
        assert!(b.is_ok());
    }
}

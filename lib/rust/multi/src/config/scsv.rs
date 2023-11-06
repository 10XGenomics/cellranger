use nom::branch::alt;
use nom::bytes::complete::{tag, take_till};
use nom::character::complete::{char, space0, space1};
use nom::combinator::peek;
use nom::error::{make_error, ErrorKind};
use nom::multi::{many0, many1, many_till};
use nom::sequence::{delimited, preceded, terminated};
use nom::{Err, IResult, InputTake, Slice};
use nom_locate::LocatedSpan;
use std::convert::AsRef;
use std::fmt::{Display, Formatter};
use std::path::{Path, PathBuf};
use std::sync::Arc;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct XtraData {
    filename: Arc<PathBuf>,
}

impl XtraData {
    pub fn new(path: impl AsRef<Path>) -> Self {
        XtraData {
            filename: Arc::new(path.as_ref().to_owned()),
        }
    }
}

impl<P: AsRef<Path>> From<P> for XtraData {
    fn from(p: P) -> XtraData {
        XtraData::new(p)
    }
}

impl Display for XtraData {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.filename.display())
    }
}

pub type Span<'a> = LocatedSpan<&'a str, XtraData>;

#[derive(Debug, Clone, Copy)]
pub struct SectionHdr<'a> {
    fragment: &'a str,
    offset: usize,
    line: u32,
    column: usize,
}

impl<'a> From<&Span<'a>> for SectionHdr<'a> {
    fn from(span: &Span<'a>) -> Self {
        SectionHdr {
            fragment: span.fragment(),
            offset: span.location_offset(),
            line: span.location_line(),
            column: span.get_utf8_column(),
        }
    }
}

impl<'a> SectionHdr<'a> {
    pub fn fragment(&self) -> &'a str {
        self.fragment
    }
    pub fn location_offset(&self) -> usize {
        self.offset
    }
    pub fn location_line(&self) -> u32 {
        self.line
    }
    pub fn get_utf8_column(&self) -> usize {
        self.column
    }
}

impl<'a> Display for SectionHdr<'a> {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(f, "[{}]", self.fragment)
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Section<'a, T: Display = SectionHdr<'a>> {
    pub name: T,
    pub rows: Vec<Vec<Span<'a>>>,
}

impl<'a, T: Display> Section<'a, T> {
    fn new<R, S>(name: T, rows: S) -> Self
    where
        R: IntoIterator<Item = Span<'a>>,
        S: IntoIterator<Item = R>,
    {
        Section {
            name,
            rows: rows
                .into_iter()
                .filter_map(|row| {
                    let row = row.into_iter().collect::<Vec<_>>();
                    if row.is_empty() {
                        None
                    } else {
                        Some(row)
                    }
                })
                .collect::<Vec<_>>(),
        }
    }
}

fn line_ending(input: Span<'_>) -> IResult<Span<'_>, Span<'_>> {
    let f = input.fragment();
    match f.get(0..=0) {
        Some("\r") => {
            if f.get(1..=1) == Some("\n") {
                Ok((input.slice(2..), input.slice(..2)))
            } else {
                Ok((input.slice(1..), input.slice(..1)))
            }
        }
        Some("\n") => Ok((input.slice(1..), input.slice(..1))),
        _ => Err(Err::Error(make_error(input, ErrorKind::CrLf))),
    }
}

fn take_until_newline(input: Span<'_>) -> IResult<Span<'_>, Span<'_>> {
    let f = input.fragment();
    for i in 0..f.len() {
        match f.get(i..=i) {
            Some("\r") | Some("\n") => return Ok((input.slice(i..), input.slice(..i))),
            _ => {}
        }
    }
    Err(Err::Error(make_error(input, ErrorKind::Eof)))
}

fn comment(input: Span<'_>) -> IResult<Span<'_>, Span<'_>> {
    preceded(
        space0,
        delimited(char('#'), take_until_newline, line_ending),
    )(input)
}

pub(crate) fn eof(input: Span<'_>) -> IResult<Span<'_>, Span<'_>> {
    let (input, _) = space0(input)?;
    if input.fragment().is_empty() {
        return Ok((input.clone(), input));
    }
    Err(Err::Error(make_error(input, ErrorKind::Eof)))
}

pub(crate) fn take_until_eof(input: Span<'_>) -> IResult<Span<'_>, Span<'_>> {
    let l = input.fragment().len();
    let s = input.slice(l..);
    Ok((s, input.slice(..l)))
}

fn line_end(input: Span<'_>) -> IResult<Span<'_>, Span<'_>> {
    alt((eof, line_ending, comment))(input)
}

fn header_chars(input: Span<'_>) -> IResult<Span<'_>, Span<'_>> {
    let f = input.fragment();
    let l = f.len();
    let mut i = 0usize;
    let mut ws = 0usize;
    loop {
        let s = &f[i..=i];
        if l <= i || s == "\r" || s == "\n" || s == "]" {
            break;
        }
        // keep track of the longest run of whitespace
        if s.chars().all(|c| c.is_ascii_whitespace()) {
            ws += 1;
        } else {
            ws = 0;
        }
        i += 1;
    }
    i -= ws;
    Ok((input.slice(i..), input.slice(..i)))
}

fn header(input: Span<'_>) -> IResult<Span<'_>, Span<'_>> {
    delimited(
        many0(alt((space1, tag(","), line_ending, comment))),
        terminated(
            terminated(
                delimited(
                    char('['),
                    delimited(space0, header_chars, space0),
                    char(']'),
                ),
                many0(alt((tag(","), space1))),
            ),
            line_end,
        ),
        space0,
    )(input)
}

fn inside_double_quotes(input: Span<'_>) -> IResult<Span<'_>, Span<'_>> {
    let f = input.fragment();
    let l = f.len();
    let mut i = 0usize;
    loop {
        if l <= i {
            break;
        }
        if i + 1 < l {
            let twochar = &f[i..=(i + 1)];
            if twochar == "\\\"" || twochar == "\"\"" {
                i += 2;
                continue;
            }
        }
        if &f[i..=i] == "\"" {
            break;
        }
        i += 1;
    }
    Ok((input.slice(i..), input.slice(..i)))
}

fn double_quoted(input: Span<'_>) -> IResult<Span<'_>, Span<'_>> {
    delimited(char('"'), inside_double_quotes, char('"'))(input)
}

fn quoted(input: Span<'_>) -> IResult<Span<'_>, Span<'_>> {
    delimited(space0, double_quoted, space0)(input)
}

// this could technically be IResult<&str, ()>, no?
fn end_cell(input: Span<'_>) -> IResult<Span<'_>, Span<'_>> {
    preceded(space0, alt((peek(tag(",")), peek(line_end), peek(comment))))(input)
}

pub(crate) fn rstrip_space(input: Span<'_>) -> Span<'_> {
    let mut l = input.fragment().len();
    for c in input.fragment().chars().rev() {
        if c == ' ' || c == '\t' {
            l -= 1;
        } else {
            break;
        }
    }
    input.slice(..l)
}

fn not_end_cell(input: Span<'_>) -> IResult<Span<'_>, Span<'_>> {
    let (s, r) = terminated(
        alt((
            take_till(|c| c == ',' || c == '\r' || c == '\n'),
            take_until_eof,
        )),
        space0,
    )(input)?;
    Ok((s, rstrip_space(r)))
}

fn cell1(input: Span<'_>) -> IResult<Span<'_>, Span<'_>> {
    delimited(space0, alt((quoted, not_end_cell)), end_cell)(input)
}

fn cell2(input: Span<'_>) -> IResult<Span<'_>, Span<'_>> {
    delimited(
        alt((delimited(space0, tag(","), space0), space0)),
        alt((quoted, not_end_cell)),
        end_cell,
    )(input)
}

fn row(input: Span<'_>) -> IResult<Span<'_>, Vec<Span<'_>>> {
    // check if row is just a newline here, before parsing cells.
    if let Ok((input, _)) = line_end(input.clone()) {
        return Ok((input, vec![]));
    }

    let (mut input, cell) = cell1(input)?;
    let mut cells = vec![cell];
    // custom many_till because I don't want to insert(0, x) into a vec
    loop {
        if let Ok((input, _)) = line_end(input.clone()) {
            if cells.iter().all(|c| c.fragment() == &"") {
                // if all cells are empty, return an empty row
                return Ok((input, vec![]));
            }
            return Ok((input, cells));
        }
        let (input_, cell) = cell2(input)?;
        cells.push(cell);
        input = input_;
    }
}

fn section(input: Span<'_>) -> IResult<Span<'_>, Section<'_, SectionHdr<'_>>> {
    let (input, name) = header(input)?;
    let (input, (rows, _)) = many_till(row, peek(alt((header, eof))))(input)?;
    Ok((input, Section::new(SectionHdr::from(&name), rows)))
}

pub fn section_csv(input: Span<'_>) -> IResult<Span<'_>, Vec<Section<'_, SectionHdr<'_>>>> {
    many1(section)(input)
}

pub fn plain_csv<T: Display>(name: T, input: Span<'_>) -> IResult<Span<'_>, Section<'_, T>> {
    let (input, _) = input.take_split(0);
    let (input, (rows, _)) = many_till(row, peek(alt((header, eof))))(input)?;
    Ok((input, Section::new(name, rows)))
}

#[cfg(test)]
mod tests {
    use super::{section_csv, Section, SectionHdr, Span, XtraData};
    use anyhow::Result;
    use std::convert::Into;

    #[derive(Debug, Clone, PartialEq, Eq)]
    pub struct SectionStr<'a> {
        name: &'a str,
        rows: Vec<Vec<&'a str>>,
    }

    impl<'a> SectionStr<'a> {
        pub fn new(name: &'a str, rows: Vec<Vec<&'a str>>) -> Self {
            SectionStr { name, rows }
        }
    }

    impl<'a> From<Section<'a, SectionHdr<'a>>> for SectionStr<'a> {
        fn from(s: Section<'a, SectionHdr<'a>>) -> Self {
            SectionStr {
                name: s.name.fragment(),
                rows: s
                    .rows
                    .into_iter()
                    .map(|row| row.into_iter().map(|v| *v.fragment()).collect::<Vec<_>>())
                    .collect::<Vec<_>>(),
            }
        }
    }

    #[test]
    fn one_section() -> Result<()> {
        let xtra = XtraData::new("tests::one_section");
        let scsv = Span::new_extra(
            r#"
            # comment
            [hi]
            yes,1, 2"#,
            xtra,
        );
        let (_input, parsed) = section_csv(scsv)?;
        let expected = SectionStr::new("hi", vec![vec!["yes", "1", "2"]]);
        assert_eq!(expected, parsed[0].clone().into());
        Ok(())
    }

    #[test]
    fn two_sections() -> Result<()> {
        let xtra = XtraData::new("tests::two_sections");
        let scsv = Span::new_extra(
            r#"
            # i am a comment
            [hi] # more comments
            yes,1,3

            #i am also a comment
            [bye] #hey, look, a comment!
            no,2,4,6
            "adios",8

            # i am also a comment!!!
            "#,
            xtra,
        );
        let (_input, parsed) = section_csv(scsv)?;
        let expected = vec![
            SectionStr::new("hi", vec![vec!["yes", "1", "3"]]),
            SectionStr::new("bye", vec![vec!["no", "2", "4", "6"], vec!["adios", "8"]]),
        ];
        assert_eq!(
            expected,
            parsed.into_iter().map(Into::into).collect::<Vec<_>>(),
        );
        Ok(())
    }

    #[test]
    #[ignore = "since # is not a valid character in any field, this is not triggered yet."]
    fn avoid_early_comment_termination() -> Result<()> {
        let xtra = XtraData::new("tests::avoid_early_comment_termination");
        let scsv = Span::new_extra(
            r#"
            # i am a comment
            [hi] #this should be a valid comment too
            col1,col2,col3
            val1,val#2,val3 # real comment"#,
            xtra,
        );
        let (_input, parsed) = section_csv(scsv)?;
        let expected = vec![SectionStr::new(
            "hi",
            vec![vec!["col1", "col2", "col3"], vec!["val1", "val#2", "val3"]],
            // Currently reporting this instead:
            //vec![vec!["col1", "col2", "col3"], vec!["val1", "val#2", "val3 # real comment"]],
        )];
        assert_eq!(
            expected,
            parsed.into_iter().map(Into::into).collect::<Vec<_>>(),
        );
        Ok(())
    }

    #[test]
    fn excel_garbage() -> Result<()> {
        let xtra = XtraData::new("tests::excel_garbage");
        let scsv = Span::new_extra(
            r#"
            ,,# comment for funsies
            ,,
            [garbage],,
            this,1,
            ,,"2, 2",
            totally,",3",
              total garbage  ,  "4, "  ,
            ,,
            ,,# comment
            "#,
            xtra,
        );
        let (_input, parsed) = section_csv(scsv)?;
        let expected = vec![SectionStr::new(
            "garbage",
            vec![
                vec!["this", "1", ""],
                vec!["", "", "2, 2", ""],
                vec!["totally", ",3", ""],
                vec!["total garbage", "4, ", ""],
                // next line should technically be vec!["", "", ""], but see the
                // avoid_early_comment_termination test for issues when parsing inline comments
                vec!["", "", "# comment"],
            ],
        )];
        assert_eq!(
            expected,
            parsed.into_iter().map(Into::into).collect::<Vec<_>>(),
        );
        Ok(())
    }

    #[test]
    fn rfc4180() -> Result<()> {
        let xtra = XtraData::new("tests::rfc4180");
        let scsv = Span::new_extra(
            r#"
[rfc4180]
A,B,C
aaa,"b
bb",ccc
aaa,"b""bb",ccc
aaa,"b\"bb",ccc


            "#,
            xtra,
        );
        let (_input, parsed) = section_csv(scsv)?;
        let expected = vec![SectionStr::new(
            "rfc4180",
            vec![
                vec!["A", "B", "C"],
                vec!["aaa", "b\nbb", "ccc"],
                vec!["aaa", "b\"\"bb", "ccc"],
                vec!["aaa", "b\\\"bb", "ccc"],
            ],
        )];
        assert_eq!(
            expected,
            parsed.into_iter().map(Into::into).collect::<Vec<_>>(),
        );
        Ok(())
    }

    #[test]
    fn spaces_in_headers() -> Result<()> {
        let xtra = XtraData::new("tests::spaces_in_headers");
        let scsv = Span::new_extra(
            r#"
            # comment
            [    hi there      ]
            yes,1, 2
            "#,
            xtra,
        );
        let (_input, parsed) = section_csv(scsv)?;
        let expected = SectionStr::new("hi there", vec![vec!["yes", "1", "2"]]);
        assert_eq!(expected, parsed[0].clone().into());
        Ok(())
    }
}

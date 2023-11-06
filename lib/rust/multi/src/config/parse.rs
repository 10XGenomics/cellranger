use super::scsv::{eof, rstrip_space, take_until_eof, Span};
use anyhow::{bail, Result};
use nom::branch::alt;
use nom::bytes::complete::{tag, take_till};
use nom::character::complete::{alpha1, digit1, one_of, space0, space1};
use nom::combinator::peek;
use nom::multi::many0;
use nom::sequence::{delimited, preceded, terminated};
use nom::{AsChar, IResult, InputTakeAtPosition};
use regex::bytes::Regex;
use std::fmt::{Display, Formatter};
use std::ops::RangeInclusive;
use std::str::FromStr;

#[derive(Clone, Copy)]
pub enum ParseCtx<'a, T: Display> {
    Hdr(T),
    HdrCol(T, &'a str),
    HdrRow(T, usize),
    HdrRowCol(T, usize, &'a str),
}

impl<'a, T: Display> Display for ParseCtx<'a, T> {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        use ParseCtx::{Hdr, HdrCol, HdrRow, HdrRowCol};
        match self {
            Hdr(hdr) | HdrCol(hdr, _) => write!(f, "{hdr}"),
            HdrRow(hdr, row) | HdrRowCol(hdr, row, _) => write!(f, "{hdr} row {row}"),
        }
    }
}

impl<'a, T: Display> ParseCtx<'a, T> {
    pub fn with_col(self, col: &'a str) -> Self {
        use ParseCtx::{Hdr, HdrCol, HdrRow, HdrRowCol};
        match self {
            Hdr(hdr) => HdrCol(hdr, col),
            HdrCol(hdr, _) => HdrCol(hdr, col),
            HdrRow(hdr, row) => HdrRowCol(hdr, row, col),
            HdrRowCol(hdr, row, _) => HdrRowCol(hdr, row, col),
        }
    }
}

pub trait Parse<T: Display> {
    fn parse<R>(&self, ctx: ParseCtx<'_, T>) -> Result<R>
    where
        R: FromStr,
        <R as FromStr>::Err: Display;
}

impl<'a, T: Display> Parse<T> for Span<'a> {
    fn parse<R>(&self, ctx: ParseCtx<'_, T>) -> Result<R>
    where
        R: FromStr,
        <R as FromStr>::Err: Display,
    {
        use ParseCtx::{HdrCol, HdrRowCol};
        match self.fragment().parse::<R>() {
            Ok(result) => Ok(result),
            Err(err) => bail!(
                "{} {} '{}' at line: {}, col: {}: {}",
                ctx,
                match ctx {
                    HdrRowCol(_, _, col) | HdrCol(_, col) => format!("has invalid {col}"),
                    _ => "failed to parse".to_string(),
                },
                self.fragment(),
                self.location_line(),
                self.get_utf8_column(),
                err,
            ),
        }
    }
}

fn end_val(input: Span<'_>) -> IResult<Span<'_>, Span<'_>> {
    preceded(space0, alt((tag("|"), peek(eof))))(input)
}

fn not_end_val(input: Span<'_>) -> IResult<Span<'_>, Span<'_>> {
    let (s, r) = terminated(alt((take_till(|c| c == '|'), take_until_eof)), space0)(input)?;
    Ok((s, rstrip_space(r)))
}

fn val1(input: Span<'_>) -> IResult<Span<'_>, Span<'_>> {
    delimited(space0, not_end_val, end_val)(input)
}

fn val2(input: Span<'_>) -> IResult<Span<'_>, Span<'_>> {
    delimited(many0(alt((space1, tag("|")))), not_end_val, end_val)(input)
}

pub fn parse_vec(input: Span<'_>) -> Result<Vec<Span<'_>>, nom::Err<nom::error::Error<Span<'_>>>> {
    let (mut input, val) = val1(input)?;
    let mut vals = vec![val];
    // custom many_till because I don't want to insert(0, x) into a vec
    loop {
        if let Ok((_, _)) = eof(input.clone()) {
            return Ok(vals);
        }
        let (input_, val) = val2(input)?;
        vals.push(val);
        input = input_;
    }
}

#[allow(dead_code)]
pub fn parse_ident<T>(input: T) -> IResult<T, T>
where
    T: InputTakeAtPosition + Clone,
    <T as InputTakeAtPosition>::Item: AsChar + Copy,
{
    let (s, _) = peek(alpha1)(input)?;
    s.split_at_position_complete(|c| !(c.is_alphanum() || c.as_char() == '_'))
}

#[allow(dead_code)]
fn number(input: Span<'_>) -> IResult<Span<'_>, Span<'_>> {
    let (s, _) = peek(one_of("123456789"))(input)?;
    digit1(s)
}

pub fn parse_range<T>(input: Span<'_>) -> Result<RangeInclusive<T>>
where
    T: FromStr + PartialOrd<T> + Clone,
    <T as FromStr>::Err: std::error::Error + Send + Sync + 'static,
{
    let r = Regex::new(r"^([1-9][0-9]*)(-[1-9][0-9]*)?$")?;
    if let Some(cap) = r.captures(input.fragment().as_bytes()) {
        let m = cap.get(1).unwrap();
        let b = unsafe { std::str::from_utf8_unchecked(m.as_bytes()) };
        let b = b.parse::<T>()?;
        return if let Some(e) = cap.get(2) {
            // skip the '-'
            let e = unsafe { std::str::from_utf8_unchecked(&e.as_bytes()[1..]) };
            let e = e.parse::<T>()?;
            Ok(b..=e)
        } else {
            Ok(b.clone()..=b)
        };
    }
    bail!(
        "invalid range '{}' at line: {}, col: {}",
        input.fragment(),
        input.location_line(),
        input.get_utf8_column()
    )
}

pub fn parse_prefixed_range(input: Span<'_>) -> Result<Vec<String>> {
    let r = Regex::new(r"([1-9][0-9]*)(-[1-9][0-9]*)?$")?;
    if let Some(cap) = r.captures(input.fragment().as_bytes()) {
        let m = cap.get(1).unwrap();
        let b = unsafe { std::str::from_utf8_unchecked(m.as_bytes()) };
        let b = b.parse::<usize>()?;
        let range = if let Some(e) = cap.get(2) {
            let e = unsafe { std::str::from_utf8_unchecked(&e.as_bytes()[1..]) };
            let e = e.parse::<usize>()?;
            b..=e
        } else {
            b..=b
        };
        let prefix = &input.fragment()[..m.start()];
        return Ok(range.map(|n| format!("{prefix}{n}")).collect::<Vec<_>>());
    }
    Ok(vec![input.to_string()])
}

#[cfg(test)]
mod tests {
    use super::super::scsv::{Span, XtraData};
    use super::*;
    use anyhow::Result;

    fn span(s: &str) -> Span<'_> {
        let xtra = XtraData::new("parse::tests");
        Span::new_extra(s, xtra)
    }

    #[test]
    fn test_parse_prefixed_range() -> Result<()> {
        let s = span("cmo_1-3");
        let r = parse_prefixed_range(s)?;
        assert_eq!(r, vec!["cmo_1", "cmo_2", "cmo_3"]);

        let s = span("cmo_1");
        let r = parse_prefixed_range(s)?;
        assert_eq!(r, vec!["cmo_1"]);

        let s = span("cmo");
        let r = parse_prefixed_range(s)?;
        assert_eq!(r, vec!["cmo"]);

        Ok(())
    }

    #[test]
    fn test_parse_ident() -> Result<()> {
        let (s, r) = parse_ident("a1")?;
        assert!(s.is_empty());
        assert_eq!(r, "a1");

        let r = parse_ident("1a");
        assert!(r.is_err());

        Ok(())
    }

    #[test]
    fn test_parse_range() -> Result<()> {
        let s = span("1-3");
        let r = parse_range::<u8>(s)?;
        assert_eq!(r, 1..=3u8);

        let s = span("1");
        let r = parse_range::<u8>(s)?;
        assert_eq!(r, 1..=1u8);

        let s = span("1-A");
        let r = parse_range::<u8>(s);
        assert!(r.is_err());

        Ok(())
    }

    #[test]
    fn test_parse_vec() -> Result<()> {
        let s = span("a|b|c");
        let r = parse_vec(s)?;
        let r = r
            .iter()
            .map(nom_locate::LocatedSpan::fragment)
            .collect::<Vec<_>>();
        assert_eq!(r, vec![&"a", &"b", &"c"]);

        let s = span("a");
        let r = parse_vec(s)?;
        let r = r
            .iter()
            .map(nom_locate::LocatedSpan::fragment)
            .collect::<Vec<_>>();
        assert_eq!(r, vec![&"a"]);

        Ok(())
    }
}

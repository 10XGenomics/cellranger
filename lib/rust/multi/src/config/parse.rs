#![deny(missing_docs)]
use super::scsv::{Span, eof, rstrip_space, space0, space1, take_until_eof};
use anyhow::{Result, bail};
use nom::IResult;
use nom::branch::alt;
use nom::bytes::complete::{tag, take_till};
use nom::combinator::peek;
use nom::multi::many0;
use nom::sequence::{delimited, preceded, terminated};
use regex::Regex;
use std::fmt::{Display, Formatter};
use std::ops::RangeInclusive;
use std::str::FromStr;

#[derive(Clone, Copy)]
pub(crate) enum ParseCtx<'a, T: Display> {
    Hdr(T),
    HdrCol(T, &'a str),
    HdrRow(T, usize),
    HdrRowCol(T, usize, &'a str),
}

impl<T: Display> Display for ParseCtx<'_, T> {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        use ParseCtx::{Hdr, HdrCol, HdrRow, HdrRowCol};
        match self {
            Hdr(hdr) | HdrCol(hdr, _) => write!(f, "{hdr}"),
            HdrRow(hdr, row) | HdrRowCol(hdr, row, _) => write!(f, "{hdr} row {row}"),
        }
    }
}

impl<'a, T: Display> ParseCtx<'a, T> {
    pub(crate) fn with_col(self, col: &'a str) -> Self {
        use ParseCtx::{Hdr, HdrCol, HdrRow, HdrRowCol};
        match self {
            Hdr(hdr) => HdrCol(hdr, col),
            HdrCol(hdr, _) => HdrCol(hdr, col),
            HdrRow(hdr, row) => HdrRowCol(hdr, row, col),
            HdrRowCol(hdr, row, _) => HdrRowCol(hdr, row, col),
        }
    }
}

pub(crate) trait Parse<T: Display> {
    fn parse<R>(&self, ctx: ParseCtx<'_, T>) -> Result<R>
    where
        R: FromStr,
        <R as FromStr>::Err: Display;
}

impl<T: Display> Parse<T> for Span<'_> {
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
    let (s, r) = terminated(
        alt((take_till(|c| c == '|'), |x| Ok(take_until_eof(x)))),
        space0,
    )(input)?;
    Ok((s, rstrip_space(r)))
}

fn val1(input: Span<'_>) -> IResult<Span<'_>, Span<'_>> {
    delimited(space0, not_end_val, end_val)(input)
}

fn val2(input: Span<'_>) -> IResult<Span<'_>, Span<'_>> {
    delimited(many0(alt((space1, tag("|")))), not_end_val, end_val)(input)
}

pub(crate) fn parse_vec(
    input: Span<'_>,
) -> Result<Vec<Span<'_>>, nom::Err<nom::error::Error<Span<'_>>>> {
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

pub(crate) fn parse_range<T>(input: Span<'_>) -> Result<RangeInclusive<T>>
where
    T: FromStr + PartialOrd<T> + Clone,
    <T as FromStr>::Err: std::error::Error + Send + Sync + 'static,
{
    let r = Regex::new(r"^([1-9][0-9]*)(-[1-9][0-9]*)?$")?;
    if let Some(cap) = r.captures(input.fragment()) {
        let m = cap.get(1).unwrap();
        let b = m.as_str().parse::<T>()?;
        return if let Some(e) = cap.get(2) {
            // skip the '-'
            let e = &e.as_str()[1..];
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

pub(crate) fn parse_prefixed_range(input: Span<'_>) -> Result<Vec<String>> {
    let r = Regex::new(r"([1-9][0-9]*)(-[1-9][0-9]*)?$")?;
    if let Some(cap) = r.captures(input.fragment()) {
        let m = cap.get(1).unwrap();
        let b = m.as_str();
        let b = b.parse::<usize>()?;
        let range = if let Some(e) = cap.get(2) {
            let e = &e.as_str()[1..];
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

#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
struct WellIndex {
    col: usize, // Specifying col first so that WellIndex is sorted in column-major order.
    row: char,
}
impl WellIndex {
    fn increment_column_major(&self, plate_layout: &PlateLayout) -> Option<WellIndex> {
        if self.row < plate_layout.well_max.row {
            Some(WellIndex {
                row: ((self.row as u8) + 1) as char,
                col: self.col,
            })
        } else if self.col < plate_layout.well_max.col {
            Some(WellIndex {
                row: plate_layout.well_min.row,
                col: self.col + 1,
            })
        } else {
            None
        }
    }
    fn new(col: usize, row: char, plate_layout: &PlateLayout) -> Result<WellIndex> {
        let well_index = WellIndex { col, row };
        if well_index >= plate_layout.well_min && well_index <= plate_layout.well_max {
            return Ok(well_index);
        }
        anyhow::bail!(
            "Invalid well coordinate. Well coordinates must be between {} and {}.",
            plate_layout.well_min,
            plate_layout.well_max
        )
    }
}
impl Display for WellIndex {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(f, "{:02}{}", self.col, self.row)
    }
}
struct PlateLayout {
    well_min: WellIndex,
    well_max: WellIndex,
}
impl PlateLayout {
    const WELL_96: PlateLayout = PlateLayout {
        well_min: WellIndex { col: 1, row: 'A' },
        well_max: WellIndex { col: 12, row: 'H' },
    };
    fn expand_range_along_plate_columns(
        &self,
        start_well: WellIndex,
        end_well: WellIndex,
    ) -> Result<Vec<WellIndex>> {
        /* Expands a range of plate coordinates by traversing columns.
        e.g. Given a range 1B to 3C, return all wells from 1B through 1H (assuming a 96 well plate),
        the entirety of column 2, and 3A through 3C. This logic is consistent with how samples are specified for Scale. */
        if start_well > end_well {
            anyhow::bail!("Invalid range. Range start must be less than or equal to range end.")
        }
        let mut wells: Vec<WellIndex> = Vec::new();
        let mut current_well = Some(start_well);

        while let Some(well) = current_well {
            wells.push(well);

            if well >= end_well {
                break;
            }
            current_well = well.increment_column_major(self);
        }

        Ok(wells)
    }
}

pub(crate) fn parse_prefixed_plate_range(input: Span<'_>) -> Result<Vec<String>> {
    let plate_layout = PlateLayout::WELL_96;
    let r = Regex::new(r"^([^0-9]*)([0-9]{2})([A-Z])(?:-([0-9]{2})([A-Z]))?$")?;
    let Some(cap) = r.captures(input.fragment()) else {
        bail!(
            "invalid range '{}' at line: {}, col: {}. Expected coordinates like <prefix><two-digit row><upper case column letter>, e.g. prefix01A",
            input.fragment(),
            input.location_line(),
            input.get_utf8_column()
        )
    };
    let prefix = cap.get(1).unwrap().as_str();
    let col_start = cap.get(2).unwrap().as_str().parse::<usize>()?;
    let row_start = cap.get(3).unwrap().as_str().chars().next().unwrap();
    let start_well = WellIndex::new(col_start, row_start, &plate_layout)?;

    let Some(m) = cap.get(4) else {
        return Ok(vec![format!("{prefix}{start_well}")]);
    };
    let col_end = m.as_str().parse::<usize>()?;
    let row_end = cap.get(5).unwrap().as_str().chars().next().unwrap();

    let end_well = WellIndex::new(col_end, row_end, &plate_layout)?;

    Ok(plate_layout
        .expand_range_along_plate_columns(start_well, end_well)?
        .into_iter()
        .map(|well_index| format!("{prefix}{well_index}"))
        .collect())
}

#[cfg(test)]
mod tests {
    use super::super::scsv::{Span, XtraData};
    use super::*;
    use anyhow::Result;
    use nom::character::complete::alpha1;
    use nom::{AsChar, InputTakeAtPosition};

    fn span(s: &str) -> Span<'_> {
        let xtra = XtraData::new("parse::tests");
        Span::new_extra(s, xtra)
    }

    pub(crate) fn parse_ident<T>(input: T) -> IResult<T, T>
    where
        T: InputTakeAtPosition + Clone,
        <T as InputTakeAtPosition>::Item: AsChar + Copy,
    {
        let (s, _) = peek(alpha1)(input)?;
        s.split_at_position_complete(|c| !(c.is_alphanum() || c.as_char() == '_'))
    }

    #[test]
    fn test_parse_prefixed_plate_range() -> Result<()> {
        let s = span("QRT01A");
        let r = parse_prefixed_plate_range(s)?;
        assert_eq!(r, ["QRT01A"]);
        let s = span("QRT01A-03B");
        let r = parse_prefixed_plate_range(s)?;
        assert_eq!(
            r,
            vec![
                "QRT01A", "QRT01B", "QRT01C", "QRT01D", "QRT01E", "QRT01F", "QRT01G", "QRT01H",
                "QRT02A", "QRT02B", "QRT02C", "QRT02D", "QRT02E", "QRT02F", "QRT02G", "QRT02H",
                "QRT03A", "QRT03B"
            ]
        );

        let s = span("QRT1A-3B");
        let err = parse_prefixed_plate_range(s).unwrap_err().to_string();
        assert_eq!(
            err,
            "invalid range 'QRT1A-3B' at line: 1, col: 1. Expected coordinates like <prefix><two-digit row><upper case column letter>, e.g. prefix01A"
        );

        let s = span("QRT01a-02B");
        let err = parse_prefixed_plate_range(s).unwrap_err().to_string();
        assert_eq!(
            err,
            "invalid range 'QRT01a-02B' at line: 1, col: 1. Expected coordinates like <prefix><two-digit row><upper case column letter>, e.g. prefix01A"
        );

        let s = span("QRT03B-01A");
        let err = parse_prefixed_plate_range(s).unwrap_err().to_string();
        assert_eq!(
            err,
            "Invalid range. Range start must be less than or equal to range end."
        );

        let s = span("QRT00A-01A");
        let err = parse_prefixed_plate_range(s).unwrap_err().to_string();
        assert_eq!(
            err,
            "Invalid well coordinate. Well coordinates must be between 01A and 12H."
        );

        let s = span("QRT03D-03A");
        let err = parse_prefixed_plate_range(s).unwrap_err().to_string();
        assert_eq!(
            err,
            "Invalid range. Range start must be less than or equal to range end."
        );
        Ok(())
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

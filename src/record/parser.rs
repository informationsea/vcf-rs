use super::VCFRecord;
use crate::U8Vec;
use nom::{
    self, branch::alt, bytes::complete::is_not, bytes::complete::tag, bytes::complete::take_while1,
    character::is_digit, combinator::opt, combinator::recognize, eof, named, sequence::tuple,
};
use std::str;

pub fn parse_separated_values<'a, U, F, G, E>(
    result: &mut Vec<U8Vec>,
    input: &'a [u8],
    data: F,
    separator: G,
    require_one_entry: bool,
) -> nom::IResult<&'a [u8], (), E>
where
    F: Fn(&'a [u8]) -> nom::IResult<&'a [u8], &'a [u8], E>,
    G: Fn(&'a [u8]) -> nom::IResult<&'a [u8], U, E>,
    E: nom::error::ParseError<&'a [u8]>,
{
    let mut index = 0;
    let mut rest = input;

    loop {
        if let Ok((r, d)) = data(rest) {
            if index < result.len() {
                result[index].clear();
                result[index].extend_from_slice(d);
            } else {
                result.push(d.to_vec());
            }
            index += 1;
            rest = r;
        }
        if let Ok((r, _)) = separator(rest) {
            rest = r;
            continue;
        }
        if index == 0 && require_one_entry {
            return Err(nom::Err::Error(nom::error::make_error(
                input,
                nom::error::ErrorKind::SeparatedNonEmptyList,
            )));
        }
        if index <= result.len() {
            result.drain(index..);
            // remove overflowed content
        }
        return Ok((rest, ()));
    }
}

pub fn parse_nested_separated_values<'a, U, V, F, G, H, E>(
    result: &mut Vec<Vec<U8Vec>>,
    input: &'a [u8],
    data: F,
    separator_inside: H,
    separator_outside: G,
    require_one_entry: bool,
) -> nom::IResult<&'a [u8], (), E>
where
    F: Fn(&'a [u8]) -> nom::IResult<&'a [u8], &'a [u8], E>,
    G: Fn(&'a [u8]) -> nom::IResult<&'a [u8], U, E>,
    H: Fn(&'a [u8]) -> nom::IResult<&'a [u8], V, E>,
    E: nom::error::ParseError<&'a [u8]>,
{
    let mut index = 0;
    let mut rest = input;

    loop {
        while result.len() <= index {
            result.push(Vec::new());
        }
        if let Ok((r, _)) =
            parse_separated_values(&mut result[index], rest, &data, &separator_inside, true)
        {
            index += 1;
            rest = r;
        }
        if let Ok((r, _)) = separator_outside(rest) {
            rest = r;
            continue;
        }
        if index == 0 && require_one_entry {
            return Err(nom::Err::Error(nom::error::make_error(
                input,
                nom::error::ErrorKind::SeparatedNonEmptyList,
            )));
        }
        if index <= result.len() {
            result.drain(index..);
            // remove overflowed content
        }
        return Ok((rest, ()));
    }
}

pub fn parse_double_nested_separated_values<'a, U, V, W, F, G, H, I, E>(
    result: &mut Vec<Vec<Vec<U8Vec>>>,
    input: &'a [u8],
    data: F,
    separator_inside: H,
    separator_outside: G,
    separator_outside2: I,
) -> nom::IResult<&'a [u8], (), E>
where
    F: Fn(&'a [u8]) -> nom::IResult<&'a [u8], &'a [u8], E>,
    G: Fn(&'a [u8]) -> nom::IResult<&'a [u8], U, E>,
    H: Fn(&'a [u8]) -> nom::IResult<&'a [u8], V, E>,
    I: Fn(&'a [u8]) -> nom::IResult<&'a [u8], W, E>,
    E: nom::error::ParseError<&'a [u8]>,
{
    let mut index = 0;
    let mut rest = input;

    loop {
        while result.len() <= index {
            result.push(Vec::new());
        }
        if let Ok((r, _)) = parse_nested_separated_values(
            &mut result[index],
            rest,
            &data,
            &separator_inside,
            &separator_outside,
            false,
        ) {
            index += 1;
            rest = r;
        }
        if let Ok((r, _)) = separator_outside2(rest) {
            rest = r;
            continue;
        }
        if index <= result.len() {
            result.drain(index..);
            // remove overflowed content
        }
        return Ok((rest, ()));
    }
}

lazy_static! {
    static ref EMPTY_INFO: Vec<(U8Vec, Vec<U8Vec>)> = vec![(b".".to_vec(), vec![])];
}

pub fn parse_info<'a, E>(
    input: &'a [u8],
    info: &mut Vec<(U8Vec, Vec<U8Vec>)>,
) -> nom::IResult<&'a [u8], (), E>
where
    E: nom::error::ParseError<&'a [u8]>,
{
    let mut index = 0;
    let mut rest = input;

    while let Ok((r, key)) = is_not::<_, _, E>(&b"\t\r\n=;"[..])(rest) {
        if info.len() <= index {
            info.push((key.to_vec(), Vec::new()));
        } else {
            info[index].0.clear();
            info[index].0.extend_from_slice(key);
        }

        if let Ok((r, _)) = tag::<_, _, E>(b"=")(r) {
            let (r, _) = parse_separated_values(
                &mut info[index].1,
                r,
                is_not(&b"\t\r\n,;"[..]),
                tag(b","),
                false,
            )?;
            rest = r;
        } else {
            info[index].1.clear();
            rest = r;
        }
        index += 1;

        if let Ok((r, _)) = tag::<_, _, E>(b";")(rest) {
            rest = r;
        } else {
            //eprintln!("No colon: {:?}", rest);
            break;
        }
    }

    if index <= info.len() {
        info.drain(index..);
        // remove overflowed content
    }

    if info == &*EMPTY_INFO {
        info.clear();
    }

    Ok((rest, ()))
}

fn parse_float<'a, E>(data: &'a [u8]) -> nom::IResult<&'a [u8], &'a [u8], E>
where
    E: nom::error::ParseError<&'a [u8]>,
{
    alt((
        tag(b"."),
        recognize(tuple((
            take_while1(is_digit),
            opt(tuple((tag(b"."), opt(take_while1(is_digit))))),
        ))),
    ))(data)
}

fn parse_record_optional_columns<'a, E>(
    rest: &'a [u8],
    record: &mut VCFRecord,
) -> nom::IResult<&'a [u8], (), E>
where
    E: nom::error::ParseError<&'a [u8]>,
{
    let rest = match tag::<_, _, E>(b"\t")(rest) {
        Ok((rest, _)) => rest,
        Err(_) => {
            record.qual = None;
            record.filter.clear();
            record.info.clear();
            record.format.clear();
            record.genotype.clear();
            return Ok((rest, ()));
        }
    };
    let (rest, qual) = parse_float(rest)?;
    if qual == b"." {
        record.qual = None;
    } else {
        record.qual = Some(str::from_utf8(qual).unwrap().parse().unwrap());
    }
    let rest = match tag::<_, _, E>(b"\t")(rest) {
        Ok((rest, _)) => rest,
        Err(_) => {
            record.filter.clear();
            record.info.clear();
            record.format.clear();
            record.genotype.clear();
            return Ok((rest, ()));
        }
    };
    let (rest, _) = parse_separated_values(
        &mut record.filter,
        rest,
        is_not(&b"\t\r\n,"[..]),
        tag(b","),
        false,
    )?;
    if record.filter == [b"."] {
        record.filter.clear();
    }
    let rest = match tag::<_, _, E>(b"\t")(rest) {
        Ok((rest, _)) => rest,
        Err(_) => {
            record.info.clear();
            record.format.clear();
            record.genotype.clear();
            return Ok((rest, ()));
        }
    };
    let (rest, _) = parse_info(rest, &mut record.info)?;
    let rest = match tag::<_, _, E>(b"\t")(rest) {
        Ok((rest, _)) => rest,
        Err(_) => {
            record.format.clear();
            record.genotype.clear();
            return Ok((rest, ()));
        }
    };
    let (rest, _) = parse_separated_values(
        &mut record.format,
        rest,
        is_not(&b"\t\r\n:"[..]),
        tag(b":"),
        false,
    )?;
    if record.format == [b"."] {
        record.format.clear();
    }
    let rest = match tag::<_, _, E>(b"\t")(rest) {
        Ok((rest, _)) => rest,
        Err(_) => {
            record.genotype.clear();
            return Ok((rest, ()));
        }
    };
    let (rest, _) = parse_double_nested_separated_values(
        &mut record.genotype,
        rest,
        is_not(&b"\t\r\n:,"[..]),
        tag(b","),
        tag(b":"),
        tag(b"\t"),
    )?;
    if record.genotype == [[[b"."]]] {
        record.genotype.clear();
    }

    Ok((rest, ()))
}

named!(eof_parser, eof!());

pub fn parse_record<'a>(line: &'a [u8], record: &mut VCFRecord) -> nom::IResult<&'a [u8], ()> {
    let (rest, chromosome) = is_not(&b"\t\r\n"[..])(line)?;
    record.chromosome.clear();
    record.chromosome.extend_from_slice(chromosome);
    let (rest, _) = tag(b"\t")(rest)?;

    let (rest, position) = take_while1(is_digit)(rest)?;
    record.position = str::from_utf8(position).unwrap().parse().unwrap();
    let (rest, _) = tag(b"\t")(rest)?;

    let (rest, _) = parse_separated_values(
        &mut record.id,
        rest,
        is_not(&b"\t\r\n,"[..]),
        tag(b","),
        false,
    )?;
    if record.id == [b"."] {
        record.id.clear();
    }
    let (rest, _) = tag(b"\t")(rest)?;
    let (rest, reference) = is_not(&b"\t\r\n"[..])(rest)?;
    record.reference.clear();
    record.reference.extend_from_slice(reference);

    let (rest, _) = tag(b"\t")(rest)?;
    let (rest, _) = parse_separated_values(
        &mut record.alternative,
        rest,
        is_not(&b"\t\r\n,"[..]),
        tag(b","),
        false,
    )?;
    if record.alternative == [b"."] {
        record.alternative.clear();
    }
    let (rest, _) = parse_record_optional_columns(rest, record)?;
    let (rest, _) = alt((tag("\r\n"), tag("\n"), eof_parser))(rest)?;

    record.recreate_info_and_genotype_index();

    Ok((rest, ()))
}

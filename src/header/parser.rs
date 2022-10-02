use super::*;
use nom::{
    self, branch::alt, bytes::complete::is_not, bytes::complete::tag, bytes::complete::take_while,
    character::is_digit, combinator::eof, combinator::map, combinator::opt, multi::separated_list0,
    sequence::separated_pair, sequence::tuple,
};
use std::str;

pub type EntryPair<'a> = (&'a [u8], &'a [u8]);

pub fn parse_header_item(header_line: &[u8]) -> VResult<&[u8], VCFHeaderLine> {
    let line = header_line.to_vec();
    let (rest, _) = tag("##")(header_line)?;
    let (rest, contents) = parse_header_content(rest)?;
    let (rest, _) = alt((tag("\r\n"), tag("\n"), eof))(rest)?;
    Ok((rest, VCFHeaderLine { line, contents }))
}

pub fn parse_header_content(header_line_without_sharp: &[u8]) -> VResult<&[u8], VCFHeaderContent> {
    alt((
        parse_vcf_file_format_header,
        parse_vcf_contig_header,
        parse_vcf_info_header,
        parse_vcf_format_header,
        parse_vcf_alt_header,
        parse_vcf_filter_header,
        parse_other_header_item,
    ))(header_line_without_sharp)
}

pub fn parse_header_entries(value: &[u8]) -> VResult<&[u8], Vec<EntryPair>> {
    separated_list0(
        tag(b","),
        separated_pair(
            is_not(&b">,= \r\n\t"[..]),
            tag(b"="),
            alt((
                map(tuple((tag(b"\""), is_not(&b"\""[..]), tag(b"\""))), |x| x.1),
                is_not(&b">, \r\n\t"[..]),
            )),
        ),
    )(value)
}

pub fn find_key<'a>(entry_pair: &[EntryPair<'a>], key: &[u8]) -> Option<&'a [u8]> {
    entry_pair.iter().find(|(k, _)| *k == key).map(|(_, v)| *v)
}

pub fn find_key_or_error<'a>(
    entry_pair: &[EntryPair<'a>],
    key: &[u8],
    error_content: &'a [u8],
    error_message: &'static str,
) -> Result<&'a [u8], nom::Err<nom::error::VerboseError<&'a [u8]>>> {
    find_key(entry_pair, key).ok_or_else(|| {
        nom::Err::Error(nom::error::VerboseError {
            errors: vec![(
                error_content,
                nom::error::VerboseErrorKind::Context(error_message),
            )],
        })
    })
}

pub fn parse_number(value: &[u8]) -> Number {
    match value {
        b"R" => Number::Reference,
        b"A" => Number::Allele,
        b"G" => Number::Genotype,
        b"0" => Number::Zero,
        b"." => Number::Unknown,
        x if x.iter().copied().all(is_digit) => {
            Number::Number(str::from_utf8(x).unwrap().parse().unwrap())
        }
        _ => Number::Other(value.to_vec()),
    }
}

pub fn parse_value_type(value: &[u8]) -> ValueType {
    match value {
        b"String" => ValueType::String,
        b"Integer" => ValueType::Integer,
        b"Flag" => ValueType::Flag,
        b"Character" => ValueType::Character,
        b"Float" => ValueType::Float,
        _ => ValueType::Other(value.to_vec()),
    }
}

pub fn parse_vcf_file_format_header(header_line: &[u8]) -> VResult<&[u8], VCFHeaderContent> {
    let (rest, _) = tag(b"fileformat=")(header_line)?;
    let (rest, version) = take_while(|x: u8| x != b'\n' && x != b'\r')(rest)?;
    let parsed_version = match version {
        b"VCFv4.3" => VCFVersion::Vcf4_3,
        b"VCFv4.2" => VCFVersion::Vcf4_2,
        b"VCFv4.1" => VCFVersion::Vcf4_1,
        b"VCFv4.0" => VCFVersion::Vcf4_0,
        _ => VCFVersion::Other(version.to_vec()),
    };
    Ok((rest, VCFHeaderContent::FileFormat(parsed_version)))
}

pub fn parse_vcf_contig_header(header_line: &[u8]) -> VResult<&[u8], VCFHeaderContent> {
    let (rest, _) = tag(b"contig=<")(header_line)?;
    let (rest, entries) = parse_header_entries(rest)?;
    let (rest, _) = tag(b">")(rest)?;
    let id = find_key_or_error(&entries, b"ID", header_line, "No ID tag")?.to_vec();
    let length = entries
        .iter()
        .find(|(k, _)| k == b"length")
        .map(|(_, v)| str::from_utf8(v).ok())
        .unwrap_or(Option::None)
        .map(|x| x.parse::<u64>().ok())
        .unwrap_or(Option::None);
    Ok((rest, VCFHeaderContent::Contig { id, length }))
}

pub fn parse_vcf_info_header(header_line: &[u8]) -> VResult<&[u8], VCFHeaderContent> {
    let (rest, _) = tag(b"INFO=<")(header_line)?;
    let (rest, entries) = parse_header_entries(rest)?;
    let (rest, _) = tag(b">")(rest)?;
    let id = find_key_or_error(&entries, b"ID", header_line, "No ID tag")?.to_vec();
    let number = parse_number(find_key_or_error(
        &entries,
        b"Number",
        header_line,
        "No Number tag",
    )?);
    let value_type = parse_value_type(find_key_or_error(
        &entries,
        b"Type",
        header_line,
        "No Type tag",
    )?);
    let description =
        find_key_or_error(&entries, b"Description", header_line, "No Description tag")?.to_vec();
    let source = find_key(&entries, b"Source").map(|x| x.to_vec());
    let version = find_key(&entries, b"Version").map(|x| x.to_vec());

    Ok((
        rest,
        VCFHeaderContent::INFO {
            id,
            number,
            value_type,
            description,
            source,
            version,
        },
    ))
}

pub fn parse_vcf_format_header(header_line: &[u8]) -> VResult<&[u8], VCFHeaderContent> {
    let (rest, _) = tag(b"FORMAT=<")(header_line)?;
    let (rest, entries) = parse_header_entries(rest)?;
    let (rest, _) = tag(b">")(rest)?;
    let id = find_key_or_error(&entries, b"ID", header_line, "No ID tag")?.to_vec();
    let number = parse_number(find_key_or_error(
        &entries,
        b"Number",
        header_line,
        "No Number tag",
    )?);
    let value_type = parse_value_type(find_key_or_error(
        &entries,
        b"Type",
        header_line,
        "No Type tag",
    )?);
    let description =
        find_key_or_error(&entries, b"Description", header_line, "No Description tag")?.to_vec();
    let source = find_key(&entries, b"Source").map(|x| x.to_vec());
    let version = find_key(&entries, b"Version").map(|x| x.to_vec());

    Ok((
        rest,
        VCFHeaderContent::FORMAT {
            id,
            number,
            value_type,
            description,
            source,
            version,
        },
    ))
}

pub fn parse_vcf_filter_header(header_line: &[u8]) -> VResult<&[u8], VCFHeaderContent> {
    let (rest, _) = tag(b"FILTER=<")(header_line)?;
    let (rest, entries) = parse_header_entries(rest)?;
    let (rest, _) = tag(b">")(rest)?;
    let id = find_key_or_error(&entries, b"ID", header_line, "No ID tag")?.to_vec();
    let description =
        find_key_or_error(&entries, b"Description", header_line, "No Description tag")?.to_vec();

    Ok((rest, VCFHeaderContent::FILTER { id, description }))
}

pub fn parse_vcf_alt_header(header_line: &[u8]) -> VResult<&[u8], VCFHeaderContent> {
    let (rest, _) = tag(b"ALT=<")(header_line)?;
    let (rest, entries) = parse_header_entries(rest)?;
    let (rest, _) = tag(b">")(rest)?;
    let id = find_key_or_error(&entries, b"ID", header_line, "No ID tag")?.to_vec();
    let description =
        find_key_or_error(&entries, b"Description", header_line, "No Description tag")?.to_vec();

    Ok((rest, VCFHeaderContent::ALT { id, description }))
}

pub fn parse_other_header_item(header_line: &[u8]) -> VResult<&[u8], VCFHeaderContent> {
    let (rest, _) = is_not(&b"\r\n"[..])(header_line)?;
    Ok((rest, VCFHeaderContent::Other))
}

pub fn parse_samples(header_line: &[u8]) -> VResult<&[u8], Vec<U8Vec>> {
    let (rest, data) = tuple((
        tag(b"#CHROM\tPOS\tID\tREF\tALT"),
        opt(tuple((
            tag("\tQUAL"),
            opt(tuple((
                tag("\tFILTER"),
                opt(tuple((
                    tag("\tINFO"),
                    opt(tuple((
                        tag("\tFORMAT"),
                        opt(tuple((
                            tag("\t"),
                            separated_list0(tag("\t"), is_not(&b"\t\r\n"[..])),
                        ))),
                    ))),
                ))),
            ))),
        ))),
        alt((tag("\r\n"), tag("\n"))),
    ))(header_line)?;

    let samples: Vec<U8Vec> = data
        .1
        .map(|(_, v)| v)
        .flatten()
        .map(|(_, v)| v)
        .flatten()
        .map(|(_, v)| v)
        .flatten()
        .map(|(_, v)| v)
        .flatten()
        .map(|(_, v)| v.iter().map(|x| x.to_vec()).collect())
        .unwrap_or_else(Vec::new);

    Ok((rest, samples))
}

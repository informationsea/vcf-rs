#[macro_use]
extern crate failure;

use std::borrow::Cow;
use std::collections::HashMap;
use std::fmt::{self, Write};
use std::io::{self, BufRead};
use std::iter::Iterator;
use std::str::FromStr;

use indexmap::IndexMap;

#[derive(Fail, Debug)]
pub enum VCFParseError {
    #[fail(display = "IO error: {}", error)]
    IoError { error: io::Error },
    #[fail(display = "Error: {}", _0)]
    Other(&'static str),
    #[fail(display = "Position is not number: {}", _0)]
    PositionIsNotNumber(String),
    #[fail(display = "Too small number of columns")]
    NotEnoughColumns,
}

#[derive(Debug, Clone, Hash, PartialEq, Eq, PartialOrd, Ord)]
pub enum Number {
    Reference,
    Allele,
    Genotype,
    Zero,
    Number(i32),
    Unknown,
    Other(String),
}

#[derive(Debug, Clone, Hash, PartialEq, Eq, PartialOrd, Ord)]
pub enum ValueType {
    String,
    Integer,
    Flag,
    Character,
    Float,
    Other(String),
}

#[derive(Debug, Clone, Hash, PartialEq, Eq, PartialOrd, Ord)]
pub enum VCFHeaderContent {
    INFO {
        id: String,
        number: Number,
        value_type: ValueType,
        description: String,
        source: Option<String>,
        version: Option<String>,
    },
    FORMAT {
        id: String,
        number: Number,
        value_type: ValueType,
        description: String,
        source: Option<String>,
        version: Option<String>,
    },
    Contig {
        id: String,
        length: Option<u64>,
    },
    FileFormat(String),
    Other,
}

fn vcf_header_line_helper<'a>(
    line: &'a str,
) -> Result<(&'a str, &'a str, HashMap<&'a str, &'a str>), VCFParseError> {
    if !line.starts_with("##") {
        return Err(VCFParseError::Other("Header line should starts with ##"));
    }
    if let Some(equal) = line.find('=') {
        let key = &line[2..equal];
        let value = &line[(equal + 1)..];

        if &line[(equal + 1)..(equal + 2)] == "<" && line.ends_with('>') {
            Ok((
                key,
                value,
                vcf_header_parse_helper(&line[(equal + 2)..(line.len() - 1)]),
            ))
        } else {
            Ok((key, value, HashMap::new()))
        }
    } else {
        Ok((&line[2..], &line[0..0], HashMap::new()))
    }
}

#[derive(Debug, Clone, Hash, PartialEq, Eq, PartialOrd, Ord)]
enum HeaderParseState {
    Key,
    Value,
    ValueInQuote,
    QuoteEnded,
}

fn vcf_header_parse_helper<'a>(info: &'a str) -> HashMap<&'a str, &'a str> {
    let mut result = HashMap::new();

    let mut key_start: usize = 0;
    let mut key_end: usize = 0;
    let mut value_start: usize = 0;
    let mut current: usize = 0;
    let mut state = HeaderParseState::Key;

    for ch in info.chars() {
        match state {
            HeaderParseState::Key => {
                if ch == '=' {
                    key_end = current;
                    value_start = current + ch.len_utf8();
                    state = HeaderParseState::Value;
                }
            }
            HeaderParseState::Value => {
                if ch == ',' {
                    result.insert(&info[key_start..key_end], &info[value_start..current]);
                    key_start = current + ch.len_utf8();
                    state = HeaderParseState::Key;
                } else if ch == '"' {
                    value_start = current + ch.len_utf8();
                    state = HeaderParseState::ValueInQuote;
                }
            }
            HeaderParseState::ValueInQuote => {
                if ch == '"' {
                    result.insert(&info[key_start..key_end], &info[value_start..current]);
                    state = HeaderParseState::QuoteEnded;
                }
            }
            HeaderParseState::QuoteEnded => {
                if ch == ',' {
                    state = HeaderParseState::Key;
                    key_start = current + ch.len_utf8();
                }
            }
        }
        current += ch.len_utf8();
    }

    if state == HeaderParseState::Value {
        result.insert(&info[key_start..key_end], &info[value_start..current]);
    }

    result
}

#[derive(Debug, Clone, Hash, PartialEq, Eq, PartialOrd, Ord)]
pub struct VCFHeaderLine {
    pub line: String,
    pub contents: VCFHeaderContent,
}

impl FromStr for VCFHeaderLine {
    type Err = VCFParseError;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let parsed_result = vcf_header_line_helper(s)?;

        let contents = match parsed_result.0 {
            "fileformat" => VCFHeaderContent::FileFormat(parsed_result.1.to_string()),
            "INFO" | "FORMAT" => {
                let id = parsed_result
                    .2
                    .get("ID")
                    .map(|x| x.to_string())
                    .unwrap_or_else(|| "###NO_ID###".to_string());
                let number = match parsed_result.2.get("Number") {
                    Some(&"A") => Number::Allele,
                    Some(&"G") => Number::Genotype,
                    Some(&"R") => Number::Reference,
                    Some(&"0") => Number::Zero,
                    Some(&".") => Number::Unknown,
                    Some(s) => s
                        .parse::<i32>()
                        .map(Number::Number)
                        .unwrap_or_else(|_| Number::Other(s.to_string())),
                    None => Number::Unknown,
                };
                let value_type = match parsed_result.2.get("Type") {
                    Some(&"Integer") => ValueType::Integer,
                    Some(&"String") => ValueType::String,
                    Some(&"Flag") => ValueType::Flag,
                    Some(&"Float") => ValueType::Float,
                    Some(&"Character") => ValueType::Character,
                    Some(s) => ValueType::Other(s.to_string()),
                    None => ValueType::Other("".to_string()),
                };
                let description = parsed_result
                    .2
                    .get("Description")
                    .map(|x| x.to_string())
                    .unwrap_or_else(|| "##NO_DESCRIPTION##".to_string());
                let source = parsed_result.2.get("Source").map(|x| x.to_string());
                let version = parsed_result.2.get("Version").map(|x| x.to_string());

                if parsed_result.0 == "INFO" {
                    VCFHeaderContent::INFO {
                        id,
                        number,
                        value_type,
                        description,
                        source,
                        version,
                    }
                } else if parsed_result.0 == "FORMAT" {
                    VCFHeaderContent::FORMAT {
                        id,
                        number,
                        value_type,
                        description,
                        source,
                        version,
                    }
                } else {
                    unreachable!()
                }
            }
            "contig" => VCFHeaderContent::Contig {
                id: parsed_result
                    .2
                    .get("ID")
                    .map(|x| x.to_string())
                    .unwrap_or_else(|| "###NO_ID###".to_string()),
                length: parsed_result
                    .2
                    .get("length")
                    .map(|x| x.parse::<u64>().ok())
                    .unwrap_or(None),
            },
            _ => VCFHeaderContent::Other,
        };

        Ok(VCFHeaderLine {
            line: s.to_string(),
            contents,
        })
    }
}

#[derive(Debug, Clone, Hash, PartialEq, Eq, PartialOrd, Ord)]
pub struct VCFHeader {
    pub items: Vec<VCFHeaderLine>,
    pub samples: Vec<String>,
}

#[derive(Debug)]
pub struct VCFReader<R: io::BufRead> {
    reader: R,
    header: VCFHeader,
}

impl<R: io::BufRead> VCFReader<R> {
    pub fn header(&self) -> &VCFHeader {
        &self.header
    }
}

impl<R: io::Read> VCFReader<io::BufReader<R>> {
    pub fn new(read: R) -> Result<Self, VCFParseError> {
        let mut reader = io::BufReader::new(read);
        let mut header = VCFHeader {
            items: Vec::new(),
            samples: Vec::new(),
        };

        let mut line = String::new();
        loop {
            line.clear();
            let read_bytes = reader
                .read_line(&mut line)
                .map_err(|e| VCFParseError::IoError { error: e })?;
            if read_bytes == 0 {
                break;
            }

            if line.starts_with("##") {
                header.items.push(line.trim().parse::<VCFHeaderLine>()?);
            } else if line.starts_with('#') {
                // TODO: check header
                let elements: Vec<_> = line.trim().split('\t').collect();
                for one in elements.iter().skip(9) {
                    header.samples.push(one.to_string());
                }
                break;
            }
        }

        Ok(VCFReader { reader, header })
    }
}

impl<R: io::BufRead> Iterator for VCFReader<R> {
    type Item = Result<VCFRecord, VCFParseError>;
    fn next(&mut self) -> Option<Self::Item> {
        let mut line = String::new();
        let result = self.reader.read_line(&mut line);
        match result {
            Ok(read_bytes) => {
                if read_bytes == 0 {
                    None
                } else {
                    match VCFRecord::parse_line(&line, &self.header.samples) {
                        Ok(record) => Some(Ok(record)),
                        Err(e) => Some(Err(e)),
                    }
                }
            }
            Err(e) => Some(Err(VCFParseError::IoError { error: e })),
        }
    }
}

fn encode_vcf_value(value: &str) -> Cow<str> {
    let mut processed_char_num = 0;
    let mut processed_byte_num = 0;

    for one in value.chars() {
        match one {
            ':' | ';' | '=' | '%' | ',' | '\r' | '\n' | '\t' => {
                break;
            }
            _ => {
                processed_byte_num += one.len_utf8();
                processed_char_num += 1;
            }
        }
    }

    if processed_byte_num == value.len() {
        return Cow::Borrowed(value);
    }

    let mut encoded = value[0..processed_byte_num].to_string();

    for one in value.chars().skip(processed_char_num) {
        match one {
            ':' => encoded.push_str("%3A"),
            ';' => encoded.push_str("%3B"),
            '=' => encoded.push_str("%3D"),
            '%' => encoded.push_str("%25"),
            ',' => encoded.push_str("%2C"),
            '\r' => encoded.push_str("%0D"),
            '\n' => encoded.push_str("%0A"),
            '\t' => encoded.push_str("%09"),
            _ => encoded.push(one),
        }
    }

    Cow::Owned(encoded)
}

fn decode_vcf_value(value: &str) -> String {
    let mut decoded = String::new();
    let mut value_chars = value.chars();

    while let Some(next) = value_chars.next() {
        match next {
            '%' => {
                let first = value_chars.next();
                let second = value_chars.next();
                match (first, second) {
                    (Some('3'), Some('A')) => decoded.push(':'),
                    (Some('3'), Some('B')) => decoded.push(';'),
                    (Some('3'), Some('D')) => decoded.push('='),
                    (Some('2'), Some('5')) => decoded.push('%'),
                    (Some('2'), Some('C')) => decoded.push(','),
                    (Some('0'), Some('D')) => decoded.push('\r'),
                    (Some('0'), Some('A')) => decoded.push('\n'),
                    (Some('0'), Some('9')) => decoded.push('\t'),
                    _ => {
                        decoded.push('%');
                        if let Some(x) = first {
                            decoded.push(x);
                        }
                        if let Some(x) = second {
                            decoded.push(x);
                        }
                    }
                }
            }
            _ => {
                decoded.push(next);
            }
        }
    }

    decoded
}

fn dot_value(value: &[&str], index: usize) -> Option<String> {
    if let Some(value) = value.get(index) {
        match *value {
            "." | "" => None,
            _ => Some(decode_vcf_value(value)),
        }
    } else {
        None
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct VCFRecord {
    pub chromosome: String,
    pub position: u64,
    pub id: Option<String>,
    pub reference: String,
    pub alternative: Vec<String>,
    pub quality: Option<String>,
    pub filter: Vec<String>,
    pub info: IndexMap<String, Vec<String>>,
    pub format: Vec<String>,
    pub call: HashMap<String, HashMap<String, Vec<String>>>,
}

impl VCFRecord {
    pub fn parse_line(line: &str, samples: &[String]) -> Result<VCFRecord, VCFParseError> {
        let elements: Vec<_> = line.trim().split('\t').collect();

        if elements.len() < 5 {
            return Err(VCFParseError::NotEnoughColumns);
        }

        let info = if elements.len() < 8 || elements[7] == "." {
            IndexMap::new()
        } else {
            elements[7]
                .split(';')
                .map(|x| {
                    let mut y = x.splitn(2, '=');
                    let first = y.next().unwrap();
                    let second = y.next();
                    (
                        first.to_string(),
                        second
                            .map(|x| {
                                x.split(',')
                                    .map(|z| decode_vcf_value(z))
                                    .collect::<Vec<_>>()
                            })
                            .unwrap_or_else(|| vec![]),
                    )
                })
                .collect()
        };

        let format = if elements.len() < 9 || elements[8] == "." {
            vec![]
        } else {
            elements[8]
                .split(':')
                .map(|x| decode_vcf_value(x))
                .collect()
        };

        let call = elements
            .iter()
            .skip(9)
            .zip(samples)
            .map(|(v, k)| {
                (
                    k.to_string(),
                    v.split(':')
                        .zip(format.iter())
                        .filter(|(x, _)| *x != ".")
                        .map(|(x, f)| {
                            (
                                f.to_string(),
                                x.split(',').map(|z| decode_vcf_value(z)).collect(),
                            )
                        })
                        .collect(),
                )
            })
            .collect();

        Ok(VCFRecord {
            chromosome: decode_vcf_value(elements[0]),
            position: elements[1]
                .parse::<u64>()
                .map_err(|_| VCFParseError::PositionIsNotNumber(elements[1].to_string()))?,
            id: dot_value(&elements, 2),
            reference: decode_vcf_value(elements[3]),
            alternative: elements[4]
                .split(',')
                .map(|x| decode_vcf_value(x))
                .collect(),
            quality: dot_value(&elements, 5),
            filter: if elements.len() < 7 || elements[6] == "." {
                vec![]
            } else {
                elements[6]
                    .split(',')
                    .map(|x| decode_vcf_value(x))
                    .collect()
            },
            info,
            format,
            call,
        })
    }

    pub fn write_line<W: Write>(&self, writer: &mut W, samples: &[String]) -> fmt::Result {
        write!(
            writer,
            "{}\t{}\t",
            encode_vcf_value(&self.chromosome),
            self.position,
        )?;
        write_str_or_dot(writer, &self.id)?;
        write!(writer, "\t{}\t", self.reference)?;
        write_encoded_vec_or_dot(writer, &self.alternative)?;
        writer.write_char('\t')?;
        write_str_or_dot(writer, &self.quality)?;
        writer.write_char('\t')?;
        write_encoded_vec_or_dot(writer, &self.filter)?;

        writer.write_char('\t')?;
        if self.info.is_empty() {
            writer.write_char('.')?;
        } else {
            for (i, (k, v)) in self.info.iter().enumerate() {
                if i != 0 {
                    writer.write_char(';')?;
                }
                writer.write_str(k)?;
                if v.is_empty() {
                    // skip
                } else {
                    writer.write_char('=')?;
                    write_encoded_vec_or_dot(writer, &v)?;
                }
            }
        }

        if !samples.is_empty() {
            writer.write_char('\t')?;
            if self.format.is_empty() {
                writer.write_char('.')?;
            } else {
                for (i, one) in self.format.iter().enumerate() {
                    if i != 0 {
                        writer.write_char(':')?;
                    }
                    writer.write_str(&encode_vcf_value(&one))?;
                }
            }
            writer.write_char('\t')?;
            for (si, one_sample) in samples.iter().enumerate() {
                if si != 0 {
                    writer.write_char('\t')?;
                }
                if let Some(call_result) = self.call.get(one_sample) {
                    for (i, one_format) in self.format.iter().enumerate() {
                        if i != 0 {
                            writer.write_char(':')?;
                        }
                        if let Some(one_call) = call_result.get(one_format) {
                            write_encoded_vec_or_dot(writer, one_call)?;
                        } else {
                            writer.write_char('.')?;
                        }
                    }
                } else {
                    writer.write_char('.')?;
                }
            }
        }

        writer.write_char('\n')?;

        Ok(())
    }
}

fn encode_str_or_dot(value: &Option<String>) -> Cow<str> {
    value
        .as_ref()
        .map(|x| encode_vcf_value(x.as_str()))
        .unwrap_or(Cow::Borrowed("."))
}

fn write_str_or_dot<W: Write>(writer: &mut W, value: &Option<String>) -> fmt::Result {
    write!(
        writer,
        "{}",
        value
            .as_ref()
            .map(|x| encode_vcf_value(x.as_str()))
            .unwrap_or(Cow::Borrowed("."))
    )?;
    Ok(())
}

fn write_encoded_vec_or_dot<W: Write>(writer: &mut W, value: &[String]) -> fmt::Result {
    if value.is_empty() {
        writer.write_str(".")?;
    } else {
        for (i, v) in value.iter().enumerate() {
            if i != 0 {
                writer.write_str(",")?;
            }
            write!(writer, "{}", encode_vcf_value(v))?;
        }
    }

    Ok(())
}

fn str_or_dot(value: &Option<String>) -> &str {
    value.as_ref().map(|x| x.as_str()).unwrap_or(".")
}

#[cfg(test)]
mod test;

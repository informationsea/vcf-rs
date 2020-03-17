use super::{U8Vec, VCFError, VCFErrorKind, VResult};
use std::collections::{hash_map::Keys, HashMap};
use std::io::BufRead;
use std::str::FromStr;
mod parser;

pub use parser::parse_header_item;

/// A number of entries of INFO or FORMAT.
#[derive(Debug, Clone, Hash, PartialEq, Eq, PartialOrd, Ord)]
pub enum Number {
    Reference,
    Allele,
    Genotype,
    Zero,
    Number(i32),
    Unknown,
    Other(U8Vec),
}

/// An entry value type of INFO or FORMAT.
#[derive(Debug, Clone, Hash, PartialEq, Eq, PartialOrd, Ord)]
pub enum ValueType {
    String,
    Integer,
    Flag,
    Character,
    Float,
    Other(U8Vec),
}

#[derive(Debug, Clone, Hash, PartialEq, Eq, PartialOrd, Ord)]
pub enum VCFVersion {
    Vcf4_3,
    Vcf4_2,
    Vcf4_1,
    Vcf4_0,
    Other(U8Vec),
}

/// A content of header line.
#[derive(Debug, Clone, Hash, PartialEq, Eq, PartialOrd, Ord)]
pub enum VCFHeaderContent {
    INFO {
        id: U8Vec,
        number: Number,
        value_type: ValueType,
        description: U8Vec,
        source: Option<U8Vec>,
        version: Option<U8Vec>,
    },
    FORMAT {
        id: U8Vec,
        number: Number,
        value_type: ValueType,
        description: U8Vec,
        source: Option<U8Vec>,
        version: Option<U8Vec>,
    },
    ALT {
        id: U8Vec,
        description: U8Vec,
    },
    FILTER {
        id: U8Vec,
        description: U8Vec,
    },
    Contig {
        id: U8Vec,
        length: Option<u64>,
    },
    FileFormat(VCFVersion),
    Other,
}

#[derive(Debug, Clone, Hash, PartialEq, Eq, PartialOrd, Ord)]
pub struct VCFHeaderInfoFormat<'a> {
    pub id: &'a [u8],
    pub number: &'a Number,
    pub value_type: &'a ValueType,
    pub description: &'a [u8],
    pub source: Option<&'a [u8]>,
    pub version: Option<&'a [u8]>,
}

#[derive(Debug, Clone, Hash, PartialEq, Eq, PartialOrd, Ord)]
pub struct VCFHeaderFilterAlt<'a> {
    pub id: &'a [u8],
    pub description: &'a [u8],
}

/// A header line.
#[derive(Debug, Clone, Hash, PartialEq, Eq, PartialOrd, Ord)]
pub struct VCFHeaderLine {
    line: U8Vec,
    contents: VCFHeaderContent,
}

impl VCFHeaderLine {
    pub fn from_bytes(line: &[u8], line_num: u64) -> Result<Self, VCFError> {
        parse_header_item(line)
            .map_err(|_| VCFErrorKind::HeaderParseError(line_num).into())
            .map(|x| x.1)
    }
    pub fn line(&self) -> &[u8] {
        &self.line
    }

    pub fn contents(&self) -> &VCFHeaderContent {
        &self.contents
    }
}

impl FromStr for VCFHeaderLine {
    type Err = VCFError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        parse_header_item(s.as_bytes())
            .map_err(|_| VCFErrorKind::HeaderParseError(0).into())
            .map(|x| x.1)
    }
}

/// VCF header struct.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct VCFHeader {
    items: Vec<VCFHeaderLine>,
    samples: Vec<U8Vec>,
    info_key: HashMap<U8Vec, usize>,
    format_key: HashMap<U8Vec, usize>,
    alt_key: HashMap<U8Vec, usize>,
    filter_key: HashMap<U8Vec, usize>,
    sample_to_index: HashMap<U8Vec, usize>,
}

impl VCFHeader {
    pub fn new(items: Vec<VCFHeaderLine>, samples: Vec<U8Vec>) -> VCFHeader {
        VCFHeader {
            info_key: create_info_key(&items),
            format_key: create_format_key(&items),
            alt_key: create_alt_key(&items),
            filter_key: create_filter_key(&items),
            sample_to_index: samples
                .iter()
                .enumerate()
                .map(|(k, v)| (v.to_vec(), k))
                .collect(),
            items,
            samples,
        }
    }
    pub fn items(&self) -> &[VCFHeaderLine] {
        &self.items
    }

    pub fn samples(&self) -> &[U8Vec] {
        &self.samples
    }

    pub fn info_list(&self) -> Keys<U8Vec, usize> {
        self.info_key.keys()
    }

    pub fn info<'a>(&'a self, key: &[u8]) -> Option<VCFHeaderInfoFormat<'a>> {
        self.info_key
            .get(key)
            .map(|x| match &self.items[*x].contents() {
                VCFHeaderContent::INFO {
                    id,
                    number,
                    value_type,
                    description,
                    source,
                    version,
                } => VCFHeaderInfoFormat {
                    id,
                    number,
                    value_type,
                    description,
                    source: source.as_ref().map(|x| -> &[u8] { &x }),
                    version: version.as_ref().map(|x| -> &[u8] { &x }),
                },
                _ => unreachable!(),
            })
    }

    pub fn format_list(&self) -> Keys<U8Vec, usize> {
        self.format_key.keys()
    }
    pub fn format<'a>(&'a self, key: &[u8]) -> Option<VCFHeaderInfoFormat<'a>> {
        self.format_key
            .get(key)
            .map(|x| match &self.items[*x].contents() {
                VCFHeaderContent::FORMAT {
                    id,
                    number,
                    value_type,
                    description,
                    source,
                    version,
                } => VCFHeaderInfoFormat {
                    id,
                    number,
                    value_type,
                    description,
                    source: source.as_ref().map(|x| -> &[u8] { &x }),
                    version: version.as_ref().map(|x| -> &[u8] { &x }),
                },
                _ => unreachable!(),
            })
    }

    pub fn alt_list(&self) -> Keys<U8Vec, usize> {
        self.alt_key.keys()
    }

    pub fn alt<'a>(&'a self, key: &[u8]) -> Option<VCFHeaderFilterAlt<'a>> {
        self.alt_key
            .get(key)
            .map(|x| match &self.items[*x].contents() {
                VCFHeaderContent::ALT { id, description } => VCFHeaderFilterAlt { id, description },
                _ => unreachable!(),
            })
    }

    pub fn filter_list(&self) -> Keys<U8Vec, usize> {
        self.filter_key.keys()
    }

    pub fn filter<'a>(&'a self, key: &[u8]) -> Option<VCFHeaderFilterAlt<'a>> {
        self.filter_key
            .get(key)
            .map(|x| match &self.items[*x].contents() {
                VCFHeaderContent::FILTER { id, description } => {
                    VCFHeaderFilterAlt { id, description }
                }
                _ => unreachable!(),
            })
    }

    pub fn sample_index(&self, sample_name: &[u8]) -> Option<usize> {
        self.sample_to_index.get(sample_name).cloned()
    }
}

fn create_info_key(header_line: &[VCFHeaderLine]) -> HashMap<U8Vec, usize> {
    header_line
        .iter()
        .enumerate()
        .filter_map(|x| match &x.1.contents {
            VCFHeaderContent::INFO { id, .. } => Some((id.to_vec(), x.0)),
            _ => None,
        })
        .collect()
}

fn create_format_key(header_line: &[VCFHeaderLine]) -> HashMap<U8Vec, usize> {
    header_line
        .iter()
        .enumerate()
        .filter_map(|x| match &x.1.contents {
            VCFHeaderContent::FORMAT { id, .. } => Some((id.to_vec(), x.0)),
            _ => None,
        })
        .collect()
}

fn create_alt_key(header_line: &[VCFHeaderLine]) -> HashMap<U8Vec, usize> {
    header_line
        .iter()
        .enumerate()
        .filter_map(|x| match &x.1.contents {
            VCFHeaderContent::ALT { id, .. } => Some((id.to_vec(), x.0)),
            _ => None,
        })
        .collect()
}

fn create_filter_key(header_line: &[VCFHeaderLine]) -> HashMap<U8Vec, usize> {
    header_line
        .iter()
        .enumerate()
        .filter_map(|x| match &x.1.contents {
            VCFHeaderContent::FILTER { id, .. } => Some((id.to_vec(), x.0)),
            _ => None,
        })
        .collect()
}

pub fn parse_header<R: BufRead>(
    reader: &mut R,
) -> Result<(u64, Option<U8Vec>, VCFHeader), VCFError> {
    let mut line_num: u64 = 0;
    let mut items = Vec::new();

    loop {
        let mut buffer = Vec::new();
        line_num += 1;
        reader.read_until(b'\n', &mut buffer)?;
        if buffer.starts_with(b"##") {
            let item = parser::parse_header_item(&buffer)
                .map_err::<VCFError, _>(|_| VCFErrorKind::HeaderParseError(line_num).into())?;
            items.push(item.1);
        } else if buffer.starts_with(b"#") {
            let samples = parser::parse_samples(&buffer)
                .map_err::<VCFError, _>(|_| VCFErrorKind::HeaderParseError(line_num).into())?
                .1;
            return Ok((line_num, None, VCFHeader::new(items, samples)));
        } else {
            return Ok((line_num, Some(buffer), VCFHeader::new(items, Vec::new())));
        }
    }
}

#[cfg(test)]
mod test;

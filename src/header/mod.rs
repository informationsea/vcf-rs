use super::*;

/// A number of entries of INFO or FORMAT.
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

/// An entry value type of INFO or FORMAT.
#[derive(Debug, Clone, Hash, PartialEq, Eq, PartialOrd, Ord)]
pub enum ValueType {
    String,
    Integer,
    Flag,
    Character,
    Float,
    Other(String),
}

/// A content of header line.
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

/// A header line.
#[derive(Debug, Clone, Hash, PartialEq, Eq, PartialOrd, Ord)]
pub struct VCFHeaderLine {
    pub line: String,
    pub contents: VCFHeaderContent,
}

impl VCFHeaderLine {
    pub fn new(s: &str) -> Result<Self, VCFParseError> {
        s.parse::<VCFHeaderLine>()
    }
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

/// VCF header struct.
#[derive(Debug, Clone, Hash, PartialEq, Eq, PartialOrd, Ord)]
pub struct VCFHeader {
    pub items: Vec<VCFHeaderLine>,
    pub samples: Vec<String>,
}

#[cfg(test)]
mod test;
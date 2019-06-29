//! Rust implementation of vcf parser/writer without strict validation.
//!
//! # Example
//!
//! ```
//! use vcf::*;
//! use flate2::read::MultiGzDecoder;
//! use std::fs::File;
//!
//! # fn main() { let _ = run(); }
//! # fn run() -> Result<(), VCFParseError> {
//! let mut vcf_reader = VCFReader::new(MultiGzDecoder::new(
//!     File::open("testfiles/ALL.chr20.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.20-34001094-34168504-subset.vcf.gz")
//!         .map_err(|e| VCFParseError::IoError{error:e})?
//! ))?;
//! assert_eq!(vcf_reader.header().items.len(), 255);
//! assert_eq!(vcf_reader.header().items[0],
//!     VCFHeaderLine{
//!         line: "##fileformat=VCFv4.1".to_string(),
//!         contents: VCFHeaderContent::FileFormat("VCFv4.1".to_string())
//!     });
//! assert_eq!(
//!     vcf_reader.header().samples,
//!     vec!["HG00096", "HG00097", "HG00099"]
//! );
//!
//! // load next record
//! let record = vcf_reader.iter().next().unwrap()?;
//! assert_eq!(record.chromosome, "20");
//! assert_eq!(record.position, 34001111);
//! assert_eq!(record.id, vec!["rs565014200"]);
//! assert_eq!(record.reference, "T");
//! assert_eq!(record.alternative, vec!["C"]);
//! assert_eq!(record.quality, Some("100".to_string()));
//! assert_eq!(record.info["AN"], vec!["6"]); // vcf-rs does not validate a number of entries and type
//! assert_eq!(record.call["HG00096"]["GT"], vec!["0|0"]);
//! #   Ok(())
//! # }
//! ```

#[macro_use]
extern crate failure;

mod header;

mod read;
mod record;
mod write;
use std::borrow::Cow;
use std::collections::HashMap;
use std::io::{self, Write};
use std::iter::Iterator;
use std::str::FromStr;

use indexmap::IndexMap;

pub use header::{Number, VCFHeader, VCFHeaderContent, VCFHeaderLine, ValueType};
pub use read::VCFReader;
pub use write::VCFWriter;
pub use record::VCFRecord;

/// Error returned when something wrong in a vcf format.
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

impl From<io::Error> for VCFParseError {
    fn from(e: io::Error) -> Self {
        VCFParseError::IoError { error: e }
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

fn encode_str_or_dot(value: &Option<String>) -> Cow<str> {
    value
        .as_ref()
        .map(|x| encode_vcf_value(x.as_str()))
        .unwrap_or(Cow::Borrowed("."))
}

fn write_str_or_dot<W: io::Write>(writer: &mut W, value: &Option<String>) -> io::Result<()> {
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

fn write_encoded_vec_or_dot<W: Write>(writer: &mut W, value: &[String]) -> io::Result<()> {
    if value.is_empty() {
        writer.write_all(b".")?;
    } else {
        for (i, v) in value.iter().enumerate() {
            if i != 0 {
                writer.write_all(b",")?;
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
#[macro_export]
macro_rules! hash {
    ( $($v:expr),* ) => {
        {
            let mut map = std::collections::HashMap::new();
            $(
                map.insert($v.0, $v.1);
            )*
            map
        }
    }
}

#[cfg(test)]
#[macro_export]
macro_rules! shash {
    ( $($v:expr),* ) => {
        {
            let mut map = std::collections::HashMap::new();
            $(
                map.insert($v.0.to_string(), $v.1);
            )*
            map
        }
    }
}

#[cfg(test)]
#[macro_export]
macro_rules! sihash {
    ( $($v:expr),* ) => {
        {
            let mut map = IndexMap::new();
            $(
                map.insert($v.0.to_string(), $v.1);
            )*
            map
        }
    }
}

#[cfg(test)]
#[macro_export]
macro_rules! sshash {
    ( $($v:expr),* ) => {
        {
            let mut map = std::collections::HashMap::new();
            $(
                map.insert($v.0.to_string(), $v.1.to_string());
            )*
            map
        }
    }
}

#[cfg(test)]
#[macro_export]
macro_rules! svec {
    ( $( $v:expr ),* ) => {
        {
            let mut tmp_vec = Vec::new();
            $(
                tmp_vec.push($v.to_string());
            )*
            tmp_vec
        }
    };
}

#[cfg(test)]
pub mod test;

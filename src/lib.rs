//! Rust implementation of vcf parser/writer
//!
//! # Example
//!
//! ```
//! use vcf::{VCFReader, U8Vec, VCFHeaderFilterAlt, VCFError, VCFRecord};
//! use flate2::read::MultiGzDecoder;
//! use std::fs::File;
//! use std::io::BufReader;
//!
//! fn main() -> Result<(), VCFError> {
//!     let mut reader = VCFReader::new(BufReader::new(MultiGzDecoder::new(File::open(
//!         "./testfiles/1kGP-subset.vcf.gz",
//!     )?)))?;
//!
//!     // access FILTER contents
//!     assert_eq!(
//!         Some(VCFHeaderFilterAlt {
//!             id: b"PASS",
//!             description: b"All filters passed"
//!         }),
//!         reader.header().filter(b"PASS")
//!     );
//!
//!     // access INFO contents
//!     assert_eq!(
//!         b"Stop position of the interval",
//!         reader.header().info(b"END").unwrap().description
//!     );
//!
//!     // prepare VCFRecord object
//!     let mut vcf_record = VCFRecord::new(reader.header());
//!
//!     // read one record
//!     reader.next_record(&mut vcf_record)?;
//!
//!    // get record attributes
//!    assert_eq!(vcf_record.chromosome, b"13");
//!    assert_eq!(vcf_record.position, 32872836);
//!    assert_eq!(vcf_record.id, Vec::<U8Vec>::new());
//!    assert_eq!(vcf_record.reference, b"A");
//!    assert_eq!(vcf_record.alternative, vec![b"C"]);
//!    assert_eq!(vcf_record.qual, Some(495.23));
//!    assert_eq!(vcf_record.info(b"AC"), Some(&vec![b"1".to_vec()]));
//!    assert_eq!(
//!        vcf_record.genotype(b"SRP150637__HG00099", b"GT"),
//!        Some(&vec![b"0/0".to_vec()])
//!    );
//!    assert_eq!(
//!        vcf_record.genotype(b"SRP150637__HG00099", b"AD"),
//!        Some(&vec![b"31".to_vec(), b"0".to_vec()])
//!    );
//!
//!     Ok(())
//! }
//!
//! ```

#[macro_use]
extern crate failure;

#[macro_use]
extern crate lazy_static;

use std::io::prelude::*;
use std::rc::Rc;

mod error;
mod header;
mod record;

pub use error::{VCFError, VCFErrorKind};
pub use header::{
    Number, VCFHeader, VCFHeaderContent, VCFHeaderFilterAlt, VCFHeaderInfoFormat, VCFHeaderLine,
    VCFVersion, ValueType,
};
pub use record::VCFRecord;
pub type U8Vec = Vec<u8>;
type VResult<I, O> = nom::IResult<I, O, nom::error::VerboseError<I>>;

pub struct VCFReader<R: BufRead> {
    buffer: Vec<u8>,
    unprocessed_line: Option<Vec<u8>>,
    current_line: u64,
    reader: R,
    vcf_header: Rc<VCFHeader>,
}

impl<R: BufRead> VCFReader<R> {
    pub fn new(mut reader: R) -> Result<Self, VCFError> {
        let (current_line, unprocessed_line, vcf_header) = header::parse_header(&mut reader)?;
        Ok(VCFReader {
            buffer: Vec::new(),
            unprocessed_line,
            current_line,
            reader,
            vcf_header: Rc::new(vcf_header),
        })
    }

    /// Read next record.
    /// Return false if no record was remained.
    pub fn next_record(&mut self, record: &mut record::VCFRecord) -> Result<bool, VCFError> {
        if let Some(val) = self.unprocessed_line.as_ref() {
            self.buffer.extend(val);
            self.unprocessed_line = None;
        } else {
            self.buffer.clear();
            self.reader.read_until(b'\n', &mut self.buffer)?;
            self.current_line += 1;
        }
        if self.buffer.is_empty() {
            return Ok(false);
        }

        record::parse_record(&self.buffer, record)
            .map_err(|_| VCFErrorKind::RecordParseError(self.current_line))?;

        Ok(true)
    }

    pub fn header(&self) -> Rc<header::VCFHeader> {
        self.vcf_header.clone()
    }
}

pub struct VCFWriter<W: Write> {
    writer: W,
    // header: VCFHeader,
}

impl<W: Write> VCFWriter<W> {
    pub fn new(mut writer: W, header: &VCFHeader) -> Result<Self, VCFError> {
        for one in header.items() {
            writer.write_all(one.line())?;
        }
        write!(writer, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")?;
        if !header.samples().is_empty() {
            write!(writer, "\tFORMAT")?;
            for one in header.samples() {
                writer.write_all(b"\t")?;
                writer.write_all(one)?;
            }
        }
        writer.write_all(b"\n")?;

        Ok(VCFWriter { writer })
    }

    pub fn write_record(&mut self, vcf_record: &VCFRecord) -> Result<(), VCFError> {
        vcf_record.write_record(&mut self.writer)?;
        Ok(())
    }
}

#[cfg(test)]
mod test;

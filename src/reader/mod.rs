use super::{VCFHeader, VCFHeaderLine, VCFParseError, VCFRecord};

use std::io;
use std::io::prelude::*;
#[derive(Debug)]
pub struct VCFReader<R: BufRead> {
    reader: R,
    header: VCFHeader,
}

impl<R: BufRead> VCFReader<R> {
    pub fn header(&self) -> &VCFHeader {
        &self.header
    }
}

impl<R: Read> VCFReader<io::BufReader<R>> {
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

impl<R: BufRead> VCFReader<R> {
    pub fn next_item(&mut self) -> Option<Result<VCFRecord, VCFParseError>> {
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

    pub fn iter(&mut self) -> Iter<'_, R> {
        Iter { vcf_reader: self }
    }
}

impl<R: BufRead> IntoIterator for VCFReader<R> {
    type Item = Result<VCFRecord, VCFParseError>;
    type IntoIter = IntoIter<R>;

    fn into_iter(self) -> IntoIter<R> {
        IntoIter { vcf_reader: self }
    }
}

impl<'a, R: BufRead> IntoIterator for &'a mut VCFReader<R> {
    type Item = Result<VCFRecord, VCFParseError>;
    type IntoIter = Iter<'a, R>;

    fn into_iter(self) -> Iter<'a, R> {
        Iter { vcf_reader: self }
    }
}

#[derive(Debug)]
pub struct Iter<'a, R: BufRead> {
    vcf_reader: &'a mut VCFReader<R>,
}

impl<'a, R: BufRead> Iterator for Iter<'a, R> {
    type Item = Result<VCFRecord, VCFParseError>;
    fn next(&mut self) -> Option<Self::Item> {
        self.vcf_reader.next_item()
    }
}

#[derive(Debug)]
pub struct IntoIter<R: BufRead> {
    vcf_reader: VCFReader<R>,
}

impl<R: BufRead> Iterator for IntoIter<R> {
    type Item = Result<VCFRecord, VCFParseError>;
    fn next(&mut self) -> Option<Self::Item> {
        self.vcf_reader.next_item()
    }
}

#[cfg(test)]
mod test;
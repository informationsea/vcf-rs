mod parser;

use crate::{U8Vec, VCFError, VCFErrorKind, VCFHeader};
pub use parser::parse_record;
use std::collections::HashMap;
use std::io::{self, Write};
use std::usize;

pub const NOT_FOUND: usize = usize::MAX;

#[derive(Debug, Clone, PartialEq)]
pub struct VCFRecord {
    header: VCFHeader,
    pub chromosome: U8Vec,
    pub position: u64,
    pub id: Vec<U8Vec>,
    pub reference: U8Vec,
    pub alternative: Vec<U8Vec>,
    pub qual: Option<f64>,
    pub filter: Vec<U8Vec>,
    pub info: Vec<(U8Vec, Vec<U8Vec>)>,
    info_index: HashMap<U8Vec, usize>,
    pub format: Vec<U8Vec>,
    format_index: HashMap<U8Vec, usize>,
    pub genotype: Vec<Vec<Vec<U8Vec>>>,
}

impl VCFRecord {
    pub fn new(header: VCFHeader) -> Self {
        VCFRecord {
            header,
            chromosome: vec![],
            position: 0,
            id: vec![],
            reference: vec![],
            alternative: vec![],
            qual: None,
            filter: vec![],
            info: vec![],
            info_index: HashMap::new(),
            format: vec![],
            format_index: HashMap::new(),
            genotype: vec![],
        }
    }

    pub fn from_bytes(line: &[u8], line_num: u64, header: VCFHeader) -> Result<Self, VCFError> {
        let mut record = VCFRecord::new(header);
        record.parse_bytes(line, line_num)?;
        Ok(record)
    }

    pub fn parse_bytes(&mut self, line: &[u8], line_num: u64) -> Result<(), VCFError> {
        parse_record::<nom::error::VerboseError<_>>(line, self)
            .map_err(|_| -> VCFError { VCFErrorKind::RecordParseError(line_num).into() })?;
        Ok(())
    }

    pub fn header(&self) -> &VCFHeader {
        &self.header
    }

    pub fn info(&self, key: &[u8]) -> Option<&Vec<U8Vec>> {
        self.info_index
            .get(key)
            .map(|x| self.info.get(*x).map(|y| &y.1))
            .flatten()
    }

    pub fn info_mut(&mut self, key: &[u8]) -> Option<&mut Vec<U8Vec>> {
        self.info_index
            .get(key)
            .cloned()
            .map(move |x| self.info.get_mut(x).map(|y| &mut y.1))
            .flatten()
    }

    pub fn insert_info(&mut self, key: &[u8], mut values: Vec<U8Vec>) -> Option<Vec<U8Vec>> {
        if let Some(x) = self.info_mut(key) {
            let mut ret = Vec::new();
            ret.append(x);
            x.append(&mut values);
            Some(ret)
        } else {
            self.info.push((key.to_vec(), values));
            self.info_index.insert(key.to_vec(), self.info.len() - 1);
            None
        }
    }

    pub fn genotype(&self, sample_name: &[u8], key: &[u8]) -> Option<&Vec<U8Vec>> {
        self.header
            .sample_index(sample_name)
            .map(|x| {
                self.format_index
                    .get(key)
                    .map(|y| self.genotype.get(x).map(|z| z.get(*y)))
            })
            .flatten()
            .flatten()
            .flatten()
    }

    pub fn genotype_mut(&mut self, sample_name: &[u8], key: &[u8]) -> Option<&mut Vec<U8Vec>> {
        self.header
            .sample_index(sample_name)
            .map(move |x| {
                self.format_index
                    .get(key)
                    .cloned()
                    .map(move |y| self.genotype.get_mut(x).map(|z| z.get_mut(y)))
            })
            .flatten()
            .flatten()
            .flatten()
    }

    pub fn insert_genotype(
        &mut self,
        sample_name: &[u8],
        key: &[u8],
        mut values: Vec<U8Vec>,
    ) -> Option<Vec<U8Vec>> {
        if let Some(sample_index) = self.header.sample_index(sample_name) {
            if let Some(x) = self.genotype_mut(sample_name, key) {
                let mut ret = Vec::new();
                ret.append(x);
                x.append(&mut values);
                Some(ret)
            } else {
                if let Some(format_index) = self.format_index.get(key) {
                    while self.genotype[sample_index].len() < *format_index {
                        self.genotype[sample_index].push(vec![b".".to_vec()]);
                    }
                } else {
                    self.format.push(key.to_vec());
                    self.format_index
                        .insert(key.to_vec(), self.format.len() - 1);
                    while self.genotype[sample_index].len() < self.format.len() - 1 {
                        self.genotype[sample_index].push(vec![b".".to_vec()]);
                    }
                }
                self.genotype[sample_index].push(values);
                None
            }
        } else {
            None
        }
    }

    /// Recreate info and genotype index cache.
    /// Please call this method if you modify info and format field manually.
    pub fn recreate_info_and_genotype_index(&mut self) {
        // create info_index
        for v in self.info_index.values_mut() {
            *v = NOT_FOUND;
        }
        for (i, k) in self.info.iter().enumerate() {
            if let Some(x) = self.info_index.get_mut(&k.0) {
                *x = i;
            } else {
                self.info_index.insert(k.0.to_vec(), i);
            }
        }

        // create format_index
        for v in self.format_index.values_mut() {
            *v = NOT_FOUND;
        }
        for (i, k) in self.format.iter().enumerate() {
            if let Some(x) = self.format_index.get_mut(k) {
                *x = i;
            } else {
                self.format_index.insert(k.to_vec(), i);
            }
        }
    }
}

fn write_array(writer: &mut impl Write, array: &[Vec<u8>], delimiter: &[u8]) -> io::Result<()> {
    if array.is_empty() {
        writer.write_all(b".")?;
    } else {
        for (i, one) in array.iter().enumerate() {
            if i != 0 {
                writer.write_all(delimiter)?;
            }
            writer.write_all(one)?;
        }
    }

    Ok(())
}

fn write_info(writer: &mut impl Write, info: &[(U8Vec, Vec<U8Vec>)]) -> io::Result<()> {
    if info.is_empty() {
        writer.write_all(b".")?;
    } else {
        for (i, (k, v)) in info.iter().enumerate() {
            if i != 0 {
                writer.write_all(b";")?;
            }
            writer.write_all(k)?;
            if !v.is_empty() {
                writer.write_all(b"=")?;
                for (j, x) in v.iter().enumerate() {
                    if j != 0 {
                        writer.write_all(b",")?;
                    }
                    writer.write_all(x)?;
                }
            }
        }
    }

    Ok(())
}

impl VCFRecord {
    pub fn write_record<W: Write>(&self, mut writer: W) -> io::Result<()> {
        writer.write_all(&self.chromosome)?;
        writer.write_all(b"\t")?;
        write!(writer, "{}\t", self.position)?;
        write_array(&mut writer, &self.id, b",")?;
        writer.write_all(b"\t")?;
        writer.write_all(&self.reference)?;
        writer.write_all(b"\t")?;
        write_array(&mut writer, &self.alternative, b",")?;
        writer.write_all(b"\t")?;
        if let Some(qual) = self.qual.as_ref() {
            if (qual.round() - *qual).abs() < 0.000_000_01 {
                write!(writer, "{:.1}", qual)?;
            } else {
                write!(writer, "{}", qual)?;
            }
        } else {
            writer.write_all(b".")?;
        }
        writer.write_all(b"\t")?;
        write_array(&mut writer, &self.filter, b",")?;
        writer.write_all(b"\t")?;
        write_info(&mut writer, &self.info)?;
        if !self.format.is_empty() {
            writer.write_all(b"\t")?;
            write_array(&mut writer, &self.format, b":")?;
            for one_genotype in self.genotype.iter() {
                writer.write_all(b"\t")?;
                for (i, v) in one_genotype.iter().enumerate() {
                    if i != 0 {
                        writer.write_all(b":")?;
                    }
                    write_array(&mut writer, v, b",")?;
                }
            }
        }

        writer.write_all(b"\n")?;
        Ok(())
    }
}

#[cfg(test)]
mod test;

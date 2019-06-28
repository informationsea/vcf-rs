use super::{VCFHeader, VCFRecord};

use std::io;


pub struct VCFWriter<W: io::Write> {
    writer: W,
    header: VCFHeader,
}

impl<W: io::Write> VCFWriter<W> {
    pub fn new(mut writer: W, header: VCFHeader) -> io::Result<VCFWriter<W>> {
        for one in &header.items {
            writeln!(writer, "{}", one.line)?;
        }
        write!(
            writer,
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT"
        )?;
        for one in &header.samples {
            write!(writer, "\t{}", one)?;
        }
        writeln!(writer)?;

        Ok(VCFWriter { writer, header })
    }

    pub fn write_record(&mut self, record: &VCFRecord) -> io::Result<()> {
        record.write_line(&mut self.writer, &self.header.samples)
    }
}

#[cfg(test)]
mod test;
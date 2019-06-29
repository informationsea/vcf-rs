use super::{VCFHeader, VCFRecord};

use std::io;

/// Write VCF file as plain text.
/// If you want to write VCF file as compressed VCF,
/// please use [bgzip-rs](https://crates.io/crates/bgzip) to compress.
/// # Examples
/// ## Filter existing VCF file
///
/// ```
/// use std::io;
/// use vcf::*;
/// use std::fs::File;
/// use flate2::read::MultiGzDecoder;
///
/// # fn test_writer_filter() -> Result<(), VCFParseError> {
/// let vcf_reader = VCFReader::new(MultiGzDecoder::new(
///     File::open("testfiles/ALL.chr20.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.20-34001094-34168504-subset.vcf.gz")?
/// ))?;
/// let mut header = vcf_reader.header().clone();
/// header.items.push(VCFHeaderLine::new("##vcffilter=filter out odd")?);
/// let mut writer = VCFWriter::new(File::create("target/test-vcf-doc.vcf")?, header)?;
/// for one in vcf_reader {
///     let one_record = one?;
///     if one_record.position % 2 == 0 {
///         writer.write_record(&one_record)?;
///     }
/// }
/// #    Ok(())
/// # }
/// ```
///
/// ## Create VCF file from scrach
/// ```
/// use std::io;
/// use vcf::*;
/// use indexmap::indexmap;
/// use std::str;
///
/// # fn test_writer_from_scrach() -> io::Result<()> {
/// let mut write_buf: Vec<u8> = Vec::new();
/// let header = VCFHeader {
///     items: vec![
///         VCFHeaderLine::new("##fileformat=VCFv4.3").unwrap(),
///         VCFHeaderLine::new("##INFO=<ID=AF,Number=A,Type=Float,Description=\"Estimated allele frequency in the range (0,1)\">").unwrap(),
///         VCFHeaderLine::new("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">").unwrap(),
///     ],
///     samples: vec!["sampleA".to_string(), "sampleB".to_string()],
/// };
/// let mut writer = VCFWriter::new(&mut write_buf, header)?;
///
/// let record1 = VCFRecord {
///     chromosome: "10".to_string(),
///     position: 1234,
///     id: vec![],
///     reference: "A".to_string(),
///     alternative: vec!["T".to_string(), "G".to_string()],
///     filter: vec![],
///     quality: None,
///     info: indexmap!("AF".to_string() => vec!["0.2".to_string()]),
///     format: vec!["GT".to_string()],
///     call: vec![(
///         "sampleA".to_string(),
///         vec![("GT".to_string(), vec!["0|1".to_string()])]
///             .into_iter()
///             .collect(),
///     )]
///     .into_iter()
///     .collect(),
/// };
///
/// writer.write_record(&record1)?;
///
/// assert_eq!(
///     r#"##fileformat=VCFv4.3
/// ##INFO=<ID=AF,Number=A,Type=Float,Description="Estimated allele frequency in the range (0,1)">
/// ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
/// #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sampleA	sampleB
/// 10	1234	.	A	T,G	.	.	AF=0.2	GT	0|1	.
/// "#,
///     str::from_utf8(&write_buf).unwrap()
/// );
///
/// #    Ok(())
/// # }
/// ```
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

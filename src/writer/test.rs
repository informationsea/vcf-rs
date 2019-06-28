use super::*;
use crate::reader::VCFReader;
use std::fs::File;
use flate2::read::MultiGzDecoder;
use std::io::prelude::*;

#[test]
fn test_writer() {
    let mut vcf_reader = VCFReader::new(MultiGzDecoder::new(
        File::open("testfiles/ALL.chr20.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.20-34001094-34168504-subset.vcf.gz")
            .unwrap()
    )).unwrap();

    let records: Vec<_> = vcf_reader.iter().collect();

    let mut raw_data = Vec::new();
    MultiGzDecoder::new(
        File::open("testfiles/ALL.chr20.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.20-34001094-34168504-subset.vcf.gz")
            .unwrap()
    ).read_to_end(&mut raw_data).unwrap();

    let mut write_buf: Vec<u8> = Vec::new();
    let mut writer = VCFWriter::new(&mut write_buf, vcf_reader.header().clone());
}
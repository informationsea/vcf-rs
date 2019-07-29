#![feature(test)]

extern crate test;
use std::fs::File;
use flate2::read::MultiGzDecoder;
use vcf::VCFReader;

#[bench]
fn load_vcf(b: &mut test::Bencher) {
    b.iter(|| {
        let mut vcf_reader = VCFReader::new(
            MultiGzDecoder::new(
                File::open("testfiles/ALL.chr20.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.20-34001094-34168504-subset.vcf.gz")
                    .unwrap())).unwrap();
        for _one in vcf_reader {
            // ignore
        }
    })
}
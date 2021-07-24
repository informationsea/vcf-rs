vcf-rs
======
[![Build](https://github.com/informationsea/vcf-rs/actions/workflows/build.yml/badge.svg)](https://github.com/informationsea/vcf-rs/actions/workflows/build.yml)
[![Creates.io](http://meritbadge.herokuapp.com/vcf)](https://crates.io/crates/vcf)
[![doc-rs](https://docs.rs/vcf/badge.svg)](https://docs.rs/vcf)
[![GitHub](https://img.shields.io/github/license/informationsea/vcf-rs)](https://github.com/informationsea/vcf-rs)
[![GitHub top language](https://img.shields.io/github/languages/top/informationsea/vcf-rs)](https://github.com/informationsea/vcf-rs)

Rust implementation of VCF parser.

License
-------
Apache 2.0


Example
-------

```rust
use vcf::*;
use flate2::read::MultiGzDecoder;
use std::fs::File;

fn usage_test() -> Result<(), VCFError> {
    let mut reader = VCFReader::new(BufReader::new(MultiGzDecoder::new(File::open(
        "./testfiles/NA12878-subset.vcf.gz",
    )?)))?;

    // access FILTER contents
    assert_eq!(
        Some(VCFHeaderFilterAlt {
            id: b"PASS",
            description: b"All filters passed"
        }),
        reader.header().filter(b"PASS")
    );

    // access INFO contents
    assert_eq!(
        b"Stop position of the interval",
        reader.header().info(b"END").unwrap().description
    );

    // prepare VCFRecord object
    let mut vcf_record = VCFRecord::new(reader.header());

    // read one record
    reader.next_record(&mut vcf_record)?;

    // get record attributes
    assert_eq!(vcf_record.chromosome, b"13");
    assert_eq!(vcf_record.position, 32889968);
    assert_eq!(vcf_record.id, Vec::<U8Vec>::new());
    assert_eq!(vcf_record.reference, b"G");
    assert_eq!(vcf_record.alternative, vec![b"A"]);
    assert_eq!(vcf_record.qual, Some(25743.5));
    assert_eq!(vcf_record.info(b"AC"), Some(&vec![b"54".to_vec()]));
    assert_eq!(
        vcf_record.genotype(b"ERP001775_HiSeq2000_SAMEA1531955-1", b"GT"),
        Some(&vec![b"1/1".to_vec()])
    );
    assert_eq!(
        vcf_record.genotype(b"ERP001775_HiSeq2000_SAMEA1531955-1", b"AD"),
        Some(&vec![b"0".to_vec(), b"14".to_vec()])
    );

    Ok(())
}
```
use super::*;

use std::io::BufReader;

#[test]
fn test_reader1() -> Result<(), VCFError> {
    let mut simple_vcf = BufReader::new(&include_bytes!("../testfiles/simple1.vcf")[..]);
    let mut vcf_reader = VCFReader::new(&mut simple_vcf).unwrap();
    let mut vcf_record = vcf_reader.empty_record();
    let mut record_count = 0;
    while vcf_reader.next_record(&mut vcf_record)? {
        record_count += 1;
    }
    assert_eq!(record_count, 3);
    Ok(())
}

#[test]
fn test_reader2() -> Result<(), VCFError> {
    let mut simple_vcf = BufReader::new(&include_bytes!("../testfiles/1kGP-subset.vcf")[..]);
    let mut vcf_reader = VCFReader::new(&mut simple_vcf).unwrap();
    let mut vcf_record = vcf_reader.empty_record();
    let mut record_count = 0;
    while vcf_reader.next_record(&mut vcf_record)? {
        record_count += 1;
    }
    assert_eq!(record_count, 306);
    Ok(())
}

#[test]
fn test_reader3() -> Result<(), VCFError> {
    let mut simple_vcf = BufReader::new(&include_bytes!("../testfiles/bad1.vcf")[..]);
    let mut vcf_reader = VCFReader::new(&mut simple_vcf).unwrap();
    let mut vcf_record = vcf_reader.empty_record();
    let mut record_count = 0;
    while vcf_reader.next_record(&mut vcf_record)? {
        record_count += 1;
    }
    assert_eq!(record_count, 2);
    Ok(())
}

#[test]
fn test_writer() -> Result<(), VCFError> {
    let vcf_bytes = include_bytes!("../testfiles/simple1.vcf");
    let mut simple_vcf = BufReader::new(&vcf_bytes[..]);
    let mut vcf_reader = VCFReader::new(&mut simple_vcf)?;
    let mut vcf_record = vcf_reader.empty_record();
    let mut buffer = Vec::new();
    let mut vcf_writer = VCFWriter::new(&mut buffer, &vcf_reader.header())?;
    while vcf_reader.next_record(&mut vcf_record)? {
        vcf_writer.write_record(&vcf_record)?;
    }
    assert_eq!(&vcf_bytes[..], &buffer[..]);

    Ok(())
}

use flate2::read::MultiGzDecoder;
use std::fs::File;

#[test]
#[allow(clippy::unreadable_literal)]
fn usage_test() -> Result<(), VCFError> {
    let mut reader = VCFReader::new(BufReader::new(MultiGzDecoder::new(File::open(
        "./testfiles/1kGP-subset.vcf.gz",
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
    let mut vcf_record = reader.empty_record();

    // read one record
    reader.next_record(&mut vcf_record)?;

    // get record attributes
    assert_eq!(vcf_record.chromosome, b"13");
    assert_eq!(vcf_record.position, 32872836);
    assert_eq!(vcf_record.id, Vec::<U8Vec>::new());
    assert_eq!(vcf_record.reference, b"A");
    assert_eq!(vcf_record.alternative, vec![b"C"]);
    assert_eq!(vcf_record.qual, Some(495.23));
    assert_eq!(vcf_record.info(b"AC"), Some(&vec![b"1".to_vec()]));
    assert_eq!(
        vcf_record.genotype(b"SRP150637__HG00099", b"GT"),
        Some(&vec![b"0/0".to_vec()])
    );
    assert_eq!(
        vcf_record.genotype(b"SRP150637__HG00099", b"AD"),
        Some(&vec![b"31".to_vec(), b"0".to_vec()])
    );

    Ok(())
}

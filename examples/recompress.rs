use clap::{App, Arg};
use flate2::read::MultiGzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression;
use std::boxed::Box;
use std::fs::File;
use std::io::{BufReader, Read};

pub fn main() -> Result<(), vcf::VCFError> {
    let matches = App::new("Recompress VCF")
        .version("0.1")
        .author("Yasunobu Okamura")
        .about("Parse and recompress VCF file")
        .arg(
            Arg::with_name("input")
                .index(1)
                .takes_value(true)
                .required(true)
                .help("Input VCF"),
        )
        .arg(
            Arg::with_name("output")
                .long("output")
                .short("o")
                .takes_value(true)
                .required(true)
                .help("Output VCF"),
        )
        .get_matches();
    let input_vcf_path = matches.value_of("input").unwrap();
    let output_vcf_path = matches.value_of("output").unwrap();

    let reader: BufReader<Box<dyn Read>> = BufReader::new(
        if input_vcf_path.ends_with(".gz") || input_vcf_path.ends_with(".bgz") {
            Box::new(MultiGzDecoder::new(File::open(input_vcf_path)?))
        } else {
            Box::new(File::open(input_vcf_path)?)
        },
    );
    let mut vcf_reader = vcf::VCFReader::new(reader)?;
    let mut record = vcf_reader.empty_record();

    let writer = GzEncoder::new(File::create(output_vcf_path)?, Compression::default());
    let mut vcf_writer = vcf::VCFWriter::new(writer, &vcf_reader.header())?;

    while let Ok(_) = vcf_reader.next_record(&mut record) {
        vcf_writer.write_record(&record)?;
    }

    Ok(())
}

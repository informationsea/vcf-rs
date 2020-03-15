use super::*;
use std::collections::HashSet;
use std::io::BufReader;

#[test]
fn test_header_item_parser() {
    assert!(parse_header_item(b"#CHROM").is_err());
    assert_eq!(
        parse_header_item(b"##fileformat=VCFv4.2\n").unwrap().1,
        VCFHeaderLine {
            contents: VCFHeaderContent::FileFormat(VCFVersion::Vcf4_2),
            line: b"##fileformat=VCFv4.2\n".to_vec()
        }
    );
}

#[test]
fn test_parse_number() {
    assert_eq!(parser::parse_number(b"A"), Number::Allele);
    assert_eq!(parser::parse_number(b"R"), Number::Reference);
    assert_eq!(parser::parse_number(b"G"), Number::Genotype);
    assert_eq!(parser::parse_number(b"0"), Number::Zero);
    assert_eq!(parser::parse_number(b"."), Number::Unknown);
    assert_eq!(parser::parse_number(b"1"), Number::Number(1));
    assert_eq!(parser::parse_number(b"10"), Number::Number(10));
    assert_eq!(parser::parse_number(b"100"), Number::Number(100));
    assert_eq!(parser::parse_number(b"X"), Number::Other(b"X".to_vec()));
}

#[test]
fn test_parse_value_type() {
    assert_eq!(parser::parse_value_type(b"String"), ValueType::String);
    assert_eq!(parser::parse_value_type(b"Integer"), ValueType::Integer);
    assert_eq!(parser::parse_value_type(b"Flag"), ValueType::Flag);
    assert_eq!(parser::parse_value_type(b"Character"), ValueType::Character);
    assert_eq!(parser::parse_value_type(b"Float"), ValueType::Float);
    assert_eq!(
        parser::parse_value_type(b"XXX"),
        ValueType::Other(b"XXX".to_vec())
    );
}

#[test]
fn test_parse_vcf_file_format_header() {
    assert_eq!(
        parser::parse_vcf_file_format_header(b"fileformat=VCFv4.0").unwrap(),
        (&b""[..], VCFHeaderContent::FileFormat(VCFVersion::Vcf4_0))
    );
    assert_eq!(
        parser::parse_vcf_file_format_header(b"fileformat=VCFv4.1").unwrap(),
        (&b""[..], VCFHeaderContent::FileFormat(VCFVersion::Vcf4_1))
    );
    assert_eq!(
        parser::parse_vcf_file_format_header(b"fileformat=VCFv4.2").unwrap(),
        (&b""[..], VCFHeaderContent::FileFormat(VCFVersion::Vcf4_2))
    );
    assert_eq!(
        parser::parse_vcf_file_format_header(b"fileformat=VCFv4.3").unwrap(),
        (&b""[..], VCFHeaderContent::FileFormat(VCFVersion::Vcf4_3))
    );
    assert_eq!(
        parser::parse_vcf_file_format_header(b"fileformat=VCFv3.0").unwrap(),
        (
            &b""[..],
            VCFHeaderContent::FileFormat(VCFVersion::Other(b"VCFv3.0".to_vec()))
        )
    );
}

#[test]
fn test_parse_header_entries() {
    assert_eq!(
        parser::parse_header_entries(b"ID=1,length=249250621>").unwrap(),
        (
            &b">"[..],
            vec![(&b"ID"[..], &b"1"[..]), (&b"length"[..], &b"249250621"[..]),]
        )
    );

    assert_eq!(
        parser::parse_header_entries(
            b"ID=MQ,Number=1,Type=Float,Description=\"RMS Mapping Quality\">"
        )
        .unwrap(),
        (
            &b">"[..],
            vec![
                (&b"ID"[..], &b"MQ"[..]),
                (&b"Number"[..], &b"1"[..]),
                (&b"Type"[..], &b"Float"[..]),
                (&b"Description"[..], &b"RMS Mapping Quality"[..])
            ]
        )
    );
}

#[test]
#[allow(clippy::unreadable_literal)]
fn test_parse_vcf_contig_header() {
    assert!(parser::parse_vcf_contig_header(b"contig=<length=1234>").is_err());
    assert_eq!(
        parser::parse_header_content(b"contig=<ID=chr1>").unwrap(),
        (
            &b""[..],
            VCFHeaderContent::Contig {
                id: b"chr1".to_vec(),
                length: None
            }
        )
    );

    assert_eq!(
        parser::parse_header_content(b"contig=<ID=chr1,length=249250621>").unwrap(),
        (
            &b""[..],
            VCFHeaderContent::Contig {
                id: b"chr1".to_vec(),
                length: Some(249250621)
            }
        )
    )
}

#[test]
#[allow(clippy::unreadable_literal)]
fn test_parse_vcf_alt_header() {
    assert!(parser::parse_vcf_alt_header(b"ALT=<ID=foo>").is_err());
    assert!(parser::parse_vcf_alt_header(b"ALT=<Description=\"bar\">").is_err());
    assert_eq!(
        parser::parse_header_content(b"ALT=<ID=foo,Description=\"bar\">").unwrap(),
        (
            &b""[..],
            VCFHeaderContent::ALT {
                id: b"foo".to_vec(),
                description: b"bar".to_vec(),
            }
        )
    );
}

#[test]
#[allow(clippy::unreadable_literal)]
fn test_parse_vcf_filter_header() {
    assert!(parser::parse_vcf_filter_header(b"ALT=<ID=foo>").is_err());
    assert!(parser::parse_vcf_filter_header(b"ALT=<Description=\"bar\">").is_err());
    assert_eq!(
        parser::parse_header_content(b"ALT=<ID=foo,Description=\"bar\">").unwrap(),
        (
            &b""[..],
            VCFHeaderContent::ALT {
                id: b"foo".to_vec(),
                description: b"bar".to_vec(),
            }
        )
    );
}

#[test]
fn test_parse_vcf_info_header() {
    assert!(parser::parse_vcf_info_header(
        b"INFO=<Number=0,Type=Flag,Description=\"Were any of the samples downsampled?\">"
    )
    .is_err());
    assert!(parser::parse_vcf_info_header(
        b"INFO=<ID=DS,Type=Flag,Description=\"Were any of the samples downsampled?\">"
    )
    .is_err());
    assert!(parser::parse_vcf_info_header(
        b"INFO=<ID=DS,Number=0,Description=\"Were any of the samples downsampled?\">"
    )
    .is_err());
    assert!(parser::parse_vcf_info_header(b"INFO=<ID=DS,Number=0,Type=Flag>").is_err());

    assert_eq!(
        parser::parse_header_content(
            b"INFO=<ID=DS,Number=0,Type=Flag,Description=\"Were any of the samples downsampled?\">"
        )
        .unwrap(),
        (
            &b""[..],
            VCFHeaderContent::INFO {
                id: b"DS".to_vec(),
                number: Number::Zero,
                value_type: ValueType::Flag,
                description: b"Were any of the samples downsampled?".to_vec(),
                source: None,
                version: None,
            }
        )
    );

    assert_eq!(
        parser::parse_header_content(
            b"INFO=<ID=DS,Number=0,Type=Flag,Description=\"Were any of the samples downsampled?\",Source=\"GATK\",Version=4>"
        )
        .unwrap(),
        (
            &b""[..],
            VCFHeaderContent::INFO {
                id: b"DS".to_vec(),
                number: Number::Zero,
                value_type: ValueType::Flag,
                description: b"Were any of the samples downsampled?".to_vec(),
                source: Some(b"GATK".to_vec()),
                version: Some(b"4".to_vec()),
            }
        )
    );
}

#[test]
fn test_parse_vcf_format_header() {
    assert!(parser::parse_vcf_file_format_header(
        b"INFO=<Number=0,Type=Flag,Description=\"Were any of the samples downsampled?\">"
    )
    .is_err());
    assert!(parser::parse_vcf_file_format_header(
        b"INFO=<ID=DS,Type=Flag,Description=\"Were any of the samples downsampled?\">"
    )
    .is_err());
    assert!(parser::parse_vcf_file_format_header(
        b"INFO=<ID=DS,Number=0,Description=\"Were any of the samples downsampled?\">"
    )
    .is_err());
    assert!(parser::parse_vcf_file_format_header(b"INFO=<ID=DS,Number=0,Type=Flag>").is_err());

    assert_eq!(
        parser::parse_header_content(
            b"FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
        )
        .unwrap(),
        (
            &b""[..],
            VCFHeaderContent::FORMAT {
                id: b"GT".to_vec(),
                number: Number::Number(1),
                value_type: ValueType::String,
                description: b"Genotype".to_vec(),
                source: None,
                version: None,
            }
        )
    );

    assert_eq!(
        parser::parse_header_content(
            b"FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\",Source=GATK,Version=3>"
        )
        .unwrap(),
        (
            &b""[..],
            VCFHeaderContent::FORMAT {
                id: b"GT".to_vec(),
                number: Number::Number(1),
                value_type: ValueType::String,
                description: b"Genotype".to_vec(),
                source: Some(b"GATK".to_vec()),
                version: Some(b"3".to_vec()),
            }
        )
    );
}

#[test]
fn test_parse_samples() {
    assert_eq!(
        Vec::<U8Vec>::new(),
        parser::parse_samples(b"#CHROM\tPOS\tID\tREF\tALT\n")
            .unwrap()
            .1
    );

    assert_eq!(
        Vec::<U8Vec>::new(),
        parser::parse_samples(b"#CHROM\tPOS\tID\tREF\tALT\tQUAL\n")
            .unwrap()
            .1
    );

    assert_eq!(
        Vec::<U8Vec>::new(),
        parser::parse_samples(b"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\n")
            .unwrap()
            .1
    );

    assert_eq!(
        Vec::<U8Vec>::new(),
        parser::parse_samples(b"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
            .unwrap()
            .1
    );

    assert_eq!(
        Vec::<U8Vec>::new(),
        parser::parse_samples(b"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\n")
            .unwrap()
            .1
    );

    assert_eq!(
        vec![b"A".to_vec()],
        parser::parse_samples(b"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\n")
            .unwrap()
            .1
    );

    assert_eq!(
        vec![b"A".to_vec(), b"BCD".to_vec()],
        parser::parse_samples(b"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tBCD\n")
            .unwrap()
            .1
    );
}

#[test]
#[allow(clippy::unreadable_literal)]
fn test_parse_header() -> Result<(), VCFError> {
    let mut simple1_vcf = BufReader::new(&include_bytes!("../../testfiles/simple1.vcf")[..]);

    let (line_count, unprocessed_line, header) = parse_header(&mut simple1_vcf)?;
    assert_eq!(line_count, 20);
    assert_eq!(unprocessed_line, None);

    assert_eq!(
        header.samples,
        vec![
            b"ERP001775_HiSeq2000_SAMEA1531955-1".to_vec(),
            b"ERP001775_HiSeq2000_SAMEA1531955-2".to_vec()
        ]
    );

    let expected_items = vec![
            VCFHeaderLine {
                line: b"##fileformat=VCFv4.2\n".to_vec(),
                contents: VCFHeaderContent::FileFormat(VCFVersion::Vcf4_2)
            },
            VCFHeaderLine {
                line: b"##FILTER=<ID=PASS,Description=\"All filters passed\">\n".to_vec(),
                contents: VCFHeaderContent::FILTER{
                    id: b"PASS".to_vec(),
                    description: b"All filters passed".to_vec()
                }
            },
            VCFHeaderLine {
                line: b"##ALT=<ID=NON_REF,Description=\"Represents any possible alternative allele at this location\">\n".to_vec(),
                contents: VCFHeaderContent::ALT{
                    id: b"NON_REF".to_vec(),
                    description: b"Represents any possible alternative allele at this location".to_vec()
                }
            },
            VCFHeaderLine {
                line: b"##FILTER=<ID=LowQual,Description=\"Low quality\">\n".to_vec(),
                contents: VCFHeaderContent::FILTER{
                    id: b"LowQual".to_vec(),
                    description: b"Low quality".to_vec()
                }
            },
            VCFHeaderLine {
                line: b"##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed\">\n".to_vec(),
                contents: VCFHeaderContent::FORMAT{
                    id: b"AD".to_vec(),
                    number: Number::Reference,
                    value_type: ValueType::Integer,
                    description: b"Allelic depths for the ref and alt alleles in the order listed".to_vec(),
                    source: None,
                    version: None,
                }
            },
            VCFHeaderLine {
                line: b"##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Approximate read depth (reads with MQ=255 or with bad mates are filtered)\">\n".to_vec(),
                contents: VCFHeaderContent::FORMAT{
                    id: b"DP".to_vec(),
                    number: Number::Number(1),
                    value_type: ValueType::Integer,
                    description: b"Approximate read depth (reads with MQ=255 or with bad mates are filtered)".to_vec(),
                    source: None,
                    version: None,
                }
            },
            VCFHeaderLine {
                line: b"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n".to_vec(),
                contents: VCFHeaderContent::FORMAT{
                    id: b"GT".to_vec(),
                    number: Number::Number(1),
                    value_type: ValueType::String,
                    description: b"Genotype".to_vec(),
                    source: None,
                    version: None,
                }
            },
            VCFHeaderLine {
                line: b"##GATKCommandLine=<ID=HaplotypeCaller,CommandLine=\"HaplotypeCaller  --contamination-fraction-to-filter 0.0 --emit-ref-confidence GVCF --output ERP001775_HiSeq2000_SAMEA1531955-1.hs37d5.g.vcf.gz --intervals /cromwell-executions/VariantCall/6f09e738-2160-49dd-8287-ff8ef9368d27/call-HaplotypeCaller/shard-20/inputs/1101158470/0020-scattered.interval_list --input /cromwell-executions/VariantCall/6f09e738-2160-49dd-8287-ff8ef9368d27/call-HaplotypeCaller/shard-20/inputs/-929370236/ERP001775_HiSeq2000_SAMEA1531955-1.hs37d5.bam --reference /cromwell-executions/VariantCall/6f09e738-2160-49dd-8287-ff8ef9368d27/call-HaplotypeCaller/shard-20/inputs/865204270/hs37d5.fa  --use-new-qual-calculator true --use-old-qual-calculator false --annotate-with-num-discovered-alleles false --heterozygosity 0.001 --indel-heterozygosity 1.25E-4 --heterozygosity-stdev 0.01 --standard-min-confidence-threshold-for-calling 30.0 --max-alternate-alleles 6 --max-genotype-count 1024 --sample-ploidy 2 --num-reference-samples-if-no-call 0 --genotyping-mode DISCOVERY --genotype-filtered-alleles false --output-mode EMIT_VARIANTS_ONLY --all-site-pls false --gvcf-gq-bands 1 --gvcf-gq-bands 2 --gvcf-gq-bands 3 --gvcf-gq-bands 4 --gvcf-gq-bands 5 --gvcf-gq-bands 6 --gvcf-gq-bands 7 --gvcf-gq-bands 8 --gvcf-gq-bands 9 --gvcf-gq-bands 10 --gvcf-gq-bands 11 --gvcf-gq-bands 12 --gvcf-gq-bands 13 --gvcf-gq-bands 14 --gvcf-gq-bands 15 --gvcf-gq-bands 16 --gvcf-gq-bands 17 --gvcf-gq-bands 18 --gvcf-gq-bands 19 --gvcf-gq-bands 20 --gvcf-gq-bands 21 --gvcf-gq-bands 22 --gvcf-gq-bands 23 --gvcf-gq-bands 24 --gvcf-gq-bands 25 --gvcf-gq-bands 26 --gvcf-gq-bands 27 --gvcf-gq-bands 28 --gvcf-gq-bands 29 --gvcf-gq-bands 30 --gvcf-gq-bands 31 --gvcf-gq-bands 32 --gvcf-gq-bands 33 --gvcf-gq-bands 34 --gvcf-gq-bands 35 --gvcf-gq-bands 36 --gvcf-gq-bands 37 --gvcf-gq-bands 38 --gvcf-gq-bands 39 --gvcf-gq-bands 40 --gvcf-gq-bands 41 --gvcf-gq-bands 42 --gvcf-gq-bands 43 --gvcf-gq-bands 44 --gvcf-gq-bands 45 --gvcf-gq-bands 46 --gvcf-gq-bands 47 --gvcf-gq-bands 48 --gvcf-gq-bands 49 --gvcf-gq-bands 50 --gvcf-gq-bands 51 --gvcf-gq-bands 52 --gvcf-gq-bands 53 --gvcf-gq-bands 54 --gvcf-gq-bands 55 --gvcf-gq-bands 56 --gvcf-gq-bands 57 --gvcf-gq-bands 58 --gvcf-gq-bands 59 --gvcf-gq-bands 60 --gvcf-gq-bands 70 --gvcf-gq-bands 80 --gvcf-gq-bands 90 --gvcf-gq-bands 99 --floor-blocks false --indel-size-to-eliminate-in-ref-model 10 --use-alleles-trigger false --disable-optimizations false --just-determine-active-regions false --dont-genotype false --do-not-run-physical-phasing false --use-filtered-reads-for-annotations false --correct-overlapping-quality false --adaptive-pruning false --do-not-recover-dangling-branches false --recover-dangling-heads false --consensus false --dont-trim-active-regions false --max-disc-ar-extension 25 --max-gga-ar-extension 300 --padding-around-indels 150 --padding-around-snps 20 --kmer-size 10 --kmer-size 25 --dont-increase-kmer-sizes-for-cycles false --allow-non-unique-kmers-in-ref false --num-pruning-samples 1 --min-dangling-branch-length 4 --recover-all-dangling-branches false --max-num-haplotypes-in-population 128 --min-pruning 2 --adaptive-pruning-initial-error-rate 0.001 --pruning-lod-threshold 2.302585092994046 --max-unpruned-variants 100 --debug-assembly false --debug-graph-transformations false --capture-assembly-failure-bam false --error-correct-reads false --kmer-length-for-read-error-correction 25 --min-observations-for-kmer-to-be-solid 20 --likelihood-calculation-engine PairHMM --base-quality-score-threshold 18 --pair-hmm-gap-continuation-penalty 10 --pair-hmm-implementation FASTEST_AVAILABLE --pcr-indel-model CONSERVATIVE --phred-scaled-global-read-mismapping-rate 45 --native-pair-hmm-threads 4 --native-pair-hmm-use-double-precision false --bam-writer-type CALLED_HAPLOTYPES --dont-use-soft-clipped-bases false --min-base-quality-score 10 --smith-waterman JAVA --max-mnp-distance 0 --min-assembly-region-size 50 --max-assembly-region-size 300 --assembly-region-padding 100 --max-reads-per-alignment-start 50 --active-probability-threshold 0.002 --max-prob-propagation-distance 50 --force-active false --interval-set-rule UNION --interval-padding 0 --interval-exclusion-padding 0 --interval-merging-rule ALL --read-validation-stringency SILENT --seconds-between-progress-updates 10.0 --disable-sequence-dictionary-validation false --create-output-bam-index true --create-output-bam-md5 false --create-output-variant-index true --create-output-variant-md5 false --lenient false --add-output-sam-program-record true --add-output-vcf-command-line true --cloud-prefetch-buffer 40 --cloud-index-prefetch-buffer -1 --disable-bam-index-caching false --sites-only-vcf-output false --help false --version false --showHidden false --verbosity INFO --QUIET false --use-jdk-deflater false --use-jdk-inflater false --gcs-max-retries 20 --gcs-project-for-requester-pays  --disable-tool-default-read-filters false --minimum-mapping-quality 20 --disable-tool-default-annotations false --enable-all-annotations false --allow-old-rms-mapping-quality-annotation-data false\",Version=\"4.1.3.0\",Date=\"October 3, 2019 5:19:41 PM UTC\">\n".to_vec(),
                contents: VCFHeaderContent::Other
            },
            VCFHeaderLine {
                line: b"##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count in genotypes, for each ALT allele, in the same order as listed\">\n".to_vec(),
                contents: VCFHeaderContent::INFO{
                    id: b"AC".to_vec(),
                    number: Number::Allele,
                    value_type: ValueType::Integer,
                    description: b"Allele count in genotypes, for each ALT allele, in the same order as listed".to_vec(),
                    source: None,
                    version: None,
                }
            },
            VCFHeaderLine {
                line: b"##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency, for each ALT allele, in the same order as listed\">\n".to_vec(),
                contents: VCFHeaderContent::INFO{
                    id: b"AF".to_vec(),
                    number: Number::Allele,
                    value_type: ValueType::Float,
                    description: b"Allele Frequency, for each ALT allele, in the same order as listed".to_vec(),
                    source: None,
                    version: None,
                }
            },
            VCFHeaderLine {
                line: b"##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles in called genotypes\">\n".to_vec(),
                contents: VCFHeaderContent::INFO{
                    id: b"AN".to_vec(),
                    number: Number::Number(1),
                    value_type: ValueType::Integer,
                    description: b"Total number of alleles in called genotypes".to_vec(),
                    source: None,
                    version: None,
                }
            },
            VCFHeaderLine {
                line: b"##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Approximate read depth; some reads may have been filtered\">\n".to_vec(),
                contents: VCFHeaderContent::INFO{
                    id: b"DP".to_vec(),
                    number: Number::Number(1),
                    value_type: ValueType::Integer,
                    description: b"Approximate read depth; some reads may have been filtered".to_vec(),
                    source: None,
                    version: None,
                }
            },
            VCFHeaderLine {
                line: b"##contig=<ID=13,length=115169878>\n".to_vec(),
                contents: VCFHeaderContent::Contig {
                    id: b"13".to_vec(),
                    length: Some(115169878),
                }
            },
            VCFHeaderLine {
                line: b"##contig=<ID=14,length=107349540>\n".to_vec(),
                contents: VCFHeaderContent::Contig {
                    id: b"14".to_vec(),
                    length: Some(107349540),
                }
            },
            VCFHeaderLine {
                line: b"##source=CombineGVCFs\n".to_vec(),
                contents: VCFHeaderContent::Other
            },
            VCFHeaderLine {
                line: b"##source=GenotypeGVCFs\n".to_vec(),
                contents: VCFHeaderContent::Other
            },
            VCFHeaderLine {
                line: b"##source=HaplotypeCaller\n".to_vec(),
                contents: VCFHeaderContent::Other
            },
            VCFHeaderLine {
                line: b"##bcftools_viewVersion=1.9+htslib-1.9\n".to_vec(),
                contents: VCFHeaderContent::Other
            },
            VCFHeaderLine {
                line: b"##bcftools_viewCommand=view -O z -o /tmp/NA12878-subset.vcf.gz -r 13:32889150-32975410,17:41194315-41277931 NA12878_comparison.genotyped.vcf.gz; Date=Sat Mar  7 20:29:56 2020\n".to_vec(),
                contents: VCFHeaderContent::Other
            }
        ];

    for (x, y) in expected_items.iter().zip(header.items.iter()) {
        assert_eq!(x, y);
    }
    assert_eq!(expected_items.len(), header.items.len());

    assert_eq!(
        header.info(b"AC"),
        Some(VCFHeaderInfoFormat {
            id: b"AC",
            number: &Number::Allele,
            value_type: &ValueType::Integer,
            description:
                b"Allele count in genotypes, for each ALT allele, in the same order as listed",
            source: None,
            version: None,
        })
    );

    assert_eq!(header.info(b"XX"), None);
    assert_eq!(header.info_key.len(), 4);

    assert_eq!(
        header.format(b"AD"),
        Some(VCFHeaderInfoFormat {
            id: b"AD",
            number: &Number::Reference,
            value_type: &ValueType::Integer,
            description: b"Allelic depths for the ref and alt alleles in the order listed",
            source: None,
            version: None,
        })
    );
    assert_eq!(header.format(b"AC"), None);
    assert_eq!(header.format_key.len(), 3);
    assert_eq!(
        header
            .format_list()
            .map(|x| -> &[u8] { &x })
            .collect::<HashSet<&[u8]>>(),
        [&b"AD"[..], &b"DP"[..], &b"GT"[..]]
            .iter()
            .cloned()
            .collect::<HashSet<&[u8]>>()
    );

    assert_eq!(
        header.filter(b"LowQual"),
        Some(VCFHeaderFilterAlt {
            id: b"LowQual",
            description: b"Low quality",
        })
    );
    assert_eq!(header.filter(b"AC"), None);
    assert_eq!(header.filter_key.len(), 2);
    assert_eq!(
        header
            .filter_list()
            .map(|x| -> &[u8] { &x })
            .collect::<HashSet<&[u8]>>(),
        [&b"PASS"[..], &b"LowQual"[..]]
            .iter()
            .cloned()
            .collect::<HashSet<&[u8]>>()
    );

    assert_eq!(header.alt(b"DEL"), None);
    assert_eq!(header.alt_key.len(), 1);
    assert_eq!(
        header
            .alt_list()
            .map(|x| -> &[u8] { &x })
            .collect::<HashSet<&[u8]>>(),
        [&b"NON_REF"[..]]
            .iter()
            .cloned()
            .collect::<HashSet<&[u8]>>()
    );
    assert_eq!(
        header.alt(b"NON_REF"),
        Some(VCFHeaderFilterAlt {
            id: b"NON_REF",
            description: b"Represents any possible alternative allele at this location",
        })
    );

    Ok(())
}

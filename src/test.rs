use super::*;

#[test]
fn test_decode_vcf_value() {
    assert_eq!(decode_vcf_value("hoge_hoge"), "hoge_hoge");
    assert_eq!(decode_vcf_value("hoge%3Ahoge"), "hoge:hoge");
    assert_eq!(
        decode_vcf_value("hoge%3A%3B%3D%25%2C%0D%0A%09%xx"),
        "hoge:;=%,\r\n\t%xx"
    );
}

#[test]
fn test_encode_vcf_value() {
    assert_eq!(
        encode_vcf_value("hoge:;=%,\r\n\t"),
        "hoge%3A%3B%3D%25%2C%0D%0A%09"
    );
}

#[test]
fn test_vcf_header_parse_helper() {
    assert_eq!(vcf_header_parse_helper("ID=EAS_AF,Number=A,Type=Float,Description=\"Allele frequency in the EAS populations calculated from AC and AN, in the range (0,1)\""), 
        vec![("ID", "EAS_AF"), ("Number", "A"), ("Type", "Float"), ("Description", "Allele frequency in the EAS populations calculated from AC and AN, in the range (0,1)")].into_iter().collect()
    );
    assert_eq!(vcf_header_parse_helper("ID=EAS_AF,Number=A,Type=Float,Description=\"Allele frequency in the EAS populations calculated from AC and AN, in the range (0,1)\",Source=1000genomes,Version=phase3"), 
        vec![("ID", "EAS_AF"), ("Number", "A"), ("Type", "Float"), ("Description", "Allele frequency in the EAS populations calculated from AC and AN, in the range (0,1)"), ("Source", "1000genomes"), ("Version", "phase3")].into_iter().collect()
    );
    assert_eq!(
        vcf_header_parse_helper("ID=CN88,Description=\"Copy number allele: 88 copies\""),
        vec![
            ("ID", "CN88"),
            ("Description", "Copy number allele: 88 copies")
        ]
        .into_iter()
        .collect()
    );
    assert_eq!(
        vcf_header_parse_helper("ID=1,assembly=b37,length=249250621"),
        vec![("ID", "1"), ("assembly", "b37"), ("length", "249250621")]
            .into_iter()
            .collect()
    );
}

#[test]
fn test_vcf_header_line_helper() {
    assert_eq!(
        vcf_header_line_helper("##INFO=<ID=CIPOS,Number=2,Type=Integer,Description=\"Confidence interval around POS for imprecise variants\">").unwrap(), 
        (
            "INFO", "<ID=CIPOS,Number=2,Type=Integer,Description=\"Confidence interval around POS for imprecise variants\">",
            vec![
            ("ID", "CIPOS"),
            ("Number", "2"),
            ("Type", "Integer"),
            ("Description", "Confidence interval around POS for imprecise variants")
        ]
        .into_iter()
        .collect()
        )
    );

    assert_eq!(
        vcf_header_line_helper("##source=1000GenomesPhase3Pipeline").unwrap(),
        ("source", "1000GenomesPhase3Pipeline", HashMap::new())
    );

    assert_eq!(
        vcf_header_line_helper("##fileformat=VCFv4.1").unwrap(),
        ("fileformat", "VCFv4.1", HashMap::new())
    );
}

#[test]
fn test_header_line_from_str() {
    assert_eq!(
        VCFHeaderLine::from_str("##INFO=<ID=hoge,Number=A,Type=String,Description=\"foo\">")
            .unwrap(),
        VCFHeaderLine {
            line: "##INFO=<ID=hoge,Number=A,Type=String,Description=\"foo\">".to_string(),
            contents: VCFHeaderContent::INFO {
                id: "hoge".to_string(),
                number: Number::Allele,
                value_type: ValueType::String,
                description: "foo".to_string(),
                source: None,
                version: None,
            }
        }
    );

    assert_eq!(
        VCFHeaderLine::from_str("##INFO=<ID=hoge,Number=R,Type=Integer,Description=\"foo\">")
            .unwrap(),
        VCFHeaderLine {
            line: "##INFO=<ID=hoge,Number=R,Type=Integer,Description=\"foo\">".to_string(),
            contents: VCFHeaderContent::INFO {
                id: "hoge".to_string(),
                number: Number::Reference,
                value_type: ValueType::Integer,
                description: "foo".to_string(),
                source: None,
                version: None,
            }
        }
    );

    assert_eq!(
        VCFHeaderLine::from_str(
            "##INFO=<ID=bar,Number=G,Type=Float,Description=\"foo\",Source=\"ss\",Version=vv>"
        )
        .unwrap(),
        VCFHeaderLine {
            line:
                "##INFO=<ID=bar,Number=G,Type=Float,Description=\"foo\",Source=\"ss\",Version=vv>"
                    .to_string(),
            contents: VCFHeaderContent::INFO {
                id: "bar".to_string(),
                number: Number::Genotype,
                value_type: ValueType::Float,
                description: "foo".to_string(),
                source: Some("ss".to_string()),
                version: Some("vv".to_string()),
            }
        }
    );

    assert_eq!(
        VCFHeaderLine::from_str("##INFO=<ID=DB,Number=0,Type=Flag,Description=\"foo\">").unwrap(),
        VCFHeaderLine {
            line: "##INFO=<ID=DB,Number=0,Type=Flag,Description=\"foo\">".to_string(),
            contents: VCFHeaderContent::INFO {
                id: "DB".to_string(),
                number: Number::Zero,
                value_type: ValueType::Flag,
                description: "foo".to_string(),
                source: None,
                version: None,
            }
        }
    );
}

#[allow(clippy::unreadable_literal)]
#[test]
fn test_write_line() -> io::Result<()> {
    {
        let mut line = Vec::new();
        VCFRecord {
            chromosome: "19".to_string(),
            position: 11472995,
            id: svec!["hoge".to_string()],
            reference: "A".to_string(),
            alternative: svec!["C", "AC"],
            quality: None,
            filter: Vec::new(),
            info: IndexMap::new(),
            format: Vec::new(),
            call: HashMap::new(),
        }
        .write_line(&mut line, &[])?;

        assert_eq!(
            &b"19\t11472995\thoge\tA\tC,AC\t.\t.\t.\n"[..],
            &line as &[u8]
        );
    }

    {
        let mut line = Vec::new();
        VCFRecord {
            chromosome: "19".to_string(),
            position: 11472995,
            id: svec!["hoge".to_string()],
            reference: "A".to_string(),
            alternative: svec!["C", "AC"],
            quality: None,
            filter: Vec::new(),
            info: IndexMap::new(),
            format: Vec::new(),
            call: HashMap::new(),
        }
        .write_line(&mut line, &svec!["foo"])?;

        assert_eq!(
            &b"19\t11472995\thoge\tA\tC,AC\t.\t.\t.\t.\t.\n"[..],
            &line as &[u8]
        );
    }

    {
        let mut line = Vec::new();
        VCFRecord {
            chromosome: "19".to_string(),
            position: 11472995,
            id: svec!["hoge".to_string()],
            reference: "A".to_string(),
            alternative: svec!["C", "AC"],
            quality: None,
            filter: Vec::new(),
            info: sihash![
                ("HOGE", Vec::<String>::new()),
                ("FOO", svec!["TEST"]),
                ("BAR", svec!["TEST1", "TEST2"])
            ],
            format: svec!["GT", "AD"],
            call: shash![(
                "foo",
                shash![("GT", svec!["0/0"]), ("AD", svec!["10", "20"])]
            )],
        }
        .write_line(&mut line, &svec!["foo"])?;

        assert_eq!(
            &b"19\t11472995\thoge\tA\tC,AC\t.\t.\tHOGE;FOO=TEST;BAR=TEST1,TEST2\tGT:AD\t0/0:10,20\n"[..],
            &line as &[u8]
        );
    }

    {
        let mut line = Vec::new();
        VCFRecord {
            chromosome: "19".to_string(),
            position: 11472995,
            id: svec!["hoge", "foo", "bar"],
            reference: "A".to_string(),
            alternative: svec!["C", "AC"],
            quality: None,
            filter: Vec::new(),
            info: sihash![
                ("HOGE", Vec::<String>::new()),
                ("FOO", svec!["TEST"]),
                ("BAR", svec!["TEST1", "TEST2"])
            ],
            format: svec!["GT", "AD"],
            call: shash![(
                "foo",
                shash![("GT", svec!["0/0"]), ("AD", svec!["10", "20"])]
            )],
        }
        .write_line(&mut line, &svec!["foo", "bar"])?;

        assert_eq!(
            &b"19\t11472995\thoge;foo;bar\tA\tC,AC\t.\t.\tHOGE;FOO=TEST;BAR=TEST1,TEST2\tGT:AD\t0/0:10,20\t.\n"[..],
            &line as &[u8]
        );
    }

    {
        let mut line = Vec::new();
        VCFRecord {
            chromosome: "19".to_string(),
            position: 11472995,
            id: svec!["hoge", "foo"],
            reference: "A".to_string(),
            alternative: svec!["C", "AC"],
            quality: None,
            filter: Vec::new(),
            info: sihash![
                ("HOGE", Vec::<String>::new()),
                ("FOO", svec!["TEST"]),
                ("BAR", svec!["TEST1", "TEST2"])
            ],
            format: svec!["GT", "AD"],
            call: shash![
                (
                    "foo",
                    shash![("GT", svec!["0/0"]), ("AD", svec!["10", "20"])]
                ),
                ("bar", shash![("GT", svec!["0/1"])])
            ],
        }
        .write_line(&mut line, &svec!["foo", "bar"])?;

        assert_eq!(
            &b"19\t11472995\thoge;foo\tA\tC,AC\t.\t.\tHOGE;FOO=TEST;BAR=TEST1,TEST2\tGT:AD\t0/0:10,20\t0/1:.\n"[..],
            &line as &[u8]
        );
    }

    Ok(())
}


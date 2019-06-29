use super::*;

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

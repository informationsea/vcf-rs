use super::*;
use crate::VCFError;
use nom::{self, bytes::complete::is_not, bytes::complete::tag, bytes::complete::take_while};

fn create_header() -> VCFHeader {
    let vcf_data = include_bytes!("../../testfiles/simple1.vcf");
    let mut vcf_reader = io::BufReader::new(&vcf_data[..]);
    let (_, _, header) = crate::header::parse_header(&mut vcf_reader).unwrap();
    header
}

#[test]
fn test_parse_nested_separated_values() {
    let mut result = Vec::new();
    let (rest, _) =
        parser::parse_nested_separated_values::<_, _, _, _, _, (&[u8], nom::error::ErrorKind)>(
            &mut result,
            b"A:X:Z,B:Y,C,D:E:F:G 1",
            is_not(&b",: "[..]),
            tag(b":"),
            tag(b","),
            false,
        )
        .unwrap();
    assert_eq!(rest, b" 1");
    assert_eq!(
        result,
        vec![
            vec![&b"A"[..], &b"X"[..], &b"Z"[..]],
            vec![&b"B"[..], &b"Y"[..]],
            vec![&b"C"[..]],
            vec![&b"D"[..], &b"E"[..], &b"F"[..], &b"G"[..]]
        ]
    );

    let (rest, _) =
        parser::parse_nested_separated_values::<_, _, _, _, _, (&[u8], nom::error::ErrorKind)>(
            &mut result,
            b" 1",
            is_not(&b",: "[..]),
            tag(b":"),
            tag(b","),
            false,
        )
        .unwrap();
    assert_eq!(rest, b" 1");
    assert_eq!(result, Vec::<Vec<&[u8]>>::new());
}

#[test]
fn test_parse_separated_values() {
    let mut result = Vec::new();
    let (rest, _) = parser::parse_separated_values::<_, _, _, (&[u8], nom::error::ErrorKind)>(
        &mut result,
        b"A X,B Y,C=Z_x",
        is_not(&b",_"[..]),
        tag(b","),
        false,
    )
    .unwrap();
    assert_eq!(rest, b"_x");
    assert_eq!(result, vec![&b"A X"[..], &b"B Y"[..], &b"C=Z"[..]]);

    let (rest, _) = parser::parse_separated_values::<_, _, _, (&[u8], nom::error::ErrorKind)>(
        &mut result,
        b"1:2",
        is_not(&b":"[..]),
        tag(b":"),
        false,
    )
    .unwrap();
    assert_eq!(rest, b"");
    assert_eq!(result, vec![&b"1"[..], &b"2"[..]]);

    let (rest, _) = parser::parse_separated_values::<_, _, _, (&[u8], nom::error::ErrorKind)>(
        &mut result,
        b"1:: 2",
        take_while(|x| x != b':' && x != b' '),
        tag(b":"),
        false,
    )
    .unwrap();
    assert_eq!(rest, b" 2");
    assert_eq!(result, vec![&b"1"[..], &b""[..], &b""[..]]);
}

#[test]
fn test_parse_info() {
    let test_info1 = b"AC=54;AF=1;AN=54;DB;DP=749\t";
    let mut parsed_info = Vec::new();
    let (rest, _) =
        parser::parse_info::<(&[u8], nom::error::ErrorKind)>(test_info1, &mut parsed_info).unwrap();
    assert_eq!(
        parsed_info,
        vec![
            (b"AC".to_vec(), vec![b"54".to_vec()]),
            (b"AF".to_vec(), vec![b"1".to_vec()]),
            (b"AN".to_vec(), vec![b"54".to_vec()]),
            (b"DB".to_vec(), vec![]),
            (b"DP".to_vec(), vec![b"749".to_vec()]),
        ]
    );
    assert_eq!(rest, b"\t");
}

#[test]
#[allow(clippy::unreadable_literal)]
#[allow(clippy::cognitive_complexity)]
fn test_parse_record1() -> Result<(), VCFError> {
    let test_record1 = &b"13\t32889968\trs206119,rs60320776\tG\tA\t25743.5\t.\tAC=54;AF=1;AN=54;DP=749\tGT:AD:DP\t1/1:0,14:14\t1/1:0,19:19\n"[..];
    let mut record = VCFRecord::new(Rc::new(create_header()));
    assert_eq!(parse_record(&test_record1, &mut record), Ok((&b""[..], ())));

    let mut record2 = record.clone();
    record2.parse_bytes(test_record1, 1)?;
    assert_eq!(record, record2);
    let record3 = VCFRecord::from_bytes(test_record1, 1, record.header.clone())?;
    assert_eq!(record, record3);

    assert_eq!(record.chromosome, b"13");
    assert_eq!(record.position, 32889968);
    assert_eq!(record.id, vec![&b"rs206119"[..], &b"rs60320776"[..]]);
    assert_eq!(record.reference, b"G");
    assert_eq!(record.alternative, vec![b"A"]);
    assert_eq!(record.qual, Some(25743.5));
    assert_eq!(record.filter, Vec::<&[u8]>::new());
    assert_eq!(
        record.info,
        vec![
            (b"AC".to_vec(), vec![b"54".to_vec()]),
            (b"AF".to_vec(), vec![b"1".to_vec()]),
            (b"AN".to_vec(), vec![b"54".to_vec()]),
            (b"DP".to_vec(), vec![b"749".to_vec()]),
        ]
    );
    assert_eq!(record.info(b"AC"), Some(&vec![b"54".to_vec()]));
    assert_eq!(record.info(b"AF"), Some(&vec![b"1".to_vec()]));
    assert_eq!(record.info(b"AN"), Some(&vec![b"54".to_vec()]));
    assert_eq!(record.info(b"DP"), Some(&vec![b"749".to_vec()]));
    assert_eq!(record.info_mut(b"DP"), Some(&mut vec![b"749".to_vec()]));
    if let Some(x) = record.info_mut(b"DP") {
        x.push(b"XXXX".to_vec());
    }
    assert_eq!(
        record.info(b"DP"),
        Some(&vec![b"749".to_vec(), b"XXXX".to_vec()])
    );
    assert_eq!(
        record.insert_info(b"AC", vec![b"X".to_vec(), b"Y".to_vec()]),
        Some(vec![b"54".to_vec()])
    );
    assert_eq!(
        record.info(b"AC"),
        Some(&vec![b"X".to_vec(), b"Y".to_vec()])
    );
    assert_eq!(record.insert_info(b"Z", vec![b"1".to_vec()]), None);
    assert_eq!(record.info(b"Z"), Some(&vec![b"1".to_vec()]));
    assert_eq!(
        record.info,
        vec![
            (b"AC".to_vec(), vec![b"X".to_vec(), b"Y".to_vec()]),
            (b"AF".to_vec(), vec![b"1".to_vec()]),
            (b"AN".to_vec(), vec![b"54".to_vec()]),
            (b"DP".to_vec(), vec![b"749".to_vec(), b"XXXX".to_vec()]),
            (b"Z".to_vec(), vec![b"1".to_vec()])
        ]
    );

    assert_eq!(record.format, vec![&b"GT"[..], &b"AD"[..], &b"DP"[..]]);
    assert_eq!(
        record.genotype,
        vec![
            vec![
                vec![&b"1/1"[..]],
                vec![&b"0"[..], &b"14"[..]],
                vec![&b"14"[..]]
            ],
            vec![
                vec![&b"1/1"[..]],
                vec![&b"0"[..], &b"19"[..]],
                vec![&b"19"[..]]
            ]
        ]
    );
    assert_eq!(
        record.genotype(b"ERP001775_HiSeq2000_SAMEA1531955-1", b"GT"),
        Some(&vec![b"1/1".to_vec()])
    );
    assert_eq!(
        record.genotype(b"ERP001775_HiSeq2000_SAMEA1531955-1", b"AD"),
        Some(&vec![b"0".to_vec(), b"14".to_vec()])
    );
    assert_eq!(
        record.genotype(b"ERP001775_HiSeq2000_SAMEA1531955-1", b"DP"),
        Some(&vec![b"14".to_vec()])
    );

    assert_eq!(
        record.genotype_mut(b"ERP001775_HiSeq2000_SAMEA1531955-1", b"DP"),
        Some(&mut vec![b"14".to_vec()])
    );
    if let Some(x) = record.genotype_mut(b"ERP001775_HiSeq2000_SAMEA1531955-1", b"DP") {
        x.push(b"XXX".to_vec());
    }
    assert_eq!(
        record.genotype(b"ERP001775_HiSeq2000_SAMEA1531955-1", b"DP"),
        Some(&vec![b"14".to_vec(), b"XXX".to_vec()])
    );

    assert_eq!(
        record.insert_genotype(
            b"ERP001775_HiSeq2000_SAMEA1531955-1",
            b"AD",
            vec![b"A".to_vec(), b"B".to_vec()]
        ),
        Some(vec![b"0".to_vec(), b"14".to_vec()])
    );
    assert_eq!(
        record.genotype(b"ERP001775_HiSeq2000_SAMEA1531955-1", b"AD"),
        Some(&vec![b"A".to_vec(), b"B".to_vec()])
    );

    assert_eq!(
        record.insert_genotype(
            b"ERP001775_HiSeq2000_SAMEA1531955-1",
            b"X",
            vec![b"100".to_vec()]
        ),
        None
    );
    assert_eq!(
        record.insert_genotype(
            b"ERP001775_HiSeq2000_SAMEA1531955-2",
            b"X",
            vec![b"200".to_vec()]
        ),
        None
    );
    assert_eq!(
        record.insert_genotype(
            b"ERP001775_HiSeq2000_SAMEA1531955-2",
            b"Y",
            vec![b"300".to_vec()]
        ),
        None
    );
    assert_eq!(
        record.insert_genotype(
            b"ERP001775_HiSeq2000_SAMEA1531955-2",
            b"Z",
            vec![b"400".to_vec()]
        ),
        None
    );
    assert_eq!(
        record.insert_genotype(
            b"ERP001775_HiSeq2000_SAMEA1531955-1",
            b"Z",
            vec![b"500".to_vec()]
        ),
        None
    );
    assert_eq!(
        record.format,
        vec![
            &b"GT"[..],
            &b"AD"[..],
            &b"DP"[..],
            &b"X"[..],
            &b"Y"[..],
            &b"Z"[..]
        ]
    );
    assert_eq!(
        record.genotype,
        vec![
            vec![
                vec![&b"1/1"[..]],
                vec![&b"A"[..], &b"B"[..]],
                vec![&b"14"[..], &b"XXX"[..]],
                vec![&b"100"[..]],
                vec![&b"."[..]],
                vec![&b"500"[..]]
            ],
            vec![
                vec![&b"1/1"[..]],
                vec![&b"0"[..], &b"19"[..]],
                vec![&b"19"[..]],
                vec![&b"200"[..]],
                vec![&b"300"[..]],
                vec![&b"400"[..]]
            ]
        ]
    );

    Ok(())
}

#[test]
#[allow(clippy::unreadable_literal)]
fn test_parse_record2() -> Result<(), VCFError> {
    let test_record1 = &b"13\t32889968\trs206119,rs60320776\tG\tA"[..];
    let mut record = VCFRecord::new(Rc::new(create_header()));
    assert_eq!(parse_record(&test_record1, &mut record), Ok((&b""[..], ())));
    assert_eq!(record.chromosome, b"13");
    assert_eq!(record.position, 32889968);
    assert_eq!(record.id, vec![&b"rs206119"[..], &b"rs60320776"[..]]);
    assert_eq!(record.reference, b"G");
    assert_eq!(record.alternative, vec![b"A"]);
    Ok(())
}

#[test]
#[allow(clippy::unreadable_literal)]
fn test_parse_record3() -> Result<(), VCFError> {
    let test_record1 = &b"13\t32889968\trs206119,rs60320776\tG\tA\t25743.5\r\n"[..];
    let mut record = VCFRecord::new(Rc::new(create_header()));

    assert_eq!(parse_record(&test_record1, &mut record), Ok((&b""[..], ())));
    assert_eq!(record.chromosome, b"13");
    assert_eq!(record.position, 32889968);
    assert_eq!(record.id, vec![&b"rs206119"[..], &b"rs60320776"[..]]);
    assert_eq!(record.reference, b"G");
    assert_eq!(record.alternative, vec![b"A"]);
    assert_eq!(record.qual, Some(25743.5));
    Ok(())
}

#[test]
#[allow(clippy::unreadable_literal)]
fn test_parse_record3_1() -> Result<(), VCFError> {
    let test_record1 = &b"13\t32889968\trs206119,rs60320776\tG\tA\t8.39728e+06\r\n"[..];
    let mut record = VCFRecord::new(Rc::new(create_header()));

    assert_eq!(parse_record(&test_record1, &mut record), Ok((&b""[..], ())));
    assert_eq!(record.chromosome, b"13");
    assert_eq!(record.position, 32889968);
    assert_eq!(record.id, vec![&b"rs206119"[..], &b"rs60320776"[..]]);
    assert_eq!(record.reference, b"G");
    assert_eq!(record.alternative, vec![b"A"]);
    assert_eq!(record.qual, Some(8.39728e+06));
    Ok(())
}

#[test]
#[allow(clippy::unreadable_literal)]
fn test_parse_record4() -> Result<(), VCFError> {
    let test_record1 = &b"13\t32889968\trs206119,rs60320776\tG\tA\t25743.5\tPASS"[..];
    let mut record = VCFRecord::new(Rc::new(create_header()));

    assert_eq!(parse_record(&test_record1, &mut record), Ok((&b""[..], ())));
    assert_eq!(record.chromosome, b"13");
    assert_eq!(record.position, 32889968);
    assert_eq!(record.id, vec![&b"rs206119"[..], &b"rs60320776"[..]]);
    assert_eq!(record.reference, b"G");
    assert_eq!(record.alternative, vec![b"A"]);
    assert_eq!(record.qual, Some(25743.5));
    assert_eq!(record.filter, vec![&b"PASS"[..]]);
    Ok(())
}

#[test]
#[allow(clippy::unreadable_literal)]
fn test_parse_record5() -> Result<(), VCFError> {
    let test_record1 =
        &b"13\t32889968\trs206119,rs60320776\tG\tA\t25743.5\t.\tAC=54;AF=1;AN=54;DP=749\n"[..];
    let mut record = VCFRecord::new(Rc::new(create_header()));

    assert_eq!(parse_record(&test_record1, &mut record), Ok((&b""[..], ())));
    assert_eq!(record.chromosome, b"13");
    assert_eq!(record.position, 32889968);
    assert_eq!(record.id, vec![&b"rs206119"[..], &b"rs60320776"[..]]);
    assert_eq!(record.reference, b"G");
    assert_eq!(record.alternative, vec![b"A"]);
    assert_eq!(record.qual, Some(25743.5));
    assert_eq!(record.filter, Vec::<&[u8]>::new());
    assert_eq!(
        record.info,
        vec![
            (b"AC".to_vec(), vec![b"54".to_vec()]),
            (b"AF".to_vec(), vec![b"1".to_vec()]),
            (b"AN".to_vec(), vec![b"54".to_vec()]),
            (b"DP".to_vec(), vec![b"749".to_vec()]),
        ]
    );
    Ok(())
}

#[test]
#[allow(clippy::unreadable_literal)]
fn test_parse_record6() -> Result<(), VCFError> {
    let test_record1 = &b"13\t32889968\t.\tG\t.\t.\t.\tAC=54;AF=1;AN=54;DP=749\tGT:AD:DP"[..];
    let mut record = VCFRecord::new(Rc::new(create_header()));

    assert_eq!(parse_record(&test_record1, &mut record), Ok((&b""[..], ())));
    assert_eq!(record.chromosome, b"13");
    assert_eq!(record.position, 32889968);
    assert_eq!(record.id, Vec::<&[u8]>::new());
    assert_eq!(record.reference, b"G");
    assert_eq!(record.alternative, Vec::<&[u8]>::new());
    assert_eq!(record.qual, None);
    assert_eq!(record.filter, Vec::<&[u8]>::new());
    assert_eq!(
        record.info,
        vec![
            (b"AC".to_vec(), vec![b"54".to_vec()]),
            (b"AF".to_vec(), vec![b"1".to_vec()]),
            (b"AN".to_vec(), vec![b"54".to_vec()]),
            (b"DP".to_vec(), vec![b"749".to_vec()]),
        ]
    );
    assert_eq!(record.format, vec![&b"GT"[..], &b"AD"[..], &b"DP"[..]]);
    Ok(())
}

#[test]
#[allow(clippy::unreadable_literal)]
#[allow(clippy::cognitive_complexity)]
fn test_parse_record7() -> Result<(), VCFError> {
    let test_record1 = &b"13\t32889968\trs206119,rs60320776\tG\tA\t25743.5\t.\tAC=54;AF=1;AN=54;DP=749\tGT:AD:DP\t1/1:0,14:14\t1/1:0,19:19\n"[..];
    let mut record = VCFRecord::new(Rc::new(create_header()));

    assert_eq!(parse_record(&test_record1, &mut record), Ok((&b""[..], ())));
    assert_eq!(record.chromosome, b"13");
    assert_eq!(record.position, 32889968);
    assert_eq!(record.id, vec![&b"rs206119"[..], &b"rs60320776"[..]]);
    assert_eq!(record.reference, b"G");
    assert_eq!(record.alternative, vec![b"A"]);
    assert_eq!(record.qual, Some(25743.5));
    assert_eq!(record.filter, Vec::<&[u8]>::new());
    assert_eq!(
        record.info,
        vec![
            (b"AC".to_vec(), vec![b"54".to_vec()]),
            (b"AF".to_vec(), vec![b"1".to_vec()]),
            (b"AN".to_vec(), vec![b"54".to_vec()]),
            (b"DP".to_vec(), vec![b"749".to_vec()]),
        ]
    );
    assert_eq!(record.format, vec![&b"GT"[..], &b"AD"[..], &b"DP"[..]]);
    assert_eq!(
        record.genotype,
        vec![
            vec![
                vec![&b"1/1"[..]],
                vec![&b"0"[..], &b"14"[..]],
                vec![&b"14"[..]]
            ],
            vec![
                vec![&b"1/1"[..]],
                vec![&b"0"[..], &b"19"[..]],
                vec![&b"19"[..]]
            ]
        ]
    );
    assert_eq!(record.info(b"AC"), Some(&vec![b"54".to_vec()]));
    assert_eq!(record.info(b"AF"), Some(&vec![b"1".to_vec()]));
    assert_eq!(record.info(b"XX"), None);

    let test_record1 = &b"13\t32889968\trs206119,rs60320776\tG\tA"[..];
    assert_eq!(parse_record(&test_record1, &mut record), Ok((&b""[..], ())));
    assert_eq!(record.chromosome, b"13");
    assert_eq!(record.position, 32889968);
    assert_eq!(record.id, vec![&b"rs206119"[..], &b"rs60320776"[..]]);
    assert_eq!(record.reference, b"G");
    assert_eq!(record.alternative, vec![b"A"]);
    assert_eq!(record.qual, None);
    assert_eq!(record.filter, Vec::<U8Vec>::new());
    assert_eq!(record.info, vec![]);
    assert_eq!(record.format, Vec::<U8Vec>::new());
    assert_eq!(record.genotype, Vec::<Vec<Vec<U8Vec>>>::new());
    assert_eq!(record.info(b"XX"), None);
    assert_eq!(record.info(b"AC"), None);
    assert_eq!(record.info(b"AF"), None);

    let test_record1 = &b"13\t32889968\trs206119,rs60320776\tG\tA\t25743.5\r\n"[..];
    assert_eq!(parse_record(&test_record1, &mut record), Ok((&b""[..], ())));
    assert_eq!(record.chromosome, b"13");
    assert_eq!(record.position, 32889968);
    assert_eq!(record.id, vec![&b"rs206119"[..], &b"rs60320776"[..]]);
    assert_eq!(record.reference, b"G");
    assert_eq!(record.alternative, vec![b"A"]);
    assert_eq!(record.qual, Some(25743.5));

    let test_record1 = &b"13\t32889968\trs206119,rs60320776\tG\tA\t25743.5\tPASS"[..];

    assert_eq!(parse_record(&test_record1, &mut record), Ok((&b""[..], ())));
    assert_eq!(record.chromosome, b"13");
    assert_eq!(record.position, 32889968);
    assert_eq!(record.id, vec![&b"rs206119"[..], &b"rs60320776"[..]]);
    assert_eq!(record.reference, b"G");
    assert_eq!(record.alternative, vec![b"A"]);
    assert_eq!(record.qual, Some(25743.5));
    assert_eq!(record.filter, vec![&b"PASS"[..]]);

    let test_record1 =
        &b"13\t32889968\trs206119,rs60320776\tG\tA\t25743.5\t.\tAC=54;AF=1;AN=54;DP=749\n"[..];

    assert_eq!(parse_record(&test_record1, &mut record), Ok((&b""[..], ())));
    assert_eq!(record.chromosome, b"13");
    assert_eq!(record.position, 32889968);
    assert_eq!(record.id, vec![&b"rs206119"[..], &b"rs60320776"[..]]);
    assert_eq!(record.reference, b"G");
    assert_eq!(record.alternative, vec![b"A"]);
    assert_eq!(record.qual, Some(25743.5));
    assert_eq!(record.filter, Vec::<&[u8]>::new());
    assert_eq!(
        record.info,
        vec![
            (b"AC".to_vec(), vec![b"54".to_vec()]),
            (b"AF".to_vec(), vec![b"1".to_vec()]),
            (b"AN".to_vec(), vec![b"54".to_vec()]),
            (b"DP".to_vec(), vec![b"749".to_vec()]),
        ]
    );
    assert_eq!(record.info(b"AC"), Some(&vec![b"54".to_vec()]));
    assert_eq!(record.info(b"AF"), Some(&vec![b"1".to_vec()]));
    assert_eq!(record.info(b"XX"), None);

    let test_record1 = &b"13\t32889968\t.\tG\t.\t.\t.\tAC=54;AF=1;AN=54;DP=749\tGT:AD:DP"[..];

    assert_eq!(parse_record(&test_record1, &mut record), Ok((&b""[..], ())));
    assert_eq!(record.chromosome, b"13");
    assert_eq!(record.position, 32889968);
    assert_eq!(record.id, Vec::<&[u8]>::new());
    assert_eq!(record.reference, b"G");
    assert_eq!(record.alternative, Vec::<&[u8]>::new());
    assert_eq!(record.qual, None);
    assert_eq!(record.filter, Vec::<&[u8]>::new());
    assert_eq!(
        record.info,
        vec![
            (b"AC".to_vec(), vec![b"54".to_vec()]),
            (b"AF".to_vec(), vec![b"1".to_vec()]),
            (b"AN".to_vec(), vec![b"54".to_vec()]),
            (b"DP".to_vec(), vec![b"749".to_vec()]),
        ]
    );
    assert_eq!(record.format, vec![&b"GT"[..], &b"AD"[..], &b"DP"[..]]);

    let test_record1 = &b"13\t32889968\t.\tG\t.\t.\t.\t.\n"[..];

    assert_eq!(parse_record(&test_record1, &mut record), Ok((&b""[..], ())));
    assert_eq!(record.chromosome, b"13");
    assert_eq!(record.position, 32889968);
    assert_eq!(record.id, Vec::<&[u8]>::new());
    assert_eq!(record.reference, b"G");
    assert_eq!(record.alternative, Vec::<&[u8]>::new());
    assert_eq!(record.qual, None);
    assert_eq!(record.filter, Vec::<&[u8]>::new());
    assert_eq!(record.info, vec![]);
    Ok(())
}

#[test]
fn test_write() {
    let mut vcf_record = VCFRecord::new(Rc::new(create_header()));
    let test_record1 = &b"13\t32889968\t.\tG\tA\t.\t.\tAC=1;AF=0.5\n"[..];
    parse_record(test_record1, &mut vcf_record).unwrap();
    let mut write_data = Vec::new();
    vcf_record.write_record(&mut write_data).unwrap();
    assert_eq!(write_data, test_record1);

    let test_record2 =
        &b"13\t32889968\t123\tG\tA\t102.0\tLOW,HIGH\tAC=1;AF=0.5\tGT:DP\t1/1:10\t0/1:20\n"[..];
    parse_record(test_record2, &mut vcf_record).unwrap();
    let mut write_data = Vec::new();
    vcf_record.write_record(&mut write_data).unwrap();
    assert_eq!(write_data, test_record2);
}

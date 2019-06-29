use super::*;

#[allow(clippy::unreadable_literal)]
#[test]
fn test_parse_vcf_record() -> Result<(), VCFParseError> {
    assert_eq!(VCFRecord::parse_line("19\t11472995\t.\tA\tC,AC\t.\tPASS\tCONTQ=93;DP=84;ECNT=1;GERMQ=93;MBQ=20,20,25;MFRL=208,245,147;MMQ=60,60,60;MPOS=11,34;NALOD=-2.009,1.16;NLOD=5.61,9.79;POPAF=6,6;SEQQ=3;STRANDQ=5;TLOD=3.84,6.63\tGT:AD:AF:DP:F1R2:F2R1:SB\t0/1/2:20,7,6:0.199,0.211:33:6,4,2:14,3,4:2,18,0,13\t0/0:35,3,0:0.084,0.023:38:13,2,0:22,1,0:0,35,0,3",
        &["SRP171128-4213-4Met".to_string(), "SRP171128-4213N".to_string()])?,
        VCFRecord {
            chromosome: "19".to_string(),
            position: 11472995,
            id: Vec::new(),
            reference: "A".to_string(),
            alternative: svec!["C", "AC"],
            quality: None,
            filter: svec!["PASS"],
            info: sihash![
                ("CONTQ", svec!["93"]),
                ("DP", svec!["84"]),
                ("ECNT", svec!["1"]),
                ("GERMQ", svec!["93"]),
                ("MBQ", svec!["20", "20", "25"]),
                ("MFRL", svec!["208", "245", "147"]),
                ("MMQ", svec!["60", "60", "60"]),
                ("MPOS", svec!["11", "34"]),
                ("NALOD", svec!["-2.009", "1.16"]),
                ("NLOD", svec!["5.61", "9.79"]),
                ("POPAF", svec!["6", "6"]),
                ("SEQQ", svec!["3"]),
                ("STRANDQ", svec!["5"]),
                ("TLOD", svec!["3.84", "6.63"])
            ],
            format: svec!["GT", "AD", "AF", "DP", "F1R2", "F2R1", "SB"],
            call: shash![
                ("SRP171128-4213-4Met", shash![
                    ("GT", svec!["0/1/2"]),
                    ("AD", svec!["20", "7", "6"]),
                    ("AF", svec!["0.199", "0.211"]),
                    ("DP", svec!["33"]),
                    ("F1R2", svec!["6", "4", "2"]),
                    ("F2R1", svec!["14", "3", "4"]),
                    ("SB", svec!["2", "18", "0", "13"])
                ]),
                ("SRP171128-4213N", shash![
                    ("GT", svec!["0/0"]),
                    ("AD", svec!["35", "3", "0"]),
                    ("AF", svec!["0.084", "0.023"]),
                    ("DP", svec!["38"]),
                    ("F1R2", svec!["13", "2", "0"]),
                    ("F2R1", svec!["22", "1", "0"]),
                    ("SB", svec!["0", "35", "0", "3"])
                ])
            ]
        }
    );

    assert_eq!(VCFRecord::parse_line("19\t11472995\t.\tA\tC,AC\t.\tPASS\tCONTQ=93;DP=84;ECNT=1;GERMQ=93;MBQ=20,20,25;MFRL=208,245,147;MMQ=60,60,60;MPOS=11,34;NALOD=-2.009,1.16;NLOD=5.61,9.79;POPAF=6,6;SEQQ=3;STRANDQ=5;TLOD=3.84,6.63\tGT:AD:AF:DP:F1R2:F2R1:SB\t0/1/2:20,7,6:0.199,0.211:33:6,4,2:14,3,4:2,18,0,13",
        &["SRP171128-4213-4Met".to_string(), "SRP171128-4213N".to_string()])?,
        VCFRecord {
            chromosome: "19".to_string(),
            position: 11472995,
            id: Vec::new(),
            reference: "A".to_string(),
            alternative: svec!["C", "AC"],
            quality: None,
            filter: svec!["PASS"],
            info: sihash![
                ("CONTQ", svec!["93"]),
                ("DP", svec!["84"]),
                ("ECNT", svec!["1"]),
                ("GERMQ", svec!["93"]),
                ("MBQ", svec!["20", "20", "25"]),
                ("MFRL", svec!["208", "245", "147"]),
                ("MMQ", svec!["60", "60", "60"]),
                ("MPOS", svec!["11", "34"]),
                ("NALOD", svec!["-2.009", "1.16"]),
                ("NLOD", svec!["5.61", "9.79"]),
                ("POPAF", svec!["6", "6"]),
                ("SEQQ", svec!["3"]),
                ("STRANDQ", svec!["5"]),
                ("TLOD", svec!["3.84", "6.63"])
            ],
            format: svec!["GT", "AD", "AF", "DP", "F1R2", "F2R1", "SB"],
            call: shash![
                ("SRP171128-4213-4Met", shash![
                    ("GT", svec!["0/1/2"]),
                    ("AD", svec!["20", "7", "6"]),
                    ("AF", svec!["0.199", "0.211"]),
                    ("DP", svec!["33"]),
                    ("F1R2", svec!["6", "4", "2"]),
                    ("F2R1", svec!["14", "3", "4"]),
                    ("SB", svec!["2", "18", "0", "13"])
                ])
            ]
        }
    );

    assert_eq!(VCFRecord::parse_line("19\t11472995\t.\tA\tC,AC\t.\tPASS\tCONTQ=93;DP=84;ECNT=1;GERMQ=93;MBQ=20,20,25;MFRL=208,245,147;MMQ=60,60,60;MPOS=11,34;NALOD=-2.009,1.16;NLOD=5.61,9.79;POPAF=6,6;SEQQ=3;STRANDQ=5;TLOD=3.84,6.63\tGT:AD:AF:DP:F1R2:F2R1:SB",
        &["SRP171128-4213-4Met".to_string(), "SRP171128-4213N".to_string()])?,
        VCFRecord {
            chromosome: "19".to_string(),
            position: 11472995,
            id: Vec::new(),
            reference: "A".to_string(),
            alternative: svec!["C", "AC"],
            quality: None,
            filter: svec!["PASS"],
            info: sihash![
                ("CONTQ", svec!["93"]),
                ("DP", svec!["84"]),
                ("ECNT", svec!["1"]),
                ("GERMQ", svec!["93"]),
                ("MBQ", svec!["20", "20", "25"]),
                ("MFRL", svec!["208", "245", "147"]),
                ("MMQ", svec!["60", "60", "60"]),
                ("MPOS", svec!["11", "34"]),
                ("NALOD", svec!["-2.009", "1.16"]),
                ("NLOD", svec!["5.61", "9.79"]),
                ("POPAF", svec!["6", "6"]),
                ("SEQQ", svec!["3"]),
                ("STRANDQ", svec!["5"]),
                ("TLOD", svec!["3.84", "6.63"])
            ],
            format: svec!["GT", "AD", "AF", "DP", "F1R2", "F2R1", "SB"],
            call: HashMap::new()
        }
    );

    assert_eq!(VCFRecord::parse_line("19\t11472995\t.\tA\tC,AC\t.\tPASS\tCONTQ=93;DP=84;ECNT=1;GERMQ=93;MBQ=20,20,25;MFRL=208,245,147;MMQ=60,60,60;MPOS=11,34;NALOD=-2.009,1.16;NLOD=5.61,9.79;POPAF=6,6;SEQQ=3;STRANDQ=5;TLOD=3.84,6.63",
        &["SRP171128-4213-4Met".to_string(), "SRP171128-4213N".to_string()])?,
        VCFRecord {
            chromosome: "19".to_string(),
            position: 11472995,
            id: Vec::new(),
            reference: "A".to_string(),
            alternative: svec!["C", "AC"],
            quality: None,
            filter: svec!["PASS"],
            info: sihash![
                ("CONTQ", svec!["93"]),
                ("DP", svec!["84"]),
                ("ECNT", svec!["1"]),
                ("GERMQ", svec!["93"]),
                ("MBQ", svec!["20", "20", "25"]),
                ("MFRL", svec!["208", "245", "147"]),
                ("MMQ", svec!["60", "60", "60"]),
                ("MPOS", svec!["11", "34"]),
                ("NALOD", svec!["-2.009", "1.16"]),
                ("NLOD", svec!["5.61", "9.79"]),
                ("POPAF", svec!["6", "6"]),
                ("SEQQ", svec!["3"]),
                ("STRANDQ", svec!["5"]),
                ("TLOD", svec!["3.84", "6.63"])
            ],
            format: Vec::new(),
            call: HashMap::new()
        }
    );

    assert_eq!(
        VCFRecord::parse_line(
            "19\t11472995\t.\tA\tC,AC\t.\tPASS",
            &[
                "SRP171128-4213-4Met".to_string(),
                "SRP171128-4213N".to_string()
            ]
        )?,
        VCFRecord {
            chromosome: "19".to_string(),
            position: 11472995,
            id: Vec::new(),
            reference: "A".to_string(),
            alternative: svec!["C", "AC"],
            quality: None,
            filter: svec!["PASS"],
            info: IndexMap::new(),
            format: Vec::new(),
            call: HashMap::new()
        }
    );

    assert_eq!(
        VCFRecord::parse_line(
            "19\t11472995\t.\tA\tC,AC\t.\t.",
            &[
                "SRP171128-4213-4Met".to_string(),
                "SRP171128-4213N".to_string()
            ]
        )?,
        VCFRecord {
            chromosome: "19".to_string(),
            position: 11472995,
            id: Vec::new(),
            reference: "A".to_string(),
            alternative: svec!["C", "AC"],
            quality: None,
            filter: Vec::new(),
            info: IndexMap::new(),
            format: Vec::new(),
            call: HashMap::new()
        }
    );

    assert_eq!(
        VCFRecord::parse_line(
            "19\t11472995\t.\tA\tC,AC\t.",
            &[
                "SRP171128-4213-4Met".to_string(),
                "SRP171128-4213N".to_string()
            ]
        )?,
        VCFRecord {
            chromosome: "19".to_string(),
            position: 11472995,
            id: Vec::new(),
            reference: "A".to_string(),
            alternative: svec!["C", "AC"],
            quality: None,
            filter: Vec::new(),
            info: IndexMap::new(),
            format: Vec::new(),
            call: HashMap::new()
        }
    );

    assert_eq!(
        VCFRecord::parse_line(
            "19\t11472995\t.\tA\tC,AC\t20",
            &[
                "SRP171128-4213-4Met".to_string(),
                "SRP171128-4213N".to_string()
            ]
        )?,
        VCFRecord {
            chromosome: "19".to_string(),
            position: 11472995,
            id: Vec::new(),
            reference: "A".to_string(),
            alternative: svec!["C", "AC"],
            quality: Some("20".to_string()),
            filter: Vec::new(),
            info: IndexMap::new(),
            format: Vec::new(),
            call: HashMap::new()
        }
    );

    assert_eq!(
        VCFRecord::parse_line(
            "19\t11472995\t.\tA\tC,AC",
            &[
                "SRP171128-4213-4Met".to_string(),
                "SRP171128-4213N".to_string()
            ]
        )?,
        VCFRecord {
            chromosome: "19".to_string(),
            position: 11472995,
            id: Vec::new(),
            reference: "A".to_string(),
            alternative: svec!["C", "AC"],
            quality: None,
            filter: Vec::new(),
            info: IndexMap::new(),
            format: Vec::new(),
            call: HashMap::new()
        }
    );

    assert_eq!(
        VCFRecord::parse_line(
            "19\t11472995\thoge\tA\tC,AC",
            &[
                "SRP171128-4213-4Met".to_string(),
                "SRP171128-4213N".to_string()
            ]
        )?,
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
            call: HashMap::new()
        }
    );

    assert!(VCFRecord::parse_line("19\t11472995\thoge\tA", &[]).is_err());

    Ok(())
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


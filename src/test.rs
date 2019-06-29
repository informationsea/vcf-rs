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


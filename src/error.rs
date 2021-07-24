use thiserror::Error;

#[derive(Error, Debug)]
pub enum VCFError {
    #[error("Failed to parse header at line: {}", _0)]
    HeaderParseError(u64),
    #[error("Failed to parse record at line: {}", _0)]
    RecordParseError(u64),
    #[error("I/O Error")]
    IoError(#[from] std::io::Error),
    #[error("Utf8 Error")]
    Utf8Error(#[from] std::str::Utf8Error),
}

use failure::{Backtrace, Context, Fail};
use std::fmt::{self, Display};

#[derive(Debug)]
pub struct VCFError {
    inner: failure::Context<VCFErrorKind>,
}

#[derive(Copy, Clone, Eq, PartialEq, Debug, Fail)]
pub enum VCFErrorKind {
    #[fail(display = "Failed to parse header at line: {}", _0)]
    HeaderParseError(u64),
    #[fail(display = "Failed to parse record at line: {}", _0)]
    RecordParseError(u64),
    #[fail(display = "I/O Error")]
    IoError,
    #[fail(display = "Utf8 Error")]
    Utf8Error,
}

impl Fail for VCFError {
    fn cause(&self) -> Option<&dyn Fail> {
        self.inner.cause()
    }

    fn backtrace(&self) -> Option<&Backtrace> {
        self.inner.backtrace()
    }
}

impl Display for VCFError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        Display::fmt(&self.inner, f)
    }
}

impl VCFError {
    pub fn kind(&self) -> VCFErrorKind {
        *self.inner.get_context()
    }
}

impl From<VCFErrorKind> for VCFError {
    fn from(kind: VCFErrorKind) -> VCFError {
        VCFError {
            inner: Context::new(kind),
        }
    }
}

impl From<Context<VCFErrorKind>> for VCFError {
    fn from(inner: Context<VCFErrorKind>) -> VCFError {
        VCFError { inner }
    }
}

impl From<std::io::Error> for VCFError {
    fn from(e: std::io::Error) -> VCFError {
        e.context(VCFErrorKind::IoError).into()
    }
}

impl From<std::str::Utf8Error> for VCFError {
    fn from(e: std::str::Utf8Error) -> VCFError {
        e.context(VCFErrorKind::Utf8Error).into()
    }
}

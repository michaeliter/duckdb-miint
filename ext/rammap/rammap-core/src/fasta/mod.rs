pub mod record;
pub mod reader;
pub mod stream;

pub use record::Record;
pub use reader::{Reader, open, FastxError, parse_fasta_bytes};
pub use stream::{FastaStreamer, FastqStreamer};

pub use reader::read_fasta;

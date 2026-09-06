//! `help_functions.readfq`, ported.
//!
//! The reference uses Heng Li's `readfq`, and it has two behaviours worth
//! naming because they are load-bearing and easy to "fix" by accident:
//!
//! 1. **Spaces in the header become underscores** (`last[1:].replace(" ", "_")`).
//!    Accessions are later split on `_`, and the score is appended with `_`, so
//!    this is not cosmetic.
//! 2. **Lines are chomped, not `rstrip`ped.** The reference drops a trailing
//!    `\n` and nothing else, so a CRLF file carries its `\r` into the name,
//!    sequence and quality -- in the reference too. Do not "fix" that here; it
//!    would be a silent divergence.
//!
//!    It used to use `l[:-1]`, removing the final character whatever it was,
//!    which made a file with no trailing newline yield a record with **no
//!    quality at all** and crash the caller (Finding 11). That is fixed in the
//!    Python now, and this matches the fixed behaviour. A quality string that
//!    is genuinely shorter than the sequence still yields `None`, and still
//!    crashes the reference -- see `run` in main.rs.

/// One record. `qual` is `None` for fasta input.
#[derive(Debug, Clone, PartialEq)]
pub struct Record {
    pub name: String,
    pub seq: String,
    pub qual: Option<String>,
}

/// Drop a trailing newline if there is one -- the reference's `_chomp`.
///
/// Not `trim_end`: only `\n` goes. A CRLF file keeps its `\r`, in the
/// reference and here alike.
fn chop(line: &str) -> &str {
    line.strip_suffix('\n').unwrap_or(line)
}

/// Parse fastq/fasta exactly as the reference's generator does.
pub fn read(text: &str) -> Vec<Record> {
    let lines: Vec<&str> = text.split_inclusive('\n').collect();
    let mut out = Vec::new();
    let mut i = 0usize;
    let mut last: Option<String> = None;

    loop {
        if last.is_none() {
            // Look for the next header.
            while i < lines.len() {
                let l = lines[i];
                i += 1;
                if l.starts_with('>') || l.starts_with('@') {
                    last = Some(chop(l).to_string());
                    break;
                }
            }
        }
        let header = match last.take() {
            Some(h) => h,
            None => break,
        };

        let name = header[1..].replace(' ', "_");
        let mut seqs: Vec<&str> = Vec::new();
        let mut next_header: Option<String> = None;
        while i < lines.len() {
            let l = lines[i];
            i += 1;
            if l.starts_with('@') || l.starts_with('+') || l.starts_with('>') {
                next_header = Some(chop(l).to_string());
                break;
            }
            seqs.push(chop(l));
        }

        let is_fastq = matches!(&next_header, Some(h) if h.starts_with('+'));
        if !is_fastq {
            // fasta record
            out.push(Record {
                name,
                seq: seqs.concat(),
                qual: None,
            });
            match next_header {
                Some(h) => last = Some(h),
                None => break,
            }
            continue;
        }

        let seq = seqs.concat();
        let mut quals: Vec<&str> = Vec::new();
        let mut leng = 0usize;
        let mut completed = false;
        while i < lines.len() {
            let l = lines[i];
            i += 1;
            let q = chop(l);
            quals.push(q);
            leng += q.chars().count();
            if leng >= seq.chars().count() {
                last = None;
                out.push(Record {
                    name: name.clone(),
                    seq: seq.clone(),
                    qual: Some(quals.concat()),
                });
                completed = true;
                break;
            }
        }
        if !completed {
            // EOF before enough quality: the reference yields a fasta record
            // and stops entirely.
            out.push(Record {
                name,
                seq,
                qual: None,
            });
            break;
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reads_a_basic_fastq() {
        let r = read("@r1\nACGT\n+\nIIII\n@r2\nTTTT\n+\nJJJJ\n");
        assert_eq!(r.len(), 2);
        assert_eq!(r[0].name, "r1");
        assert_eq!(r[0].seq, "ACGT");
        assert_eq!(r[0].qual.as_deref(), Some("IIII"));
        assert_eq!(r[1].name, "r2");
    }

    /// Spaces become underscores, which matters because accessions are split on
    /// `_` to recover the score.
    #[test]
    fn spaces_in_the_header_become_underscores() {
        let r = read("@read 1 strand=+\nACGT\n+\nIIII\n");
        assert_eq!(r[0].name, "read_1_strand=+");
    }

    /// Finding 11, now fixed in the Python: a file whose last line has no
    /// newline parses normally. Before the fix the whole quality string was
    /// dropped and the caller died with
    /// `TypeError: 'NoneType' object is not iterable`.
    #[test]
    fn a_missing_final_newline_parses_normally() {
        let with = read("@r1\nACGT\n+\nIIII\n");
        let without = read("@r1\nACGT\n+\nIIII");
        assert_eq!(with[0].qual.as_deref(), Some("IIII"));
        assert_eq!(without[0].qual.as_deref(), Some("IIII"));
        assert_eq!(with, without);
    }

    #[test]
    fn a_missing_final_newline_does_not_lose_the_last_read() {
        let r = read("@r1\nACGT\n+\nIIII\n@r2\nTTTT\n+\nJJJJ");
        assert_eq!(r.len(), 2);
        assert_eq!(r[1].name, "r2");
        assert_eq!(r[1].qual.as_deref(), Some("JJJJ"));
    }

    /// CRLF is deliberately NOT handled: the reference carries the `\r` into
    /// the name, sequence and quality, so this does too.
    #[test]
    fn crlf_carries_the_carriage_return_through() {
        let r = read("@r1\r\nACGT\r\n+\r\nIIII\r\n");
        assert_eq!(r[0].name, "r1\r");
        assert_eq!(r[0].seq, "ACGT\r");
        assert_eq!(r[0].qual.as_deref(), Some("IIII\r"));
    }

    /// Quality genuinely shorter than the sequence still yields `None`, with or
    /// without a trailing newline. That path is untouched by the fix and still
    /// crashes the reference, so the port reports it.
    #[test]
    fn quality_shorter_than_sequence_yields_none() {
        let r = read("@r1\nACGTA\n+\nIII\n");
        assert_eq!(r[0].seq, "ACGTA");
        assert_eq!(r[0].qual, None);
    }

    #[test]
    fn multiline_sequence_and_quality_are_joined() {
        let r = read("@r1\nAC\nGT\n+\nII\nII\n");
        assert_eq!(r[0].seq, "ACGT");
        assert_eq!(r[0].qual.as_deref(), Some("IIII"));
    }

    #[test]
    fn fasta_records_have_no_quality() {
        let r = read(">r1\nACGT\n>r2\nTTTT\n");
        assert_eq!(r.len(), 2);
        assert_eq!(r[0].qual, None);
        assert_eq!(r[1].seq, "TTTT");
    }

    #[test]
    fn empty_input_yields_nothing() {
        assert!(read("").is_empty());
        assert!(read("\n\n").is_empty());
    }
}

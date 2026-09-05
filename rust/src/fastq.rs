//! `help_functions.readfq`, ported.
//!
//! The reference uses Heng Li's `readfq`, and it has two behaviours worth
//! naming because they are load-bearing and easy to "fix" by accident:
//!
//! 1. **Spaces in the header become underscores** (`last[1:].replace(" ", "_")`).
//!    Accessions are later split on `_`, and the score is appended with `_`, so
//!    this is not cosmetic.
//! 2. **Every line is truncated with `l[:-1]`, not `rstrip`.** That removes the
//!    final character unconditionally, assuming it is a newline. A file whose
//!    last line has no trailing newline therefore loses its last character --
//!    for a fastq, one quality value. Reproduced deliberately; see the test.

/// One record. `qual` is `None` for fasta input.
#[derive(Debug, Clone, PartialEq)]
pub struct Record {
    pub name: String,
    pub seq: String,
    pub qual: Option<String>,
}

/// Strip exactly one trailing character, as Python's `l[:-1]` does.
///
/// Note this is not `trim_end`: it drops the last character whatever it is, and
/// on a `\r\n` file it leaves the `\r` in place -- which the reference also
/// does, so a CRLF fastq carries `\r` into the sequence in both.
fn chop(line: &str) -> &str {
    let mut chars = line.chars();
    chars.next_back();
    chars.as_str()
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

    /// Without a trailing newline the final record comes back with **no
    /// quality at all**, not a truncated one: `l[:-1]` shortens the line, the
    /// running length never reaches `len(seq)`, the loop hits EOF, and the
    /// reference falls through to `yield name, (seq, None)`.
    ///
    /// Asked of the reference rather than reasoned about -- the first version of
    /// this test asserted `Some("III")` and was wrong. Downstream this is
    /// Finding 11: the caller iterates the quality string and dies with
    /// `TypeError: 'NoneType' object is not iterable`.
    #[test]
    fn a_missing_final_newline_drops_the_quality_entirely() {
        let with = read("@r1\nACGT\n+\nIIII\n");
        let without = read("@r1\nACGT\n+\nIIII");
        assert_eq!(with[0].qual.as_deref(), Some("IIII"));
        assert_eq!(without[0].qual, None);
    }

    /// Quality shorter than the sequence: same path, same outcome.
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

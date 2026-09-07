//! 2-bit packed nucleotide sequences.
//!
//! Sequences are the largest single thing the port keeps resident: 865 MB of a
//! 2238 MB accounted heap on SIRV_real_full, one raw ASCII byte per base. Two
//! bits per base takes that to 229 MB.
//!
//! # This breaks byte-identity on non-ACGT input, deliberately
//!
//! `A`, `C`, `G` and `T` round-trip exactly. **Anything else -- `N`, lowercase,
//! IUPAC ambiguity codes -- is packed as `A` and comes back as `A`**, because two
//! bits cannot represent a fifth symbol. That is observable in at least two
//! places: `get_kmer_minimizers` picks minimizers by lexicographic order on the
//! k-mer string, and `'N'` (78) sorts between `'G'` (71) and `'T'` (84), so an
//! `N` changes which k-mer wins a window; and parasail's matrix is built for
//! `"ACGT"`, scoring anything else as 0 rather than as a match, which changes the
//! alignment path.
//!
//! This is an accepted divergence, not an oversight. Every corpus in
//! `bench/corpora.tsv` is pure ACGT -- checked, zero non-ACGT bases -- so the
//! equivalence harness cannot see it, which is *Finding 5*'s lesson and the
//! reason `from_bytes` counts the substitutions and the loader reports them. An
//! earlier design kept a sparse table of exception positions to stay exact; it
//! was dropped as deliberate scope, on the grounds that ONT and PacBio basecalls
//! are ACGT.
use std::sync::Arc;

/// A sequence, four bases to the byte. `len` is the base count, which the packed
/// length cannot imply because it rounds up to a whole byte.
#[derive(Clone)]
pub struct PackedSeq {
    bits: Arc<[u8]>,
    len: u32,
}

#[inline]
fn code(b: u8) -> u8 {
    match b {
        b'A' => 0,
        b'C' => 1,
        b'G' => 2,
        b'T' => 3,
        // Everything else becomes `A`; see the module note.
        _ => 0,
    }
}

const DECODE: [u8; 4] = *b"ACGT";

impl PackedSeq {
    /// Pack `seq`, returning the sequence and how many bases were not ACGT.
    pub fn from_bytes(seq: &[u8]) -> (Self, usize) {
        let mut bits = vec![0u8; seq.len().div_ceil(4)];
        let mut substituted = 0usize;
        for (i, &b) in seq.iter().enumerate() {
            if !matches!(b, b'A' | b'C' | b'G' | b'T') {
                substituted += 1;
            }
            // Base 0 in the high bits, so the packed order reads left to right.
            bits[i / 4] |= code(b) << (6 - 2 * (i % 4));
        }
        (
            Self {
                bits: bits.into(),
                len: u32::try_from(seq.len()).expect("a read shorter than 4 Gbp"),
            },
            substituted,
        )
    }

    pub fn len(&self) -> usize {
        self.len as usize
    }

    #[allow(dead_code)]
    pub fn is_empty(&self) -> bool {
        self.len == 0
    }

    /// Unpack into `out`, replacing its contents.
    ///
    /// Callers reuse one buffer per loop rather than allocating: the whole point
    /// of packing is to stop paying per-read allocation, and this runs once per
    /// read and once per alignment candidate.
    pub fn unpack_into(&self, out: &mut Vec<u8>) {
        out.clear();
        out.reserve(self.len());
        for i in 0..self.len() {
            let byte = self.bits[i / 4];
            out.push(DECODE[((byte >> (6 - 2 * (i % 4))) & 0b11) as usize]);
        }
    }

    /// Unpack to a fresh `Vec`. For the output writers, which run once per
    /// surviving cluster rather than once per read.
    pub fn to_bytes(&self) -> Vec<u8> {
        let mut v = Vec::new();
        self.unpack_into(&mut v);
        v
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn acgt_round_trips_at_every_length_offset() {
        // Lengths 0..=32 cover every position within the 4-base byte, including
        // the partial final byte.
        let alphabet = b"ACGT";
        for n in 0..=32usize {
            let seq: Vec<u8> = (0..n).map(|i| alphabet[i % 4]).collect();
            let (p, sub) = PackedSeq::from_bytes(&seq);
            assert_eq!(sub, 0);
            assert_eq!(p.len(), n);
            assert_eq!(p.to_bytes(), seq, "length {n}");
        }
    }

    #[test]
    fn packing_is_a_quarter_of_the_bytes() {
        let (p, _) = PackedSeq::from_bytes(&vec![b'C'; 1000]);
        assert_eq!(p.bits.len(), 250);
    }

    /// The documented divergence, pinned so it cannot change silently.
    #[test]
    fn non_acgt_becomes_a_and_is_counted() {
        let (p, sub) = PackedSeq::from_bytes(b"ACGNTacgt");
        assert_eq!(sub, 5, "N and the four lowercase bases");
        assert_eq!(p.to_bytes(), b"ACGATAAAA");
    }

    #[test]
    fn unpack_into_reuses_the_buffer() {
        let (a, _) = PackedSeq::from_bytes(b"ACGTACGT");
        let (b, _) = PackedSeq::from_bytes(b"TTTT");
        let mut buf = Vec::new();
        a.unpack_into(&mut buf);
        assert_eq!(buf, b"ACGTACGT");
        b.unpack_into(&mut buf);
        assert_eq!(buf, b"TTTT", "the previous contents must not linger");
    }
}

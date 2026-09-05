//! The phred probability tables, frozen from the reference environment.
//!
//! GENERATED -- do not hand-edit. Regenerate with the snippet in bench/README.md.
//!
//! WHY THESE ARE CONSTANTS RATHER THAN COMPUTED
//! --------------------------------------------
//! The reference builds them as `10 ** (-(ord(c) - 33) / 10.0)`, and computing
//! the same expression in Rust does NOT give the same bits. Measured on this
//! machine, for `%` (phred 4, i.e. `10 ** -0.4`):
//!
//! | how                                   | bits               |
//! |---------------------------------------|--------------------|
//! | CPython `10 ** -0.4`                  | `3fd97a967f7524b2` |
//! | Rust `10f64.powf(x)` at runtime       | `3fd97a967f7524b3` |
//! | C `pow(10.0, -0.4)` on the same libm  | `3fd97a967f7524b3` |
//! | Rust `powf` on a *constant* (folded)  | `3fd97a967f7524b2` |
//!
//! One value of the 95 printable ones differs, by one ULP -- and it is `%`,
//! a perfectly ordinary ONT quality character. Simulated reads and PacBio CCS
//! reads do not contain it, so every corpus except real ONT agreed; real SIRV
//! and Drosophila reads diverged in the 13th significant digit of the score.
//!
//! The deeper point is that the *reference* is not portable here either: its
//! value comes from whatever libm its interpreter was built against, and this
//! interpreter disagrees with this machine's own C library. Freezing the table
//! makes the port reproducible everywhere and identical to the pinned
//! reference; computing it would make the port agree with neither reliably.
//! `rust/tests/phred_oracle.rs` re-checks the frozen values against the live
//! reference, so a changed environment is a test failure rather than a silent
//! drift.
//!
//! See PORTING.md, Finding 12.
//!
//! Stored as raw `u64` bit patterns rather than `f64` literals: `f64::from_bits`
//! only became usable in a `const` context in Rust 1.83, and there is no reason
//! for this port to require a toolchain that new. Converting on access is free,
//! and the bits make it obvious that these are exact values copied from the
//! reference rather than something re-derived here.

/// `D` in the reference: capped at 0.79433. Used by the rolling no-error
/// product in `expected_number_of_erroneous_kmers`, and by `cluster.py`'s
/// `phred_char_to_p`.
static CAPPED_BITS: [u64; 128] = [
    0x3fe96b26bf8769ec,
    0x3fe96b26bf8769ec,
    0x3fe96b26bf8769ec,
    0x3fe96b26bf8769ec,
    0x3fe96b26bf8769ec,
    0x3fe96b26bf8769ec,
    0x3fe96b26bf8769ec,
    0x3fe96b26bf8769ec,
    0x3fe96b26bf8769ec,
    0x3fe96b26bf8769ec,
    0x3fe96b26bf8769ec,
    0x3fe96b26bf8769ec,
    0x3fe96b26bf8769ec,
    0x3fe96b26bf8769ec,
    0x3fe96b26bf8769ec,
    0x3fe96b26bf8769ec,
    0x3fe96b26bf8769ec,
    0x3fe96b26bf8769ec,
    0x3fe96b26bf8769ec,
    0x3fe96b26bf8769ec,
    0x3fe96b26bf8769ec,
    0x3fe96b26bf8769ec,
    0x3fe96b26bf8769ec,
    0x3fe96b26bf8769ec,
    0x3fe96b26bf8769ec,
    0x3fe96b26bf8769ec,
    0x3fe96b26bf8769ec,
    0x3fe96b26bf8769ec,
    0x3fe96b26bf8769ec,
    0x3fe96b26bf8769ec,
    0x3fe96b26bf8769ec,
    0x3fe96b26bf8769ec,
    0x3fe96b26bf8769ec,
    0x3fe96b26bf8769ec,
    0x3fe96b230bcdc434,
    0x3fe430cd74f6d478,
    0x3fe009b9cf334252,
    0x3fd97a967f7524b2,
    0x3fd43d136248490f,
    0x3fd0137987dd704c,
    0x3fc98a13577c93c0,
    0x3fc44960c576b375,
    0x3fc01d3f2d9684d0,
    0x3fb999999999999a,
    0x3fb455b5a30b035c,
    0x3fb0270ac3f8a9fa,
    0x3fa9a9294b8536e9,
    0x3fa46211ff90ea2a,
    0x3fa030dc4ea03a72,
    0x3f99b8c272fbe6de,
    0x3f946e75df96dc9a,
    0x3f903ab3d12bc2c4,
    0x3f89c86515bda14e,
    0x3f847ae147ae147b,
    0x3f8044914f3c02b0,
    0x3f79d811398ddcc0,
    0x3f7487543c6a9257,
    0x3f704e74cc73ee88,
    0x3f69e7c6e43390b7,
    0x3f6493cec2631f18,
    0x3f60585e4c78b079,
    0x3f59f7861b7937a3,
    0x3f54a050de314dd8,
    0x3f50624dd2f1a9fc,
    0x3f4a074ee52cd119,
    0x3f44acda94717d66,
    0x3f406c4363887513,
    0x3f3a1721471fe40d,
    0x3f34b96be9c2da2c,
    0x3f30763f01e8e5ad,
    0x3f2a26fd472780c1,
    0x3f24c604e2c75fb6,
    0x3f208040b1c10b13,
    0x3f1a36e2eb1c432d,
    0x3f14d2a58423da81,
    0x3f108a4876c1311e,
    0x3f0a46d238da54eb,
    0x3f04df4dd27fe99e,
    0x3f009456549be1bd,
    0x3efa56cb36416f83,
    0x3ef4ebfdd286009a,
    0x3ef09e6a4f05e62b,
    0x3eea66cde934de7d,
    0x3ee4f8b588e368f1,
    0x3ee0a88469b64867,
    0x3eda76da579b81ca,
    0x3ed50574fa4843ef,
    0x3ed0b2a4a866547e,
    0x3eca86f0875fcf94,
    0x3ec5123c2b678c69,
    0x3ec0bccb0ed19a15,
    0x3eba97107e6fd6ab,
    0x3eb51f0b20f71864,
    0x3eb0c6f7a0b5ed8d,
    0x3eaaa73a42bd40a5,
    0x3ea52be1dfaf9b08,
    0x3ea0d12a61d3698c,
    0x3e9ab76dda3d53fd,
    0x3e9538c06c4ca610,
    0x3e90db6355ec7054,
    0x3e8ac7ab4ae8f688,
    0x3e8545a6cb8cabbc,
    0x3e80e5a280c5ad1d,
    0x3e7ad7f29abcaf48,
    0x3e75529502310084,
    0x3e70efe7e62615a0,
    0x3e6ae843cfb8a8e0,
    0x3e655f8b14fddcca,
    0x3e60fa3389d6eb40,
    0x3e5af89eefe0b3b9,
    0x3e556c8908ba5ed3,
    0x3e5104856fa3bc97,
    0x3e4b0904013c482e,
    0x3e45798ee2308c3a,
    0x3e410edd9b5a66d0,
    0x3e3b197309d68910,
    0x3e35869ca62d53da,
    0x3e31193c10cb1708,
    0x3e2b29ec0fbe4534,
    0x3e2593b259808fc7,
    0x3e2123a0d3c84be6,
    0x3e1b3a6f1905fa7b,
    0x3e15a0d000fd068c,
    0x3e112e0be826d695,
    0x3e0b4afc2bc3d7b3,
    0x3e05adf5a1786da6,
    0x3e01387d51bddcae,
    0x3dfb5b934e11be74,
];

/// `D_no_min` in the reference: no cap. Used for `error_rate`, which gates
/// the `--q` filter and is written into `final_cluster_origins.tsv`.
static UNCAPPED_BITS: [u64; 128] = [
    0x409f2d0c9c4b9258,
    0x4098c392a10b6613,
    0x4093abb39f263d21,
    0x408f400000000000,
    0x4088d2a03986f199,
    0x4083b7a8a4390b7b,
    0x407f52fee8b01d8c,
    0x4078e1b6f87865d8,
    0x4073c3a4edfa9759,
    0x406f66095d5c7f53,
    0x4068f0d6e36fa846,
    0x4063cfa880d5eb43,
    0x405f791f6509fb68,
    0x4059000000000000,
    0x4053dbb36138c148,
    0x404f8c4106c1abfc,
    0x40490f3253c017a0,
    0x4043e7c5939384ad,
    0x403f9f6e4990f227,
    0x40391e6de449ff75,
    0x4033f3df1c59536e,
    0x402fb2a734897866,
    0x40292db2b73b2f86,
    0x4024000000000000,
    0x401fc5ebcec13542,
    0x40193d00d2348997,
    0x40140c28430012e6,
    0x400fd93c1f526ddf,
    0x40094c583ada5b53,
    0x40041857e9d4cc5f,
    0x3fffec982d5bb8af,
    0x3ff95bb8f6d46053,
    0x3ff4248ef8fc2604,
    0x3ff0000000000000,
    0x3fe96b230bcdc434,
    0x3fe430cd74f6d478,
    0x3fe009b9cf334252,
    0x3fd97a967f7524b2,
    0x3fd43d136248490f,
    0x3fd0137987dd704c,
    0x3fc98a13577c93c0,
    0x3fc44960c576b375,
    0x3fc01d3f2d9684d0,
    0x3fb999999999999a,
    0x3fb455b5a30b035c,
    0x3fb0270ac3f8a9fa,
    0x3fa9a9294b8536e9,
    0x3fa46211ff90ea2a,
    0x3fa030dc4ea03a72,
    0x3f99b8c272fbe6de,
    0x3f946e75df96dc9a,
    0x3f903ab3d12bc2c4,
    0x3f89c86515bda14e,
    0x3f847ae147ae147b,
    0x3f8044914f3c02b0,
    0x3f79d811398ddcc0,
    0x3f7487543c6a9257,
    0x3f704e74cc73ee88,
    0x3f69e7c6e43390b7,
    0x3f6493cec2631f18,
    0x3f60585e4c78b079,
    0x3f59f7861b7937a3,
    0x3f54a050de314dd8,
    0x3f50624dd2f1a9fc,
    0x3f4a074ee52cd119,
    0x3f44acda94717d66,
    0x3f406c4363887513,
    0x3f3a1721471fe40d,
    0x3f34b96be9c2da2c,
    0x3f30763f01e8e5ad,
    0x3f2a26fd472780c1,
    0x3f24c604e2c75fb6,
    0x3f208040b1c10b13,
    0x3f1a36e2eb1c432d,
    0x3f14d2a58423da81,
    0x3f108a4876c1311e,
    0x3f0a46d238da54eb,
    0x3f04df4dd27fe99e,
    0x3f009456549be1bd,
    0x3efa56cb36416f83,
    0x3ef4ebfdd286009a,
    0x3ef09e6a4f05e62b,
    0x3eea66cde934de7d,
    0x3ee4f8b588e368f1,
    0x3ee0a88469b64867,
    0x3eda76da579b81ca,
    0x3ed50574fa4843ef,
    0x3ed0b2a4a866547e,
    0x3eca86f0875fcf94,
    0x3ec5123c2b678c69,
    0x3ec0bccb0ed19a15,
    0x3eba97107e6fd6ab,
    0x3eb51f0b20f71864,
    0x3eb0c6f7a0b5ed8d,
    0x3eaaa73a42bd40a5,
    0x3ea52be1dfaf9b08,
    0x3ea0d12a61d3698c,
    0x3e9ab76dda3d53fd,
    0x3e9538c06c4ca610,
    0x3e90db6355ec7054,
    0x3e8ac7ab4ae8f688,
    0x3e8545a6cb8cabbc,
    0x3e80e5a280c5ad1d,
    0x3e7ad7f29abcaf48,
    0x3e75529502310084,
    0x3e70efe7e62615a0,
    0x3e6ae843cfb8a8e0,
    0x3e655f8b14fddcca,
    0x3e60fa3389d6eb40,
    0x3e5af89eefe0b3b9,
    0x3e556c8908ba5ed3,
    0x3e5104856fa3bc97,
    0x3e4b0904013c482e,
    0x3e45798ee2308c3a,
    0x3e410edd9b5a66d0,
    0x3e3b197309d68910,
    0x3e35869ca62d53da,
    0x3e31193c10cb1708,
    0x3e2b29ec0fbe4534,
    0x3e2593b259808fc7,
    0x3e2123a0d3c84be6,
    0x3e1b3a6f1905fa7b,
    0x3e15a0d000fd068c,
    0x3e112e0be826d695,
    0x3e0b4afc2bc3d7b3,
    0x3e05adf5a1786da6,
    0x3e01387d51bddcae,
    0x3dfb5b934e11be74,
];

/// `D[c]` -- capped.
#[inline]
pub fn capped(c: u8) -> f64 {
    f64::from_bits(CAPPED_BITS[c as usize])
}

/// `D_no_min[c]` -- uncapped.
#[inline]
pub fn uncapped(c: u8) -> f64 {
    f64::from_bits(UNCAPPED_BITS[c as usize])
}

/// Every entry, for the oracle in `tests/phred_oracle.rs`. Not used by the
/// binary itself, which only ever looks up single characters.
#[allow(dead_code)]
pub fn all() -> impl Iterator<Item = (usize, f64, f64)> {
    (0..128).map(|i| {
        (
            i,
            f64::from_bits(CAPPED_BITS[i]),
            f64::from_bits(UNCAPPED_BITS[i]),
        )
    })
}

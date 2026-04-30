extern crate libc;
extern crate serde;

use c2rust_bitfields::BitfieldStruct;
use serde::{Deserialize, Serialize};
use std::cmp::Ordering;
use std::convert::From;

#[repr(C, align(1))]
#[derive(BitfieldStruct, Clone, Copy, Default, Debug, PartialEq, PartialOrd)]
pub struct Var32 {
    /* enc: 8
     * alen: 2
     * rlen: 2
     * position: 20 */
    #[bitfield(name = "enc", ty = "libc::uint32_t", bits = "0..=7")]
    #[bitfield(name = "alen", ty = "libc::uint32_t", bits = "8..=9")]
    #[bitfield(name = "rlen", ty = "libc::uint32_t", bits = "10..=11")]
    #[bitfield(name = "position", ty = "libc::uint32_t", bits = "12..=31")]
    data: [u8; 4],
}

const fn init_bitset() -> u128 {
    (1 << 'A' as usize) | (1 << 'C' as usize) | (1 << 'G' as usize) | (1 << 'T' as usize)
}

// use a bitset to check if incoming letters are ACGT and issue warning if not.
const DNA_BITS: u128 = init_bitset();

// TODO: since ref[0] == alt[0], we can get one more base.
#[repr(C, align(1))]
#[derive(Clone, Copy, Default, Debug, PartialEq, PartialOrd)]
pub struct PRA {
    position: u32,
    reference: [char; 3],
    alternate: [char; 3],
}

#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct LongVariant {
    pub position: u32,
    pub idx: u32,
    pub sequence: Vec<u32>,
}

// implement this as we need to exclude idx from the eq.
impl PartialEq for LongVariant {
    fn eq(&self, other: &Self) -> bool {
        self.position == other.position && self.sequence == other.sequence
    }
}

impl Eq for LongVariant {}

// implement this as we need to exclude idx from the eq.
impl PartialOrd for LongVariant {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        if self.position != other.position {
            return self.position.partial_cmp(&other.position);
        }
        self.sequence.partial_cmp(&other.sequence)
    }
}

impl Ord for LongVariant {
    fn cmp(&self, other: &Self) -> Ordering {
        self.partial_cmp(other).unwrap()
    }
}

pub const MAX_COMBINED_LEN: usize = 4;

pub(crate) const LOOKUP: [u32; 128] = [
    3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
    3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
    3, 0, 3, 1, 3, 3, 3, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
    3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
];

pub(crate) const RLOOKUP: [char; 4] = ['A', 'C', 'G', 'T'];

impl From<u32> for Var32 {
    #[inline]
    fn from(u: u32) -> Self {
        unsafe { std::mem::transmute::<u32, Var32>(u) }
    }
}
impl From<Var32> for u32 {
    #[inline]
    fn from(v: Var32) -> Self {
        unsafe { std::mem::transmute::<Var32, u32>(v) }
    }
}

#[inline(always)]
fn bit_test(set: u128, bit: usize) -> bool {
    set & (1 << bit as u32) != 0
}

#[inline]
pub fn encode(pos: u32, ref_allele: &[u8], alt_allele: &[u8], warn: &mut i32) -> u32 {
    let mut v: Var32 = Var32::default();
    v.set_position(pos);

    if ref_allele.len() + alt_allele.len() > MAX_COMBINED_LEN {
        // too large to encode. but we signal that by setting both lengths to 3.
        v.set_alen(3);
        v.set_rlen(3);
        return u32::from(v);
    }
    v.set_alen(alt_allele.len() as u32);
    v.set_rlen(ref_allele.len() as u32);
    let mut ra: u32 = 0;

    for a in ref_allele.iter() {
        if !bit_test(DNA_BITS, *a as usize) && *warn < 10 {
            *warn += 1;
            eprintln!(
                "[warning] found non ACGT REF character '{}', encoding as 'T' for (1-based) position: {}",
                *a as char, pos + 1
            );
        }
        ra *= 4;
        ra += LOOKUP[*a as usize];
    }

    for a in alt_allele.iter() {
        if !bit_test(DNA_BITS, *a as usize) && *warn < 10 {
            *warn += 1;
            eprintln!(
                "[warning] found non ACGT ALT character '{}', encoding as 'T' for (1-based) position: {}",
                *a as char, pos + 1
            );
        }
        ra *= 4;
        ra += LOOKUP[*a as usize];
    }

    v.set_enc(ra);

    u32::from(v)
}

/// Decodes a packed `Var32` encoding back into reference and alternate allele sequences.
///
/// Used by **BED/tab annotation**: when annotating BED or tabular input in position-scan
/// mode (no REF/ALT columns), the annotation path calls this to turn stored Var32 entries
/// into `(ref, alt)` pairs so it can enumerate and annotate every variant at each position.
///
/// The low 8 bits of the encoding store bases as 2 bits each (A=0, C=1, G=2, T=3), with
/// reference bases first, then alternate bases. Reference and alternate lengths come from
/// the `Var32` layout (2 bits each for `rlen` and `alen`).
///
/// Returns `None` when both lengths are 3, which is the sentinel used for variants that
/// are too long to fit in 32 bits and are stored as [`LongVariant`] instead.
///
/// Returns `Some((reference, alt))` where each sequence is a `Vec<u8>` of ASCII bytes
/// `b'A'`, `b'C'`, `b'G'`, or `b'T'`.
pub fn decode_to_alleles(enc: u32) -> Option<(Vec<u8>, Vec<u8>)> {
    let v: Var32 = Var32::from(enc);
    let rlen = v.rlen() as usize;
    let alen = v.alen() as usize;
    if rlen == 3 && alen == 3 {
        return None;
    } // sentinel for LongVariant

    let mut e = v.enc();
    let mut alt = vec![0u8; alen];
    for i in (0..alen).rev() {
        alt[i] = RLOOKUP[(e & 3) as usize] as u8;
        e >>= 2;
    }
    let mut reference = vec![0u8; rlen];
    for i in (0..rlen).rev() {
        reference[i] = RLOOKUP[(e & 3) as usize] as u8;
        e >>= 2;
    }
    Some((reference, alt))
}

#[inline]
pub fn decode(enc: u32) -> PRA {
    let v: Var32 = unsafe { std::mem::transmute::<u32, Var32>(enc) };

    let mut result = PRA {
        position: v.position(),
        ..Default::default()
    };

    let mut e = v.enc();
    let h = (v.alen() - 1) as usize;
    for i in 0..=h {
        result.alternate[h - i] = RLOOKUP[(e & 3) as usize];
        e >>= 2;
    }

    let h = (v.rlen() - 1) as usize;
    for i in 0..=h {
        result.reference[h - i] = RLOOKUP[(e & 3) as usize];
        e >>= 2;
    }
    result
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_encode() {
        let mut w = 0;
        let e = encode(423432, b"A", b"ACA", &mut w);
        assert_eq!(e, 1734379268);
    }

    #[test]
    fn test_decode() {
        let d = decode(1734379268);
        assert_eq!(d.position, 423432);

        eprintln!("{:?}", d);

        assert_eq!(d.reference, ['A', '\0', '\0']);
        assert_eq!(d.alternate, ['A', 'C', 'A']);
    }

    #[test]
    fn test_size() {
        assert_eq!(std::mem::size_of::<Var32>(), 4);
    }

    #[test]
    fn test_ordering() {
        let mut a = Var32 {
            ..Default::default()
        };
        a.set_position(12334);
        a.set_alen(3);
        a.set_rlen(2);
        a.set_enc(63);

        let mut b = Var32 {
            ..Default::default()
        };
        b.set_position(12333);
        b.set_alen(3);
        b.set_rlen(2);
        b.set_enc(63);

        assert_eq!(true, a > b);
        b.set_position(a.position());
        assert_eq!(a, b);

        b.set_alen(2);
        assert_eq!(true, a > b);

        b.set_alen(a.alen());

        b.set_enc(a.enc() - 1);
        assert_eq!(true, a > b);
    }

    #[test]
    fn test_decode_to_alleles_roundtrip() {
        let mut w = 0;
        let enc = encode(423432, b"A", b"ACA", &mut w);
        let result = decode_to_alleles(enc);
        assert!(result.is_some());
        let (r, a) = result.unwrap();
        assert_eq!(r, b"A");
        assert_eq!(a, b"ACA");
    }

    #[test]
    fn test_decode_to_alleles_snv() {
        let mut w = 0;
        let enc = encode(100, b"A", b"T", &mut w);
        let result = decode_to_alleles(enc);
        assert!(result.is_some());
        let (r, a) = result.unwrap();
        assert_eq!(r, b"A");
        assert_eq!(a, b"T");
    }

    #[test]
    fn test_decode_to_alleles_sentinel() {
        let mut w = 0;
        // Combined length > 4 triggers sentinel (rlen=3, alen=3)
        let enc = encode(100, b"AAA", b"CCC", &mut w);
        let result = decode_to_alleles(enc);
        assert!(result.is_none());
    }

    #[test]
    fn test_var32() {
        let mut v = Var32 {
            ..Default::default()
        };
        v.set_position(12334);
        v.set_alen(3);
        v.set_rlen(2);
        v.set_enc(63);

        assert_eq!(v.position(), 12334);
        assert_eq!(v.alen(), 3);
        assert_eq!(v.rlen(), 2);
        assert_eq!(v.enc(), 63);

        let w: u32 = v.into();

        assert_eq!(unsafe { std::mem::transmute::<Var32, u32>(v) }, 50522943);
        assert_eq!(w, 50522943);

        let vv: Var32 = w.into();
        v.set_alen(1);
        assert_ne!(v, vv);
    }
}

extern crate libc;
extern crate serde;
use c2rust_bitfields::BitfieldStruct;
use serde::{Deserialize, Serialize};
use std::cmp::Ordering;
use std::convert::From;

#[repr(C, align(4))]
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

pub const MAX_COMBINED_LEN: usize = 4;
pub const MAX_COMBINED_LEN_64: usize = 17;

#[repr(C, align(8))]
#[derive(BitfieldStruct, Clone, Copy, Default, Debug, PartialEq, Serialize, Deserialize)]
// TODO: store u64 as variant separate from u32 indexes.
pub struct Var64 {

    /* 64 bits 
     - 20 for position == 44 bits
     - 5 each for rlen, alen (for max of 31 bases each) == 34
     leaves 34 for enc.
    enc: 34
     
    */
    #[bitfield(name = "enc", ty = "libc::uint64_t", bits = "0..=33")]
    #[bitfield(name = "alen", ty = "libc::uint64_t", bits = "34..=38")]
    #[bitfield(name = "rlen", ty = "libc::uint64_t", bits = "39..=43")]
    #[bitfield(name = "position", ty = "libc::uint64_t", bits = "44..=63")]
    data: [u8; 8],
}

#[repr(C, align(8))]
#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct LongVariant {
    pub position: u32,
    pub idx: u32,
    pub sequence: Vec<u32>,
}

// TODO: since ref[0] == alt[0], we can get one more base.

#[repr(C, align(1))]
#[derive(Clone, Copy, Default, Debug, PartialEq, PartialOrd)]
pub struct PRA {
    position: u32,
    reference: [char; 3],
    alternate: [char; 3],
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
        return self.sequence.partial_cmp(&other.sequence);
    }
}

impl Ord for LongVariant {
    fn cmp(&self, other: &Self) -> Ordering {
        self.partial_cmp(&other).unwrap()
    }
}


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

impl From<u64> for Var64 {
    #[inline]
    fn from(u: u64) -> Self {
        unsafe { std::mem::transmute::<u64, Var64>(u) }
    }
}

impl From<Var64> for u64 {
    #[inline]
    fn from(v: Var64) -> Self {
        unsafe { std::mem::transmute::<Var64, u64>(v) }
    }
}

#[inline]
pub fn encode(pos: u32, ref_allele: &[u8], alt_allele: &[u8]) -> u32 {
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
        ra *= 4;
        ra += LOOKUP[*a as usize];
    }

    for a in alt_allele.iter() {
        ra *= 4;
        ra += LOOKUP[*a as usize];
    }

    v.set_enc(ra);

    return u32::from(v);
}

#[inline]
pub fn encode64(pos: u32, ref_allele: &[u8], alt_allele: &[u8]) -> u64 {
    let mut v: Var64 = Var64::default();
    v.set_position(pos as u64);

    if ref_allele.len() + alt_allele.len() > MAX_COMBINED_LEN_64 {
        // too large to encode. but we signal that by setting both lengths to 3.
        v.set_alen(0);
        v.set_rlen(0);
        return u64::from(v);
    }
    v.set_alen(alt_allele.len() as u64);
    v.set_rlen(ref_allele.len() as u64);
    let mut ra: u64 = 0;

    for a in ref_allele.iter() {
        ra *= 4;
        ra += LOOKUP[*a as usize] as u64;
    }

    for a in alt_allele.iter() {
        ra *= 4;
        ra += LOOKUP[*a as usize] as u64;
    }

    v.set_enc(ra);

    return u64::from(v);
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
        let e = encode(423432, b"A", b"ACA");
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
        assert_eq!(std::mem::size_of::<Var64>(), 8);
        assert_eq!(std::mem::size_of::<LongVariant>(), 32);
        /*
        let s = Vec::new();
        let l = LongVariant{idx: 0, position: 0, sequence: s};
        assert_eq!(bincode::DefaultOptions::new().serialized_size(&l).expect("WTF"), 32);
        */
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

    #[test]
    fn test_var64() {
        let mut v = Var64 {
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

        /*
        let w: u64 = v.into();
        assert_eq!(unsafe { std::mem::transmute::<Var64, u64>(v) }, w);
        */

    }
}

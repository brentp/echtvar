extern crate libc;
use std::convert::From;
use c2rust_bitfields::BitfieldStruct;

#[repr(C, align(1))]
#[derive(BitfieldStruct, Clone, Copy, Default, Debug, PartialEq)]
pub struct Var32 {
    /* enc: 8
     * alen: 2
     * rlen: 2
     * position: 22 */
    #[bitfield(name = "enc", ty = "libc::uint32_t", bits = "0..=7")]
    #[bitfield(name = "alen", ty = "libc::uint32_t", bits = "8..=9")]
    #[bitfield(name = "rlen", ty = "libc::uint32_t", bits = "10..=11")]
    #[bitfield(name = "position", ty = "libc::uint32_t", bits = "12..=31")]
    data: [u8; 4]
}

impl From<u32> for Var32 {
	fn from(u:u32) -> Self {
		return unsafe { std::mem::transmute::<u32, Var32>(u) };
	}
}
impl From<Var32> for u32 {
	fn from(v: Var32) -> Self {
		return unsafe { std::mem::transmute::<Var32, u32>(v) };
	}
}


#[cfg(test)]
mod tests {
	use super::*;

    #[test]
    fn test_size() {
        assert_eq!(std::mem::size_of::<Var32>(), 4);
    }

	#[test]
	fn test_var32() {

		let mut v = Var32{ ..Default::default() };
		v.set_position(12334);
		v.set_alen(3);
		v.set_rlen(2);
		v.set_enc(63);

		assert_eq!(v.position(), 12334);
		assert_eq!(v.alen(), 3);
		assert_eq!(v.rlen(), 2);
		assert_eq!(v.enc(), 63);

		let w:u32 = v.into();

		assert_eq!(unsafe { std::mem::transmute::<Var32, u32>(v) }, 50522943);
		assert_eq!(w, 50522943);

		let vv:Var32 = w.into();
		v.set_alen(1);
		assert_ne!(v, vv);

	}
}



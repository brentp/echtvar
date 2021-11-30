
#[inline]
pub fn encode(i: i32) -> u32 {
   return ((i << 1) ^ (i >> 31)) as u32
}

#[inline]
pub fn decode(e: u32) -> i32 {
   return ((e >> 1) as i32) ^ -((e & 1) as i32)
}


#[cfg(test)]
mod tests {
	use super::*;
	#[test]
	fn test_roundtrip() {
        let values : [i32; 6] = [0, -10, 33, -123, -123213, i32::MAX];

        for v in values {
            let e = encode(v);
            let d = decode(e);
            //eprintln!("{}, {}, {}", e, e == u32::MAX, u32::MAX);
            assert_eq!(d, v);
        }

    }
}
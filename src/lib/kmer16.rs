
const LOOKUP: [u32; 128] = [
    3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
    3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
    3, 0, 3, 1, 3, 3, 3, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
    3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
];

pub type K16 = u32;

pub type K16s = Vec<K16>;

pub fn encode(dna: &[u8]) -> K16s {

	let mut result: Vec<K16> = vec![0; 2 + ((dna.len() - 1) >> 3)];
	result[0] = dna.len() as u32;

    for (i, a) in dna.iter().enumerate() {
		let idx = 1 + (i>>3);
		result[idx] *= 4;
        result[idx] += LOOKUP[*a as usize];
    }
	result
}



#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_encode() {
        let dna = b"ACAAAA";
        let e = encode(dna);
        assert_eq!(e.len(), 2);
        assert_eq!(e[0], dna.len() as u32)
    }

    #[test]
    fn test_equal() {
        let dna = b"ACAAAA";
        let a = encode(dna);
        let dnb = b"ACAAAA";
        let b = encode(dnb);
		assert_eq!(a, b);
	}

    #[test]
    fn test_ne() {
        let dna = b"ACAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
        let a = encode(dna);
        let dnb = b"ACAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAG";
        let b = encode(dnb);
		assert_ne!(a, b);
	}

	#[test]
    fn test_single_base() {
        let dna = b"A";
        let a = encode(dna);
        assert_eq!(a.len(), 2);
        let dnb = b"T";
        let b = encode(dnb);
		assert_ne!(a, b)
    }

}


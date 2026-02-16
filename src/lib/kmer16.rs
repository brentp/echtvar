const LOOKUP: [u32; 128] = [
    3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
    3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
    3, 0, 3, 1, 3, 3, 3, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
    3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
];

pub const RLOOKUP: [u8; 4] = [b'A', b'C', b'G', b'T'];

pub fn decode_var(sequence: &[u32]) -> (Vec<u8>, Vec<u8>) {
    let ref_len = sequence[0] as usize;
    let alt_len = sequence[1] as usize;
    let total = ref_len + alt_len;
    let mut bases = Vec::with_capacity(total);
    for i in 0..total {
        let word_idx = 2 + (i >> 4);
        let bases_in_word = std::cmp::min(16, total - (i & !15));
        let bit_in_word = i & 15;
        let shift = (bases_in_word - 1 - bit_in_word) * 2;
        let code = (sequence[word_idx] >> shift) & 3;
        bases.push(RLOOKUP[code as usize]);
    }
    let ref_allele = bases[..ref_len].to_vec();
    let alt_allele = bases[ref_len..].to_vec();
    (ref_allele, alt_allele)
}

pub fn encode_var(ref_allele: &[u8], alt_allele: &[u8]) -> Vec<u32> {
    let mut result: Vec<u32> = vec![0; 3 + ((alt_allele.len() + ref_allele.len() - 1) >> 4)];
    result[0] = ref_allele.len() as u32;
    result[1] = alt_allele.len() as u32;

    let mut i = 0;

    for a in ref_allele.iter() {
        let idx = 2 + (i >> 4);
        result[idx] <<= 2;
        result[idx] += LOOKUP[*a as usize];
        i += 1;
    }
    for a in alt_allele.iter() {
        let idx = 2 + (i >> 4);
        result[idx] <<= 2;
        result[idx] += LOOKUP[*a as usize];
        i += 1;
    }

    result
}

pub fn encode(dna: &[u8]) -> Vec<u32> {
    let mut result: Vec<u32> = vec![0; 2 + ((dna.len() - 1) >> 4)];
    result[0] = dna.len() as u32;

    for (i, a) in dna.iter().enumerate() {
        let idx = 1 + (i >> 4);
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

    #[test]
    fn test_pair() {
        let ra = b"AAAAC";
        let aa = b"A";

        let a = encode_var(ra, aa);
        assert_eq!(a.len(), 3);

        let b = encode_var(aa, ra);
        assert_eq!(b.len(), 3);

        assert_ne!(a, b)
    }

    #[test]
    fn test_decode_var_roundtrip_single() {
        let ref_a = b"A";
        let alt_a = b"C";
        let enc = encode_var(ref_a, alt_a);
        let (r, a) = decode_var(&enc);
        assert_eq!(r, ref_a);
        assert_eq!(a, alt_a);
    }

    #[test]
    fn test_decode_var_roundtrip_16_bases() {
        let ref_a = b"ACGTACGT";
        let alt_a = b"TGCATGCA";
        let enc = encode_var(ref_a, alt_a);
        let (r, a) = decode_var(&enc);
        assert_eq!(r, ref_a);
        assert_eq!(a, alt_a);
    }

    #[test]
    fn test_decode_var_roundtrip_multi_word() {
        let ref_a = b"ACGTACGTACGTACGTA";
        let alt_a = b"G";
        let enc = encode_var(ref_a, alt_a);
        let (r, a) = decode_var(&enc);
        assert_eq!(r, ref_a);
        assert_eq!(a, alt_a);
    }

    #[test]
    fn test_decode_var_roundtrip_both_long() {
        let ref_a = b"ACGTACGTACGTACGTAC";
        let alt_a = b"TGCATGCATGCATGCATG";
        let enc = encode_var(ref_a, alt_a);
        let (r, a) = decode_var(&enc);
        assert_eq!(r, ref_a);
        assert_eq!(a, alt_a);
    }
}

use crate::fields;
use crate::kmer16;
use crate::var32;
use crate::zigzag;
use bincode::Options;
use rust_htslib::bcf;
use std::io::prelude::*;
use std::{fs, io, str};

use byteorder::{LittleEndian, ReadBytesExt};
use std::io::BufReader;

use stream_vbyte::{decode::decode, x86::Ssse3};

use ieee754::Ieee754;

#[derive(Debug, Clone, Copy)]
pub enum Value {
    Int(i32),
    Float(f32),
}

impl Value {
    pub fn value(self) -> f64 {
        match self {
            Value::Int(i) => i as f64,
            Value::Float(f) => f as f64,
        }
    }
}

#[derive(Debug)]
pub struct EchtVars {
    pub zip: zip::ZipArchive<std::fs::File>,
    pub chrom: String,
    last_rid: i32,
    pub start: u32,
    pub var32s: Vec<u32>,
    pub longs: Vec<var32::LongVariant>,
    // the values for a chunk are stored in values.
    pub values: Vec<Vec<u32>>,

    // for storing values used by fasteval
    pub evalues: Vec<Value>,
    // values.len() == fields.len() and fields[i] indicates how we
    // handle values[i]
    pub fields: Vec<fields::Field>,
    buffer: Vec<u8>,

	warn: i32,

    // lookup for categorical fields.
    // these will be empty for non-categorical.
    pub strings: Vec<Vec<std::string::String>>,
}

pub trait Variant {
    fn chrom(&self) -> std::string::String;
    fn rid(&self) -> i32;
    fn position(&self) -> u32;
    fn alleles(&self) -> Vec<&[u8]>;
}

#[inline]
pub fn strip_chr(chrom: std::string::String) -> std::string::String {
    if chrom.len() < 4 {
        return chrom;
    }
    let bchrom = chrom.as_bytes();
    if bchrom[0] as char == 'c' && bchrom[1] as char == 'h' && bchrom[2] as char == 'r' {
        return chrom[3..].to_string();
    }
    return chrom;
}

#[inline]
pub fn bstrip_chr(chrom: &str) -> &str {
    if chrom.len() < 4 {
        return chrom;
    }
    let bchrom = chrom.as_bytes();
    if bchrom[0] as char == 'c' && bchrom[1] as char == 'h' && bchrom[2] as char == 'r' {
        return &chrom[3..];
    }
    return chrom;
}

impl Variant for bcf::record::Record {
    fn chrom(&self) -> std::string::String {
        let rid = self.rid().unwrap();
        let n: &[u8] = self.header().rid2name(rid as u32).unwrap();
        str::from_utf8(n).unwrap().to_string()
    }

    fn rid(&self) -> i32 {
        self.rid().unwrap() as i32
    }

    fn position(&self) -> u32 {
        self.pos() as u32
    }

    fn alleles(&self) -> Vec<&[u8]> {
        self.alleles()
    }
}

impl EchtVars {
    pub fn open(path: &str) -> Self {
        let ep = std::path::Path::new(&*path);
        let file = fs::File::open(ep).unwrap_or_else(|e| panic!("error accessing zip file: \"{}\" ({})", path, e));
        let mut result = EchtVars {
            zip: zip::ZipArchive::new(file).unwrap_or_else(|e| panic!("error opening: \"{}\" as zip file ({})", path, e)),
            chrom: "".to_string(),
            last_rid: -1,
            start: u32::MAX,
            var32s: vec![],
            longs: vec![],
            values: vec![],
            evalues: vec![],
            fields: vec![],
            buffer: vec![],
            strings: vec![],
			warn: 0,
        };

        {
            let mut fc = result
                .zip
                .by_name("echtvar/config.json")
                .expect("unable to open echtvar/config.json");
            let mut contents = String::new();
            fc.read_to_string(&mut contents)
                .expect("eror reading config.json");
            drop(fc);
            let mut flds: Vec<fields::Field> = json5::from_str(&contents).unwrap();
            for fld in flds.iter_mut() {
                fld.values_i = result.fields.len();
                if fld.ftype == fields::FieldType::Categorical {
                    // read in the strings for this field. replace ';' with ',' to handle the filter field.
                    let fname = format!("echtvar/strings/{}.txt", fld.alias);
                    let fh = result
                        .zip
                        .by_name(&fname)
                        .expect("error opening strings file");
                    result.strings.push(
                        BufReader::new(fh)
                            .lines()
                            .map(|l| l.unwrap().replace(";", ","))
                            .collect(),
                    );
                    // update missing value to be the index of the missing_string
                    let strings_len = result.strings[result.strings.len() - 1].len();
                    fld.missing_value = result.strings[result.strings.len() - 1]
                        .iter()
                        .position(|s| s == &fld.missing_string)
                        .unwrap_or(strings_len) as i32;
                    // if it wasn't in the list, add it.
                    if fld.missing_value == strings_len as i32 {
                        let rl = result.strings.len() - 1;
                        result.strings[rl].push(fld.missing_string.clone());
                    }
                } else {
                    result.strings.push(Vec::new());
                }
                let f = fld.clone();
                result.fields.push(f);
            }
            result.values.resize(result.fields.len(), vec![]);
            result.evalues.resize(result.fields.len(), Value::Int(0));
        }
        eprintln!("fields: {:?}", result.fields);
        result
    }

    pub fn update_header(self: &mut EchtVars, header: &mut bcf::header::Header, path: &str) {
        for e in &self.fields {
            header.push_record(
                format!(
                    "##INFO=<ID={},Number={},Type={},Description=\"{}\">",
                    e.alias,
                    e.number,
                    if e.ftype == fields::FieldType::Integer {
                        "Integer"
                    } else if e.ftype == fields::FieldType::Categorical {
                        "String"
                    } else {
                        "Float"
                    },
                    if &e.description.to_string() == "added by echtvar"{
                        format!("added by echtvar from {}", path)
                    } else {
                        e.description.to_string()
                    }
                )
                .as_bytes(),
            );
        }
    }

    #[inline(always)]
    pub fn set_position(
        self: &mut EchtVars,
        rid: i32,
        chromosome: String,
        position: u32,
    ) -> io::Result<()> {
        if rid == self.last_rid && position >> 20 == self.start >> 20 {
            return Ok(());
        }
        self.last_rid = rid;
        self.start = position >> 20 << 20; // round to 20 bits.
        self.chrom = strip_chr(chromosome);
        let base_path = format!("echtvar/{}/{}", self.chrom, position >> 20);

        for fi in self.fields.iter_mut() {
            // RUST-TODO: use .fill function. problems with double borrow.
            let path = format!("{}/{}.bin", base_path, fi.alias);
            let rzip = self.zip.by_name(&path);
            match rzip {
                Ok(mut iz) => {
                    let n = iz.read_u32::<LittleEndian>()? as usize;
                    self.buffer
                        .resize(iz.size() as usize - std::mem::size_of::<u32>(), 0x0);
                    iz.read_exact(&mut self.buffer)?;
                    self.values[fi.values_i].resize(n, 0x0);
                    // TODO: use skip to first position.
                    let bytes_decoded =
                        decode::<Ssse3>(&self.buffer, n, &mut self.values[fi.values_i]);

                    if bytes_decoded != self.buffer.len() {
                        return Err(std::io::Error::new(
                            std::io::ErrorKind::Other,
                            "didn't read expected number of values from zip",
                        ));
                    }
                }
                _ => {
                    // if we don't find this file we set everything to empty.
                    self.values[fi.values_i].clear();
                    self.var32s.clear();
                }
            }
        }

        if self.values[0].len() > 0 {
            let path = format!("{}/var32.bin", base_path);
            let mut iz = self.zip.by_name(&path)?;
            let n = iz.read_u32::<LittleEndian>()? as usize;
            //eprintln!("n:{}", n);
            self.buffer
                .resize(iz.size() as usize - std::mem::size_of::<u32>(), 0x0);
            iz.read_exact(&mut self.buffer)?;

            self.var32s.resize(n, 0x0);
            let bytes_decoded = decode::<Ssse3>(&self.buffer, n, &mut self.var32s);
            // cumsum https://users.rust-lang.org/t/inplace-cumulative-sum-using-iterator/56532/3
            self.var32s.iter_mut().fold(0, |acc, x| {
                *x += acc;
                *x
            });

            if bytes_decoded != self.buffer.len() {
                return Err(std::io::Error::new(
                    std::io::ErrorKind::Other,
                    "didn't read expected number of values from zip",
                ));
            }
        }

        if self.var32s.len() > 0 {
            let long_path = format!("{}/too-long-for-var32.enc", base_path);
            let mut iz = self.zip.by_name(&long_path)?;
            self.buffer.clear();
            iz.read_to_end(&mut self.buffer)?;
            self.longs = bincode::DefaultOptions::new()
                .deserialize(&self.buffer)
                .expect("error decoding long variants");
        } else {
            self.longs.clear();
        }

        Ok(())
    }

    #[inline]
    fn get_int_value(self: &EchtVars, fld: &fields::Field, idx: usize) -> i32 {
        let v: u32 = self.values[fld.values_i][idx];
        return if v == u32::MAX {
            fld.missing_value as i32
        } else {
            if fld.zigzag {
                zigzag::decode(v) as i32
            } else {
                v as i32
            }
        };
    }

    #[inline]
    fn get_float_value(self: &EchtVars, fld: &fields::Field, idx: usize) -> f32 {
        let v: u32 = self.values[fld.values_i][idx];
        return if v == u32::MAX {
            if fld.missing_value == 0x7F800001 { Ieee754::from_bits(0x7F800001) } else { fld.missing_value as f32 }
        } else {
            if fld.zigzag {
                (zigzag::decode(v) as f32) / (fld.multiplier as f32)
            } else {
                (v as f32) / (fld.multiplier as f32)
            }
        };
    }

    pub fn update_expr_values<T: Variant>(
        self: &mut EchtVars,
        variant: &mut T,
        expr_values: &mut Vec<f64>,
    ) {
        let pos = variant.position();
        let rid = variant.rid();
        if rid != self.last_rid || pos >> 20 != self.start >> 20 {
            let chrom = variant.chrom();
            let _ = self.set_position(rid, chrom, pos);
        }

        let alleles = variant.alleles();
        if alleles.len() != 2 {
            panic!(
                "[echtvar] variants must be decomposed before running. got variant with {} alleles at {}:{} ({:?}). see: https://github.com/brentp/echtvar/wiki/decompose",
                alleles.len() - 1,
                variant.chrom(),
                variant.position() + 1,
                variant.alleles()
            );
        }
        let eidx = if alleles[0].len() + alleles[1].len() <= crate::var32::MAX_COMBINED_LEN {
            let enc = var32::encode(pos, alleles[0], alleles[1], &mut self.warn);
            self.var32s.binary_search(&enc)
        } else {
            let l = var32::LongVariant {
                position: pos,
                sequence: kmer16::encode_var(alleles[0], alleles[1]),
                idx: 0,
            };
            let r = self.longs.binary_search(&l);
            match r {
                Ok(idx) => Ok(self.longs[idx].idx as usize),
                Err(_) => Err(0),
            }
        };
        match eidx {
            Ok(idx) => {
                for fld in &self.fields {
                    if fld.ftype == fields::FieldType::Integer
                        || fld.ftype == fields::FieldType::Categorical
                    {
                        let val = self.get_int_value(fld, idx);
                        self.evalues[fld.values_i] = Value::Int(val);
                        expr_values[fld.values_i] = val as f64
                    } else if fld.ftype == fields::FieldType::Float {
                        let val = self.get_float_value(fld, idx);
                        self.evalues[fld.values_i] = Value::Float(val);
                        expr_values[fld.values_i] = val as f64
                    } else {
                        panic!("not implemented");
                    }
                }
            }
            Err(_) => {
                for fld in &self.fields {
                    if fld.ftype == fields::FieldType::Integer
                        || fld.ftype == fields::FieldType::Categorical
                    {
                        // for Categorical missing_value has been set to the index of missing_string
                        let val = fld.missing_value as i32;
                        self.evalues[fld.values_i] = Value::Int(val);
                        expr_values[fld.values_i] = val as f64
                    } else if fld.ftype == fields::FieldType::Float {
                        let val = if fld.missing_value == 0x7F800001 { Ieee754::from_bits(0x7F800001) } else { fld.missing_value as f32 };
                        self.evalues[fld.values_i] = Value::Float(val);
                        expr_values[fld.values_i] = val as f64
                    }
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    //#[test]
    fn test_read() {
        let mut e = EchtVars::open("ec.zip");
        e.set_position(22, "chr21".to_string(), 5030088).ok();

        assert_eq!(e.fields.len(), 3);
        assert_eq!(e.values[0].len(), 46881);
        assert_eq!(e.values[1].len(), e.var32s.len());

        assert_eq!(e.longs[0].position, 5030185);
    }

    //#[test]
    fn test_search() {
        let mut e = EchtVars::open("ec.zip");
        e.set_position(22, "chr21".to_string(), 5030088).ok();

        let mut vals = vec![];
        vals.resize(3, 0.0);

        pub struct Var<'a> {
            chrom: std::string::String, //b"chr21"
            pos: u32,                   // 5030087,
            alleles: Vec<&'a [u8]>,     //vec!["C", "T"],
        }

        impl<'a> Variant for Var<'a> {
            fn chrom(&self) -> std::string::String {
                self.chrom.clone()
            }
            fn position(&self) -> u32 {
                self.pos
            }
            fn rid(&self) -> i32 {
                1
            }
            fn alleles(&self) -> Vec<&[u8]> {
                self.alleles.clone()
            }
        }

        let mut variant = Var {
            chrom: "chr21".to_string(),
            pos: 5030087,
            alleles: vec![b"C", b"T"],
        };

        let idx = e.update_expr_values(&mut variant, &mut vals);
        eprintln!("vals:{:?} {:?}", vals, idx);
        assert_eq!(vals[1], 2.0);
    }
}

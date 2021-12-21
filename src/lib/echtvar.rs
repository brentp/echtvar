use crate::fields;
use crate::var32;
use crate::zigzag;
use rust_htslib::bcf;
use std::io::prelude::*;
use std::{fs, io, str};

use byteorder::{LittleEndian, ReadBytesExt};

use stream_vbyte::{decode::decode, x86::Ssse3};

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
    pub evalues: Vec<f64>,
    // values.len() == fields.len() and fields[i] indicates how we
    // handle values[i]
    pub fields: Vec<fields::Field>,
    buffer: Vec<u8>,
}

impl EchtVars {
    pub fn open(path: &str) -> Self {
        let ep = std::path::Path::new(&*path);
        let file = fs::File::open(ep).expect("error accessing zip file");
        let mut result = EchtVars {
            zip: zip::ZipArchive::new(file).expect("error opening zip file"),
            chrom: "".to_string(),
            last_rid: -1,
            start: u32::MAX,
            var32s: vec![],
            longs: vec![],
            values: vec![],
            evalues: vec![],
            fields: vec![],
            buffer: vec![],
        };

        {
            let mut f = result
                .zip
                .by_name("echtvar/config.json")
                .expect("unable to open echtvar/config.json");
            let mut contents = String::new();
            f.read_to_string(&mut contents)
                .expect("eror reading config.json");
            let flds: Vec<fields::Field> = json5::from_str(&contents).unwrap();
            eprintln!("fields: {:?}", flds);
            for fld in flds {
                let mut f = fld.clone();
                f.values_i = result.fields.len();
                result.fields.push(f);
            }
            result.values.resize(result.fields.len(), vec![]);
            result.evalues.resize(result.fields.len(), 0.0);
        }
        result
    }

    pub fn update_header(self: &mut EchtVars, header: &mut bcf::header::Header) {
        for e in &self.fields {
            header.push_record(
                format!(
                    "##INFO=<ID={},Type={},Number=1,Description=\"{}\">",
                    e.alias,
                    if e.ftype == fields::FieldType::Integer {
                        "Integer"
                    } else {
                        "Float"
                    },
                    "added by echtvar from ..."
                )
                .as_bytes(),
            );
        }
    }

    /*
    pub fn fill(self: &EchtVars, fi: &mut EchtVar<u32>, path: String) -> io::Result<()> {
        //eprintln!("path:{}", path);
        let mut iz = self.zip.by_name(&path)?;
        let n = iz.read_u32::<LittleEndian>()? as usize;
        //eprintln!("n:{}", n);
        self.buffer
            .resize(iz.size() as usize - std::mem::size_of::<u32>(), 0x0);
        iz.read_exact(&mut self.buffer)?;
        fi.values.resize(n, 0x0);
        // TODO: use skip to first position.
        let bytes_decoded = decode::<Ssse3>(&self.buffer, n, &mut fi.values);

        if bytes_decoded != self.buffer.len() {
            return Err(std::io::Error::new(
                std::io::ErrorKind::Other,
                "didn't read expected number of values from zip",
            ));
        }
        Ok(())
    }
    */

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
        self.chrom = chromosome;
        let base_path = format!("echtvar/{}/{}", self.chrom, position >> 20);

        for fi in self.fields.iter_mut() {
            // RUST-TODO: use .fill function. problems with double borrow.
            let path = format!("{}/{}.bin", base_path, fi.alias);
            //self.fill(fi, path)?;
            let rzip = self.zip.by_name(&path);
            match rzip {
                Ok(mut iz) => {
                    let n = iz.read_u32::<LittleEndian>()? as usize;
                    //eprintln!("n:{}", n);
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
            let long_path = format!("{}/too-long-for-var32.txt", base_path);
            let mut iz = self.zip.by_name(&long_path)?;
            self.buffer.clear();
            iz.read_to_end(&mut self.buffer)?;
            self.longs = serde_json::from_slice(&self.buffer)?;
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
    fn get_float_value(self: &EchtVars, fld: &fields::Field, idx: usize) -> f32 {
        let v: u32 = self.values[fld.values_i][idx];
        return if v == u32::MAX {
            fld.missing_value as f32
        } else {
            if fld.zigzag {
                (zigzag::decode(v) as f32) / fld.multiplier as f32
            } else {
                v as f32 / fld.multiplier as f32
            }
        };
    }

    pub fn check_and_update_variant(
        self: &mut EchtVars,
        variant: &mut bcf::record::Record,
        header: &bcf::header::HeaderView,
    ) -> bool {
        let pos = variant.pos() as u32;
        let rid = variant.rid().unwrap() as i32;
        if rid != self.last_rid || pos >> 20 != self.start >> 20 {
            let n: &[u8] = header.rid2name(rid as u32).unwrap();
            let chrom = str::from_utf8(n).unwrap().to_string();
            let _ = self.set_position(rid, chrom, pos);
        }

        // TODO prototype fasteval stuff here.
        let alleles = variant.alleles();
        let eidx = if alleles[0].len() + alleles[1].len() <= crate::var32::MAX_COMBINED_LEN {
            let enc = var32::encode(pos, alleles[0], alleles[1]);
            self.var32s.binary_search(&enc)
        } else {
            let l = var32::LongVariant {
                position: pos,
                reference: unsafe { str::from_utf8_unchecked(alleles[0]) }.to_string(),
                alternate: unsafe { str::from_utf8_unchecked(alleles[1]) }.to_string(),
                idx: 0,
            };
            let r = self.longs.binary_search(&l);
            match r {
                Ok(idx) => Ok(self.longs[idx].idx as usize),
                Err(_) => Err(0),
            }
            //Ok(0)
            //Err(0)
        };
        match eidx {
            Ok(idx) => {
                for fld in &self.fields {
                    if fld.ftype == fields::FieldType::Integer {
                        let val = [self.get_int_value(fld, idx)];
                        variant
                            .push_info_integer(fld.alias.as_bytes(), &val)
                            .expect(&format!("error adding integer {}", fld.alias).to_string());
                    } else if fld.ftype == fields::FieldType::Float {
                        let val = [self.get_float_value(fld, idx)];
                        variant
                            .push_info_float(fld.alias.as_bytes(), &val)
                            .expect(&format!("error adding float {}", fld.alias).to_string());
                    } else {
                        panic!("not implemented");
                    }
                }
            }
            Err(_) => {
                for fld in &self.fields {
                    if fld.ftype == fields::FieldType::Integer {
                        let val = [fld.missing_value as i32];
                        variant
                            .push_info_integer(fld.alias.as_bytes(), &val)
                            .expect(&format!("error adding integer {}", fld.alias).to_string());
                    } else if fld.ftype == fields::FieldType::Float {
                        let val = [fld.missing_value as f32];
                        variant
                            .push_info_float(fld.alias.as_bytes(), &val)
                            .expect(&format!("error adding float {}", fld.alias).to_string());
                    }
                }
            }
        }

        true
    }

    pub fn values(
        self: &mut EchtVars,
        rid: i32,
        chromosome: &[u8],
        position: u32,
        reference: &[u8],
        alternate: &[u8],
        values: &mut Vec<i32>,
    ) -> io::Result<()> {
        self.set_position(
            rid,
            unsafe { std::str::from_utf8_unchecked(chromosome).to_string() },
            position,
        )?;

        let e = var32::encode(position, reference, alternate);
        values.clear();

        let eidx = self.var32s.binary_search(&e);
        match eidx {
            Ok(idx) => {
                for e in &self.fields {
                    let v: u32 = self.values[e.values_i][idx];
                    if v == u32::MAX {
                        values.push(e.missing_value as i32);
                    } else {
                        if e.zigzag {
                            values.push(zigzag::decode(v) as i32 / e.multiplier as i32);
                        } else {
                            values.push(v as i32 / e.multiplier as i32);
                        }
                    }
                }
            }
            Err(_) => {
                // variant nDDDDDot found. fill with missing values.
                for e in &self.fields {
                    values.push(e.missing_value as i32);
                }
            }
        };

        eprintln!("r:{:?}, {}, {:?}", eidx, e, &self.var32s[..10]);

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_read() {
        let mut e = EchtVars::open("ec.zip");
        e.set_position(22, "chr21".to_string(), 5030088).ok();

        assert_eq!(e.fields.len(), 2);
        assert_eq!(e.values[0].len(), 46881);
        assert_eq!(e.values[1].len(), e.var32s.len());

        assert_eq!(e.longs[0].position, 5030185);
    }

    #[test]
    fn test_search() {
        let mut e = EchtVars::open("ec.zip");
        e.set_position(22, "chr21".to_string(), 5030088).ok();

        let mut vals = vec![];

        let idx = e.values(1, b"chr21", 5030087, b"C", b"T", &mut vals).ok();
        eprintln!("{:?} {:?}", vals, idx);
    }
}

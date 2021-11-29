extern crate bitpacking;
extern crate byteorder;
extern crate rust_htslib;
use crate::rust_htslib::bcf::record::{Buffer, Record};
use crate::rust_htslib::bcf::{Read, Reader};
use bitpacking::{BitPacker, BitPacker4x};
use echtvar_lib;
use echtvar_lib::var32;
use echtvar_lib::zigzag;
use std::borrow::{Borrow, BorrowMut};
use std::io::{self, Write};

use byteorder::{LittleEndian, ReadBytesExt, WriteBytesExt};
use std::str;

use zip::write::FileOptions;

#[inline]
fn get_int_field<'a, B: BorrowMut<Buffer> + Borrow<Buffer> + 'a>(
    rec: &Record,
    field: &[u8],
    buffer: B,
) -> i32 {
    return match rec
        .info_shared_buffer(field, buffer)
        .integer()
        .expect("error reading info")
    {
        Some(v) => v[0],
        None => -1,
    };
}

fn write_bits(
    values: &Vec<u32>,
    sorted: bool,
    zipf: &mut zip::ZipWriter<std::fs::File>,
    bitpacker: &BitPacker4x,
    compressed: &mut [u8],
) {
    // TODO: append zeros until multiple of 128.
    // TODO: allow sorted
    // TODO: errors
    zipf.write_u32::<LittleEndian>(values.len() as u32).ok();
    for i in (0..values.len() - 128).step_by(128) {
        let num_bits;
        let compressed_len;
        if sorted {
            num_bits = bitpacker.num_bits_sorted(values[i], &values[i..i + 128]);
            compressed_len = bitpacker.compress_sorted(
                values[i],
                &values[i..i + 128],
                &mut compressed[..],
                num_bits,
            );
        } else {
            num_bits = bitpacker.num_bits(&values[i..i + 128]);
            compressed_len = bitpacker.compress(&values[i..i + 128], &mut compressed[..], num_bits);
        }
        eprintln!("sorted: {}, vars:{}:{} num_bits:{}", sorted, i, i+128, num_bits);
        zipf.write_all(&compressed[..compressed_len]).ok();
    }
	// bitpacker writes blocks of 128 u32's. 
	// we write any leftovers as-is.
    let remaining = values.len() % 128;
	eprintln!("remaining: {} start: {}", remaining, values.len() - remaining);
	for i in (values.len() - remaining) .. values.len() {
		zipf.write_u32::<LittleEndian>(values[i]).ok();
	}
}

fn main() {
    if std::env::args().len() < 3 {
        println!("expecting arguments: <vcf> <zip>")
    }
    let args: Vec<String> = std::env::args().collect();
    let path = &*args[1];
    let zpath = std::path::Path::new(&*args[2]);

    eprintln!("{}", path);
    let mut vcf = Reader::from_path(path).ok().expect("Error opening vcf.");
    vcf.set_threads(2).ok();
    let header = vcf.header().clone();
    let mut buffer = Buffer::new();

    let mut acs: Vec<u32> = Vec::new();
    let bitpacker = BitPacker4x::new();

    let zfile = std::fs::File::create(&zpath).unwrap();

    let mut zipf = zip::ZipWriter::new(zfile);

    zipf.add_directory("echtvar/", FileOptions::default().large_file(true))
        .expect("error writing zip");

    let options = FileOptions::default()
        .compression_method(zip::CompressionMethod::Deflated)
        .unix_permissions(0o755);

    let mut compressed = vec![0u8; 4 * BitPacker4x::BLOCK_LEN];
    let mut last_rid: i32 = -1;
    let mut last_mod: i64 = 0;

    let mut evars: Vec<u32> = Vec::new();

    for r in vcf.records() {
        let rec = r.expect("error getting record");
        // if we hit a new chrom or a new chunk we write the last chunk and start a new one.
        if rec.rid().expect("no rid found") as i32 != last_rid || rec.pos() >> 20 != last_mod {
            if last_rid != -1 {
                if acs.len() != 0 {
                    let n: &[u8] = header.rid2name(last_rid as u32).unwrap();
                    let chrom = str::from_utf8(n).unwrap();
                    let fname = format!("echtvar/{}/{}/AC.bin", chrom, last_mod);
                    zipf.start_file(fname, options)
                        .expect("error starting file");
                    write_bits(&mut acs, false, &mut zipf, &bitpacker, &mut compressed);

                    let fname = format!("echtvar/{}/{}/var32.bin", chrom, last_mod);
                    zipf.start_file(fname, options)
                        .expect("error starting file");
                    write_bits(&mut evars, true, &mut zipf, &bitpacker, &mut compressed);
					acs.clear();
					evars.clear();
                }
            }
            last_rid = rec.rid().unwrap() as i32;
            last_mod = rec.pos() >> 20;
        }
        let ac = get_int_field(&rec, b"AC", &mut buffer);
        acs.push(ac as u32);

        let alleles = rec.alleles();

        evars.push(var32::encode(rec.pos() as u32, alleles[0], alleles[1]));
    }
    if acs.len() != 0 {
        let n: &[u8] = header.rid2name(last_rid as u32).unwrap();
        let chrom = str::from_utf8(n).unwrap();
        let fname = format!("echtvar/{}/{}/AC.bin", chrom, last_mod);
        zipf.start_file(fname, options)
            .expect("error starting file");
        write_bits(&mut acs, false, &mut zipf, &bitpacker, &mut compressed);

        let fname = format!("echtvar/{}/{}/var32.bin", chrom, last_mod);
        zipf.start_file(fname, options)
            .expect("error starting file");
        write_bits(&mut evars, true, &mut zipf, &bitpacker, &mut compressed);
    }
    zipf.finish().expect("error closing zip file");
}

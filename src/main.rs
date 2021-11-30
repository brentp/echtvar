extern crate bitpacking;
extern crate byteorder;
extern crate rust_htslib;
use crate::rust_htslib::bcf::record::{Buffer, Record};
use crate::rust_htslib::bcf::{Read, Reader};
use bitpacking::{BitPacker, BitPacker4x as BitPackerImpl};
use echtvar_lib;
use echtvar_lib::var32;
//use echtvar_lib::zigzag;
use std::borrow::{Borrow, BorrowMut};
use std::io::{Write};

use byteorder::{LittleEndian, WriteBytesExt};
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

fn write_long(zipf: &mut zip::ZipWriter<std::fs::File>, long_vars: &Vec<var32::LongVariant>) {
    eprintln!("writing {} longs", long_vars.len());
    serde_json::to_writer(zipf, &long_vars).expect("error writing long variables");
}

fn write_bits(
    values: &Vec<u32>,
    sorted: bool,
    zipf: &mut zip::ZipWriter<std::fs::File>,
    bitpacker: &BitPackerImpl,
    compressed: &mut [u8],
) {
    // TODO: errors
    zipf.write_u32::<LittleEndian>(values.len() as u32).ok();
    let mut last_value: u32 = 0;
    for i in (0..values.len() - BitPackerImpl::BLOCK_LEN).step_by(BitPackerImpl::BLOCK_LEN) {
        let num_bits;
        let compressed_len;
        if sorted {
            num_bits =
                bitpacker.num_bits_sorted(last_value, &values[i..i + BitPackerImpl::BLOCK_LEN]);
            compressed_len = bitpacker.compress_sorted(
                last_value,
                &values[i..i + BitPackerImpl::BLOCK_LEN],
                &mut compressed[..],
                num_bits,
            );
            last_value = values[i + BitPackerImpl::BLOCK_LEN - 1];
        } else {
            num_bits = bitpacker.num_bits(&values[i..i + BitPackerImpl::BLOCK_LEN]);
            compressed_len = bitpacker.compress(
                &values[i..i + BitPackerImpl::BLOCK_LEN],
                &mut compressed[..],
                num_bits,
            );
        }
        zipf.write_all(&compressed[..compressed_len]).ok();
    }
    // bitpacker writes blocks of BitPackerImpl::BLOCK_LEN u32's.
    // we write any leftovers as-is.
    let remaining = values.len() % BitPackerImpl::BLOCK_LEN;
    eprintln!(
        "total: {}, remaining: {} start: {}",
        values.len(),
        remaining,
        values.len() - remaining
    );
    for i in (values.len() - remaining)..values.len() {
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

    eprintln!("{} {}", path, BitPackerImpl::BLOCK_LEN);
    let mut vcf = Reader::from_path(path).ok().expect("Error opening vcf.");
    vcf.set_threads(2).ok();
    let header = vcf.header().clone();
    let mut buffer = Buffer::new();

    let bitpacker = BitPackerImpl::new();

    let zfile = std::fs::File::create(&zpath).unwrap();

    let mut zipf = zip::ZipWriter::new(zfile);

    zipf.add_directory("echtvar/", FileOptions::default().large_file(true))
        .expect("error writing zip");

    let options = FileOptions::default()
        .compression_method(zip::CompressionMethod::Deflated)
        .unix_permissions(0o755);

    let mut compressed = vec![0u8; 4 * BitPackerImpl::BLOCK_LEN];
    let mut last_rid: i32 = -1;
    let mut last_mod: i64 = 0;

    let mut long_vars: Vec<var32::LongVariant> = Vec::new();
    let mut var32s: Vec<u32> = Vec::new();

    let fields: Vec<&[u8]> = vec![b"AC", b"AN"];
    let mut values_vv: Vec<Vec<u32>> = fields.iter().map(|_| Vec::new()).collect();

    for r in vcf.records() {
        let rec = r.expect("error getting record");
        // if we hit a new chrom or a new chunk we write the last chunk and start a new one.
        if rec.rid().expect("no rid found") as i32 != last_rid || rec.pos() >> 20 != last_mod {
            if last_rid != -1 {
                if values_vv[0].len() != 0 {
                    let n: &[u8] = header.rid2name(last_rid as u32).unwrap();
                    let chrom = str::from_utf8(n).unwrap();

                    for (i, values) in values_vv.iter_mut().enumerate() {
                        let fname = format!("echtvar/{}/{}/{}.bin", chrom, last_mod, std::str::from_utf8(fields[i]).unwrap());
                        zipf.start_file(fname, options)
                            .expect("error starting file");
                        write_bits(values, false, &mut zipf, &bitpacker, &mut compressed);
                        values.clear();
                    }

                    let fname = format!("echtvar/{}/{}/var32.bin", chrom, last_mod);
                    zipf.start_file(fname, options)
                        .expect("error starting file");
                    write_bits(&mut var32s, true, &mut zipf, &bitpacker, &mut compressed);
                    var32s.clear();

                    let fname = format!("echtvar/{}/{}/too-long-for-var32.txt", chrom, last_mod);
                    zipf.start_file(fname, options)
                        .expect("error starting file");
                    write_long(&mut zipf, &long_vars);
                    long_vars.clear()
                }
            }
            last_rid = rec.rid().unwrap() as i32;
            last_mod = rec.pos() >> 20;
        }

        for (i, fld) in fields.iter().enumerate() {
            let v = get_int_field(&rec, fld, &mut buffer);
            values_vv[i].push(v as u32);
        }

        let alleles = rec.alleles();
        var32s.push(var32::encode(rec.pos() as u32, alleles[0], alleles[1]));

        if alleles[0].len() + alleles[1].len() > var32::MAX_COMBINED_LEN {
            long_vars.push(var32::LongVariant {
                position: rec.pos() as u32,
                reference: unsafe { str::from_utf8_unchecked(alleles[0]).to_string() },
                alternate: unsafe { str::from_utf8_unchecked(alleles[1]).to_string() },
                idx: (var32s.len() - 1) as u32,
            });
        }
    }
    if values_vv[0].len() != 0 {
        let n: &[u8] = header.rid2name(last_rid as u32).unwrap();
        let chrom = str::from_utf8(n).unwrap();

        for (i, values) in values_vv.iter_mut().enumerate() {
            let fname = format!("echtvar/{}/{}/{}.bin", chrom, last_mod, std::str::from_utf8(fields[i]).unwrap());
            zipf.start_file(fname, options)
                .expect("error starting file");
            write_bits(values, false, &mut zipf, &bitpacker, &mut compressed);
            values.clear();
        }

        let fname = format!("echtvar/{}/{}/var32.bin", chrom, last_mod);
        zipf.start_file(fname, options)
            .expect("error starting file");
        write_bits(&mut var32s, true, &mut zipf, &bitpacker, &mut compressed);

        let fname = format!("echtvar/{}/{}/too-long-for-var32.txt", chrom, last_mod);
        zipf.start_file(fname, options)
            .expect("error starting file");
        write_long(&mut zipf, &long_vars);
    }
    zipf.finish().expect("error closing zip file");
}

extern crate rust_htslib;
extern crate bitpacking;
use crate::rust_htslib::bcf::{Reader, Read};
use crate::rust_htslib::bcf::record::{Record, Buffer};
use std::borrow::{Borrow, BorrowMut};
use std::io::Write;
use std::str;
use echtvar_lib::zigzag;
use echtvar_lib::var32;
use bitpacking::{BitPacker4x, BitPacker};
use echtvar_lib;

use zip::write::FileOptions;

#[cfg(not(target_env = "msvc"))]
use jemallocator::Jemalloc;

#[cfg(not(target_env = "msvc"))]
#[global_allocator]
static GLOBAL: Jemalloc = Jemalloc;


//fn get_float_field<'a, B: BorrowMut<Buffer> + Borrow<Buffer> + 'a>(rec: &Record, field: &[u8], buffer: B) -> f32 {
fn get_float_field<'a, B: BorrowMut<Buffer> + 'a>(rec: &Record, field: &[u8], buffer: B) -> f32 {
	return match rec.info_shared_buffer(field, buffer).float().expect("error reading info") {
		Some(v) => v[0],
		None => -1.0,
	};
}

#[inline]
fn get_int_field<'a, B: BorrowMut<Buffer> + Borrow<Buffer> + 'a>(rec: &Record, field: &[u8], buffer: B) -> i32 {
	return match rec.info_shared_buffer(field, buffer).integer().expect("error reading info") {
		Some(v) => v[0],
		None => -1,
	};
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

	let mut acs:Vec<u32> = Vec::new();
    let bitpacker = BitPacker4x::new();

	let zfile = std::fs::File::create(&zpath).unwrap();

	let mut zip = zip::ZipWriter::new(zfile);

	zip.add_directory("echtvar/", FileOptions::default().large_file(true)).expect("error writing zip");

	let options = FileOptions::default()
	    .compression_method(zip::CompressionMethod::Deflated)
		.unix_permissions(0o755);

	let mut compressed = vec![0u8; 4 * BitPacker4x::BLOCK_LEN];
	let mut last_rid:i32 = -1;

	let mut evars : Vec<u32> = Vec::new();

	for r in vcf.records() {
		let rec = r.expect("error getting record");
		if rec.rid().expect("no rid found") as i32 != last_rid {
  		  let n: &[u8] = header.rid2name(rec.rid().expect("bad rid")).unwrap();
		  let chrom = str::from_utf8(n).unwrap();
		  let fname = format!("echtvar/{}/AC.bin", chrom);
	      zip.start_file(fname, options).expect("error starting file");
		  last_rid = rec.rid().unwrap() as i32;
		}
		let ac = get_int_field(&rec, b"AC", & mut buffer);
		acs.push(ac as u32);

		let alleles = rec.alleles();

		evars.push(var32::encode(rec.pos() as u32, alleles[0], alleles[1]));

        if acs.len() == 128 {
			/*
			let mut ics = vec![0u32; BitPacker4x::BLOCK_LEN];
			for (i, v) in acs.iter().enumerate() {
				ics[i] = zigzag::encode(*v as i32 - acs[0] as i32);
			}
			//eprintln!("{:?}", ics);
			*/

            let num_bits = bitpacker.num_bits(&acs[..]);
            let compressed_len = bitpacker.compress(&acs[..], &mut compressed[..], num_bits);
			//let compressed_len = (num_bits * 16) as usize;
            eprintln!("uncompressed bits: {}, compressed bits: {}, num_bits: {}", acs.len() * std::mem::size_of::<u32>(), compressed_len * std::mem::size_of::<u8>(), num_bits);
			//zip.write_all(&compressed[..compressed_len]).expect("OK");

			/*
			let mut decompressed = vec![0u32; BitPacker4x::BLOCK_LEN];
			bitpacker.decompress(&compressed[..compressed_len], &mut decompressed[..], num_bits);
			if decompressed != acs {
				eprintln!("expected same values!");

			}
			*/


            acs.clear();
        }

		//eprintln!("{}:{} AF: {} AC: {}", std::str::from_utf8(n).unwrap(), rec.pos(), af, ac).expect("error writing to stderr");

	}

	zip.start_file("echtvar/21/variants.bin", options).expect("error starting file");
	for i in (0..evars.len()-128).step_by(128) {
        let num_bits = bitpacker.num_bits_sorted(evars[i], &evars[i..i+128]);
		eprintln!("vars:{} num_bits:{}", i, num_bits);
        let compressed_len = bitpacker.compress_sorted(evars[i], &evars[i..i+128], &mut compressed[..], num_bits);
	    zip.write_all(&compressed[..compressed_len]).expect("OK");

	}

	// https://docs.rs/bitpacking/0.8.4/bitpacking/trait.BitPacker.html


	eprintln!("len: {}", acs.len());
	eprintln!("{:?}", acs);
	zip.finish().expect("error closing zip file");

}

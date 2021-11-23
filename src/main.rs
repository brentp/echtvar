extern crate rust_htslib;
extern crate bitpacking;
use crate::rust_htslib::bcf::{Reader, Read};
use crate::rust_htslib::bcf::record::{Record, Buffer};
use std::borrow::{Borrow, BorrowMut};
use std::io::Write;
use bitpacking::{BitPacker4x, BitPacker};
use echtvar_lib;

//fn get_float_field<'a, B: BorrowMut<Buffer> + Borrow<Buffer> + 'a>(rec: &Record, field: &[u8], buffer: B) -> f32 {
fn get_float_field<'a, B: BorrowMut<Buffer> + 'a>(rec: &Record, field: &[u8], buffer: B) -> f32 {
	return match rec.info_shared_buffer(field, buffer).float().expect("error reading info") {
		Some(v) => v[0],
		None => -1.0,
	};
}
fn get_int_field<'a, B: BorrowMut<Buffer> + Borrow<Buffer> + 'a>(rec: &Record, field: &[u8], buffer: B) -> i32 {
	return match rec.info_shared_buffer(field, buffer).integer().expect("error reading info") {
		Some(v) => v[0],
		None => -1,
	};
}

fn main() {
    let mut stderr = std::io::stderr();

    if std::env::args().len() < 2 {
        println!("expecting arguments: <vcf>")
    }
	//writeln!(stderr, "{}", std::mem::size_of::<Var32>()).expect("error writing to stderr");
    let args: Vec<String> = std::env::args().collect();
    let path = &*args[1];
	writeln!(stderr, "{}", path).expect("error writing to stderr");
    let mut vcf = Reader::from_path(path).ok().expect("Error opening vcf.");
	let header = vcf.header().clone();
	let mut buffer = Buffer::new();

	let mut acs:Vec<u32> = Vec::new();
    let bitpacker = BitPacker4x::new();


	for r in vcf.records() {
		let mut rec = r.expect("error getting record");
		let n: &[u8] = header.rid2name(rec.rid().expect("bad rid")).unwrap();
		let nac: [i32; 1] = [333];

		let af = get_float_field(&rec, b"AF", & mut buffer);
		//rec.push_info_integer(b"AC", &nac).ok();
		let ac = get_int_field(&rec, b"AC", & mut buffer);
        //if ac > 3 { continue; }
		acs.push(ac as u32);


        //if acs.len() * std::mem::size_of::<u32>() == 128 {
        if acs.len() == 128 {
            let num_bits = bitpacker.num_bits(&acs[..]);
            let mut compressed = vec![0u8; 4 * BitPacker4x::BLOCK_LEN];
            let compressed_len = bitpacker.compress(&acs[..], &mut compressed[..], num_bits);
            writeln!(stderr, "uncompressed bits: {}, compressed bits: {}", acs.len() * std::mem::size_of::<u32>(), compressed_len).expect("error writing to stderr");
            acs.clear()

        }

		//writeln!(stderr, "{}:{} AF: {} AC: {}", std::str::from_utf8(n).unwrap(), rec.pos(), af, ac).expect("error writing to stderr");

	}

	// https://docs.rs/bitpacking/0.8.4/bitpacking/trait.BitPacker.html


	writeln!(stderr, "len: {}", acs.len()).expect("error writing to stderr");



}

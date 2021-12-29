use echtvar_lib::{fields, var32, zigzag, kmer16};
use rust_htslib::bcf::record::{Buffer, Record};
use rust_htslib::bcf::{Read as BCFRead, Reader};
use stream_vbyte::{encode::encode, x86::Sse41};

use std::borrow::{Borrow, BorrowMut};
use std::fs::File;
use std::io::Read;
use std::io::Write;

use byteorder::{LittleEndian, WriteBytesExt};
use std::str;

use zip::write::FileOptions;

#[inline]
fn get_int_field<'a, B: BorrowMut<Buffer> + Borrow<Buffer> + 'a>(
    rec: &Record,
    field: &[u8],
    buffer: B,
    default: i32,
) -> i32 {
    return match rec
        .info_shared_buffer(field, buffer)
        .integer()
        .unwrap_or(None) // this becomes default below.
    {
        Some(v) => v[0],
        None => default,
    };
}

#[inline]
fn get_float_field<'a, B: BorrowMut<Buffer> + Borrow<Buffer> + 'a>(
    rec: &Record,
    field: &[u8],
    buffer: B,
    default: f32,
) -> f32 {
    return match rec
        .info_shared_buffer(field, buffer)
        .float()
        .unwrap_or(None) // this becomes default below.
    {
        Some(v) => v[0],
        None => default,
    };
}

fn write_long(
    zipf: &mut zip::ZipWriter<std::io::BufWriter<std::fs::File>>,
    long_vars: &mut Vec<var32::LongVariant>,
    indexes: Vec<usize>,
) {
    let rev_index = argsort(&indexes);
    for l in long_vars.iter_mut() {
        l.idx = rev_index[l.idx as usize] as u32;
    }
    long_vars.sort();
    serde_json::to_writer(zipf, &long_vars).expect("error writing long variables");
}

fn write_bits(
    values: &mut Vec<u32>,
    sorted: bool,
    zipf: &mut zip::ZipWriter<std::io::BufWriter<std::fs::File>>,
    compressed: &mut [u8],
) {
    zipf.write_u32::<LittleEndian>(values.len() as u32).ok();

    if sorted {
        // delta coding
        let mut last = values[0];
        for i in 1..values.len() {
            if values[i] < last {
                panic!("variants out of order at index: {}", i);
            }
            let tmp = values[i];
            values[i] -= last;
            last = tmp;
        }
    }

    let encoded_len = encode::<Sse41>(&values, compressed);
    zipf.write_all(&compressed[..encoded_len]).ok();
}

// https://stackoverflow.com/a/69764256
pub fn argsort<T: Ord>(data: &[T]) -> Vec<usize> {
    let mut indices = (0..data.len()).collect::<Vec<_>>();
    indices.sort_by_key(|&i| &data[i]);
    indices
}

// https://stackoverflow.com/questions/69764803/how-to-sort-a-vector-by-indices-in-rust/69774341#69774341
pub fn sort_by_indices<T: Ord>(data: &mut [T], mut indices: Vec<usize>) {
    for idx in 0..data.len() {
        if indices[idx] != idx {
            let mut current_idx = idx;
            loop {
                let target_idx = indices[current_idx];
                indices[current_idx] = current_idx;
                if indices[target_idx] == target_idx {
                    break;
                }
                data.swap(current_idx, target_idx);
                current_idx = target_idx;
            }
        }
    }
}

fn is_sorted<T: std::cmp::PartialOrd>(data: &Vec<T>) -> bool {
    for i in 1..data.len() {
        if data[i] < data[i - 1] {
            return false;
        }
    }
    return true;
}

pub fn encoder_main(vpath: &str, opath: &str, jpath: &str) {
    let zpath = std::path::Path::new(opath);
    let jpath = std::path::Path::new(jpath);

    let mut json = String::new();
    File::open(jpath)
        .expect("error opening json file")
        .read_to_string(&mut json)
        .expect("error parsing json file");
    let fields: Vec<fields::Field> =
        json5::from_str(&json).expect("error reading json into fields");

    let mut vcf = Reader::from_path(vpath).ok().expect("Error opening vcf.");
    vcf.set_threads(2).ok();
    let header = vcf.header().clone();
    let mut buffer = Buffer::new();

    let zfile = std::fs::File::create(&zpath).unwrap();
    let fbuffer = std::io::BufWriter::with_capacity(65536, zfile);
    let mut zipf = zip::ZipWriter::new(fbuffer);

    zipf.add_directory("echtvar/", FileOptions::default().large_file(true))
        .expect("error writing zip");

    let options = FileOptions::default()
        .compression_method(zip::CompressionMethod::Deflated)
        .unix_permissions(0o755);

    zipf.start_file("echtvar/config.json", options)
        .expect("error starting json file");
    zipf.write_all(
        serde_json::to_string_pretty(&fields)
            .expect("error serializing json config")
            .as_bytes(),
    )
    .ok();

    let mut compressed = vec![0u8; 50_000_000]; // TODO: set this based on values.len
    let mut last_rid: i32 = -1;
    let mut last_mod: i64 = 0;

    let mut long_vars: Vec<var32::LongVariant> = Vec::new();
    let mut var32s: Vec<u32> = Vec::new();

    let mut values_vv: Vec<Vec<u32>> = fields.iter().map(|_| Vec::new()).collect();

    for r in vcf.records() {
        let rec = r.expect("error getting record");
        // if we hit a new chrom or a new chunk we write the last chunk and start a new one.
        if rec.rid().expect("no rid found") as i32 != last_rid || rec.pos() >> 20 != last_mod {
            if last_rid != -1 {
                if values_vv[0].len() != 0 {
                    let n: &[u8] = header.rid2name(last_rid as u32).unwrap();
                    let chrom = str::from_utf8(n).unwrap();

                    // we just assume it's unsorted and apply the permutation
                    let indexes = argsort(&var32s);

                    for (i, values) in values_vv.iter_mut().enumerate() {
                        let fname =
                            format!("echtvar/{}/{}/{}.bin", chrom, last_mod, fields[i].alias);
                        zipf.start_file(fname, options)
                            .expect("error starting file");
                        sort_by_indices(values, indexes.clone());
                        write_bits(values, false, &mut zipf, &mut compressed);
                        values.clear();
                    }

                    let fname = format!("echtvar/{}/{}/var32.bin", chrom, last_mod);
                    zipf.start_file(fname, options)
                        .expect("error starting file");
                    sort_by_indices(&mut var32s, indexes.clone());
                    if !is_sorted(&var32s) {
                        eprintln!("BAD\nBAD\nBAD\nBAD");
                    }
                    write_bits(&mut var32s, true, &mut zipf, &mut compressed);
                    var32s.clear();

                    let fname = format!("echtvar/{}/{}/too-long-for-var32.txt", chrom, last_mod);
                    zipf.start_file(fname, options)
                        .expect("error starting file");
                    write_long(&mut zipf, &mut long_vars, indexes);
                    long_vars.clear()
                }
            }
            last_rid = rec.rid().unwrap() as i32;
            last_mod = rec.pos() >> 20;
        }

        for (i, fld) in fields.iter().enumerate() {
            // NOTE that internally, we always use u32::MAX as missing and then replace this with `missing_value` when annotating.
            let v = match fld.ftype {
                fields::FieldType::Integer => {
                    let val =
                        get_int_field(&rec, fld.field.as_bytes(), &mut buffer, fld.missing_value);
                    if val == fld.missing_value {
                        u32::MAX
                    } else {
                        if fld.zigzag {
                            zigzag::encode(val)
                        } else {
                            val as u32
                        }
                    }
                }
                fields::FieldType::Float => {
                    let val = get_float_field(
                        &rec,
                        fld.field.as_bytes(),
                        &mut buffer,
                        fld.missing_value as f32,
                    );
                    if (val - fld.missing_value as f32).abs() <= f32::EPSILON {
                        u32::MAX
                    } else {
                        let v = (val * fld.multiplier as f32) as i32;
                        if fld.zigzag {
                            zigzag::encode(v)
                        } else {
                            v as u32
                        }
                    }
                }
                fields::FieldType::Categorical => panic!("not implemented"),
            };

            values_vv[i].push(v);
        }

        let alleles = rec.alleles();
        var32s.push(var32::encode(rec.pos() as u32, alleles[0], alleles[1]));

        if alleles[0].len() + alleles[1].len() > var32::MAX_COMBINED_LEN {
            long_vars.push(var32::LongVariant {
                position: rec.pos() as u32,
                reference: kmer16::encode(alleles[0]),
                alternate: kmer16::encode(alleles[1]),
                idx: (var32s.len() - 1) as u32,
            });
        }
    }
    if values_vv[0].len() != 0 {
        let indexes = argsort(&var32s);
        let n: &[u8] = header.rid2name(last_rid as u32).unwrap();
        let chrom = str::from_utf8(n).unwrap();

        for (i, values) in values_vv.iter_mut().enumerate() {
            let fname = format!("echtvar/{}/{}/{}.bin", chrom, last_mod, fields[i].alias);
            zipf.start_file(fname, options)
                .expect("error starting file");
            sort_by_indices(values, indexes.clone());
            write_bits(values, false, &mut zipf, &mut compressed);
            values.clear();
        }

        let fname = format!("echtvar/{}/{}/var32.bin", chrom, last_mod);
        zipf.start_file(fname, options)
            .expect("error starting file");
        sort_by_indices(&mut var32s, indexes.clone());
        write_bits(&mut var32s, true, &mut zipf, &mut compressed);

        let fname = format!("echtvar/{}/{}/too-long-for-var32.txt", chrom, last_mod);
        zipf.start_file(fname, options)
            .expect("error starting file");
        write_long(&mut zipf, &mut long_vars, indexes);
    }
    zipf.finish().expect("error closing zip file");
}

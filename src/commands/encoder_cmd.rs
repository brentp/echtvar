use bincode::Options;
use echtvar_lib::{echtvar::bstrip_chr, fields, kmer16, var32, zigzag};
use rust_htslib::bcf::header::{TagLength, TagType, HeaderRecord};
use rust_htslib::bcf::record::{Buffer, Record};
use rust_htslib::bcf::{Read as BCFRead, Reader};
use stream_vbyte::{encode::encode, x86::Sse41};

use std::borrow::{Borrow, BorrowMut};
use std::collections::HashMap;
use std::fs::File;
use std::io::{Read, Write};

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

fn get_string_field<'a, B: BorrowMut<Buffer> + Borrow<Buffer> + 'a>(
    rec: &Record,
    field: &[u8],
    buffer: B,
    default: &String,
    lookup: &mut HashMap<String, u32>,
) -> u32 {
    let s = if field == b"FILTER" {
        let hdr = rec.header();
        let f = rec
            .filters()
            .map(|f| unsafe { String::from_utf8_unchecked(hdr.id_to_name(f)) })
            .collect::<Vec<String>>()
            .join(";");
        if f.len() == 0 {
            default.clone()
        } else {
            f
        }
    } else {
        match rec
            .info_shared_buffer(field, buffer)
            .string()
            .unwrap_or(None)
        {
            Some(v) => unsafe { String::from_utf8_unchecked(v[0].to_vec()) },
            None => default.to_string(),
        }
    };
    // lookup from string -> idx so we can get the reverse when decoding.
    let l = lookup.len() as u32;
    *lookup.entry(s).or_insert(l)
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

    let bc = bincode::DefaultOptions::new()
        .serialize(long_vars)
        .expect("error serializing long vars");
    zipf.write_all(&bc).expect("error writing long variants");
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
    assert!(data.len() == indices.len());
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

fn hdr_info_id2description(
    mut hrecs: Vec<HeaderRecord>,
    id: &String,
    default: &std::string::String,
) -> std::string::String {
    hrecs.retain(|rec| match rec {
        HeaderRecord::Info {key: _, values: v} => &v["ID"] == id,
        _ => false}
    );
    if hrecs.len() != 1 {
        panic!("Field {} is either not present in the header or present multiple times!", id);
    };
    let description = match hrecs.first().unwrap() {
        HeaderRecord::Info {key: _, values: v} => if v.contains_key("Description") { &v["Description"] } else { default },
        _ => default,
    };
    return description.to_string();
}

pub fn encoder_main(vpaths: Vec<&str>, opath: &str, jpath: &str) {
    let zpath = std::path::Path::new(opath);
    let jpath = std::path::Path::new(jpath);

    let mut json = String::new();
    File::open(jpath)
        .expect("error opening json file")
        .read_to_string(&mut json)
        .expect("error parsing json file");
    let mut fields: Vec<fields::Field> =
        json5::from_str(&json).expect("error reading json into fields");

    let mut vcf = if !(*vpaths[0]).eq("/dev/stdin") && !(*vpaths[0]).eq("-") {
        Reader::from_path(vpaths[0])
            .ok()
            .expect("Error opening vcf.")
    } else {
        Reader::from_stdin().ok().expect("Error opening stdin vcf.")
    };

    vcf.set_threads(2).ok();
    let header = vcf.header().clone();
    let mut buffer = Buffer::new();
    // hashmap of hashmaps e.g. {'alias': {'A': 0, 'B': 1}, 'othercol': {'PASS': 0, 'FAIL': 1}}
    let mut lookups = HashMap::new();

    for f in fields.iter_mut() {
        let (tt, tl) = if f.field == "FILTER" {
            (TagType::String, TagLength::Variable)
        } else {
            header
                .info_type(&(f.field.as_bytes()))
                .expect(&format!("unable to find header type for {}", f.field).to_string())
        };
        match tt {
            TagType::Integer => f.ftype = fields::FieldType::Integer,
            TagType::Float => {
                f.ftype = fields::FieldType::Float;
                if f.multiplier == 1 {
                    eprintln!(
                        "[echtvar] warning! using a multiplier of 1 for float field {}.",
                        f.field
                    );
                    eprintln!("\tIf this field contains values less than 1, use a multiplier so the values can be stored as integers.");
                    eprintln!("\tLarger multipliers result in higher precision.");
                }
            },
            TagType::String /* TagType::Flag */ => {
              f.ftype = fields::FieldType::Categorical;
              // and a new table into lookups for this field
              lookups.entry(f.alias.clone()).or_insert(HashMap::new());
            },
            _ => panic!(
                "[echtvar] unsupported field type: {:?} for field {}",
                tt, f.field
            ),
        };
        match tl {
            TagLength::Fixed(value) => f.number = value.to_string(),
            TagLength::AltAlleles => f.number = "A".to_string(),
            TagLength::Alleles => f.number = "R".to_string(),
            TagLength::Genotypes => f.number = "G".to_string(),
            TagLength::Variable => f.number = ".".to_string(),
            _ => panic!(
                "[echtvar] unsupported field length: {:?} for field {}",
                tl, f.field
            ),
        };
        println!("Old description for field {}: {}", f.field, f.description);
        f.description = hdr_info_id2description(header.header_records(), &f.field, &f.description);
        println!("New description for field {}: {}", f.field, f.description);
    }

    let zfile = std::fs::File::create(&zpath).unwrap();
    let fbuffer = std::io::BufWriter::with_capacity(65536, zfile);
    let mut zipf = zip::ZipWriter::new(fbuffer);

    zipf.add_directory("echtvar/", FileOptions::default().large_file(true))
        .expect("error writing zip");

    let options = FileOptions::default()
        .compression_method(zip::CompressionMethod::Deflated)
        .large_file(true)
        .unix_permissions(0o755);

    zipf.start_file("echtvar/config.json", options)
        .expect("error starting json file");
    zipf.write_all(
        serde_json::to_string_pretty(&fields)
            .expect("error serializing json config")
            .as_bytes(),
    )
    .ok();

    let mut compressed = vec![0u8; 50_000_000];
    let mut last_rid: i32 = -1;
    let mut last_mod: i64 = 0;

    let mut long_vars: Vec<var32::LongVariant> = Vec::new();
    let mut var32s: Vec<u32> = Vec::new();
    let mut n_long_vars = 0;
    let mut n_vars = 0u64;

    let mut values_vv: Vec<Vec<u32>> = fields.iter().map(|_| Vec::new()).collect();

    for (i, vpath) in vpaths.iter().enumerate() {
        if i > 0 {
            vcf = if !(*vpath).eq("/dev/stdin") && !(*vpath).eq("-") {
                Reader::from_path(vpath).ok().expect("Error opening vcf.")
            } else {
                match Reader::from_stdin() {
                    Ok(file) => file,
                    Err(error) => panic!("problem opening from stdin: {:?}", error),
                }
            };
            vcf.set_threads(2).ok();
        }
        eprintln!("[echtvar] adding VCF:{}", vpath);
        let mut warn = 0;

        for r in vcf.records() {
            let rec = r.expect("error getting record");
            n_vars += 1;
            // if we hit a new chrom or a new chunk we write the last chunk and start a new one.
            if rec.rid().expect("no rid found") as i32 != last_rid || rec.pos() >> 20 != last_mod {
                if last_rid != -1 {
                    if values_vv[0].len() != 0 {
                        let n: &[u8] = header.rid2name(last_rid as u32).unwrap();
                        let chrom = bstrip_chr(str::from_utf8(n).unwrap());

                        // we just assume it's unsorted and apply the permutation
                        let indexes = argsort(&var32s);

                        for (i, values) in values_vv.iter_mut().enumerate() {
                            let fname =
                                format!("echtvar/{}/{}/{}.bin", chrom, last_mod, fields[i].alias);
                            zipf.start_file(fname, options)
                                .expect("error starting file");
                            sort_by_indices(values, indexes.clone());
                            if compressed.len() < 5 * values.len() {
                                eprintln!(
                                    "[echtvar] resizing compressed array to: {}",
                                    5 * values.len()
                                );
                                compressed.resize(5 * values.len(), 0x0);
                            }
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

                        let fname =
                            format!("echtvar/{}/{}/too-long-for-var32.enc", chrom, last_mod);
                        zipf.start_file(fname, options)
                            .expect("error starting file");
                        write_long(&mut zipf, &mut long_vars, indexes);
                        n_long_vars += long_vars.len();
                        long_vars.clear();
                    }
                }
                if last_rid != rec.rid().unwrap() as i32 {
                    last_rid = rec.rid().unwrap() as i32;
                    let n: &[u8] = header.rid2name(last_rid as u32).unwrap();
                    let chrom = bstrip_chr(str::from_utf8(n).unwrap());
                    eprintln!("[echtvar] on chromosome {:?}", chrom);
                }
                last_mod = rec.pos() >> 20;
            }

            for (i, fld) in fields.iter().enumerate() {
                // NOTE that internally, we always use u32::MAX as missing and then replace this with `missing_value` when annotating.
                let v = match fld.ftype {
                    fields::FieldType::Integer => {
                        let val = get_int_field(
                            &rec,
                            fld.field.as_bytes(),
                            &mut buffer,
                            fld.missing_value,
                        );
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
                    fields::FieldType::Categorical => {
                        let val = get_string_field(
                            &rec,
                            fld.field.as_bytes(),
                            &mut buffer,
                            &fld.missing_string,
                            lookups.get_mut(&fld.alias).unwrap(),
                        );
                        val
                    }
                };

                values_vv[i].push(v);
            }

            let mut alleles = rec.alleles();
            if alleles.len() == 1 {
                alleles.push(&alleles[0]);
            } else if alleles.len() != 2 {
                last_rid = rec.rid().unwrap() as i32;
                let n: &[u8] = header.rid2name(last_rid as u32).unwrap();
                let chrom = bstrip_chr(str::from_utf8(n).unwrap());
                panic!(
                    "[echtvar] variants must be decomposed before running {}:{} {:?}. see: https://github.com/brentp/echtvar/wiki/decompose",
                    chrom,
                    rec.pos(),
                    rec.alleles()
                        .iter()
                        .map(|a| String::from_utf8_lossy(a))
                        .collect::<Vec<_>>()
                );
            }
            var32s.push(var32::encode(
                rec.pos() as u32,
                alleles[0],
                alleles[1],
                &mut warn,
            ));

            if alleles[0].len() + alleles[1].len() > var32::MAX_COMBINED_LEN {
                long_vars.push(var32::LongVariant {
                    position: rec.pos() as u32,
                    sequence: kmer16::encode_var(alleles[0], alleles[1]),
                    idx: (var32s.len() - 1) as u32,
                });
            }
        }
        if values_vv[0].len() != 0 {
            let indexes = argsort(&var32s);
            let n: &[u8] = header.rid2name(last_rid as u32).unwrap();
            let chrom = bstrip_chr(str::from_utf8(n).unwrap());

            for (i, values) in values_vv.iter_mut().enumerate() {
                let fname = format!("echtvar/{}/{}/{}.bin", chrom, last_mod, fields[i].alias);
                zipf.start_file(fname, options)
                    .expect("error starting file");
                sort_by_indices(values, indexes.clone());
                if compressed.len() < 5 * values.len() {
                    eprintln!(
                        "[echtvar] resizing compressed array to: {}",
                        5 * values.len()
                    );
                    compressed.resize(5 * values.len(), 0x0);
                }
                write_bits(values, false, &mut zipf, &mut compressed);
                values.clear();
            }

            let fname = format!("echtvar/{}/{}/var32.bin", chrom, last_mod);
            zipf.start_file(fname, options)
                .expect("error starting file");
            sort_by_indices(&mut var32s, indexes.clone());
            write_bits(&mut var32s, true, &mut zipf, &mut compressed);

            let fname = format!("echtvar/{}/{}/too-long-for-var32.enc", chrom, last_mod);
            zipf.start_file(fname, options)
                .expect("error starting file");
            n_long_vars += long_vars.len();
            write_long(&mut zipf, &mut long_vars, indexes);
            long_vars.clear();
            var32s.clear();
        }
    }
    for (afld, lookup) in lookups.iter() {
        let fname = format!("echtvar/strings/{}.txt", afld);
        zipf.start_file(fname, options)
            .expect("error starting file");
        eprintln!(
            "[echtvar] found {} unique values for {})",
            lookup.len(),
            afld
        );
        // use array to make sure we save in sorted order.
        let mut arr = vec![""; lookup.len()];
        for (name, idx) in lookup.iter() {
            arr[*idx as usize] = name;
        }
        // write the string entries it order to the file.
        // when decoding, we can index into this array to get the string value from
        // the encoded integer.
        for v in arr.iter() {
            zipf.write(v.as_bytes()).expect("error writing to zip file");
            zipf.write(b"\n").expect("error writing to zip file");
        }
    }

    zipf.finish().expect("error closing zip file");
    let pct = 100.0 * (n_long_vars as f64) / (n_vars as f64);
    eprintln!(
        "[echtvar] wrote {n_vars} total variants and {n_long_vars} long variants ({pct:.2}%)"
    );
}

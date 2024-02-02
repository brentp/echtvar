use std::error::Error;
use std::time;

use rust_htslib::bcf::{header::Header, Format, Writer};
use rust_htslib::bcf::{Read as BCFRead, Reader};

use echtvar_lib::echtvar::{EchtVars, Value};
use echtvar_lib::fields;

use fasteval::Compiler;
use fasteval::Evaler;

pub fn annotate_main(
    vpath: &str,
    opath: &str,
    include_expr: Option<&str>,
    epaths: Vec<&str>,
) -> Result<(), Box<dyn Error>> {
    let mut ipath = vpath;
    if ipath == "-" || ipath == "stdin" {
        ipath = "/dev/stdin";
    }
    let mut vcf = Reader::from_path(ipath).expect("Error opening vcf.");
    vcf.set_threads(2).expect("error setting threads");
    let header_view = vcf.header();
    let mut header = Header::from_template(header_view);

    // each given echtvar file updates the header.
    let mut echts: Vec<EchtVars> = epaths
        .iter()
        .map(|p| {
            let mut e = EchtVars::open(p);
            e.update_header(&mut header, p);
            e
        })
        .collect();
    EchtVars::add_cmd_header(&mut header, ipath, opath, &include_expr, epaths);

    let parser = fasteval::Parser::new();
    let mut slab = fasteval::Slab::new();
    let mut ns = fasteval::EmptyNamespace;
    let mut expr_values = Vec::with_capacity(echts.len());

    for (i, e) in echts.iter().enumerate() {
        // a vector within expr_values for each echtvar file.
        expr_values.push(Vec::with_capacity(e.fields.len()));
        // handle the expression stuff.
        for (j, fld) in e.fields.iter().enumerate() {
            expr_values[i].push(0.0_f64);
            unsafe {
                slab.ps
                    .add_unsafe_var(fld.alias.clone(), &expr_values[i][j])
            }
            // this sets, e.g missense = 1
            // so user can do: impact == missense or any numerical operation.
            if !e.strings[j].is_empty() {
                if e.strings[j].len() > 256 {
                    eprintln!("[echtvar] not exposing field '{}' for expressions as it has cardinality > 256", fld.alias)
                } else {
                    for (k, s) in e.strings[j].iter().enumerate() {
                        // TODO: use another method, don't need unsafe_var here as these will be
                        // constants.
                        unsafe { slab.ps.add_unsafe_var(s.clone(), &(k as f64)) }
                    }
                }
            }
        }
    }
    let mut is_compiled = false;
    let compiled = if let Some(uexpr) = include_expr {
        is_compiled = true;
        parser
            .parse(uexpr, &mut slab.ps)?
            .from(&slab.ps)
            .compile(&slab.ps, &mut slab.cs)
    } else {
        // just compile an empty expression that never is evaluated.
        parser
            .parse("true", &mut slab.ps)?
            .from(&slab.ps)
            .compile(&slab.ps, &mut slab.cs)
    };

    // TODO: handle stdout
    let mut ovcf = Writer::from_path(
        opath,
        &header,
        false,
        if opath.ends_with("bcf") {
            Format::Bcf
        } else {
            Format::Vcf
        },
    )
    .expect("error opening bcf for output");
    ovcf.set_threads(2).expect("error setting threads");
    let oheader_view = ovcf.header().clone();

    let start = time::Instant::now();
    let mut n = 0u64;
    let mut n_written = 0u64;
    let mut modu = 10000u64;

    for r in vcf.records() {
        let mut record = r.expect("error reading record");
        ovcf.translate(&mut record);

        if n > 0 && n % modu == 0 {
            let rid = record.rid().unwrap();
            let chrom = std::str::from_utf8(oheader_view.rid2name(rid).unwrap()).unwrap();
            let mili = time::Instant::now().duration_since(start).as_millis();
            if n >= 3 * modu && modu < 1000000 {
                modu *= 10;
            }
            if n >= 100 * modu && modu < 300000000 {
                modu *= 10;
            }
            eprintln!(
                "[echtvar] {}:{} evaluated {} variants ({} / second). wrote {} variants.",
                chrom,
                record.pos(),
                n,
                1000 * (n as u128) / mili,
                n_written,
            );
        }
        n += 1;
        // First check if the variant is *, skip those
        if String::from_utf8_lossy(record.alleles()[1]) == "*" {
            let rid = record.rid().unwrap();
            let chrom = std::str::from_utf8(oheader_view.rid2name(rid).unwrap()).unwrap();
            eprintln!(
                "contig {} pos {} alt has * value, skipping annotation, outputting entry as-is",
                &chrom,
                record.pos() + 1
            );
            ovcf.write(&record).expect("failed to write record");
            continue;
        }
        // this updates evalues and fills expr values

        for (i, e) in echts.iter_mut().enumerate() {
            e.update_expr_values(&mut record, &mut expr_values[i]);
        }

        if is_compiled && fasteval::eval_compiled!(compiled, &slab, &mut ns) == 0.0 {
            continue;
        }
        n_written += 1;

        for e in echts.iter() {
            for fld in &e.fields {
                let v = e.evalues[fld.values_i];

                match v {
                    Value::Int(i) => match fld.ftype {
                        fields::FieldType::Categorical => {
                            // categorical missing_value must be set to the index of the missing_string
                            assert!(i >= 0, "can't have missing value for categorical!");
                            let val = [e.strings[fld.values_i][i as usize].as_bytes()];
                            record
                                .push_info_string(fld.alias.as_bytes(), &val)
                                .unwrap_or_else(|_| {
                                    panic!("{}", format!("error adding string for {}", fld.alias))
                                });
                        }
                        _ => {
                            let val = [i];
                            record
                                .push_info_integer(fld.alias.as_bytes(), &val)
                                .unwrap_or_else(|_| {
                                    panic!("{}", format!("error adding integer {}", fld.alias))
                                });
                        }
                    },
                    Value::Float(f) => {
                        let val = [f];
                        record
                            .push_info_float(fld.alias.as_bytes(), &val)
                            .unwrap_or_else(|_| {
                                panic!("{}", format!("error adding float {}", fld.alias))
                            });
                    }
                }
            }
        }
        ovcf.write(&record).expect("failed to write record");
    }
    let mili = std::cmp::max(1, time::Instant::now().duration_since(start).as_millis());
    eprintln!(
        "[echtvar] evaluated {} variants ({} / second). wrote {} variants.",
        n,
        1000 * (n as u128) / mili,
        n_written,
    );

    /*
    //let ep = std::path::Path::new(&*epaths[0]);
    let file = fs::File::open(ep).expect("error accessing zip file");
    let mut archive = zip::ZipArchive::new(&file).expect("error opening zip file");


    // encoded 46881 u32s into 60515 bytes
    let mut iz = archive.by_name("echtvar/chr21/4/gnomad_AN.bin").expect("unable to open file");

    let n = iz.read_u32::<LittleEndian>().ok().expect("error reading number of values from zip file") as usize;
    let mut comr = vec![0 as u8; iz.size() as usize - 4];
    iz.read_exact(&mut comr)?;

    eprintln!("number of values: {} bytes:{}", n, iz.size());

    let mut nums: Vec<u32> = Vec::new();
    nums.resize(n, 0);

    let n_d = decode::<Ssse3>(&comr, n, &mut nums);
    eprintln!("{} {} {:?}", n_d, iz.size(), &nums[..1000]);
    */

    Ok(())
}

use std::error::Error;
use std::time;

use rust_htslib::bcf::{header::Header, Format, Writer, record::Record};
use rust_htslib::bcf::{Read as BCFRead, Reader};

use echtvar_lib::echtvar::EchtVars;
use echtvar_lib::echtvar::Value;

use fasteval::Compiler;
use fasteval::Evaler;

use std::thread;
use std::sync::mpsc::sync_channel;

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
    let mut vcf = Reader::from_path(ipath).ok().expect("Error opening vcf.");
    vcf.set_threads(2).expect("error setting threads");
    let header_view = vcf.header();
    let mut header = Header::from_template(&header_view);

    let mut e = EchtVars::open(&*epaths[0]);
    e.update_header(&mut header);

    let parser = fasteval::Parser::new();
    let mut slab = fasteval::Slab::new();
    let mut ns = fasteval::EmptyNamespace;
    let mut expr_values = vec![];

    for (i, fld) in e.fields.iter().enumerate() {
        expr_values.push(0.0 as f64);
        unsafe { slab.ps.add_unsafe_var(fld.alias.clone(), &expr_values[i]) }
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

    let mut fake_ovcf = Writer::from_path(
        "/dev/null",
        &header,
        false,
        if opath.ends_with("bcf") {
            Format::Bcf
        } else {
            Format::Vcf
        },
    )
    .ok()
    .expect("error opening bcf for output");

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
    .ok()
    .expect("error opening bcf for output");
    ovcf.set_threads(2).expect("error setting threads");


    let oheader_view = ovcf.header().clone();

    let start = time::Instant::now();
    let mut n = 0;
    let mut n_written = 0;
    let mut modu = 10000;

    
    let (tx, rx) = sync_channel::<Record>(100);
    
    thread::spawn(move|| {
       while let Ok(record) = rx.recv() {
          ovcf.write(&record).expect("failed to write record");
          drop(record)
       }
    });

    for r in vcf.records() {
        let mut record = r.expect("error reading record");
        fake_ovcf.translate(&mut record);

        if n > 0 && n % modu == 0 {
            let rid = record.rid().unwrap();
            let chrom = std::str::from_utf8(oheader_view.rid2name(rid).unwrap()).unwrap();
            let mili = time::Instant::now().duration_since(start).as_millis();
            if n >= 3 * modu && modu < 1000000 {
                modu *= 10;
            }

            eprintln!(
                "[echtvar] {}:{} evaluated {} variants ({} / second). wrote {} variants.",
                chrom,
                record.pos(),
                n,
                1000 * n / mili,
                n_written,
            );
        }
        n += 1;
        // this updates evalues and fills expr values
        e.update_expr_values(&mut record, &mut expr_values);

        if is_compiled {
            if fasteval::eval_compiled!(compiled, &slab, &mut ns) == 0.0 {
                continue;
            }
        }
        n_written += 1;

        for fld in &e.fields {
            let v = e.evalues[fld.values_i];

            match v {
                Value::Int(i) => {
                    let val = [i];
                    record
                        .push_info_integer(fld.alias.as_bytes(), &val)
                        .expect(&format!("error adding integer {}", fld.alias).to_string());
                }
                Value::Float(f) => {
                    let val = [f];
                    record
                        .push_info_float(fld.alias.as_bytes(), &val)
                        .expect(&format!("error adding float {}", fld.alias).to_string());
                }
            }
        }
        tx.send(record).expect("error sending record");
        //ovcf.write(&record).expect("failed to write record");
    }
    drop(tx);
    let mili = time::Instant::now().duration_since(start).as_millis();
    eprintln!(
        "[echtvar] evaluated {} variants ({} / second). wrote {} variants.",
        n,
        1000 * n / mili,
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

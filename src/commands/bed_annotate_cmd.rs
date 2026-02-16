#![cfg(feature = "bed")]

use std::collections::HashMap;
use std::error::Error;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::time;

use echtvar_lib::echtvar::{EchtVars, Value, Variant};
use echtvar_lib::fields;

use fasteval::Compiler;
use fasteval::Evaler;

struct BedVariant {
    chrom: String,
    chrom_id: i32,
    pos: u32,
    ref_allele: Vec<u8>,
    alt_allele: Vec<u8>,
}

impl Variant for BedVariant {
    fn chrom(&self) -> String {
        self.chrom.clone()
    }

    fn rid(&self) -> i32 {
        self.chrom_id
    }
    
    fn position(&self) -> u32 {
        self.pos
    }
    
    fn alleles(&self) -> Vec<&[u8]> {
        vec![&self.ref_allele, &self.alt_allele]
    }
}

fn format_value(e: &EchtVars, fld: &fields::Field) -> String {
    let v = e.evalues[fld.values_i];
    match v {
        Value::Int(i) => match fld.ftype {
            fields::FieldType::Flag => {
                if i != 0 {
                    "true".to_string()
                } else {
                    "false".to_string()
                }
            }
            fields::FieldType::Categorical => e.strings[fld.values_i][i as usize].clone(),
            _ => i.to_string(),
        },
        Value::Float(f) => {
            if f.is_nan() || f.is_infinite() {
                ".".to_string()
            } else {
                f.to_string()
            }
        }
    }
}

pub fn bed_annotate_main(
    input_path: &str,
    output_path: &str,
    include_expr: Option<&str>,
    epaths: Vec<&str>,
    ref_col: Option<usize>,
    alt_col: Option<usize>,
) -> Result<(), Box<dyn Error>> {
    let allele_specific = ref_col.is_some() && alt_col.is_some();
    let ref_col_idx = ref_col.unwrap_or(0);
    let alt_col_idx = alt_col.unwrap_or(0);

    let mut echts: Vec<EchtVars> = epaths.iter().map(|p| EchtVars::open(p)).collect();

    // Set up fasteval expression
    let parser = fasteval::Parser::new();
    let mut slab = fasteval::Slab::new();
    let mut ns = fasteval::EmptyNamespace;
    let mut expr_values = Vec::with_capacity(echts.len());

    for (i, e) in echts.iter().enumerate() {
        expr_values.push(Vec::with_capacity(e.fields.len()));
        for (j, fld) in e.fields.iter().enumerate() {
            expr_values[i].push(0.0_f64);
            unsafe {
                slab.ps
                    .add_unsafe_var(fld.alias.clone(), &expr_values[i][j])
            }
            if !e.strings[j].is_empty() {
                if e.strings[j].len() > 256 {
                    eprintln!(
                        "[echtvar] not exposing field '{}' for expressions as it has cardinality > 256",
                        fld.alias
                    );
                } else {
                    for (k, s) in e.strings[j].iter().enumerate() {
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
        parser
            .parse("true", &mut slab.ps)?
            .from(&slab.ps)
            .compile(&slab.ps, &mut slab.cs)
    };

    // Open input BED file
    let ipath = match input_path {
        "-" | "stdin" => "/dev/stdin",
        p => p,
    };
    let reader: Box<dyn BufRead> = if ipath == "/dev/stdin" {
        Box::new(BufReader::new(std::io::stdin()))
    } else if ipath.ends_with(".gz") {
        Box::new(BufReader::new(flate2::read::GzDecoder::new(
            std::fs::File::open(ipath)?,
        )))
    } else {
        Box::new(BufReader::new(std::fs::File::open(ipath)?))
    };

    // Open output file
    let writer: Box<dyn Write> = if output_path == "-" || output_path == "/dev/stdout" {
        Box::new(BufWriter::new(std::io::stdout()))
    } else if output_path.ends_with(".gz") {
        Box::new(BufWriter::new(flate2::write::GzEncoder::new(
            std::fs::File::create(output_path)?,
            flate2::Compression::default(),
        )))
    } else {
        Box::new(BufWriter::new(std::fs::File::create(output_path)?))
    };
    let mut writer = writer;

    // Chromosome ID tracking
    let mut chrom_ids: HashMap<String, i32> = HashMap::new();
    let mut next_chrom_id: i32 = -2;

    let mut header_written = false;
    let mut input_header_cols: Vec<String> = Vec::new();

    // Collect annotation field names
    let anno_fields: Vec<String> = echts
        .iter()
        .flat_map(|e| e.fields.iter().map(|f| f.alias.clone()))
        .collect();

    let start = time::Instant::now();
    let mut n = 0u64;
    let mut n_written = 0u64;

    for line_result in reader.lines() {
        let line = line_result?;
        if line.is_empty() {
            continue;
        }

        // Handle comment/header lines
        if line.starts_with('#') {
            if !header_written {
                let cols: Vec<&str> = line[1..].split('\t').collect();
                input_header_cols = cols.iter().map(|s| s.to_string()).collect();
            }
            continue;
        }

        let fields_vec: Vec<&str> = line.split('\t').collect();
        if fields_vec.len() < 3 {
            eprintln!(
                "[echtvar] skipping line with fewer than 3 columns: {}",
                line
            );
            continue;
        }

        // Write header on first data line
        if !header_written {
            write!(writer, "#chrom\tstart\tend")?;
            if input_header_cols.len() > 3 {
                for col in &input_header_cols[3..] {
                    write!(writer, "\t{}", col)?;
                }
            } else if fields_vec.len() > 3 {
                for i in 3..fields_vec.len() {
                    write!(writer, "\tcol{}", i + 1)?;
                }
            }
            if !allele_specific {
                write!(writer, "\tREF\tALT")?;
            }
            for name in &anno_fields {
                write!(writer, "\t{}", name)?;
            }
            writeln!(writer)?;
            header_written = true;
        }

        let chrom = fields_vec[0];
        let start_pos: u32 = fields_vec[1]
            .parse()
            .expect("invalid start position in BED");

        let extra_cols: Vec<&str> = if fields_vec.len() > 3 {
            fields_vec[3..].to_vec()
        } else {
            Vec::new()
        };

        // Get or assign chromosome ID
        let chrom_id = *chrom_ids
            .entry(chrom.to_string())
            .or_insert_with(|| {
                let id = next_chrom_id;
                next_chrom_id -= 1;
                id
            });

        n += 1;

        if allele_specific {
            // Mode A: Allele-specific annotation
            let ri = ref_col_idx - 1; // convert to 0-based
            let ai = alt_col_idx - 1;
            if ri >= fields_vec.len() || ai >= fields_vec.len() {
                eprintln!(
                    "[echtvar] --ref-col or --alt-col out of range for line {}",
                    n
                );
                continue;
            }
            let ref_allele = fields_vec[ri].as_bytes().to_vec();
            let alt_allele = fields_vec[ai].as_bytes().to_vec();

            let mut variant = BedVariant {
                chrom: chrom.to_string(),
                chrom_id,
                pos: start_pos,
                ref_allele,
                alt_allele,
            };

            for (i, e) in echts.iter_mut().enumerate() {
                e.update_expr_values(&mut variant, &mut expr_values[i]);
            }

            if is_compiled && fasteval::eval_compiled!(compiled, &slab, &mut ns) == 0.0 {
                continue;
            }

            write!(
                writer,
                "{}\t{}\t{}",
                fields_vec[0], fields_vec[1], fields_vec[2]
            )?;
            for col in &extra_cols {
                write!(writer, "\t{}", col)?;
            }
            for e in echts.iter() {
                for fld in &e.fields {
                    write!(writer, "\t{}", format_value(e, fld))?;
                }
            }
            writeln!(writer)?;
            n_written += 1;
        } else {
            // Mode B: Position-scan annotation
            for e in echts.iter_mut() {
                e.set_position_by_name(chrom, start_pos)?;
            }

            // Collect all unique variants at this position across all echtvar files
            let mut all_variants: Vec<(Vec<u8>, Vec<u8>)> = Vec::new();
            for e in echts.iter() {
                let variants = e.variants_at_position(start_pos);
                for v in variants {
                    if !all_variants.contains(&v) {
                        all_variants.push(v);
                    }
                }
            }

            if all_variants.is_empty() {
                // No variants found: output with missing values
                write!(
                    writer,
                    "{}\t{}\t{}",
                    fields_vec[0], fields_vec[1], fields_vec[2]
                )?;
                for col in &extra_cols {
                    write!(writer, "\t{}", col)?;
                }
                write!(writer, "\t.\t.")?;
                for e in echts.iter() {
                    for fld in &e.fields {
                        let missing = match fld.ftype {
                            fields::FieldType::Float => ".".to_string(),
                            fields::FieldType::Flag => "false".to_string(),
                            fields::FieldType::Categorical => fld.missing_string.clone(),
                            _ => fld.missing_value.to_string(),
                        };
                        write!(writer, "\t{}", missing)?;
                    }
                }
                writeln!(writer)?;
                n_written += 1;
            } else {
                for (ref_allele, alt_allele) in &all_variants {
                    let mut variant = BedVariant {
                        chrom: chrom.to_string(),
                        chrom_id,
                        pos: start_pos,
                        ref_allele: ref_allele.clone(),
                        alt_allele: alt_allele.clone(),
                    };

                    for (i, e) in echts.iter_mut().enumerate() {
                        e.update_expr_values(&mut variant, &mut expr_values[i]);
                    }

                    if is_compiled
                        && fasteval::eval_compiled!(compiled, &slab, &mut ns) == 0.0
                    {
                        continue;
                    }

                    write!(
                        writer,
                        "{}\t{}\t{}",
                        fields_vec[0], fields_vec[1], fields_vec[2]
                    )?;
                    for col in &extra_cols {
                        write!(writer, "\t{}", col)?;
                    }
                    write!(
                        writer,
                        "\t{}\t{}",
                        std::str::from_utf8(ref_allele).unwrap_or("."),
                        std::str::from_utf8(alt_allele).unwrap_or(".")
                    )?;
                    for e in echts.iter() {
                        for fld in &e.fields {
                            write!(writer, "\t{}", format_value(e, fld))?;
                        }
                    }
                    writeln!(writer)?;
                    n_written += 1;
                }
            }
        }
    }

    let mili = std::cmp::max(1, time::Instant::now().duration_since(start).as_millis());
    eprintln!(
        "[echtvar] processed {} BED entries ({} / second). wrote {} output rows.",
        n,
        1000 * (n as u128) / mili,
        n_written,
    );

    Ok(())
}

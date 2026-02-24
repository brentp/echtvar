#![cfg(any(feature = "bed", feature = "tab"))]

use std::collections::HashMap;
use std::error::Error;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::time;

use echtvar_lib::echtvar::{EchtVars, Value, Variant};
use echtvar_lib::fields;

use fasteval::Compiler;
use fasteval::Evaler;

/// Tabular format: BED (chrom, start 0-based, end) or Tab (chrom, start 1-based, ref, alt).
/// ref_col/alt_col are passed separately for both (see tabular_annotate_main).
#[derive(Clone, Copy)]
pub enum TabularFormat {
    /// BED: columns chrom, start (0-based), end (1-based exclusive).
    Bed,
    /// Tab: columns chrom, start (1-based), ref, alt. has_header: first line is header (preserved in output) even without # prefix.
    Tab { has_header: bool },
}

struct TabularVariant {
    chrom: String,
    chrom_id: i32,
    pos: u32, // 0-based for library
    ref_allele: Vec<u8>,
    alt_allele: Vec<u8>,
}

impl Variant for TabularVariant {
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
            _ => {
                if i == fld.missing_value {
                    ".".to_string()
                } else {
                    i.to_string()
                }
            }
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

fn set_evalues_missing(echts: &mut [EchtVars]) {
    for e in echts.iter_mut() {
        for (j, fld) in e.fields.iter().enumerate() {
            e.evalues[j] = match fld.ftype {
                fields::FieldType::Float => Value::Float(f32::NAN),
                _ => Value::Int(fld.missing_value),
            };
        }
    }
}

/// Single entry point for annotating BED or tab-format files.
/// ref_col/alt_col (1-based) apply to both formats: tab default 3,4; BED allele-specific when both set.
pub fn tabular_annotate_main(
    input_path: &str,
    compressed: bool,
    output_path: &str,
    include_expr: Option<&str>,
    epaths: Vec<&str>,
    format: TabularFormat,
    ref_col: Option<usize>,
    alt_col: Option<usize>,
) -> Result<(), Box<dyn Error>> {
    let (is_tab, allele_specific, ref_col_idx, alt_col_idx, tab_has_header) = match format {
        TabularFormat::Tab { has_header } => {
            let r = ref_col.unwrap_or(3);
            let a = alt_col.unwrap_or(4);
            (true, true, r, a, has_header)
        }
        TabularFormat::Bed => {
            let as_ = ref_col.is_some() && alt_col.is_some();
            (false, as_, ref_col.unwrap_or(0), alt_col.unwrap_or(0), false)
        }
    };

    let min_cols = if is_tab {
        ref_col_idx.max(alt_col_idx)
    } else {
        3
    };

    let mut echts: Vec<EchtVars> = epaths.iter().map(|p| EchtVars::open(p)).collect();

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

    let ipath = match input_path {
        "-" | "stdin" => "/dev/stdin",
        p => p,
    };
    let reader: Box<dyn BufRead> = if ipath == "/dev/stdin" {
        Box::new(BufReader::new(std::io::stdin()))
    } else if compressed {
        Box::new(BufReader::new(flate2::read::GzDecoder::new(
            std::fs::File::open(ipath)?,
        )))
    } else {
        Box::new(BufReader::new(std::fs::File::open(ipath)?))
    };

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

    let mut chrom_ids: HashMap<String, i32> = HashMap::new();
    let mut next_chrom_id: i32 = -2;

    let mut header_written = false;
    let mut input_header_cols: Vec<String> = Vec::new();

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

        // Tab with --has-header: first line (even without #) is the header, preserved in output
        if is_tab && tab_has_header && input_header_cols.is_empty() && !line.starts_with('#') {
            let cols: Vec<&str> = line.split('\t').collect();
            input_header_cols = cols.iter().map(|s| s.to_string()).collect();
            continue;
        }

        if line.starts_with('#') {
            if !header_written {
                let cols: Vec<&str> = line[1..].split('\t').collect();
                input_header_cols = cols.iter().map(|s| s.to_string()).collect();
            }
            continue;
        }

        let fields_vec: Vec<&str> = line.split('\t').collect();
        if fields_vec.len() < min_cols {
            eprintln!(
                "[echtvar] skipping line with fewer than {} columns: {}",
                min_cols, line
            );
            continue;
        }

        let chrom = fields_vec[0];
        let (pos_0based, extra_cols, ref_allele, alt_allele, output_prefix): (u32, Vec<&str>, Vec<u8>, Vec<u8>, Vec<&str>) = if is_tab {
            // Tab: chrom, start (1-based), ref, alt (ref/alt at ref_col_idx, alt_col_idx 1-based)
            let start_1based: u32 = match fields_vec[1].parse() {
                Ok(v) => v,
                Err(_) => {
                    eprintln!("[echtvar] skipping line with invalid start position: {}", line);
                    continue;
                }
            };
            if start_1based == 0 {
                eprintln!(
                    "[echtvar] skipping line with 0 start (tab format uses 1-based start): {}",
                    line
                );
                continue;
            }
            let pos_0based = start_1based - 1;
            let ri = ref_col_idx - 1;
            let ai = alt_col_idx - 1;
            if ri >= fields_vec.len() || ai >= fields_vec.len() {
                eprintln!(
                    "[echtvar] --ref-col or --alt-col out of range for line {}",
                    n + 1
                );
                continue;
            }
            let ref_allele = fields_vec[ri].as_bytes().to_vec();
            let alt_allele = fields_vec[ai].as_bytes().to_vec();
            // Output prefix: chrom, start, ref, alt. Extra = other columns in order (excluding chrom, start, ref, alt).
            let prefix = vec![
                fields_vec[0],
                fields_vec[1],
                fields_vec[ri],
                fields_vec[ai],
            ];
            let mut extra = Vec::new();
            for (i, &f) in fields_vec.iter().enumerate() {
                let col_1based = i + 1;
                if col_1based != 1 && col_1based != 2 && col_1based != ref_col_idx && col_1based != alt_col_idx {
                    extra.push(f);
                }
            }
            (pos_0based, extra, ref_allele, alt_allele, prefix)
        } else {
            // BED: chrom, start (0-based), end
            let start_pos: u32 = match fields_vec[1].parse() {
                Ok(v) => v,
                Err(_) => {
                    eprintln!("[echtvar] skipping line with invalid start position: {}", line);
                    continue;
                }
            };
            let extra = if fields_vec.len() > 3 {
                fields_vec[3..].to_vec()
            } else {
                Vec::new()
            };
            let prefix = fields_vec[..3].to_vec();
            let (ref_allele, alt_allele) = if allele_specific {
                let ri = ref_col_idx - 1;
                let ai = alt_col_idx - 1;
                if ri >= fields_vec.len() || ai >= fields_vec.len() {
                    eprintln!(
                        "[echtvar] --ref-col or --alt-col out of range for line {}",
                        n + 1
                    );
                    continue;
                }
                (
                    fields_vec[ri].as_bytes().to_vec(),
                    fields_vec[ai].as_bytes().to_vec(),
                )
            } else {
                (Vec::new(), Vec::new())
            };
            (start_pos, extra, ref_allele, alt_allele, prefix)
        };

        let chrom_id = *chrom_ids
            .entry(chrom.to_string())
            .or_insert_with(|| {
                let id = next_chrom_id;
                next_chrom_id -= 1;
                id
            });

        n += 1;

        if allele_specific {
            // Single-variant path: reset to missing then lookup
            set_evalues_missing(&mut echts);

            let mut variant = TabularVariant {
                chrom: chrom.to_string(),
                chrom_id,
                pos: pos_0based,
                ref_allele: ref_allele.clone(),
                alt_allele: alt_allele.clone(),
            };

            for (i, e) in echts.iter_mut().enumerate() {
                e.update_expr_values(&mut variant, &mut expr_values[i]);
            }

            if is_compiled && fasteval::eval_compiled!(compiled, &slab, &mut ns) == 0.0 {
                continue;
            }

            if !header_written {
                if is_tab {
                    write!(writer, "#chrom\tstart\tref\talt")?;
                    // Extra columns: all except chrom(1), start(2), ref_col, alt_col
                    if !input_header_cols.is_empty() && input_header_cols.len() >= fields_vec.len() {
                        for (i, col) in input_header_cols.iter().enumerate() {
                            let col_1based = i + 1;
                            if col_1based != 1 && col_1based != 2 && col_1based != ref_col_idx && col_1based != alt_col_idx {
                                write!(writer, "\t{}", col)?;
                            }
                        }
                    } else if fields_vec.len() > 4 {
                        for i in 0..fields_vec.len() {
                            let col_1based = i + 1;
                            if col_1based != 1 && col_1based != 2 && col_1based != ref_col_idx && col_1based != alt_col_idx {
                                write!(writer, "\tcol{}", col_1based)?;
                            }
                        }
                    }
                } else {
                    write!(writer, "#chrom\tstart\tend")?;
                    if input_header_cols.len() > 3 {
                        for col in &input_header_cols[3..] {
                            write!(writer, "\t{}", col)?;
                        }
                    } else if fields_vec.len() > 3 {
                        for i in 3..fields_vec.len() {
                            let col_1based = i + 1;
                            let name: String = if col_1based == ref_col_idx {
                                "REF".into()
                            } else if col_1based == alt_col_idx {
                                "ALT".into()
                            } else {
                                format!("col{}", col_1based)
                            };
                            write!(writer, "\t{}", name)?;
                        }
                    }
                    if !allele_specific {
                        write!(writer, "\tREF\tALT")?;
                    }
                }
                for name in &anno_fields {
                    write!(writer, "\t{}", name)?;
                }
                writeln!(writer)?;
                header_written = true;
            }

            write!(writer, "{}", output_prefix[0])?;
            for p in output_prefix.iter().skip(1) {
                write!(writer, "\t{}", p)?;
            }
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
            // BED position-scan only
            for e in echts.iter_mut() {
                e.set_position_by_name(chrom, pos_0based)?;
            }

            let mut all_variants: Vec<(Vec<u8>, Vec<u8>)> = Vec::new();
            for e in echts.iter() {
                for v in e.variants_at_position(pos_0based) {
                    if !all_variants.contains(&v) {
                        all_variants.push(v);
                    }
                }
            }

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
                write!(writer, "\tREF\tALT")?;
                for name in &anno_fields {
                    write!(writer, "\t{}", name)?;
                }
                writeln!(writer)?;
                header_written = true;
            }

            if all_variants.is_empty() {
                write!(writer, "{}\t{}\t{}", fields_vec[0], fields_vec[1], fields_vec[2])?;
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
                            _ => ".".to_string(),
                        };
                        write!(writer, "\t{}", missing)?;
                    }
                }
                writeln!(writer)?;
                n_written += 1;
            } else {
                for (ref_a, alt_a) in &all_variants {
                    let mut variant = TabularVariant {
                        chrom: chrom.to_string(),
                        chrom_id,
                        pos: pos_0based,
                        ref_allele: ref_a.clone(),
                        alt_allele: alt_a.clone(),
                    };

                    for (i, e) in echts.iter_mut().enumerate() {
                        e.update_expr_values(&mut variant, &mut expr_values[i]);
                    }

                    if is_compiled
                        && fasteval::eval_compiled!(compiled, &slab, &mut ns) == 0.0
                    {
                        continue;
                    }

                    write!(writer, "{}\t{}\t{}", fields_vec[0], fields_vec[1], fields_vec[2])?;
                    for col in &extra_cols {
                        write!(writer, "\t{}", col)?;
                    }
                    write!(
                        writer,
                        "\t{}\t{}",
                        std::str::from_utf8(ref_a).unwrap_or("."),
                        std::str::from_utf8(alt_a).unwrap_or(".")
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
    let format_name = if is_tab { "tab" } else { "BED" };
    eprintln!(
        "[echtvar] {} annotate: {} rows, {} written in {} ms ({}/s)",
        format_name,
        n,
        n_written,
        mili,
        (n_written as u128 * 1000) / mili,
    );

    Ok(())
}

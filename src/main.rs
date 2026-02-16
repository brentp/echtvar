extern crate echtvar_lib;
pub mod commands;

#[macro_use]
extern crate clap;
extern crate fasteval;

use commands::{annotate_cmd, encoder_cmd};
#[cfg(feature = "bed")]
use commands::bed_annotate_cmd;
use std::error::Error;

const VERSION: &str = env!("CARGO_PKG_VERSION");

fn detect_format(path: &str, format_override: Option<&str>) -> &'static str {
    if let Some(fmt) = format_override {
        match fmt.to_lowercase().as_str() {
            "vcf" => return "vcf",
            #[cfg(feature = "bed")]
            "bed" => return "bed",
            "auto" => {}
            _ => {
                #[cfg(feature = "bed")]
                eprintln!(
                    "error: unknown format '{}'. Use 'auto', 'vcf', or 'bed'.",
                    fmt
                );
                #[cfg(not(feature = "bed"))]
                eprintln!(
                    "error: unknown format '{}'. Use 'auto' or 'vcf'.",
                    fmt
                );
                std::process::exit(1);
            }
        }
    }
    let lower = path.to_lowercase();
    #[cfg(feature = "bed")]
    if lower.ends_with(".bed") || lower.ends_with(".bed.gz") {
        return "bed";
    }
    if lower.ends_with(".vcf") || lower.ends_with(".vcf.gz") || lower.ends_with(".bcf") {
        "vcf"
    } else if lower.ends_with(".bed") || lower.ends_with(".bed.gz") {
        #[cfg(not(feature = "bed"))]
        {
            eprintln!(
                "error: BED format not available (build with default features). Use -f vcf for VCF/BCF.",
            );
            std::process::exit(1);
        }
        #[cfg(feature = "bed")]
        "bed"
    } else {
        eprintln!(
            "error: cannot detect format from extension '{}'. Use -f to specify.",
            path
        );
        std::process::exit(1);
    }
}

fn main() -> Result<(), Box<dyn Error>> {
    let mut app = clap_app!(echtvar =>
        (version: VERSION)
        (author: "Brent Pedersen <bpederse@gmail.com>")
        (about: "variant encoding and annotation")
        (@setting DeriveDisplayOrder)
        (@subcommand encode =>
            (about: "create an echtvar file from a population VCF/BCF")
            (@setting DeriveDisplayOrder)
            (@arg OUTPUT: +required "output zip file")
            (@arg JSON: +required "(human)-json conf file")
            (@arg VCFS: +required ... "population vcf(s) can be split by chrom")
        )
        (@subcommand anno =>
            (about: "annotate a VCF/BCF or BED file with one or more echtvar files")
            (@setting DeriveDisplayOrder)
            (@arg echtvar: -e + takes_value number_of_values(1) ... "echtvar files to annotate with. can be specified many times")
            (@arg include: -i  +takes_value number_of_values(1) "expression that determines which variants to keep in output")
            (@arg format: -f --format +takes_value "input format: auto, vcf, or bed (default: auto)")
            (@arg ref_col: --("ref-col") +takes_value "1-based column in BED for REF allele")
            (@arg alt_col: --("alt-col") +takes_value "1-based column in BED for ALT allele")
            (@arg INPUT: +required "input vcf/bcf/bed file")
            (@arg OUTPUT: +required "path to output file (format determined by extension)")
        )
    );
    let matches = app.clone().get_matches();

    if let Some(matches) = matches.subcommand_matches("encode") {
        encoder_cmd::encoder_main(
            matches.values_of("VCFS").unwrap().collect(),
            matches.value_of("OUTPUT").unwrap(),
            matches.value_of("JSON").unwrap(),
        );
    } else if let Some(matches) = matches.subcommand_matches("anno") {
        let echt_files: Vec<_> = matches.values_of("echtvar").unwrap().collect();
        let expr = matches.value_of("include");
        let input_path = matches.value_of("INPUT").unwrap();
        let output_path = matches.value_of("OUTPUT").unwrap();

        let format_override = matches.value_of("format");
        let format = detect_format(input_path, format_override);

        let ref_col = matches
            .value_of("ref_col")
            .map(|s| s.parse::<usize>().expect("--ref-col must be a positive integer"));
        let alt_col = matches
            .value_of("alt_col")
            .map(|s| s.parse::<usize>().expect("--alt-col must be a positive integer"));

        // Validation
        if ref_col.is_some() != alt_col.is_some() {
            eprintln!("error: --ref-col and --alt-col must be specified together or not at all");
            std::process::exit(1);
        }
        if format == "vcf" && (ref_col.is_some() || alt_col.is_some()) {
            eprintln!("error: --ref-col and --alt-col are only valid with BED input");
            std::process::exit(1);
        }

        match format {
            #[cfg(feature = "bed")]
            "bed" => {
                bed_annotate_cmd::bed_annotate_main(
                    input_path,
                    output_path,
                    expr,
                    echt_files,
                    ref_col,
                    alt_col,
                )?;
            }
            "vcf" => {
                annotate_cmd::annotate_main(input_path, output_path, expr, echt_files)?;
            }
            _ => unreachable!(),
        }
    } else {
        app.print_help().ok();
        println!();
    }
    Ok(())
}

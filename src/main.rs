extern crate echtvar_lib;
pub mod commands;

extern crate fasteval;

use clap::{Parser, Subcommand};
use commands::{annotate_cmd, encoder_cmd};
#[cfg(any(feature = "bed", feature = "tab"))]
use commands::tab_annotate_cmd;
use std::error::Error;

#[derive(Parser)]
#[command(name = "echtvar")]
#[command(version = env!("CARGO_PKG_VERSION"))]
#[command(author = "Brent Pedersen <bpederse@gmail.com>")]
#[command(about = "variant encoding and annotation")]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Create an echtvar file from a population VCF/BCF
    Encode(EncodeArgs),
    /// Annotate a VCF/BCF, BED, or tabular file with one or more echtvar files
    Anno(AnnoArgs),
}

#[derive(clap::Args)]
struct EncodeArgs {
    /// Output zip file
    #[arg(required = true)]
    output: String,
    /// (Human)-readable JSON config file
    #[arg(required = true)]
    json: String,
    /// Population VCF(s); can be split by chrom
    #[arg(required = true, num_args = 1..)]
    vcfs: Vec<String>,
}

#[derive(clap::Args)]
struct AnnoArgs {
    /// Echtvar files to annotate with (can be specified many times)
    #[arg(short = 'e', long, num_args = 1.., required = true)]
    echtvar: Vec<String>,
    /// Expression that determines which variants to keep in output
    #[arg(short = 'i', long)]
    include: Option<String>,
    /// Input format: auto, vcf, bed, or tab (default: auto)
    #[arg(short = 'f', long = "format")]
    #[cfg(any(feature = "bed", feature = "tab"))]
    format: Option<String>,
    /// 1-based column in BED for REF allele
    #[arg(long = "ref-col")]
    #[cfg(feature = "bed")]
    ref_col: Option<usize>,
    /// 1-based column in BED for ALT allele
    #[arg(long = "alt-col")]
    #[cfg(feature = "bed")]
    alt_col: Option<usize>,
    /// Input vcf/bcf/bed/tab file
    #[arg(required = true)]
    input: String,
    /// Path to output file (format determined by extension)
    #[arg(required = true)]
    output: String,
}

#[cfg(any(feature = "bed", feature = "tab"))]
#[derive(Clone, Copy, PartialEq)]
enum Format {
    Vcf,
    #[cfg(feature = "bed")]
    Bed,
    #[cfg(feature = "tab")]
    Tab,
}

#[cfg(any(feature = "bed", feature = "tab"))]
fn detect_format(path: &str, format_override: Option<&str>) -> (Format, bool) {
    let lower = path.to_lowercase();
    let mut parts = lower.split('.').rev();
    format_override
        .and_then(|fmt| {
            match fmt.to_lowercase().as_str() {
                "vcf" => Some(Format::Vcf),
                #[cfg(feature = "bed")]
                "bed" => Some(Format::Bed),
                #[cfg(feature = "tab")]
                "tab" => Some(Format::Tab),
                "auto" => None,
                _ => {
                    eprintln!("error: unknown format '{}'", fmt);
                    std::process::exit(1);
                }
            }
        })
        .map(|fmt| (fmt, parts.next().unwrap_or("") == "gz"))
        .unwrap_or_else(|| {
            let ext = parts.next().unwrap_or("");
            let (ext, compressed) = if ext == "gz" {
                (parts.next().unwrap_or(""), true)
            } else {
                (ext, false)
            };
            let fmt = match ext {
                "vcf" | "bcf" => Format::Vcf,
                #[cfg(feature = "bed")]
                "bed" => Format::Bed,
                #[cfg(feature = "tab")]
                "tsv" | "txt" => Format::Tab,
                _ => {
                    eprintln!(
                        "error: Cannot detect format from extension '{}'. Use -f to specify.",
                        ext
                    );
                    std::process::exit(1);
                }
            };
            (fmt, compressed)
        })
}

#[cfg(all(test, any(feature = "bed", feature = "tab")))]
mod tests {
    use super::detect_format;
    use super::Format;

    #[test]
    fn test_detect_format_override_vcf() {
        let (fmt, compressed) = detect_format("anything.xyz", Some("vcf"));
        assert_eq!(fmt, Format::Vcf);
        assert!(!compressed);
    }

    #[test]
    fn test_detect_format_override_auto_uses_extension() {
        let (fmt, _) = detect_format("file.vcf", Some("auto"));
        assert_eq!(fmt, Format::Vcf);
    }

    #[test]
    fn test_detect_format_extension_vcf() {
        let (fmt, compressed) = detect_format("file.vcf", None);
        assert_eq!(fmt, Format::Vcf);
        assert!(!compressed);
    }

    #[test]
    fn test_detect_format_extension_vcf_gz() {
        let (fmt, compressed) = detect_format("file.vcf.gz", None);
        assert_eq!(fmt, Format::Vcf);
        assert!(compressed);
    }

    #[test]
    fn test_detect_format_extension_bcf() {
        let (fmt, _) = detect_format("file.bcf", None);
        assert_eq!(fmt, Format::Vcf);
    }

    #[cfg(feature = "bed")]
    #[test]
    fn test_detect_format_override_bed() {
        let (fmt, _) = detect_format("x", Some("bed"));
        assert_eq!(fmt, Format::Bed);
    }

    #[cfg(feature = "bed")]
    #[test]
    fn test_detect_format_extension_bed() {
        let (fmt, compressed) = detect_format("file.bed", None);
        assert_eq!(fmt, Format::Bed);
        assert!(!compressed);
    }

    #[cfg(feature = "bed")]
    #[test]
    fn test_detect_format_extension_bed_gz() {
        let (fmt, compressed) = detect_format("file.bed.gz", None);
        assert_eq!(fmt, Format::Bed);
        assert!(compressed);
    }

    #[cfg(feature = "tab")]
    #[test]
    fn test_detect_format_override_tab() {
        let (fmt, _) = detect_format("x", Some("tab"));
        assert_eq!(fmt, Format::Tab);
    }

    #[cfg(feature = "tab")]
    #[test]
    fn test_detect_format_extension_tsv() {
        let (fmt, _) = detect_format("file.tsv", None);
        assert_eq!(fmt, Format::Tab);
    }

    #[cfg(feature = "tab")]
    #[test]
    fn test_detect_format_extension_txt_gz() {
        let (fmt, compressed) = detect_format("file.txt.gz", None);
        assert_eq!(fmt, Format::Tab);
        assert!(compressed);
    }
}

fn main() -> Result<(), Box<dyn Error>> {
    let cli = Cli::parse();

    match &cli.command {
        Commands::Encode(args) => {
            let vcfs: Vec<&str> = args.vcfs.iter().map(String::as_str).collect();
            encoder_cmd::encoder_main(vcfs, &args.output, &args.json);
        }
        Commands::Anno(args) => {
            let echt_files: Vec<&str> = args.echtvar.iter().map(String::as_str).collect();
            let expr = args.include.as_deref();
            let input_path = args.input.as_str();
            let output_path = args.output.as_str();

            #[cfg(any(feature = "bed", feature = "tab"))]
            let format_override = args.format.as_deref();
            #[cfg(any(feature = "bed", feature = "tab"))]
            let (format, compressed) = detect_format(input_path, format_override);

            #[cfg(feature = "bed")]
            let ref_col = args.ref_col;
            #[cfg(feature = "bed")]
            let alt_col = args.alt_col;

            #[cfg(any(feature = "bed", feature = "tab"))]
            {
                #[cfg(feature = "bed")]
                {
                    if ref_col.is_some() != alt_col.is_some() {
                        eprintln!("error: --ref-col and --alt-col must be specified together or not at all");
                        std::process::exit(1);
                    }
                    if format == Format::Vcf && (ref_col.is_some() || alt_col.is_some()) {
                        eprintln!("error: --ref-col and --alt-col are only valid with BED input");
                        std::process::exit(1);
                    }
                }

                match format {
                    #[cfg(feature = "bed")]
                    Format::Bed => {
                        tab_annotate_cmd::tabular_annotate_main(
                            input_path,
                            compressed,
                            output_path,
                            expr,
                            echt_files.clone(),
                            tab_annotate_cmd::TabularFormat::Bed {
                                ref_col,
                                alt_col,
                            },
                        )?;
                    }
                    #[cfg(feature = "tab")]
                    Format::Tab => {
                        tab_annotate_cmd::tabular_annotate_main(
                            input_path,
                            compressed,
                            output_path,
                            expr,
                            echt_files.clone(),
                            tab_annotate_cmd::TabularFormat::Tab,
                        )?;
                    }
                    Format::Vcf => {
                        annotate_cmd::annotate_main(
                            input_path,
                            output_path,
                            expr,
                            echt_files,
                        )?;
                    }
                    #[allow(unreachable_patterns)]
                    _ => unreachable!(),
                }
            }

            #[cfg(not(any(feature = "bed", feature = "tab")))]
            {
                annotate_cmd::annotate_main(
                    input_path,
                    output_path,
                    expr,
                    echt_files,
                )?;
            }
        }
    }
    Ok(())
}

extern crate echtvar_lib;
pub mod commands;

use clap::{Parser, Subcommand};
use commands::{annotate_cmd, encoder_cmd};
use std::error::Error;

const VERSION: &str = env!("CARGO_PKG_VERSION");

#[derive(Parser)]
#[command(name = "echtvar")]
#[command(version = VERSION)]
#[command(author = "Brent Pedersen <bpederse@gmail.com>")]
#[command(about = "variant encoding and annotation", long_about = None)]
struct Args {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Create an echtvar file from a population VCF/BCF
    Encode {
        /// Output zip file
        output: String,
        /// (human)-json conf file
        json: String,
        /// Population vcf(s) can be split by chrom
        #[arg(required = true)]
        vcfs: Vec<String>,
    },
    /// Annotate a VCF/BCF with one or more echtvar files
    Anno {
        /// Echtvar files to annotate with (can be specified many times)
        #[arg(short, long, required = true)]
        echtvar: Vec<String>,
        /// Expression that determines which variants to keep in output
        #[arg(short, long)]
        include: Option<String>,
        /// Input vcf or bcf
        input_vcf: String,
        /// Path to bcf/vcf.gz output file (determined by extension)
        output_vcf: String,
    },
}

fn main() -> Result<(), Box<dyn Error>> {
    let args = Args::parse();

    match args.command {
        Commands::Encode {
            output,
            json,
            vcfs,
        } => {
            let vcfs_refs: Vec<&str> = vcfs.iter().map(|s| s.as_str()).collect();
            encoder_cmd::encoder_main(vcfs_refs, &output, &json);
        }
        Commands::Anno {
            echtvar,
            include,
            input_vcf,
            output_vcf,
        } => {
            let echtvar_refs: Vec<&str> = echtvar.iter().map(|s| s.as_str()).collect();
            annotate_cmd::annotate_main(&input_vcf, &output_vcf, include.as_deref(), echtvar_refs)?;
        }
    }
    Ok(())
}

extern crate echtvar_lib;
pub mod commands;

#[macro_use]
extern crate clap;
use commands::{annotate_cmd, encoder_cmd};
use std::error::Error;

const VERSION: &str = "0.0.1";

fn main() -> Result<(), Box<dyn Error>> {
    let mut app = clap_app!(echtvar =>
        (version: VERSION)
        (author: "Brent Pedersen <bpederse@gmail.com")
        (about:"variant encoding and annotation")
        (@subcommand encode =>
            (about: "create an echtvar file from a population VCF/BCF")
            (@arg VCF: +required "population vcf")
            (@arg OUTPUT: +required "output zip file")
            (@arg JSON: +required "(human)-json conf file")
        )
        (@subcommand anno =>
            (about: "annotate a VCF/BCF with one or more echtvar files")
            (@arg echtvar: -e ... +takes_value number_of_values(1) "echtvar files to annotate with. can be specified many times")
            (@arg include: -i  +takes_value number_of_values(1) "expression that determines which variants to keep in output")
            (@arg INPUT_VCF: +required "vcf")
            (@arg OUTPUT_VCF: +required "path to bcf output file")
        )
    );
    let matches = app.clone().get_matches();

    if let Some(matches) = matches.subcommand_matches("encode") {
        encoder_cmd::encoder_main(
            matches.value_of("VCF").unwrap(),
            matches.value_of("OUTPUT").unwrap(),
            matches.value_of("JSON").unwrap(),
        );
    } else if let Some(matches) = matches.subcommand_matches("anno") {
        let echt_files: Vec<_> = matches.values_of("echtvar").unwrap().collect();
        let expr = matches.value_of("include");

        annotate_cmd::annotate_main(
            matches.value_of("INPUT_VCF").unwrap(),
            matches.value_of("OUTPUT_VCF").unwrap(),
            expr,
            echt_files,
        )?;
    } else {
        app.print_help().ok();
        print!("\n");
    }
    Ok(())
}

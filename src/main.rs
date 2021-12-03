extern crate echtvar_lib;
pub mod commands;

#[macro_use]
extern crate clap;
use commands::{annotate_cmd, encoder_cmd};

const VERSION: &str = "0.0.1";

fn main() {
    let matches = clap_app!(echtvar =>
        (version: VERSION)
        (author: "Brent Pedersen <bpederse@gmail.com")
        (about:"variant encoding and annotation")
        (@subcommand encode =>
            (about: "encode an echtvar file from a population VCF/BCF")
            (@arg VCF: +required "population vcf")
            (@arg OUTPUT: +required "output zip file")
            (@arg JSON: +required "(human)-json conf file")
        )
        (@subcommand anno =>
            (about: "annotate a VCF/BCF with one or more echtvar files")
            (@arg echtvar: -e ... +takes_value number_of_values(1) "echtvar files to annotate with. can be specified many times")
            (@arg INPUT_VCF: +required "vcf")
            (@arg OUTPUT_VCF: +required "path to bcf output file")
        )
    ).get_matches();

    if let Some(matches) = matches.subcommand_matches("encode") {
        encoder_cmd::encoder_main(
            matches.value_of("VCF").unwrap(),
            matches.value_of("OUTPUT").unwrap(),
            matches.value_of("JSON").unwrap(),
        );
    } else if let Some(matches) = matches.subcommand_matches("anno") {
        let echt_files: Vec<_> = matches.values_of("echtvar").unwrap().collect();
        annotate_cmd::annotate_main(
            matches.value_of("INPUT_VCF").unwrap(),
            matches.value_of("OUTPUT_VCF").unwrap(),
            echt_files,
        );
    }
}

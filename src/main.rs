extern crate echtvar_lib;
pub mod commands;

#[macro_use]
extern crate clap;
use commands::encoder_cmd;

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
    ).get_matches();

    if let Some(matches) = matches.subcommand_matches("encode") {
        eprintln!("{} {} {}", matches.value_of("VCF").unwrap(), matches.value_of("OUTPUT").unwrap(), matches.value_of("JSON").unwrap());
        encoder_cmd::encoder_main(matches.value_of("VCF").unwrap(), matches.value_of("OUTPUT").unwrap(), matches.value_of("JSON").unwrap());
    }



}

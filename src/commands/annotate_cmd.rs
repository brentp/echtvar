use std::fs;
use std::io;
use std::io::prelude::*;

use rust_htslib::bcf::record::{Buffer, Record};
use rust_htslib::bcf::{Read as BCFRead, Reader};
use rust_htslib::bcf::{Format, Writer, header::Header};

use byteorder::{LittleEndian, ReadBytesExt};


pub fn annotate_main(vpath: &str, opath: &str, epaths: Vec<&str>) -> io::Result<()> {

    let mut vcf = Reader::from_path(vpath).ok().expect("Error opening vcf.");
    vcf.set_threads(2).ok();
    let header_view = vcf.header().clone();
    let mut buffer = Buffer::new();
    let header = Header::from_template(&header_view);

    // TODO: handle stdout
    let mut ovcf = Writer::from_path(opath, &header, false, Format::Bcf).ok().expect("error opening bcf for output");

    let ep = std::path::Path::new(&*epaths[0]);
    let file = fs::File::open(ep).expect("error accessing zip file");
    let mut archive = zip::ZipArchive::new(&file).expect("error opening zip file");

    
    let mut iz = archive.by_name("echtvar/chr21/7/var32.bin").expect("unable to open file");


    let n = iz.read_u32::<LittleEndian>().ok().expect("error reading number of values from zip file");
    let mut comr = vec![0 as u8; iz.size() as usize - 4];
    iz.read_to_end(&mut comr)?;

    eprintln!("{} {}", n, iz.size());


    Ok(())
}

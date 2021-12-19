use std::fs;
use std::io;
use std::io::prelude::*;

use rust_htslib::bcf::record::{Buffer, Record};
use rust_htslib::bcf::{Read as BCFRead, Reader};
use rust_htslib::bcf::{Format, Writer, header::Header};

use echtvar_lib::echtvar::EchtVars;

use byteorder::{LittleEndian, ReadBytesExt};

use stream_vbyte::{
    decode::decode,
    x86::Ssse3
};


pub fn annotate_main(vpath: &str, opath: &str, epaths: Vec<&str>) -> io::Result<()> {

    let mut vcf = Reader::from_path(vpath).ok().expect("Error opening vcf.");
    vcf.set_threads(2).ok();
    let header_view = vcf.header().clone();
    let mut buffer = Buffer::new();
    let mut header = Header::from_template(&header_view);

    let mut e = EchtVars::open(&*epaths[0]);
    e.update_header(&mut header);

    // TODO: handle stdout
    let mut ovcf = Writer::from_path(opath, &header, false, Format::Bcf).ok().expect("error opening bcf for output");
    ovcf.set_threads(2);
    let oheader_view = ovcf.header().clone();

    for r in vcf.records() {
        let mut record = r.expect("failed to read record");
        // TODO:
        if e.check_and_update_variant(&mut record, &oheader_view) {
            ovcf.write(&record).expect("failed to write record");
        }
    }

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

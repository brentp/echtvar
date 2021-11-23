# TODO
make a bitpacking crate that chunks the data for Bitpacker4x


4 bits * 3 billion == 1.5GB

https://immunant.com/blog/2020/01/bitfields/

usage
=====

make a new afenc file 
```
echtvar \
   encode \
   # set a multiplier as internally, everything is stored as int
   --field AF:gnomad_af:1000000 \ 
   # integers will have a multiplier of 1.
   --field AC:gnomac_ac \
   -o gnomad.afenc \ # output file name
   $input_vcf

```

annotate a VCF with an afen file

```
echtvar annotate \
   -o $cohort.afenc.bcf \
   -a gnomad.afenc \
   -a ukbiobank.afenc \
   -f 'gnomad_af < 0.01' \
   $cohort.bcf
```

# rust

```rust
fn annotate(&self, pos:uint64, ref:string, alt:string) -> Result<Vec<int32>, Error> 

extern crate libc;

#[repr(C), align(1)]
#[derive(BitfieldStruct, Clone, Copy)]
pub struct Var32 {
    #[bitfield(name = "enc", ty = "libc::c_uint32", bits = "0..4")]
    #[bitfield(name = "alen", ty = "libc::c_uint32", bits = "4..8")]
    #[bitfield(name = "rlen", ty = "libc::c_uint32", bits = "8..12")]
    #[bitfield(name = "position", ty = "libc::c_uint32", bits = "12..32")]
    datata: [u8, 4]
}
```

type Var64 = struct {
  



bits
====

  01: AF < 1e-5
  10: AF < 1e-3
  11: AF < 1e-2
 1  : AN < 70%


  01: AF < 1e-2
  10: AF < 1e-3
  11: AN >k


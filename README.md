# TODO
make a bitpacking crate that chunks the data for Bitpacker4x


4 bits * 3 billion == 1.5GB

https://immunant.com/blog/2020/01/bitfields/

usage
=====

make a new echtvar file 
```
echtvar \
   encode \
   $input_vcf \
   echtvar-gnomad-v3.zip \
   conf.json

```

annotate a VCF with an echtvar file

```
echtvar annotate \
   -o $cohort.echtvar.bcf \
   -a gnomad.echtvar \
   -a ukbiobank.echtvar \
   -f 'gnomad_af < 0.01' \
   $cohort.bcf
```

# rust


bits
====

  01: AF < 1e-5
  10: AF < 1e-3
  11: AF < 1e-2
 1  : AN < 70%


  01: AF < 1e-2
  10: AF < 1e-3
  11: AN >k


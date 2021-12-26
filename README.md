# Echtvar
## Really, truly rapid variant annotation and filtering 

Echtvar enables rapid annotation of variants with huge pupulation datasets and
it supports filtering on those values. It chunks the genome into 1<<20 (~1 million
base) chunks, encodes each variant into a 32 bit integer (with a supplemental table
for those that can't fit). It uses [delta
encoding](https://en.wikipedia.org/wiki/Delta_encoding)
and [integer compression
](https://lemire.me/blog/2017/09/27/stream-vbyte-breaking-new-speed-records-for-integer-compression/)
to create a compact and searchable format of any integer or float columns
selected from the population file.

Once created, an echtvar file can be used to annotate variants in a VCF (or
BCF) file at a rate of >1 million variants per second.

A filter expression can be applied so that only variants that meet the
expression are written.


### usage

make (`encode`) a new echtvar file 

```
echtvar \
   encode \
   $input_vcf \
   echtvar-gnomad-v3.zip \
   conf.json

```

annotate a VCF with an echtvar file and only output variants where `gnomad_af`
from the echtvar file is < 0.01.

```
echtvar annotate \
   -o $cohort.echtvar-annotated.bcf \
   -a gnomad.echtvar \
   -a ukbiobank.echtvar \
   -i 'gnomad_af < 0.01' \
   $cohort.input.bcf
```


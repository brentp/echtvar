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
expression are written. Since `echtvar` is so fast, writing the output is a bottleneck
so filtering can actually *increase* the speed.

### usage

make (`encode`) a new echtvar file 

```
echtvar \
   encode \
   $input_population_vcf \
   echtvar-gnomad-v3.zip \
   conf.json # this defines the columns to pull from $input_vcf, and how to
name and encode them

```

See below for a description of the json file that defines which columns are
pulled from the population VCF.

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

### Configuration File for Encode

When running `echtvar encode`, a [json5](https://json5.org/) (json with
comments and other nice features) determines which columns are pulled from the
input VCF and how they are stored.

A simple example is to pull a single integer field and give it a new name (`alias`):

```
[{"field": "AC", "alias": "gnomad_AC"}]
```

This will extract the "AC" field from the INFO and labeled as "gnomad_AC" when
later used to annotate a VCF.

We can add more fields like this:

```
[
    {"field": "AC", "alias": "gnomad_AC"},
    // this JSON file is json 5 and so can have comments
    // the missing value will default to -1, but the value: -2147483648 will
    // result in '.' as it is the missing value for VCF.
    {"field": "AN", "alias":, gnomad_AN", missing_value: -2147483648},
    {
           field: "AF",
           alias: "gnomad_AF",
           missing_value: -1,
           // for floats, upon annotation, the score is divided by multiplier and stored as an integer.
           multiplier: 2000000,
          ftype: "Float"
   }
]
```

This will exctract 3 fields, the user can chooose as many as they like when encoding.
All fields in an `echtvar` file will be added (with the given alias) to any VCF it is used to annotate.

# References and Acknowledgements

Without these (and other) critical libraries, `echtvar` would not exist.

+ [rust-htslib](https://github.com/rust-bio/rust-htslib) is used for reading and writing BCF and VCF.
+ [stream-vbyte](https://lemire.me/blog/2017/09/27/stream-vbyte-breaking-new-speed-records-for-integer-compression/) is used for integer compression via the [excellent rust bindings](https://bitbucket.org/marshallpierce/stream-vbyte-rust/src/master/)
+ [fasteval](https://github.com/likebike/fasteval) is used for the expressions. It is fast and simple.

# Echtvar
## Really, truly rapid variant annotation and filtering 

[![Rust](https://github.com/brentp/echtvar/actions/workflows/ci.yml/badge.svg)](https://github.com/brentp/echtvar/actions/workflows/ci.yml)

Echtvar enables rapid annotation of variants with huge pupulation datasets and
it supports filtering on those values. It chunks the genome into 1<<20 (~1 million
) bases, encodes each variant into a 32 bit integer (with a supplemental table
for those that can't fit due to large REF and/or ALT alleles). It uses [delta
encoding](https://en.wikipedia.org/wiki/Delta_encoding)
and [integer compression
](https://lemire.me/blog/2017/09/27/stream-vbyte-breaking-new-speed-records-for-integer-compression/)
to create a compact and searchable format of any integer or float columns
selected from the population file.

Once created, an echtvar file can be used to annotate variants in a VCF (or
BCF) file at a rate of >1 million variants per second (most of the time is spent
reading and writing VCF/BCF, so this number depends on the particular file).

A filter expression can be applied so that only variants that meet the
expression are written. Since `echtvar` is so fast, writing the output is a bottleneck
so filtering can actually *increase* the speed.

### usage

make (`encode`) a new echtvar file. this is usually done once  (or download from those provided in the Release pages) 
and then the file can be re-used for the `annotate` step with each new query file.

```
echtvar \
   encode \
   echtvar-gnomad-v3.zip \
   conf.json # this defines the columns to pull from $input_vcf, and how to
   $input_population_vcf[s] \ can be split by chromosome or all in a single file.
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
   -i 'gnomad_af < 0.01' \
   $cohort.input.bcf
```

#### Configuration File for Encode

When running `echtvar encode`, a [json5](https://json5.org/) (json with
comments and other nice features) determines which columns are pulled from the
input VCF and how they are stored.

A simple example is to pull a single integer field and give it a new name (`alias`):

```
[{"field": "AC", "alias": "gnomad_AC"}]
```

This will extract the "AC" field from the INFO and labeled as "gnomad_AC" when
later used to annotate a VCF. Note that it's important to give a description/unique prefix lke "`gnomad_`" so
as not to collide with fields already in the query VCF.

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
           // higher values give better precision and worse compression.
           multiplier: 2000000,
   }
]
```

The above file will extract 3 fields, but the user can chooose as many as they like when encoding.
All fields in an `echtvar` file will be added (with the given alias) to any VCF it is used to annotate.

#### Expressions

An optional expression will determine which variants are written. It can utilize any (and only) fields present in the
echtvar file (not those present in the query VCF). An example could be:

```
-i 'gnomad_af < 0.01 && gnomad_nhomalts < 10'
```

The expressions are enabled by [fasteval](https://github.com/likebike/fasteval) with supported syntax detailed [here](https://docs.rs/fasteval/latest/fasteval/). 

In brief, the normal operators: (`&&, ||, +, -, *, /, <, <=, >, >=` and groupings `(, )`, etc) are supported and can be used to
craft an expression that returns true or false as above.

# References and Acknowledgements

Without these (and other) critical libraries, `echtvar` would not exist.

+ [rust-htslib](https://github.com/rust-bio/rust-htslib) is used for reading and writing BCF and VCF.
+ [stream-vbyte](https://lemire.me/blog/2017/09/27/stream-vbyte-breaking-new-speed-records-for-integer-compression/) is used for integer compression via the [excellent rust bindings](https://bitbucket.org/marshallpierce/stream-vbyte-rust/src/master/)
+ [fasteval](https://github.com/likebike/fasteval) is used for the expressions. It is fast and simple and awesome.
+ [bincode](https://docs.rs/bincode/latest/bincode/) is used for rapid serialization of large variants.


`echtvar` is developed in the [Jeroen De Ridder lab](https://www.umcutrecht.nl/en/research/researchers/de-ridder-jeroen-j)

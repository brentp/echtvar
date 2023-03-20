## Echtvar: Really, truly rapid variant annotation and filtering 
[![Rust](https://github.com/brentp/echtvar/actions/workflows/ci.yml/badge.svg)](https://github.com/brentp/echtvar/actions/workflows/ci.yml)

Echtvar efficiently encodes variant allele frequency and other information from huge population datasets to enable rapid (1M variants/second) annotation of genetic variants.
It chunks the genome into 1<<20 (~1 million) bases,
[encodes each variant into a 32 bit integer](https://github.com/brentp/echtvar/blob/02774b8d1cd3703b65bd2c8d7aab93af05b7940f/src/lib/var32.rs#L9-L21) (with a [supplemental table](https://github.com/brentp/echtvar/blob/02774b8d1cd3703b65bd2c8d7aab93af05b7940f/src/lib/var32.rs#L33-L38)
for those that can't fit due to large REF and/or ALT alleles). It uses the zip format, [delta
encoding](https://en.wikipedia.org/wiki/Delta_encoding)
and [integer compression
](https://lemire.me/blog/2017/09/27/stream-vbyte-breaking-new-speed-records-for-integer-compression/)
to create a compact and searchable format of any integer, float, or low-cardinality string columns
selected from the population file.

read more at the [why of echtvar](https://github.com/brentp/echtvar/wiki/why)

### Getting started.

Get a static binary and pre-encoded echtvar files for gnomad v3.1.2 (hg38) here: https://github.com/brentp/echtvar/releases/latest
That page contains exact instructions to get started with the static binary.

<details>
  <summary>:arrow_down:Download or Build instructions for linux</summary>

The linux binary is available via:

```
wget -O ~/bin/echtvar https://github.com/brentp/echtvar/releases/latest/download/echtvar \
    && chmod +x ~/bin/echtvar \
    && ~/bin/echtvar # show help
 ```

Users can make their own *echtvar* archives with `echtvar encode`, and pre-made archives for
gnomAD version 3.1.2 are [here](https://github.com/brentp/echtvar/release)

Rust users can build on linux with:

```
cargo build --release --target x86_64-unknown-linux-gnu
```

</details>

To run echtvar with an existing archive (we have several available in [releases](https://github.com/brentp/echtvar/releases/latest)) is as simple as
```
echtvar anno -e gnomad.echtvar.zip -e other.echtvar.zip input.vcf output.annotated.bcf
```

an optional filter that utilizes fields available any of the zip files can be added like:
```
-i "gnomad_popmax_af < 0.01"
```

echtvar can also accept input from stdin using "-" or "/dev/stdin" for the input argument.

### usage

##### encode 

make (`encode`) a new echtvar file. This is usually done once  (or download from those provided in the [Release pages](https://github.com/brentp/echtvar/releases/latest)) 
and then the file can be re-used for the annotation (`echtvar anno`) step with each new query file.
Note that input VCFs must be [decomposed](https://github.com/brentp/echtvar/wiki/decompose).

```
echtvar \
   encode \
   gnomad.v3.1.2.echtvar.zip \
   conf.json # this defines the columns to pull from $input_vcf, and how to
   $input_population_vcf[s] \ can be split by chromosome or all in a single file.
name and encode them

```

See below for a description of the json file that defines which columns are
pulled from the population VCF.

##### annotate 

Annotate a [**decomposed** (and normalized)](https://github.com/brentp/echtvar/wiki/decompose) VCF with an echtvar file and only output variants where `gnomad_af`
from the echtvar file is < 0.01. Note that multiple echtvar files can be specified
and the `-i` expression is optional and can be elided to output all variants.

```
echtvar anno \
   -e gnomad.v3.1.2.echtvar.v2.zip \
   -e dbsnp.echtvar.zip \
   -i 'gnomad_popmax_af < 0.01' \
   $cohort.input.bcf \
   $cohort.echtvar-annotated.filtered.bcf
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

<details>
  <summary>:arrow_down:Expand this section for detail on additional fields, including float and string types</summary>

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
           // since all values (including floats) are stored as integers, echtvar internally converts
           // any float to an integer by multiplying by `multiplier`.
           // higher values give better precision and worse compression.
           // upon annotation, the score is divided by multiplier to give a number close to the original float.
           multiplier: 2000000,
           // set zigzag to true if your data has negative values
           zigzag: true,
   }
    // echtvar will save strings as integers along with a lookup. this can work for fields with a low cardinality.
    {"field": "string_field", "alias":, gnomad_string_field", missing_string: "UNKNOWN"},
    // "FILTER" is a special case that indicates that echtvar should extract the FILTER column from the annotation vcf.
    {"field": "FILTER", "alias": "gnomad_filter"},
]
```

The above file will extract 5 fields, but the user can chooose as many as they like when encoding.
All fields in an `echtvar` file will be added (with the given alias) to any VCF it is used to annotate.

</details>

Other examples are available [here](https://github.com/brentp/echtvar/tree/main/examples)

And full examples are in the [wiki](https://github.com/brentp/echtvar/wiki)

#### Expressions

An optional expression will determine which variants are written. It can utilize any (and only) integer or float fields present in the
echtvar file (not those present in the query VCF). An example could be:

```
-i 'gnomad_af < 0.01 && gnomad_nhomalts < 10'
```

The expressions are enabled by [fasteval](https://github.com/likebike/fasteval) with supported syntax detailed [here](https://docs.rs/fasteval/latest/fasteval/). 

In brief, the normal operators: (`&&, ||, +, -, *, /, <, <=, >, >=` and groupings `(, )`, etc) are supported and can be used to
craft an expression that returns true or false as above.

# References and Acknowledgements

Without these (and other) critical libraries, `echtvar` would not exist.

+ [htslib](https://github.com/samtools/htslib) is used for reading and writing BCF and VCF via [rust-htslib](https://github.com/rust-bio/rust-htslib)
+ [stream-vbyte](https://lemire.me/blog/2017/09/27/stream-vbyte-breaking-new-speed-records-for-integer-compression/) is used for integer compression via the [excellent rust bindings](https://bitbucket.org/marshallpierce/stream-vbyte-rust/src/master/)
+ [fasteval](https://github.com/likebike/fasteval) is used for the expressions. It is fast and simple and awesome.
+ [bincode](https://docs.rs/bincode/latest/bincode/) is used for rapid serialization of large variants.


`echtvar` is developed in the [Jeroen De Ridder lab](https://www.deridderlab.nl/)

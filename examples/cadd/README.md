This directory contains scripts to convert [CADD](https://cadd.gs.washington.edu/download) to an `echtvar` archive.

**Note that a license is required for commercial use of CADD scores.**

Users can create their own archive like this:
```
# See https://cadd.gs.washington.edu/download for links to older versions, hg19, and US mirrors
wget -c https://kircherlab.bihealth.org/download/CADD/v1.7/GRCh38/whole_genome_SNVs.tsv.gz
python cadd2vcf.py whole_genome_SNVs.tsv.gz | bcftools view -O b -o cadd.v1.7.bcf
echtvar encode cadd.v1.7.hg38.zip cadd.json cadd.v1.7.bcf
```

with these `cadd2vcf.py` from this directory and contents for `cadd.json`:
```
[{
        "field": "raw",
        "alias": "cadd_raw",
        multiplier: 100,
        zigzag: true,
        missing_value: -2,
}, {
        field: "phred",
        alias: "cadd_phred",
        multiplier: 100,
        zigzag: true,
        missing_value: -2,
}]
```


The `run.sh` file contains code to recreate analyses and timings in the
`echtvar` paper.

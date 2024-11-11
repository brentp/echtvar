
v0.2.1
======
+ bump stream-vbyte for performance annotating
+ fix #43 where output vcf would have Number=, unless number was specified in config.json

v0.2.0
======
+ skip non ACGT alts (@migbro #39)
+ grab description from input vcf and propagate to annotated VCFs (@dmiller15, @migbro #34, #35)

v0.1.9
======
+ report missing float as . (thanks @dmiller15 #32)

v0.1.8
======
+ better help message indicating how to output bcf/vcf.gz (thanks @dvg-p4)
+ fix missing value in cadd example (thanks @dvg-p4)
+ doc types (thanks @enormandeau)
+ handle big files (#30)

v0.1.7 
======
+ add wiki page and more information about decomposing variants
+ clarify warning on non ACGT (#23)
+ better error messages on missing/bad zip file

v0.1.6
======
+ change VCF header (#19)
+ output fewer warnings for '*' and other alternate alleles (#22)

v0.1.5
======
+ fix build on OSX (#18 from @sstadick)
+ handle empty alts (#16 from @sstadick)
+ add script for converting dbNSFP -> VCF -> echtvar encode

v0.1.4
======
+ fix bug with filter not removing (not filtering) variants correctly.

v0.1.3
======
+ exit with error on multi-allelics (previously, only first allele was used) ( #11)
+ more docs on conf file (thanks @m-pauper)
+ fix compression error for very dense regions (#12)
+ fix error when using hts_set_threads more than once and after reading header (#12)

v0.1.2
======
+ fix #8. bug introduced for string values.

v0.1.1
======
+ support multiple echtvar files to `echtvar anno`. This makes it possible to annoate with, for example gnomad *and* dbsnp.
+ support string fields. This works for fields with a lowish number (max of thousands) of unique values.
+ support extracting FILTER field from annotation VCF. So now can annotate with, e.g. gnomad FILTER field via:
  ```
  {"field": "FILTER", "alias": "gnomad_filter"},
  ```
  in the json file for `echtvar encode`. "FILTER" is a special-case, all other `field`s are from the INFO.
  

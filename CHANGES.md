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
  

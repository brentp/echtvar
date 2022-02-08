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
  

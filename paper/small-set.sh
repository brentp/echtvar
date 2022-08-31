set -euo pipefail
ref=/hpc/cog_bioinf/GENOMES.old/1KP_GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.fa
bcftools csq -l -O u -o - --threads 3 -s - -f $ref -g Homo_sapiens.GRCh38.104.chr.gff3.gz -c CSQ HG001_GRCh38_1_22_v4.2.1_benchmark.bcf \
| slivar expr --info 'INFO.highest_impact_order < ImpactOrder.missense' -v - --pass-only -o HG001.bad.bcf

echtvar anno -e cadd.v1.6.hg38.zip HG001.bad.bcf HG001.bad.cadd.bcf

set -euo pipefail

ref=/hpc/cog_bioinf/GENOMES.old/1KP_GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.fa
python generate-all-snps.py $ref \
    | bcftools csq -l -O u -o - --threads 3 -s - -f $ref -g Homo_sapiens.GRCh38.104.chr.gff3.gz -c CSQ - \
    | bcftools +split-vep /dev/stdin -O u -c gene,Consequence -s worst \
    | bcftools annotate --threads 3 -x INFO/CSQ -O b -o hg38.csq.bcf

echtvar encode echtvar.impact.hg38.zip impact.json hg38.csq.bcf
echtvar encode echtvar.gene.hg38.zip gene.json hg38.csq.bcf

echtvar anno -e echtvar.impact.hg38.zip hg38.csq.bcf hg38.self-impact.bcf
# [echtvar] evaluated 8771197662 variants (549113 / second). wrote 8771197662 variants.
python test-self-impact.py hg38.self-impact.bcf


<<DONE
wget -O final_consensus_snv_indel_passonly_icgc.public.tgz https://dcc.icgc.org/api/v1/download?fn=/PCAWG/consensus_snv_indel/final_consensus_snv_indel_passonly_icgc.public.tgz
tar xzf final_consensus_snv_indel_passonly_icgc.public.tgz
wget -q https://data.broadinstitute.org/snowman/ofer/ICGC/data/merged_lite1.6.2.txt
awk 'BEGIN{FS=OFS="\t"}{ print $2, $20,$21,$22,$23 }' merged_lite1.6.2.txt | uniq > cancer-types.txt
DONE

<<SETUP_ARCHIVE

python ~/dbnsfp.py dbNSFP4.3a.zip \
 -f SIFT_converted_rankscore \
 -f DANN_rankscore \
 -f GERP++_RS_rankscore \
--json > dbNSFP.json

python ~/dbnsfp.py dbNSFP4.3a.zip \
 -f SIFT_converted_rankscore \
 -f DANN_rankscore \
 -f GERP++_RS_rankscore \
  | echtvar encode dbNSFP.echtvar.zip dbNSFP.json /dev/stdin

SETUP_ARCHIVE

set -uo pipefail

echo "STARTING"
for f in snv_mnv/*.bcf; do
    /usr/bin/time -v echtvar anno -e dbNSFP.echtvar.zip $f /dev/null -i 'dbsnfp_SIFT_converted_rankscore > 0.2 || dbsnfp_DANN_rankscore > 0.2 || dbsnfp_GERPpp_RS_rankscore > 0.2 ' 2> $f.time.txt
    echo "OK $f"
done

python plot-icgc.py snv_mnv/*.txt

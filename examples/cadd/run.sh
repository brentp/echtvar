set -eu
#wget -c https://kircherlab.bihealth.org/download/CADD/v1.6/GRCh38/whole_genome_SNVs.tsv.gz
#wget -c https://krishna.gs.washington.edu/download/CADD/v1.6/GRCh38/gnomad.genomes.r3.0.indel.tsv.gz

#python cadd2vcf.py whole_genome_SNVs.tsv.gz | bcftools view -O b -o cadd.v1.6.bcf 
#python cadd2vcf.py gnomad.genomes.r3.0.indel.tsv.gz | bcftools view -O b -o cadd.v1.6.indels.bcf 

#bcftools index --threads 3 cadd.v1.6.bcf &
#bcftools index --threads 3 cadd.v1.6.indels.bcf &
#echtvar encode cadd.v1.6.hg38.indels.zip cadd.json cadd.v1.6.indels.bcf  &
#echtvar encode cadd.v1.6.hg38.zip cadd.json cadd.v1.6.bcf 
wait



<<DONE
(/usr/bin/time -v echtvar anno -e cadd.v1.6.hg38.zip cadd.v1.6.bcf cadd.self-anno.bcf) | tee cadd-self-timing.txt
(/usr/bin/time -v echtvar anno -e cadd.v1.6.hg38.indels.zip cadd.v1.6.indels.bcf cadd.self-anno.indels.bcf ) | tee cadd-self-timing.indels.txt

python check.py cadd.self-anno.indels.bcf 
python check.py cadd.self-anno.bcf
DONE

(/usr/bin/time -v bcftools annotate \
        -a cadd.v1.6.bcf \
        --threads 4 \
        --columns "INFO/cadd_phred:=INFO/phred,INFO/cadd_raw:=raw" \
        -O b -o cadd.v1.6.bcftanno.bcf \
        cadd.v1.6.bcf) | tee cadd-self.bcftools.timing.txt


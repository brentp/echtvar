set -exo pipefail

qvcf=HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz
qbcf=${qvcf/.vcf.gz/.bcf}
avcf=gnomad.v3.1.2.concat.vcf.gz
abcf=gnomad.v3.1.2.concat.bcf
fasta=../human/GRCh38_full_analysis_set_plus_decoy_hla.fa

wget -qO - https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/NISTv4.2.1/GRCh38/HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz \
  | bcftools norm -m - -f $fasta -O v \
  | awk 'BEGIN{FS=OFS="\t"} ($5 != "*")' \
  | bcftools view -O z -o $qvcf --threads 3

bcftools index --force $qvcf --threads 3
bcftools view -O b -o $qbcf --threads 3 $qvcf
bcftools index --force $qbcf --threads 3
#DONE

<<DONE
bcftools concat $avcf \
    | slivar make-gnotate --prefix gnomad.slivar.v3.1.2 \
	-f AC:gnomad_ac \
	-f AN:gnomad_an \
	-f nhomalt:gnomad_nhomalt \
	-f AF:gnomad_af \
	-f AC_popmax:gnomad_popmax_ac \
	-f AN_popmax:gnomad_popmax_an \
	-f nhomalt_popmax:gnomad_popmax_nhomalt \
	-f AF_popmax:gnomad_popmax_af \
	-f AF_controls_and_biobanks:gnomad_controls_and_biobanks_af \
	-f nhomalt_controls_and_biobanks:gnomad_controls_and_biobanks_nhomalt \
      /dev/stdin

DONE

#<<SLIVAR
/usr/bin/time -v slivar expr \
    -v $qbcf \
    -g gnomad.slivar.v3.1.2.zip -o slivar-anno.bcf 2> slivar-time.txt
#SLIVAR



#<<ECHTVAR
/usr/bin/time -v \
  ~/src/echtvar/target/release/echtvar anno \
  -e gnomad.v3.1.2.echtvar.v2.zip \
  $qbcf annotated.bcf 2> echtvar-time.txt

#ECHTVAR


## VarNote prep

<<VARNOTE
java -Xmx5G -jar ~/Downloads/VarNote-1.2.0.jar Index \
  -I $avcf -O varnote-index 

echo "
@gnomAD
fields=[INFO]
info_fields=[AC, AF]
out_names=[AC:gnomAD_AC, AF:gnomAD_AF]
" > annoc
VARNOTE

#<<VARNOTE

export TMPDIR=/tmp/
mkdir -p null/temp # varnote needs this but won't create it.

/usr/bin/time -v \
java -Xmx5G -jar ~/Downloads/VarNote-1.2.0.jar Annotation \
  --query-file $qvcf \
  --d-files:db,tag=gnomAD,mode=1 $avcf  \
  --out-format VCF \
  -loj true \
  --allowLargeVariants \
  --maxVariantLength 500 \
  --thread 4 \
  --force-overlap true \
  -A annoc \
  --out-file varnote.anno.vcf.gz 2> varnote-time.txt

#VARNOTE


#<<BCFTOOLS
/usr/bin/time -v \
bcftools annotate \
  -a $abcf \
  --threads 4 \
  --columns "INFO/gnomad_AC:=INFO/AC,INFO/gnomad_AF:=INFO/AF" \
  -O b \
  -o bcft.anno.bcf \
  $qbcf 2> bcftools-time.txt

#BCFTOOLS

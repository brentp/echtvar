cd /work
mkdir -p compare
cd compare

set -exo pipefail


qvcf=HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz
qbcf=${qvcf/.vcf.gz/.bcf}
avcf=/data/gnomad/gnomad.v3.1.2.concat.vcf.gz
abcf=/data/gnomad/gnomad.v3.1.2.concat.bcf
fasta=/data/human/GRCh38_full_analysis_set_plus_decoy_hla.fa

if [[ ! -f $qvcf ]]; then
wget -qO - https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/NISTv4.2.1/GRCh38/HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz \
  | bcftools norm -m - -f $fasta -O v \
  | awk 'BEGIN{FS=OFS="\t"} ($5 != "*")' \
  | bcftools view -O z -o $qvcf --threads 3

bcftools index --force $qvcf --threads 3
bcftools view -O b -o $qbcf --threads 3 $qvcf
bcftools index --force $qbcf --threads 3

fi

if [[ ! -f gnomad.slivar.v3.1.2.zip ]]; then

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

fi

#<<SLIVAR
/usr/bin/time -v slivar expr \
    -v $qbcf \
    -g gnomad.slivar.v3.1.2.zip -o slivar-anno.bcf 2> slivar-time.txt
ssize=$(du -cm gnomad.slivar.v3.1.2.zip | grep total)
echo "#SIZE: $ssize" >> slivar-time.txt

#SLIVAR



#<<ECHTVAR
/usr/bin/time -v \
  echtvar anno \
  -e /data/gnomad/gnomad.v3.1.2.echtvar.v2.zip \
  $qbcf annotated.bcf 2> echtvar-time.txt

esize=$(du -cm /data/gnomad/gnomad.v3.1.2.echtvar.v2.zip | grep total)
echo "#SIZE: $esize" >> echtvar-time.txt

#ECHTVAR


## VarNote prep

<<DONE
varnote Index \
  -I $avcf -O varnote-index 
DONE

echo "
@gnomAD
fields=[INFO]
info_fields=[AC, AN, nhomalt, AF, AC_popmax, AN_popmax, nhomalt_popmax, AF_popmax, AF_controls_and_biobanks, nhomalt_controls_and_biobanks]
out_names=[AC:gnomad_ac, AN:gnomad_an, nhomalt:gnomad_nhomalt, AF:gnomad_af, AC_popmax:gnomad_popmax_ac, AN_popmax:gnomad_popmax_an, nhomalt_popmax:gnomad_popmax_homalt, AF_popmax:gnomad_popmax_af, AF_controls_and_biobanks:gnomad_controls_and_biobanks_af, nhomalt_controls_and_biobanks:gnomad_controls_and_biobanks_nhomalt]
" > annoc

#<<VARNOTE

export TMPDIR=/tmp/
mkdir -p null/temp # varnote needs this but won't create it.

/usr/bin/time -v \
varnote Annotation \
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

vsize=$(du -cm ${avcf}* | grep total)
echo "#SIZE: $vsize" >> varnote-time.txt
#VARNOTE


#<<BCFTOOLS
/usr/bin/time -v \
bcftools annotate \
  -a $abcf \
  --threads 4 \
  --columns "INFO/gnomad_ac:=INFO/AC,INFO/gnomad_an:=INFO/AN,INFO/gnomad_nhomalt:=INFO/nhomalt,INFO/gnomad_af:=INFO/AF,INFO/gnomad_popmax_ac:=INFO/AC_popmax,INFO/gnomad_popmax_an:=INFO/AN_popmax,INFO/gnomad_popmax_nhomalt:=INFO/nhomalt_popmax,INFO/gnomad_popmax_af:=INFO/AF_popmax,INFO/gnomad_controls_and_biobanks_af:=INFO/AF_controls_and_biobanks,INFO/gnomad_controls_and_biobanks_nhomalt:=INFO/nhomalt_controls_and_biobanks" \
  -O b \
  -o bcft.anno.bcf \
  $qbcf 2> bcftools-time.txt
bsize=$(du -cm ${bvcf} | grep total)
echo "#SIZE: $bsize" >> bcftools-time.txt

#BCFTOOLS

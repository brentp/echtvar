# make *lot* of variants
python3 make-vcf.py

# encode some of them into an echtvar file
../target/release/echtvar encode generated-exclude.vcf test.echtvar test.hjson
../target/release/echtvar anno generated-all.vcf  -e test.echtvar anno.vcf.gz

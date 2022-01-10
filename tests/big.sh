# make *lot* of variants

for mod in 1 2 3 4 5; do
	echo "mod: $mod"
	python3 make-vcf.py $mod

	# encode some of them into an echtvar file
	../target/release/echtvar encode generated-exclude.vcf test.echtvar test.hjson
	../target/release/echtvar anno generated-all.vcf  -e test.echtvar anno.vcf.gz

	# check that some variants remain unannotated
	python3 check.py anno.vcf.gz $mod

	if [[ "$mod" -eq "1" ]]; then
		echo "all check"
	# check that all variants are annotated.
	../target/release/echtvar encode generated-all.vcf test.echtvar test.hjson
	../target/release/echtvar anno generated-all.vcf  -e test.echtvar anno.vcf.gz
	python3 check.py anno.vcf.gz 1
	fi


done

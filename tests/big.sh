# make *lot* of variants

target=x86_64-unknown-linux-gnu
cargo build --target $target
echtvar=../target/$target/debug/echtvar

set -e
for mod in 1 2 3 4 5; do
	echo "mod: $mod"
	python3 make-vcf.py $mod

	# encode some of them into an echtvar file
	$echtvar encode test.echtvar test.hjson generated-exclude.vcf
	$echtvar anno generated-all.vcf  -e test.echtvar anno.vcf.gz

	# check that some variants remain unannotated
	python3 check.py anno.vcf.gz $mod

	if [[ "$mod" -eq "1" ]]; then
		echo "all check"
		# check that all variants are annotated.
		$echtvar encode test.echtvar test.hjson generated-all.vcf
		$echtvar anno generated-all.vcf  -e test.echtvar anno.vcf.gz
		python3 check.py anno.vcf.gz 1

		echo "check prefixes"
		# check 'chr' prefix mixing
		# 1. database was made with chr prefix and anno without.
		sed -e 's/^chr//' generated-all.vcf > generated-all.no-prefix.vcf
		$echtvar anno generated-all.no-prefix.vcf  -e test.echtvar anno.vcf.gz
		python3 check.py anno.vcf.gz 1

		# 2. database made without chr prefix and anno with.
		$echtvar encode test.echtvar test.hjson generated-all.no-prefix.vcf
		$echtvar anno generated-all.vcf  -e test.echtvar anno.vcf.gz
		python3 check.py anno.vcf.gz 1

                
                cat generated-all.no-prefix.vcf | perl -pe 's/^(chr1|1)|chr1/chr21/' > chr21.vcf
                $echtvar anno chr21.vcf -e test.echtvar anno.vcf.gz
                n=$(zgrep -v ^# anno.vcf.gz | grep -cv 'aval=-1')
                if [[ "$n" -ne "0" ]]; then
                    echo "annotated variants from another chrom!!!"
                    exit 1;
		fi

                rm chr21.vcf
		rm generated-all.no-prefix.vcf
	fi


done
rm generated-all.vcf generated-exclude.vcf anno.vcf.gz
echo "SUCCESS"

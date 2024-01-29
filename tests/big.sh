# make *lot* of variants

target=x86_64-unknown-linux-gnu
cargo build --target $target
echtvar=../target/$target/debug/echtvar

set -e
for mod in 2 3 4 5; do
	echo "mod: $mod"
	python3 make-vcf.py $mod

	# encode some of them into an echtvar file
	$echtvar encode test.echtvar0 test0.hjson generated-subset0.vcf

	$echtvar anno generated-all.vcf  -e test.echtvar0 anno.vcf.gz

	# check that some variants remain unannotated
	python3 check.py anno.vcf.gz $mod
	# check custom INFO Description used from config
	python3 check-string-for-issue33.py anno.vcf.gz aval "added by echtvar TEST description field"


        if [[ "mod" -ne "1" ]]; then
            # now annotate with both and all values shouldl be annotated.
            echo "two echtvars check"
	    $echtvar encode test.echtvar1 test1.hjson generated-subset1.vcf
	    $echtvar anno generated-all.vcf -e test.echtvar0 -e test.echtvar1 anno.vcf.gz
  	    python3 check.py anno.vcf.gz 1
		# check default Description used
		python3 check-string-for-issue33.py anno.vcf.gz aval1 "added by echtvar from test.echtvar1"
        fi


	if [[ "$mod" -eq "1" ]]; then
		echo "all check"
		# check that all variants are annotated.
		$echtvar encode test.echtvar0 test.hjson generated-all.vcf
		$echtvar anno generated-all.vcf  -e test.echtvar0 anno.vcf.gz
		python3 check.py anno.vcf.gz 1

		echo "check prefixes"
		# check 'chr' prefix mixing
		# 1. database was made with chr prefix and anno without.
		sed -e 's/^chr//' generated-all.vcf > generated-all.no-prefix.vcf
		$echtvar anno generated-all.no-prefix.vcf  -e test.echtvar0 anno.vcf.gz
		python3 check.py anno.vcf.gz 1

		# 2. database made without chr prefix and anno with.
		$echtvar encode test.echtvar0 test.hjson generated-all.no-prefix.vcf
		$echtvar anno generated-all.vcf  -e test.echtvar0 anno.vcf.gz
		python3 check.py anno.vcf.gz 1

                
                cat generated-all.no-prefix.vcf | perl -pe 's/^(chr1|1)|chr1/chr21/' > chr21.vcf
                $echtvar anno chr21.vcf -e test.echtvar0 anno.vcf.gz
                n=$(zgrep -v ^# anno.vcf.gz | grep -cv 'aval=-1')
                if [[ "$n" -ne "0" ]]; then
                    echo "annotated variants from another chrom!!!"
                    exit 1;
		fi

                rm chr21.vcf
		rm generated-all.no-prefix.vcf
	fi


done
rm generated-all.vcf generated-subset0.vcf anno.vcf.gz test.echtvar0 test.echtvar1
bash string.sh
echo "SUCCESS"


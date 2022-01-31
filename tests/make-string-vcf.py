import itertools
import random
import sys


#random.seed(42)

fh = open("string.vcf", "w")

header = ("""##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##FILTER=<ID=FAIL,Description="All filters passed">
##FILTER=<ID=OTHER,Description="All filters passed">
##INFO=<ID=num,Number=1,Type=Integer,Description="random string value">
##INFO=<ID=val,Number=.,Type=String,Description="random string value">
##contig=<ID=chr1,length=248956422>
##contig=<ID=1,length=248956422>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO""")

FILTERS = ["PASS", "FAIL", "OTHER"]

print(header, file=fh)

for switch in [1, 2, 3, 4, 5, 1132, 1133, 1134]:
    switch = switch<<20

    for i in range(switch - 12, switch + 12):
        for rlen in range(1, 3):
            for ref in itertools.permutations("ACGT", rlen):
                ref = "".join(ref)
                for alen in range(0, 3):
                    for balt in itertools.permutations("ACGT", alen):
                        val = random.randint(0, 100)
                        flt = random.choice(FILTERS)
                        alt = ref[0] + "".join(balt)
                        print(f"chr1\t{i}\t.\t{ref}\t{alt}\t1\t{flt}\tval=s{val};num=3", file=fh)

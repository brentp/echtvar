import itertools
import random
import sys

mod = int(sys.argv[1])

random.seed(42)

all_fh = open("generated-all.vcf", "w")
exclude_fh = open("generated-exclude.vcf", "w")

header = ("""##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed",IDX=0>
##INFO=<ID=AC,Number=A,Type=Integer,Description="Alternate allele count">
##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles">
##INFO=<ID=val,Number=1,Type=Integer,Description="random value">
##INFO=<ID=nvar,Number=1,Type=Integer,Description="variant index">
##INFO=<ID=AF,Number=A,Type=Float,Description="Alternate allele frequency">
##contig=<ID=chr1,length=248956422>
##contig=<ID=1,length=248956422>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO""")

print(header, file=all_fh)
print(header, file=exclude_fh)

nvar = 0

for switch in [1, 2, 3, 4, 5, 1132, 1133, 1134]:
    switch = switch<<20

    for i in range(switch - 32, switch + 32):
        for rlen in range(1, 5):
            for ref in itertools.permutations("ACGT", rlen):
                ref = "".join(ref)
                for alen in range(0, 5):
                    for balt in itertools.permutations("ACGT", alen):
                        val = random.randint(0, 10000000)
                        alt = ref[0] + "".join(balt)
                        print(f"chr1\t{i}\t.\t{ref}\t{alt}\t1\tPASS\tval={val};nvar={nvar}", file=all_fh)

                        if nvar % mod == 0:
                            print(f"chr1\t{i}\t.\t{ref}\t{alt}\t1\tPASS\tval={val};nvar={nvar}", file=exclude_fh)
                        nvar += 1

        for ref in ["ACCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC",
                   "ACCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCT"]:
            for alt in ["ACCCCCCCCCCCCCCCCC", "A", "ACCCCCCCCCCCCCCC"]:
                print(f"chr1\t{i}\t.\t{ref}\t{alt}\t1\tPASS\tval={val};nvar={nvar}", file=exclude_fh)
                print(f"chr1\t{i}\t.\t{ref}\t{alt}\t1\tPASS\tval={val};nvar={nvar}", file=all_fh)
                nvar += 1

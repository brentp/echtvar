import itertools
import random
import sys


#random.seed(42)

fhdb = open("string-issue-8.db.vcf", "w")
fhq = open("string-issue-8.query.vcf", "w")

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

for fh in (fhdb, fhq):
    print(header, file=fh)

for switch in [1]:
    switch = switch<<20

    ref = "A"
    alt = "C"
    i = 12345

    val = random.randint(0, 100)
    flt = random.choice(FILTERS)
    print(f"chr1\t{i}\t.\t{ref}\t{alt}\t1\t{flt}\tval=s{val};num=3", file=fhdb)
    alt = "T"
    print(f"chr1\t{i}\t.\t{ref}\t{alt}\t1\t{flt}\tval=s{val};num=3", file=fhq)


with open("string-issue-8.json", "w") as fh:
    fh.write("""[{"field": "FILTER", "alias":"test_filter", "missing_string": "OHNO"}]""")

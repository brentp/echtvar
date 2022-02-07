import sys
import gzip


for line in gzip.open(sys.argv[1], 'rt'):
    if line[0] == '#': continue

    toks = line.strip().split("\t")
    assert ";test_filter=OHNO" in toks[7]

import sys
import gzip

fh = gzip.open(sys.argv[1], 'rt')

mod = int(sys.argv[2])

i = -1
for line in fh:
    if line[0] == '#': continue
    i += 1
    if i % mod != 0: continue
    info = line.split("\t")[7].strip().split(';')
    val = [x for x in info if x.startswith('val=')][0]
    aval = [x for x in info if x.startswith('aval=')][0]
    assert val == aval[1:], (val, aval[1:], i)

import sys
import gzip


lines = 0
for line in gzip.open(sys.argv[1], 'rt'):
    lines += 1
    if line[0] == '#': continue

    toks = line.strip().split("\t")
    filt = toks[6]
    info = toks[7].split(';')
    assert f'filter={filt}' in info
    # get val== and string_val= and then make them the same.
    vs = [x.replace('string_', '') for x in info if 'val=' in x]
    assert vs[0] == vs[1], vs

assert lines > 0

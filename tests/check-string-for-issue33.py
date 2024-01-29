"""
Small test to ensure expected INFO fields generated
"""
import sys
import gzip

fpath = sys.argv[1]
field = sys.argv[2]
desc = sys.argv[3]

ok = False
for line in gzip.open(sys.argv[1], 'rt'):
    if line.startswith("##INFO=<ID={},".format(field)):
        toks = line.rstrip("\n").split(",")
        to_check = "Description=\"{}\"".format(desc)
        assert toks[3].startswith(to_check), f"Expected Desciption=\"{desc}\" for {field}, got {to_check} in {toks[3]}"
        ok = True
        break
    # only want to look through header
    elif line.startswith("#CHROM"):
        break
if not ok:
    print(f"Expected field {field} not found in INFO", file=sys.stderr)
    sys.exit(1)
else:
    sys.exit(0)
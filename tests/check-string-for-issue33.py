"""
Small test to ensure expected INFO fields generated
"""
import sys
import gzip

fpath = sys.argv[1]
field = sys.argv[2]
num = sys.argv[3]
desc = sys.argv[4]

ok = False
for line in gzip.open(sys.argv[1], 'rt'):
    if line.startswith("##INFO=<ID={},".format(field)):
        toks = line.rstrip("\n").split(",")
        num_to_check = "Number={}".format(num)
        assert toks[1].startswith(num_to_check), f"Expected Desciption=\"{num}\" for {field}, got {num_to_check} in {toks[1]}"
        desc_to_check = "Description=\"{}\"".format(desc)
        assert toks[3].startswith(desc_to_check), f"Expected Desciption=\"{desc}\" for {field}, got {desc_to_check} in {toks[3]}"
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
import sys
from cyvcf2 import VCF

keys = ["phred", "raw"]
check_keys = ["cadd_" + k for k in keys]

for v in VCF(sys.argv[1]):

    for i, k in enumerate(keys):
        a = v.INFO.get(k)
        b = v.INFO.get(check_keys[i])
        assert abs(a - b) < 0.011, (a, b, v, k, check_keys[i])


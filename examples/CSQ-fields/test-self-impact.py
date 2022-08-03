import cyvcf2
import sys

#" chr1    11876   .       T       G       100     PASS gene=DDX11L1;Consequence=non_coding;impact=non_coding

for v in cyvcf2.VCF(sys.argv[1]):
    a = v.INFO.get("Consequence", "not_found")
    b = v.INFO["impact"]
    assert a == b

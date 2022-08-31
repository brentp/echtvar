import sys
from pyfaidx import Fasta

print("""##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
""")

fa = Fasta(sys.argv[1], as_raw=True, read_ahead=1000)

alts = {'A': "CGT", 'C': 'AGT', 'G': 'ACT', 'T': 'ACG'}
for chrom in fa.keys():
    if chrom.startswith(("HLA")) or chrom.endswith(("alt", "EBV")) or "Un" in chrom or "random" in chrom:
        continue
    print(f'##contig=<ID={chrom}>')
print('#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO')
for chrom in fa.keys():
    if chrom.startswith(("HLA")) or chrom.endswith(("alt", "EBV")) or "Un" in chrom or "random" in chrom:
        continue
    seq = fa[chrom][:]
    #print(chrom, len(seq), type(seq))
    for i, c in enumerate(seq, start=1):
        if c == 'N': continue
        for alt in alts.get(c, []):
            print(f'{chrom}\t{i}\t.\t{c}\t{alt}\t100\tPASS\t.')

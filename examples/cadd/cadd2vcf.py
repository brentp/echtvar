from __future__ import print_function
from xopen import xopen
import sys

def main(precision, path):
    header = None


    hdr = """\
##fileformat=VCFv4.2
##INFO=<ID=raw,Number=1,Type=Float,Description="raw cadd score">
##INFO=<ID=phred,Number=1,Type=Float,Description="phred-scaled cadd score">
{contigs}
##CADDCOMMENT=<ID=comment,comment="{comment}">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"""

    contigs = []
    for i in range(1, 23):
        contigs.append(f"##contig=<ID={i}>")
    contigs.append('##contig=<ID=X>')
    contigs.append('##contig=<ID=Y>')


    for i, line in enumerate(xopen(path, mode='rt', threads=1)):
        if i == 0:
            print(hdr.format(comment=line.strip("# ").strip(), contigs='\n'.join(contigs)))
            continue
        elif i % 10_000_000 == 0:
            print(f'on variant: {i}: {line[0:20]}', file=sys.stderr)
        if line[0] == '#': continue
        #Chrom  Pos     Ref     Alt     RawScore        PHRED

        chrom, pos, ref, alt, raw, phred = line.split("\t", 5)
        raw = float(raw)
        phred = float(phred[:-1])
        print(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t1\tPASS\traw={raw:.2f};phred={phred:.2f}")


if __name__ == "__main__":

    import argparse
    p = argparse.ArgumentParser()
    p.add_argument("--precision", type=int, default=2, help="amount of precision to keep")
    p.add_argument("path")
    a = p.parse_args()

    main(a.precision, a.path)


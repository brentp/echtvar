"""
see: https://github.com/brentp/echtvar/wiki/dbscSNV-example
"""
import sys
import argparse

import zipfile

position_lookup = {
        "hg19": (0, 1),
        "hg38": (4, 5),
}

vcf_header = """##fileformat=VCFv4.3
##source=echtvar-dbscSNV
##dbscSNV=https://sites.google.com/site/jpopgen/dbNSFP
##reference=%s
"""

def get_float(score):
    if score == ".":
        return "."
    else:
        return "%.4g" % float(score)

def main(path, genome_build, out_fh=sys.stdout):
    z = zipfile.ZipFile(path, mode='r')

    print(vcf_header % genome_build, end="", file=out_fh)
    for c in list(range(1, 23)) + ["X", "Y"]:
        print(f'##contig=<ID=chr{c}>', file=out_fh)
    fields = ["ada_score", "rf_score"]
    for field in fields:
      print(f'##INFO=<ID={field.replace("-", "_").upper()},Number=1,Type=Float,Description="splice site prediction. Score >=0.6; SNV likely affects splicing">', file=out_fh)
    print("#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO", file=out_fh)

    chrom_col, pos_col = position_lookup[genome_build]
    chr_files = [f for f in z.namelist()]
    for c in list(range(1, 23)) + ["X", "Y"]:
        with z.open(f"dbscSNV1.1.chr{c}", 'r') as f:
            for bline in f:
                line = bline.decode('UTF-8')
                if not line.startswith("chr"):
                    arr = line.strip().split('\t')
                    chrom = arr[chrom_col]
                    vpos = arr[pos_col]
                    if not chrom == "." and not vpos == ".": 
                        if not "_alt" in chrom and not "_random" in chrom:
                            info = []
                            info.append(f'ADA_SCORE={get_float(arr[16])}')
                            info.append(f'RF_SCORE={get_float(arr[17])}')
                            
                            print(f'chr{chrom}\t{vpos}\t.\t{arr[2]}\t{arr[3]}\t32\tPASS\t{";".join(info)}', file=out_fh)

if __name__ == "__main__":

    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    p.add_argument("-g", "--genome-build", default="hg38", choices=("hg38",
        "hg19"))
    p.add_argument("zip", help="dbscSNV zip file")

    a = p.parse_args()
    main(a.zip, a.genome_build)

r"""
Example:
    First create json file:

$ python dbnsfp.py dbNSFP4.3a.zip -f SIFT4G_converted_rankscore -f SIFT_score \
   -f Polyphen2_HDIV_score -f Polyphen2_HDIV_pred --json \
   > conf.json

Then use that conf file (with same command except without --json) and pipe results  to echtvar encode.

$ python dbnsfp.py \
        dbNSFP4.3a.zip \
        -f SIFT4G_converted_rankscore \
        -f SIFT_score \
        -f Polyphen2_HDIV_score \
        -f Polyphen2_HDIV_pred \
        | echtvar encode dbsnfp.zip conf.json -

"""
import sys
import argparse
import re

import zipfile
import gzip

position_lookup = {
        "hg18": ("hg18_chr", "hg18_pos(1-based)"),
        "hg19": ("hg19_chr", "hg19_pos(1-based)"),
        "hg38": ("#chr", "pos(1-based)"),
}

lower_is_more_damaging = [
        "SIFT_score",
        "SIFT4G_score",
        "FATHHMM_score",
        "PROVEAN_score",
]

vcf_header = """##fileformat=VCFv4.3
##source=echtvar-dbnsfp
##dbNSFP=https://sites.google.com/site/jpopgen/dbNSFP
##reference=%s
"""

def get_field_lookup(z):
    readme_path = [f for f in z.namelist() if "readme" in f.lower() and not f.endswith("pdf")][0]
    fh = z.open(readme_path, 'r')
    txt = fh.read().decode()
    lookup = {}
    for m in re.finditer(r'(^\d+\t)(.*?)(?=^\d)', txt, re.MULTILINE | re.DOTALL):
      col, desc = m.groups(1)[1].replace('\n', ' ').replace('\r', '').replace('"', "'").replace('\t', '').split(':', 1)
      lookup[col] = desc.strip()
    return lookup

def getfloat(f, reducer):
    if f == ".": return '.'
    r = [x for x in f.split(';') if x != '.']
    if len(r) == 0: return '.'
    return "%.4g" % reducer(map(float, r))

def getstring(f):
    if f == '.': return '.'
    r = [x for x in f.split(';') if x != '.']
    if len(r) == 0: return '.'
    return ",".join(r)

def write_json(field_lookup, fields):
    d = []
    for field in fields:
      field = field.replace("-", "_").replace('+', 'p') # gerp++
      f = {"field": field, "alias": "dbsnfp_" + field}
      if not field.endswith("pred"):
        f["multiplier"] = 1000000
      if field in lower_is_more_damaging:
        f["missing_value"] = 999
      d.append(f)
    import json
    print(json.dumps(d, indent=4))

def main(path, genome_build, fields, json, out_fh=sys.stdout):
    z = zipfile.ZipFile(path, mode='r')

    field_lookup = get_field_lookup(z)
    if json:
      return write_json(field_lookup, fields)


    print(vcf_header % genome_build, end="", file=out_fh)
    for c in list(range(1, 23)) + ["X", "Y", "M"]:
        print(f'##contig=<ID={c}>', file=out_fh)
    header = None
    for field in fields:
      if not field in field_lookup:
         raise newException(OsError, f"field {field} not found. see readme in zip file for list of available fields")
      print(f'##INFO=<ID={field.replace("-", "_")},Number=1,Type={"String" if field.endswith("pred") else "Float"},Description="{field_lookup[field]}">', file=out_fh)
    print("#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO", file=out_fh)

    chrom_col, pos_col = position_lookup[genome_build]
    for chrom in list(range(1, 23)) + ["X", "Y", "M"]:
        fh = gzip.open(z.open(f'dbNSFP4.3a_variant.chr{chrom}.gz'), mode='rt')
        for i, line in enumerate(fh):
            toks = line.rstrip().split("\t")
            if i == 0:
                if header is None: header = toks
                assert header == toks
                continue
            d = dict(zip(header, toks))
            info = []
            for field in fields:
                if field.endswith("pred"):
                    info.append(f'{field.replace("-", "_")}={getstring(d[field])}')
                else:
                    info.append(f'{field.replace("-", "_")}={getfloat(d[field], min if field in lower_is_more_damaging else max)}')

            print(f'{d[chrom_col]}\t{d[pos_col]}\t.\t{d["ref"]}\t{d["alt"]}\t32\tPASS\t{";".join(info)}', file=out_fh)

if __name__ == "__main__":

    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    p.add_argument("-g", "--genome-build", default="hg38", choices=("hg38",
        "hg19", "hg18"))
    p.add_argument("-f", "--field", action="append", help="dbNSFP fields to extract. can be specified multiple times")
    p.add_argument("--json", action="store_true", help="output only a JSON stub for echtvar and don't convert to VCF (run this first)")
    p.add_argument("zip", help="dbnsfp zip file")

    a = p.parse_args()
    main(a.zip, a.genome_build, a.field, a.json)

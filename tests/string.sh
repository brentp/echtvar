
set -euo pipefail

python make-string-vcf.py
target=x86_64-unknown-linux-gnu
cargo build --target $target
echtvar=../target/$target/debug/echtvar

$echtvar encode string.echtvar string.hjson string.vcf
$echtvar anno -i "anno_num == 3" -e string.echtvar string.vcf string-anno.vcf.gz

python check-string.py string-anno.vcf.gz


# this outputs no variants
$echtvar anno -i "anno_num != 3" -e string.echtvar string.vcf string-anno.vcf.gz

rm -f string.vcf  string.echtvar # string-anno.vcf.gz

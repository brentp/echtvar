This example shows how one could create echtvar archives that will indicate the
gene and impact of every possible SNP in a reference genome.

The steps are:

1. Generate every possible SNP in a reference genome (3 possible alternate
   alleles per base).
2. Annotate with `bcftools csq` to get consequences and genes.
3. Use the `split-vep` plugin to get the individual CSQ fields into their own
   INFO fields.
4. encode the impact (missense/stop-gained, etc) into an echtvar archive
5. encode the gene name into an echtvar archive

For the more than 8 billion possible SNPs genome-wide, these archives compress
to 1.4GB!

We then have an additional step using `test-self-impact.py` to test that the self-annotation works. This annotation occurs at more than 500,000 variants per second on our system.

This could be used to annotate all SNPs and filter to higher impact variants,
thereby greatly reducing the computational load of annotation with consequence.
Since a smaller percentage of variants in indels, those could be sent
separately to `bcftools csq` or `vep` for annotation and still maintain a
performance increase.

#CEU Trio for Varify
## Intro
This is a Snakemake workflow using the GATK Haplotype Caller and SnpEff-CBMi to generate demo data for Varify.

The BAM files are directly from the ubiquitous CEU(Utah) trio.

* ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/pilot2_high_cov_GRCh37_bams/data/NA12878/alignment/NA12878.chrom22.ILLUMINA.bwa.CEU.high_coverage.20100311.bam
* ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/pilot2_high_cov_GRCh37_bams/data/NA12891/alignment/NA12891.chrom22.ILLUMINA.bwa.CEU.high_coverage.20100517.bam
* ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/pilot2_high_cov_GRCh37_bams/data/NA12892/alignment/NA12892.chrom22.ILLUMINA.bwa.CEU.high_coverage.20100517.bam

## Usage

```
cd /nas/is1/leipzig/CEU_trio_for_varify
. bin/activate
snakemake -j 16 vcfs/all.vcf
```

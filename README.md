#CEU Trio for Varify
## Intro
This is a Snakemake workflow using the GATK Haplotype Caller and SnpEff-CBMi to generate demo data for Varify.

The BAM files are directly from the ubiquitous CEU(Utah) trio.

* ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/pilot2_high_cov_GRCh37_bams/data/NA12878/alignment/NA12878.chrom20.ILLUMINA.bwa.CEU.high_coverage.20100311.bam
* ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/pilot2_high_cov_GRCh37_bams/data/NA12891/alignment/NA12891.chrom20.ILLUMINA.bwa.CEU.high_coverage.20100517.bam
* ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/pilot2_high_cov_GRCh37_bams/data/NA12892/alignment/NA12892.chrom20.ILLUMINA.bwa.CEU.high_coverage.20100517.bam

## Resources
* ftp://gsapubftp-anonymous@ftp.broadinstitute.org//bundle/2.8/b37/hapmap_3.3.b37.vcf.gz &
* ftp://gsapubftp-anonymous@ftp.broadinstitute.org//bundle/2.8/b37/1000G_phase1.snps.high_confidence.b37.vcf.gz
* ftp://gsapubftp-anonymous@ftp.broadinstitute.org//bundle/2.8/b37/1000G_omni2.5.b37.vcf.gz
* ftp://gsapubftp-anonymous@ftp.broadinstitute.org//bundle/2.8/b37/Mills_and_1000G_gold_standard.indels.b37.vcf.gz

## Usage
`config.json` should be altered as needed.

`filter`: can be set to `hard` or `vqsr`

```
cd /nas/is1/leipzig/CEU_trio_for_varify
mkdir vcfs gvcfs log bams gatk_resources
. bin/activate
snakemake -j 16 vcfs/all.vcf
```

import json

with open("config.json") as f:
    CONFIG = json.load(f)

rule get_bams:
     input: ["bams/"+f+".bam" for f in CONFIG["samples"]]

rule get_a_bam:
     params: ncbi_ftp = 'ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp', subdir = 'pilot3_exon_targetted_GRCh37_bams', type='exon_targetted', chrom='chrom20', date='20100311'
     output: bam="bams/{bamfile}.bam",bai="bams/{bamfile}.bam.bai"
     shell: """
            curl {params.ncbi_ftp}/technical/{params.subdir}/data/{wildcards.bamfile}/alignment/{wildcards.bamfile}.{params.chrom}.ILLUMINA.bwa.CEU.{params.type}.{params.date}.bam > {output.bam}
	    curl {params.ncbi_ftp}/technical/{params.subdir}/data/{wildcards.bamfile}/alignment/{wildcards.bamfile}.{params.chrom}.ILLUMINA.bwa.CEU.{params.type}.{params.date}.bam.bai > {output.bai}
            """

#     shell: """
#        curl ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/pilot2_high_cov_GRCh37_bams/data/NA12878/alignment/NA12878.chrom20.ILLUMINA.bwa.CEU.high_coverage.20100311.bam > bams/NA12878.bam
#	"""

#exon-targeted bams (~173MB)
#ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/pilot3_exon_targetted_GRCh37_bams/data/NA12891/alignment/NA12891.chrom22.ILLUMINA.bwa.CEU.exon_targetted.20100311.bam

rule clean:
     shell: """
     rm vcfs/*
     rm gvcfs/*
     """

def _gatk_multi_arg(flag, files):
    flag += " "
    return flag + flag.join(files)


rule gatk_haplotype_caller:
    input:
        bams="bams/{sample}.bam"
    output:
        gvcf="gvcfs/{sample}.gvcf",
        idx="gvcfs/{sample}.gvcf.idx"
    params:
        gatk_path=CONFIG.get("gatk_path", ""),
        custom=CONFIG.get("params_gatk", ""),
        ref=CONFIG.get("references").get("genome"),
        range=CONFIG.get("range")
    log:
        "log/{sample}.genotype_info.log"
    threads: 8
    run:
        bams = _gatk_multi_arg("-I", input.bams)
        shell(
            "{params.gatk_path} -T HaplotypeCaller -R {params.ref} -I {input.bams} {params.custom} "
	    "-L {params.range} "
            "--emitRefConfidence GVCF --variant_index_type LINEAR "
            "--heterozygosity {CONFIG[heterozygosity]} "
            "--indel_heterozygosity {CONFIG[indel_heterozygosity]} "
            "--dbsnp {CONFIG[known_variants][dbsnp]} -nct {threads} "
            "--variant_index_parameter 128000 -o {output.gvcf} >& {log}")


rule gatk_genotyping:
    input:
        gvcfs=expand(
            "gvcfs/{sample}.gvcf",
            sample=CONFIG["samples"])
    output:
        "vcfs/all.vcf"
    params:
        gatk_path=CONFIG.get("gatk_path", ""),
        ref=CONFIG.get("references").get("genome"),
        custom=CONFIG.get("params_gatk", "")
    log:
        "log/all.genotype.log"
    threads: 8
    run:
        gvcfs = _gatk_multi_arg(" --variant", input.gvcfs)
        shell(
            "{params.gatk_path} -T GenotypeGVCFs {gvcfs} -nt {threads} {params.custom} "
            "-R {params.ref} "
            "--dbsnp {CONFIG[known_variants][dbsnp]} -o {output} >& {log}")

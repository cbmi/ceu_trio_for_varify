import json
import hashlib

with open("config.json") as f:
    CONFIG = json.load(f)

rule get_bams:
    input:
        ["bams/" + f + ".bam" for f in CONFIG["samples"]]

# exon-targeted bams (~173MB)
# ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/pilot3_exon_targetted_GRCh37_bams/data/NA12891/alignment/NA12891.chrom22.ILLUMINA.bwa.CEU.exon_targetted.20100311.bam
rule get_a_bam:
    params:
        ncbi_ftp = 'ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp', subdir = 'pilot3_exon_targetted_GRCh37_bams', type = 'exon_targetted', chrom = 'chrom20', date = '20100311'
    output:
        bam = "bams/{bamfile}.bam", bai = "bams/{bamfile}.bam.bai"
    shell:
        """
            curl {params.ncbi_ftp}/technical/{params.subdir}/data/{wildcards.bamfile}/alignment/{wildcards.bamfile}.{params.chrom}.ILLUMINA.bwa.CEU.{params.type}.{params.date}.bam > {output.bam}
            curl {params.ncbi_ftp}/technical/{params.subdir}/data/{wildcards.bamfile}/alignment/{wildcards.bamfile}.{params.chrom}.ILLUMINA.bwa.CEU.{params.type}.{params.date}.bam.bai > {output.bai}
            """
rule all:
    input:
        "vcfs/all.snpeff.vcf",
        "vcfs/all.phased.vcf",
        "all_varify/all.snpeff.vcf"

rule clean:
    shell:
        """
     rm vcfs/*
     rm gvcfs/*
     """


def _gatk_multi_arg(flag, files):
    flag += " "
    return flag + flag.join(files)

rule get_ranges:
    output: "ranges.list"
    shell:
        """
        ./genesToRanges.py > {output}
        """

# http://bit.ly/1Dc4mXy
rule gatk_haplotype_caller:
    input:
        bams = "bams/{sample}.bam",
        ranges = "ranges.list"
    params:
        java_cmd = CONFIG.get("java_cmd", ""),
        gatk_path = CONFIG.get("gatk_path", ""),
        custom = CONFIG.get("params_gatk", ""),
        ref = CONFIG.get("references").get("genome"),
        range = CONFIG.get("range")
    output:
        gvcf = "gvcfs/{sample}.gvcf",
        idx = "gvcfs/{sample}.gvcf.idx"
    log:
        "log/{sample}.genotype_info.log"
    threads:
        8
    run:
        bams = _gatk_multi_arg("-I", input.bams)
        shell(
            "{params.java_cmd} "
            "{params.gatk_path} "
            "-T HaplotypeCaller "
            "-R {params.ref} "
            "-I {input.bams} "
            "{params.custom} "
            "-L {input.ranges} "
            "--emitRefConfidence GVCF --variant_index_type LINEAR "
            "--heterozygosity {CONFIG[heterozygosity]} "
            "--indel_heterozygosity {CONFIG[indel_heterozygosity]} "
            "--dbsnp {CONFIG[known_variants][dbsnp]} "
            "-nct {threads} "
            "--variant_index_parameter 128000 "
            "-o {output.gvcf} "
            ">& {log}"
        )

# http://bit.ly/1rcXu3u
rule gatk_genotyping:
    input:
        gvcfs = expand(
            "gvcfs/{sample}.gvcf",
            sample=CONFIG["samples"])
    params:
        java_cmd = CONFIG.get("java_cmd", ""),
        gatk_path = CONFIG.get("gatk_path", ""),
        ref = CONFIG.get("references").get("genome"),
        custom = CONFIG.get("params_gatk", "")
    output:
        "vcfs/all.vcf"
    log:
        "log/all.genotype.log"
    threads:
        8
    run:
        gvcfs = _gatk_multi_arg(" --variant", input.gvcfs)
        shell(
            "{params.java_cmd} "
            "{params.gatk_path} "
            "-R {params.ref} "
            "-T GenotypeGVCFs {gvcfs} "
            "-nt {threads} {params.custom} "
            "--dbsnp {CONFIG[known_variants][dbsnp]} "
            "-o {output} "
            ">& {log}"
        )


def _get_recal_params(wildcards):
    if wildcards.type == "snp":
        return (
            "-mode SNP -an DP -an QD -an FS -an MQRankSum -an ReadPosRankSum "
            "-resource:hapmap,known=false,training=true,truth=true,prior=15.0 {hapmap} "
            "-resource:omni,known=false,training=true,truth=true,prior=12.0 {omni} "
            "-resource:1000G,known=false,training=true,truth=false,prior=10.0 {g1k} "
            "-resource:dbsnp,known=true,training=false,truth=false,prior=2.0 {dbsnp}"
        ).format(**CONFIG["known_variants"])
    else:
        return (
            "-mode INDEL -an DP -an FS -an MQRankSum -an ReadPosRankSum "
            "-resource:mills,known=true,training=true,truth=true,prior=12.0 {mills}"
        ).format(**CONFIG["known_variants"])


# hard filtration
# this "filters out, not filters for" filterExpression
rule gatk_hard_filtration:
    input:
        vcf = "vcfs/{filename}.vcf",
        ref = CONFIG.get("references").get("genome")
    params:
        java_cmd = CONFIG.get("java_cmd", ""),
        gatk_path = CONFIG.get("gatk_path", "")
    output:
        "vcfs/{filename}.hard.vcf"
    log:
        "log/{filename}.gatk_hard_filtration.log"
    shell:
        "{params.java_cmd} "
        "{params.gatk_path} "
        "-R {input.ref} "
        "-T VariantFiltration "
        "-o {output} "
        "--variant {input.vcf} "
        "--filterExpression \"QD < 2.0 || MQ < 40.0 || FS > 60.0 || HaplotypeScore > 13.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0\" "
        "--filterName \"GATK3.0-hard-filter\" "
        ">& {log}"

rule select_passing:
    input:
        vcf = "vcfs/{filename}." + CONFIG.get("filter") + ".vcf",
        ref = CONFIG.get("references").get("genome")
    params:
        java_cmd = CONFIG.get("java_cmd", ""),
        gatk_path = CONFIG.get("gatk_path", "")
    output:
        "vcfs/{filename}.filtered.vcf"
    log:
        "log/{filename}.select_passing_variants.log"
    shell:
        "{params.java_cmd} "
        "{params.gatk_path} "
        "-R {input.ref} "
        " -T SelectVariants "
        "-o {output} "
        "--variant {input.vcf} "
        "--excludeFiltered "
        ">& {log}"

# VQSR based filtration
# requires sufficient number of samples and variants YMMV
# http://bit.ly/1w2I72P
rule gatk_variant_recalibration:
    input:
        CONFIG["known_variants"].values(),
        ref = CONFIG.get("references").get("genome"),
        vcf = "vcfs/{filename}.vcf"
    params:
        java_cmd = CONFIG.get("java_cmd", ""),
        gatk_path = CONFIG.get("gatk_path", ""),
        recal = _get_recal_params,
        custom = CONFIG.get("params_gatk", "")
    output:
        recal = temp("vcfs/variant_recal/{filename}.{type,(snp|indel)}.recal"),
        tranches = temp(
            "vcfs/variant_recal/{filename}.{type,(snp|indel)}.tranches"),
        plotting = temp(
            "vcfs/variant_recal/{filename}.{type,(snp|indel)}.plotting.R")
    log:
        "log/{filename}.{type}_recalibrate_info.log"
    threads:
        8
    shell:
        "{params.java_cmd} "
        "{params.gatk_path} "
        "-T VariantRecalibrator "
        "-R {input.ref} "
        "-input {input.vcf} "
        "{params.recal} "
        "{params.custom} "
        "-nt {threads} "
        "-tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 "
        "-recalFile {output.recal} "
        "-tranchesFile {output.tranches} "
        "-rscriptFile {output.plotting} "
        ">& {log}"

# give the recal file a pretty name
rule vqsr:
    input:
        "vcfs/{filename}.snp_recalibrated.indel_recalibrated.vcf"
    output:
        "vcfs/{filename}.vqsr.vcf"
    shell:
        "mv {input} {output}"

# this rule is smart enough to accept
# vcfs/all.snp_recalibrated.indel_recalibrated.vcf as a target
rule gatk_apply_variant_recalibration:
    input:
        ref = CONFIG.get("references").get("genome"),
        vcf = "vcfs/{filename}.vcf",
        recal = "vcfs/variant_recal/{filename}.{type}.recal",
        tranches = "vcfs/variant_recal/{filename}.{type}.tranches"
    params:
        java_cmd = CONFIG.get("java_cmd", ""),
        gatk_path = CONFIG.get("gatk_path", ""),
        mode = lambda wildcards: wildcards.type.upper(),
        custom = CONFIG.get("params_gatk", "")
    output:
        "vcfs/{filename}.{type,(snp|indel)}_recalibrated.vcf"
    log:
        "log/{filename}.{type}_recalibrate.log"
    threads:
        8
    shell:
        "{params.java_cmd} "
        "{params.gatk_path} "
        "-T ApplyRecalibration "
        "-R {input.ref} "
        "-nt {threads} "
        "-input {input.vcf} "
        "-mode {params.mode} "
        "{params.custom} "
        "-recalFile {input.recal} "
        "--ts_filter_level 99.9 "
        "-tranchesFile {input.tranches} "
        "-o {output} "
        ">& {log}"

# http://bit.ly/1EYgFZk
rule phase_by_transmission:
    input:
        vcf = "vcfs/{filename}.filtered.vcf",
        ref = CONFIG.get("references").get("genome"),
        ped = CONFIG.get("ped")
    params:
        java_cmd = CONFIG.get("java_cmd", ""),
        gatk_path = CONFIG.get("gatk_path", "")
    output:
        vcf = "vcfs/{filename}.phased.vcf",
        mvf = "vcfs/{filename}_mendelian_violations.txt"
    log:
        "log/{filename}.phase_by_transmission.log"
    shell:
        "{params.java_cmd} "
        "{params.gatk_path} "
        "-R {input.ref} "
        "-T PhaseByTransmission "
        "--variant {input.vcf} "
        "-ped {input.ped} "
        "-mvf {output.mvf} "
        "-o {output.vcf} "
        ">& {log}"

rule annotate_dbsnp:
    input:
        vcf = "vcfs/{filename}.filtered.vcf"
    params:
        java_cmd = CONFIG.get("java_cmd", ""),
        path = CONFIG.get("snpeff").get("path"),
        config = CONFIG.get("snpeff").get("config")
    output:
        vcf = "vcfs/{filename}.snpeff.vcf"
    log:
        "log/{filename}.snpeff.log"
    shell:
        "{params.java_cmd} "
        "{params.path} "
        "-c {params.config} "
        "-t hg19 "
        "-ud 10 "
        "-i vcf "
        "-o vcf {input} "
        "> {output} "
        "2> {log}"

rule varify_manifest:
    input:
        vcf = "vcfs/{filename}.snpeff.vcf"
    output:
        dirname = "{filename}_varify",
        manifest = "{filename}_varify/MANIFEST",
        vcf = "{filename}_varify/{filename}.snpeff.vcf"
    run:
        with open(input.vcf, 'rb') as vcffilehandle:
            md5_string = file_md5(vcffilehandle)
        shell("mkdir -p {0}".format(output.dirname))
        shell("ln -s ../{0} {1}".format(input.vcf, output.vcf))
        with open(output.manifest, 'w') as outfile:
            outfile.write("""
[general]
load = true

[sample]
project = CEU
batch = CEU
version = 0

[vcf]
""")
            outfile.write("file = {0}.snpeff.vcf\n".format(wildcards.filename))
            outfile.write("md5 = {0}\n".format(md5_string))


def file_md5(f, size=8192):
    "Calculates the MD5 of a file. http://stackoverflow.com/a/1131255/264696"
    md5 = hashlib.md5()
    while True:
        data = f.read(size)
        if not data:
            break
        md5.update(data)
    return md5.hexdigest()

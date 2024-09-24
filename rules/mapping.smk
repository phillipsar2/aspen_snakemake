# (1) Evaluate quality of raw reads with fastqc
# scripts/fastqc.sh

# (2) Trim reads sequenced at UCD with fastp
# Minimum length is 36 (-l 36)
# Don't filter for quality (-Q)
# Adapter trimming is enabled by default -- don't need to specify adapter seq
# Default detects and trims polyG tails for NovaSeq data
# Trim first 9 bp on each read (--trim_front1 9 --trim_front2 9)

rule fastp_trim:
    input:
        fastq = "/global/scratch/users/arphillips/raw/jgi_wgs/{sample}_R1_001.fastq.gz",
    output:
        trim = "data/raw/trimmed/{sample}.trim_1.fastq.gz",
        report = "reports/fastp/{sample}.json"
    run:
        shell("fastp -w 2 \
        -l 36 -Q \
        -i {input.fastq} \
        -o {output.trim} \
        --trim_front1 9 --trim_front2 9 \
        -j {output.report}")

# (3a) Prepare reference file
rule bwa_prep:
    input: 
        config["data"]["reference-genome"]
    output:
        index = "/global/scratch/projects/fc_moilab/projects/aspen/genome/mex_genome/genome.1MX.fasta.fai"
    shell:
        """
        ~/toolz/bwa-mem2-2.2.1_x64-linux/bwa-mem2 index {input}
        gatk
        """


# (3b) Align reads to the reference genome
# Paired-end reads in single file
rule bwa_map:
    input:
        ref = config.ref,
        index = "/global/scratch/projects/fc_moilab/projects/aspen/genome/mex_genome/genome.1MX.fasta.fai"
        fastq = "/group/jrigrp10/andropogon_shortreads/{sample}.merge.R2.fastq.gz"
    output:
        temp("data/interm/mapped_bam/{sample}.mapped.bam")
#    log:
#        "logs/bwa_mem/{sample}.log",
    conda: "bwa_map.yaml"
    threads: 8
    shell:
        "~/toolz/bwa-mem2-2.2.1_x64-linux/bwa-mem2 mem -t {threads} {input.ref} {input.r1} {input.r2} |"
        "samtools view -Sb > {output}"

# (4) Sort bams
rule samtools_sort:
    input:
        "data/interm/mapped_bam/{sample}.mapped.bam"
    output:
        temp("data/interm/sorted_bam/{sample}.sorted.bam"),
    params:
        tmp = "/scratch/aphillip/sort_bam/{sample}"
    threads: 8
    run:
        shell("mkdir -p {params.tmp}")
        shell("samtools sort -T {params.tmp} -@ {threads} {input} > {output}")
        shell("rm -rf {params.tmp}")

# (5) Add read groups
rule add_rg:
    input:
        "data/interm/sorted_bam/{sample}.sorted.bam"
    output:
        bam = temp(touch("data/interm/addrg/{sample}.rg.bam"))
    params:
        tmp = "/scratch/aphillip/addrg/{sample}",
        sample = "{sample}"
    run:
        shell("mkdir -p {params.tmp}")
        shell("gatk --java-options ""-Xmx4G"" AddOrReplaceReadGroups \
        -I {input} \
        -O {output.bam} \
        -RGID 4 \
        -RGLB lib1 \
        -RGPL illumina \
        -RGPU unit1 \
        -RGSM {params.sample} \
        --TMP_DIR {params.tmp} \
        --CREATE_INDEX true")
        shell("rm -rf {params.tmp}")

# (6) Mark duplicates
rule mark_dups:
    input:
        "data/interm/addrg/{sample}.rg.bam"
    output:
        bam = "data/interm/mark_dups/{sample}.dedup.bam",
        metrics = "qc/mark_dup/{sample}_metrics.txt"
    params:
        tmp = "/scratch/aphillip/mark_dups/{sample}"
    run:
        # Create a scratch directory
        shell("mkdir -p {params.tmp}")
        # Input bam file to output marked records. Assume bam file has been sorted. Direct to a temporary storage file (scratch).
        shell("gatk --java-options ""-Xmx10G"" MarkDuplicates \
        -I {input} \
        -O {output.bam} \
        --METRICS_FILE {output.metrics} \
        --CREATE_INDEX true \
        -MAX_FILE_HANDLES 1000 \
        --ASSUME_SORT_ORDER coordinate \
        --TMP_DIR {params.tmp}")
        # Remove scratch directory
        shell("rm -rf {params.tmp}")

# (7) Assess alignment quality metrics with qualimap
# nr is normally 100000 and -nt is normally 8, java mem size = 48
# nw is normally 400
# for higher cov, make nr 1000 and -nt 12, java mem size = 64
# tested with nr = 10000 and nw = 400, failed
rule bamqc:
    input:
        "data/interm/mark_dups/{bam}.dedup.bam"
###        "data/interm/mark_dups/I{bam}.dedup.bam"
    output:
        "reports/bamqc/{bam}_stats/qualimapReport.html"
    params:
        dir = "reports/bamqc/{bam}_stats"
    run: 
        shell("qualimap bamqc \
        -bam {input} \
        -nt 12 \
        -nr 1000 \
        -nw 400 \
        -outdir {params.dir} \
        -outformat HTML \
        --skip-duplicated \
        --java-mem-size=20G")

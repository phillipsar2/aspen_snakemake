# (1) Evaluate quality of raw reads with fastqc
# scripts/fastqc.sh

rule fastqc:
    input:
        "/global/scratch/users/arphillips/raw/jgi_wgs/{sample}.fastq.gz"
    output:
        "/global/scratch/users/arphillips/qc/fastqc/{sample}_fastqc.zip"
    params:
        tmp = "/global/scratch/users/arphillips/tmp/fastqc/{sample}",
        outdir = "/global/scratch/users/arphillips/qc/fastqc"
    conda:
        "/global/home/users/arphillips/aspen/aspen_snakemake/envs/fastqc.yaml"
    resources:
        cpus_per_task=2
   # threads: config["params"]["fastqc-threads"]
   #  resources:
    shell:
        """
        mkdir -p {params.tmp}
        fastqc -o {params.outdir} -t {resources.cpus_per_task} -d {params.tmp} -f fastq {input}
        rm -rf {params.tmp}
        """

# (2) Trim reads sequenced at UCD with fastp
# Minimum length is 36 (-l 36)
# Don't filter for quality (-Q)
# Adapter trimming is enabled by default -- don't need to specify adapter seq
# Default detects and trims polyG tails for NovaSeq data
# --cut_front is sliding window trimming from 5' ot 3'
rule fastp_trim:
    input:
        fastq = "/global/scratch/users/arphillips/raw/jgi_wgs/{sample}.fastq.gz",
    output:
        trim = "/global/scratch/users/arphillips/data/trimmed/{sample}.trim.fastq.gz"
    conda:
        "/global/home/users/arphillips/aspen/aspen_snakemake/envs/fastp.yaml"
    shell:
        """
        fastp -w 2 \
        -l 36 -Q \
        -i {input.fastq} \
        -o {output.trim} \
        --cut_front --cut_front_window_size 4 --cut_front_mean_quality 15 
        """

# (3a) Prepare reference file
# 9.6 GB per core available (core == threads)
rule bwa_prep:
    input: 
        config["data"]["reference"]["genome"]
    output:
        index = "/global/scratch/projects/fc_moilab/projects/aspen/genome/mex_genome/genome.1MX.fasta.gz.0123"
    conda: "/global/home/users/arphillips/aspen/aspen_snakemake/envs/bwa-mem2.yaml"
    shell:
        "~/toolz/bwa-mem2-2.2.1_x64-linux/bwa-mem2 index {input}"


# (3b) Align reads to the reference genome
# Paired-end reads in single file
rule bwa_map:
    input:
        ref = config["data"]["reference"]["genome"],
        index = "/global/scratch/projects/fc_moilab/projects/aspen/genome/mex_genome/genome.1MX.fasta.gz.0123",
        trim = "/global/scratch/users/arphillips/data/trimmed/{sample}.trim.fastq.gz"
    output:
        temp("/global/scratch/users/arphillips/data/interm/mapped_bam/{sample}.mapped.bam")
    conda: "/global/home/users/arphillips/aspen/aspen_snakemake/envs/envs/bwa_map.yaml"
    threads: 8
    shell:
        "~/toolz/bwa-mem2-2.2.1_x64-linux/bwa-mem2 mem -t 8 {input.ref} {input.trim} |"
        "samtools view -Sb > {output}"

# (4) Sort bams
rule samtools_sort:
    input:
        "/global/scratch/users/arphillips/data/interm/mapped_bam/{sample}.mapped.bam"
    output:
        temp("/global/scratch/users/arphillips/data/interm/sorted_bam/{sample}.sorted.bam"),
    params:
        tmp = "/global/scratch/users/arphillips/temp/sort_bam/{sample}"
    conda: "/global/home/users/arphillips/aspen/aspen_snakemake/envs/envs/samtools.yaml"
    #threads: 8
    shell:
        "mkdir -p {params.tmp}"
        "samtools sort -T {params.tmp} -@ 8 {input} > {output}"
        "rm -rf {params.tmp}"

# (5) Add read groups
rule add_rg:
    input:
        "/global/scratch/users/arphillips/data/interm/sorted_bam/{sample}.sorted.bam"
    output:
        bam = temp(touch("/global/scratch/users/arphillips/data/interm/addrg/{sample}.rg.bam"))
    params:
        tmp = "/global/scratch/users/arphillips/temp/addrg/{sample}",
        sample = "{sample}"
    conda: "/global/home/users/arphillips/aspen/aspen_snakemake/envs/envs/gatk.yaml"
    shell:
        """
        mkdir -p {params.tmp}
        gatk --java-options ""-Xmx4G"" AddOrReplaceReadGroups \
        -I {input} \
        -O {output.bam} \
        -RGID 4 \
        -RGLB lib1 \
        -RGPL illumina \
        -RGPU unit1 \
        -RGSM {params.sample} \
        --TMP_DIR {params.tmp} \
        --CREATE_INDEX true")
        rm -rf {params.tmp}
        """

# (6) Mark duplicates
rule mark_dups:
    input:
        "/global/scratch/users/arphillips/data/interm/addrg/{sample}.rg.bam"
    output:
        bam = "/global/scratch/users/arphillips/data/interm/mark_dups/{sample}.dedup.bam",
        metrics = "/global/scratch/users/arphillips/qc/mark_dup/{sample}_metrics.txt"
    params:
        tmp = "/global/scratch/users/arphillips/temp/mark_dups/{sample}"
    conda: "/global/home/users/arphillips/aspen/aspen_snakemake/envs/envs/gatk.yaml"
    shell:
        """
        # Create a scratch directory
        mkdir -p {params.tmp}
        # Input bam file to output marked records. Assume bam file has been sorted. Direct to a temporary storage file (scratch).
        gatk --java-options ""-Xmx10G"" MarkDuplicates \
        -I {input} \
        -O {output.bam} \
        --METRICS_FILE {output.metrics} \
        --CREATE_INDEX true \
        -MAX_FILE_HANDLES 1000 \
        --ASSUME_SORT_ORDER coordinate \
        --TMP_DIR {params.tmp}
        # Remove scratch directory
        rm -rf {params.tmp}
        """

# (7) Assess alignment quality metrics with qualimap
# nr is normally 100000 and -nt is normally 8, java mem size = 48
# nw is normally 400
# for higher cov, make nr 1000 and -nt 12, java mem size = 64
rule bamqc:
    input:
        "/global/scratch/users/arphillips/data/interm/mark_dups/{bam}.dedup.bam"
    output:
        "/global/scratch/users/arphillips/reports/bamqc/{bam}_stats/qualimapReport.html"
    params:
        dir = "/global/scratch/users/arphillips/reports/bamqc/{bam}_stats"
    shell:
        """
        qualimap bamqc \
        -bam {input} \
        -nt 12 \
        -nr 1000 \
        -nw 400 \
        -outdir {params.dir} \
        -outformat PDF \
        --skip-duplicated \
        --java-mem-size=64G
        """

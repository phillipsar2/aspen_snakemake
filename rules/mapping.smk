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
        "/global/scratch/projects/fc_moilab/aphillips/aspen_snakemake/envs/fastqc.yaml"
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
        trim = temp("/global/scratch/users/arphillips/data/trimmed/{sample}.trim.fastq.gz")
    conda:
        "/global/scratch/projects/fc_moilab/aphillips/aspen_snakemake/envs/fastp.yaml"
#    benchmark:
#        "/global/scratch/users/arphillips/benchmarks/{sample}.trim.benchmark.txt"
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
        index = "/global/scratch/projects/fc_moilab/projects/aspen/genome/CAM1604/Populus_tremuloides_var_CAM1604-4_HAP1_V2_release/Populus_tremuloides_var_CAM1604-4/sequences/Populus_tremuloides_var_CAM1604-4_HAP1.mainGenome.fasta.0123"
    conda: "/global/scratch/projects/fc_moilab/aphillips/aspen_snakemake/envs/bwa-mem2.yaml"
    shell:
        """
        /global/home/users/arphillips/toolz/bwa-mem2-2.2.1_x64-linux/bwa-mem2 index {input}
        samtools faidx {input}
        gatk CreateSequenceDictionary -R {input}
        """

# (3b) Align reads to the reference genome
# Paired-end reads in single file
rule bwa_map:
    input:
        ref = config["data"]["reference"]["genome"],
        index = "/global/scratch/projects/fc_moilab/projects/aspen/genome/CAM1604/Populus_tremuloides_var_CAM1604-4_HAP1_V2_release/Populus_tremuloides_var_CAM1604-4/sequences/Populus_tremuloides_var_CAM1604-4_HAP1.mainGenome.fasta.0123",
        trim = "/global/scratch/users/arphillips/data/trimmed/{sample}.trim.fastq.gz"
    output:
        temp("/global/scratch/users/arphillips/data/interm/mapped_bam/{sample}.mapped.bam")
    conda: "/global/scratch/projects/fc_moilab/aphillips/aspen_snakemake/envs/bwa_map.yaml"
#    benchmark:
#        "/global/scratch/users/arphillips/benchmarks/{sample}.bwa.benchmark.txt"
    shell:
        "/global/scratch/users/arphillips/toolz/bwa-mem2-2.2.1_x64-linux/bwa-mem2 mem -p -t 10 {input.ref} {input.trim} |"
        "samtools view -Sb > {output}"

# (4) Sort bams
rule samtools_sort:
    input:
        "/global/scratch/users/arphillips/data/interm/mapped_bam/{sample}.mapped.bam"
    output:
        temp("/global/scratch/users/arphillips/data/interm/sorted_bam/{sample}.sorted.bam"),
    params:
        tmp = "/global/scratch/users/arphillips/temp/sort_bam/{sample}"
    conda: "/global/scratch/projects/fc_moilab/aphillips/aspen_snakemake/envs/samtools.yaml"
#    benchmark:
#        "/global/scratch/users/arphillips/benchmarks/{sample}.sort.benchmark.txt"
    shell:
        """ 
        mkdir -p {params.tmp}
        samtools sort -T {params.tmp} -@ 1 {input} > {output}
        rm -rf {params.tmp}
        """

# (5) Add read groups
rule add_rg:
    input:
        "/global/scratch/users/arphillips/data/interm/sorted_bam/{sample}.sorted.bam"
    output:
        bam = temp(touch("/global/scratch/users/arphillips/data/interm/addrg/{sample}.rg.bam"))
    params:
        tmp = "/global/scratch/users/arphillips/temp/addrg/{sample}",
        sample = "{sample}",
        rg = randint(1,1000)
    conda: "/global/scratch/projects/fc_moilab/aphillips/aspen_snakemake/envs/gatk.yaml"
#    benchmark:
#        "/global/scratch/users/arphillips/benchmarks/{sample}.rg.benchmark.txt"
    shell:
        """
        mkdir -p {params.tmp}
        gatk --java-options ""-Xmx4G"" AddOrReplaceReadGroups \
        -I {input} \
        -O {output.bam} \
        -RGID {params.rg} \
        -RGLB lib1 \
        -RGPL illumina \
        -RGPU unit1 \
        -RGSM {params.sample} \
        --TMP_DIR {params.tmp} \
        --CREATE_INDEX true
        rm -rf {params.tmp}
        """

# (6) Mark duplicates
rule mark_dups:
    input:
        "/global/scratch/users/arphillips/data/interm/addrg/{sample}.rg.bam"
    output:
#        bam = "/global/scratch/users/arphillips/data/interm/mark_dups/{sample}.dedup.bam"
        bam = "/global/scratch/projects/fc_moilab/aphillips/aspen_snakemake/data/bams/{sample}.dedup.bam",
    params:
        tmp = "/global/scratch/users/arphillips/temp/mark_dups/{sample}"
    conda: "/global/scratch/projects/fc_moilab/aphillips/aspen_snakemake/envs/gatk.yaml"
#    benchmark:
#        "/global/scratch/users/arphillips/benchmarks/{sample}.dups.benchmark.txt"
    shell:
        """
        # Create a scratch directory
        mkdir -p {params.tmp}
        # Input bam file to output marked records. Assume bam file has been sorted. Direct to a temporary storage file (scratch).
        gatk --java-options ""-Xmx10G"" MarkDuplicates \
        -I {input} \
        -O {output.bam} \
        --CREATE_INDEX true \
        -MAX_FILE_HANDLES 1000 \
        --ASSUME_SORT_ORDER coordinate \
        --TMP_DIR {params.tmp} \
        --METRICS_FILE {params.tmp}.metrics
        # Remove scratch directory
        rm -rf {params.tmp}
        """

# (7) Assess alignment quality metrics with qualimap
# nr is normally 100000 and -nt is normally 8, java mem size = 48
# nw is normally 400
# for higher cov, make nr 1000 and -nt 12, java mem size = 64
# then remove all the excess files we don't need
rule bamqc:
    input:
        "/global/scratch/projects/fc_moilab/aphillips/aspen_snakemake/data/bams/{sample}.dedup.bam",
#        "/global/scratch/users/arphillips/data/interm/mark_dups/{sample}.dedup.bam"
    output:
        "/global/scratch/users/arphillips/reports/bamqc/{sample}_stats/genome_results.txt"
    params:
        dir = "/global/scratch/users/arphillips/reports/bamqc/{sample}_stats"
    conda: "/global/scratch/projects/fc_moilab/aphillips/aspen_snakemake/envs/qualimap.yaml"
#    benchmark:
#        "/global/scratch/users/arphillips/benchmarks/{sample}.bamqc.benchmark.txt"
    shell:
        """
        qualimap bamqc \
        -bam {input} \
        -nt 12 \
        -nr 1000 \
        -nw 400 \
        -outdir {params.dir} \
        -outformat HTML \
        --skip-duplicated \
        --java-mem-size=24G
        rm -r {params.dir}/css  {params.dir}/qualimapReport.html {params.dir}/images_qualimapReport  {params.dir}/raw_data_qualimapReport
        """

# (8) Assess DNA damage with mapDamage
rule mapdamage:
    input:
        bam = "/global/scratch/projects/fc_moilab/aphillips/aspen_snakemake/data/bams/{sample}.dedup.bam",
#        bam = "/global/scratch/users/arphillips/data/interm/mark_dups/{bam}.dedup.bam",
        ref = config["data"]["reference"]["genome"]
    output:
        "/global/scratch/users/arphillips/reports/mapdamage/{bam}/5pCtoT_freq.txt"
    params:
        outdir = "/global/scratch/users/arphillips/reports/mapdamage/{bam}"
    conda: "/global/home/users/arphillips/.conda/envs/mapdamage"
    benchmark:
        "/global/scratch/users/arphillips/benchmarks/{bam}.mapdamage.benchmark.txt"
    shell:
         "mapDamage -i {input.bam} -r {input.ref} -d {params.outdir}"

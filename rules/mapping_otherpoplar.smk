# (2) Trim reads sequenced at UCD with fastp
# Minimum length is 36 (-l 36)
# Don't filter for quality (-Q)
# Adapter trimming is enabled by default -- don't need to specify adapter seq
# Default detects and trims polyG tails for NovaSeq data
# --cut_front is sliding window trimming from 5' to 3'
rule fastp_trim:
    input:
        f1 = "/global/scratch/users/arphillips/raw/other_poplars/{other_pop}*1.f*q.gz",
        r2 = "/global/scratch/users/arphillips/raw/other_poplars/{other_pop}*2.f*q.gz"
    output:
        trim1 = temp("/global/scratch/users/arphillips/data/trimmed/{other_pop}.trim_1.fastq.gz"),
        trim2 = temp("/global/scratch/users/arphillips/data/trimmed/{other_pop}.trim_2.fastq.gz")
#    conda:
#        "/global/scratch/projects/fc_moilab/aphillips/aspen_snakemake/envs/fastp.yaml"
    shell:
        """
        fastp -w 2 \
        -l 36 -Q \
        -i {input.f1} \
        -I {input.r2} \
        -o {output.trim1} \
        -O {output.trim2} \
        --cut_front --cut_front_window_size 4 --cut_front_mean_quality 15 
        """

# (3a) Prepare reference file
# 9.6 GB per core available (core == threads)
rule bwa_prep:
    input: 
        config["data"]["reference"]["genome"]
    output:
        index = "/global/scratch/projects/fc_moilab/projects/aspen/genome/CAM1604/Populus_tremuloides_var_CAM1604-4_HAP1_V2_release/Populus_tremuloides_var_CAM1604-4/sequences/Populus_tremuloides_var_CAM1604-4_HAP1.mainGenome.fasta.0123"
#    conda: "/global/scratch/projects/fc_moilab/aphillips/aspen_snakemake/envs/bwa-mem2.yaml"
    shell:
        """
        /global/home/users/arphillips/toolz/bwa-mem2-2.2.1_x64-linux/bwa-mem2 index {input}
        samtools faidx {input}
        gatk CreateSequenceDictionary -R {input}
        """

# (3b) Align reads to the reference genome
# Paired-end reads
rule bwa_map:
    input:
        ref = config["data"]["reference"]["genome"],
        index = "/global/scratch/projects/fc_moilab/projects/aspen/genome/CAM1604/Populus_tremuloides_var_CAM1604-4_HAP1_V2_release/Populus_tremuloides_var_CAM1604-4/sequences/Populus_tremuloides_var_CAM1604-4_HAP1.mainGenome.fasta.0123",
        trim1 = "/global/scratch/users/arphillips/data/trimmed/{other_pop}.trim_1.fastq.gz",
        trim2 = "/global/scratch/users/arphillips/data/trimmed/{other_pop}.trim_2.fastq.gz"
    output:
        temp("/global/scratch/users/arphillips/data/interm/mapped_bam/{other_pop}.mapped.bam")
#    conda: "/global/scratch/projects/fc_moilab/aphillips/aspen_snakemake/envs/bwa_map.yaml"
    shell:
        "/global/scratch/users/arphillips/toolz/bwa-mem2-2.2.1_x64-linux/bwa-mem2 mem -t 10 {input.ref} {input.trim1} {input.trim2} |"
        "samtools view -Sb > {output}"

# (4) Sort bams
rule samtools_sort:
    input:
        "/global/scratch/users/arphillips/data/interm/mapped_bam/{other_pop}.mapped.bam"
    output:
        temp("/global/scratch/users/arphillips/data/interm/sorted_bam/{other_pop}.sorted.bam"),
    params:
        tmp = "/global/scratch/users/arphillips/temp/sort_bam/{other_pop}"
#    conda: "/global/scratch/projects/fc_moilab/aphillips/aspen_snakemake/envs/samtools.yaml"
    shell:
        """ 
        mkdir -p {params.tmp}
        samtools sort -T {params.tmp} -@ 1 {input} > {output}
        rm -rf {params.tmp}
        """

# (5) Add read groups
rule add_rg:
    input:
        "/global/scratch/users/arphillips/data/interm/sorted_bam/{other_pop}.sorted.bam"
    output:
        bam = temp(touch("/global/scratch/users/arphillips/data/interm/addrg/{other_pop}.rg.bam"))
    params:
        tmp = "/global/scratch/users/arphillips/temp/addrg/{other_pop}",
        sample = "{other_pop}",
        rg = randint(1,1000)
#    conda: "/global/scratch/projects/fc_moilab/aphillips/aspen_snakemake/envs/gatk.yaml"
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
        "/global/scratch/users/arphillips/data/interm/addrg/{other_pop}.rg.bam"
    output:
        bam = "/global/scratch/projects/fc_moilab/aphillips/aspen_snakemake/data/bams/{other_pop}.dedup.bam",
    params:
        tmp = "/global/scratch/users/arphillips/temp/mark_dups/{other_pop}"
#    conda: "/global/scratch/projects/fc_moilab/aphillips/aspen_snakemake/envs/gatk.yaml"
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
        "/global/scratch/projects/fc_moilab/aphillips/aspen_snakemake/data/bams/{other_pop}.dedup.bam",
    output:
        "/global/scratch/users/arphillips/reports/bamqc/{other_pop}_stats/genome_results.txt"
    params:
        dir = "/global/scratch/users/arphillips/reports/bamqc/{other_pop}_stats"
#    conda: "/global/scratch/projects/fc_moilab/aphillips/aspen_snakemake/envs/qualimap.yaml"
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


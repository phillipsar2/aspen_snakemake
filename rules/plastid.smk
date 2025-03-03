# Rules for assessing the plastid genome evolution of aspen

# (1) Index plastid reference
rule prep_plas:
    input: 
        ref = "/global/scratch/projects/fc_moilab/aphillips/aspen_snakemake/ref/plastid/MW376839.1.fasta"
    output:
        index = "/global/scratch/projects/fc_moilab/aphillips/aspen_snakemake/ref/plastid/MW376839.1.fasta.0123"
    conda: "/global/scratch/projects/fc_moilab/aphillips/aspen_snakemake/envs/bwa-mem2.yaml"
    shell:
        """
        /global/home/users/arphillips/toolz/bwa-mem2-2.2.1_x64-linux/bwa-mem2 index {input}
        samtools faidx {input}
        gatk CreateSequenceDictionary -R {input}
        """

# (2) Align reads to the plastid reference
# Paired-end reads in single file
rule bwa_plastid:
    input:
        ref = "/global/scratch/projects/fc_moilab/aphillips/aspen_snakemake/ref/plastid/MW376839.1.fasta",
        index = "/global/scratch/projects/fc_moilab/aphillips/aspen_snakemake/ref/plastid/MW376839.1.fasta.0123",
        trim = "/global/scratch/users/arphillips/data/trimmed/{sample}.trim.fastq.gz"
    output:
        temp("/global/scratch/users/arphillips/data/interm/mapped_chl/{sample}.mapped_chl.bam")
    conda: "/global/scratch/projects/fc_moilab/aphillips/aspen_snakemake/envs/bwa_map.yaml"
    shell:
        """
        /global/scratch/users/arphillips/toolz/bwa-mem2-2.2.1_x64-linux/bwa-mem2 mem -p -t 10 {input.ref} {input.trim} |
        samtools view -Sb > {output}
        """

# (3) Sort bams by read name (-n)
rule plas_sort:
    input:
        "/global/scratch/users/arphillips/data/interm/mapped_chl/{sample}.mapped_chl.bam"
    output:
        temp("/global/scratch/users/arphillips/data/interm/sorted_chl/{sample}.sorted_chl.bam")
    params:
        tmp = "/global/scratch/users/arphillips/temp/sort_bam/{sample}.sorted_chl"
    conda: "/global/scratch/projects/fc_moilab/aphillips/aspen_snakemake/envs/samtools.yaml"
    shell:
        """ 
        mkdir -p {params.tmp}
        samtools sort -T {params.tmp} -@ 1 -n {input} > {output}
        rm -rf {params.tmp}
        """

# (4) Extract the relevant reads
# remove singletons -s /dev/null
# -n keep read names as they are
# --exclude-flags 4 remove unmapped reads
rule extract_plas:
    input:
        "/global/scratch/users/arphillips/data/interm/sorted_chl/{sample}.sorted_chl.bam"
    output:
        r1 = "/global/scratch/users/arphillips/data/plastid/fastq/{sample}.R1.fastq.gz",
        r2 = "/global/scratch/users/arphillips/data/plastid/fastq/{sample}.R2.fastq.gz"
#    params:
#        r1 = "/global/scratch/users/arphillips/data/plastid/fastq/{sample}.R1.fastq",
#        r2 = "/global/scratch/users/arphillips/data/plastid/fastq/{sample}.R2.fastq"
    conda: "/global/scratch/projects/fc_moilab/aphillips/aspen_snakemake/envs/samtools.yaml"
    shell:
        """
        samtools fastq -1 {output.r1} -2 {output.r2} -0 /dev/null -s /dev/null -n --exclude-flags 4 {input} 
        """

configfile: "profiles/config.yaml"

SAMPLE = glob_wildcards("/global/scratch/users/arphillips/raw/jgi_wgs/{sample}.fastq.gz").sample

# =================================================================================================
#     Target Rules
# =================================================================================================
rule all:
    input:
        #fastqc = expand("/global/scratch/users/arphillips/qc/fastqc/{sample}_fastqc.zip", sample = SAMPLE),
        fastp = expand("/global/scratch/users/arphillips/data/trimmed/{sample}.trim.fastq.gz" , sample = SAMPLE),
        bwa_prep = "/global/scratch/projects/fc_moilab/projects/aspen/genome/mex_genome/genome.1MX.fasta.fai"


# =================================================================================================
#     Rule Modules
# =================================================================================================
include: "rules/mapping.smk"

configfile: "config.yaml"

SAMPLE = glob_wildcards("/global/scratch/users/arphillips/raw/jgi_wgs/{sample}.fastq.gz").sample

# =================================================================================================
#     Target Rules
# =================================================================================================
rule all:
    input:
        fastqc = expand("/global/scratch/users/arphillips/qc/fastqc/{sample}/{sample}_fastqc.zip", sample = SAMPLE)


# =================================================================================================
#     Rule Modules
# =================================================================================================
include: "rules/mapping.smk"

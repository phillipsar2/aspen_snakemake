configfile: "profiles/config.yaml"

from random import randint

# Sample names (JGI filenames)
SAMPLE = glob_wildcards("/global/scratch/users/arphillips/raw/jgi_wgs/{sample}.fastq.gz").sample

# BAMs to process
BAM = glob_wildcards("/global/scratch/users/arphillips/data/interm/mark_dups/{bam}.dedup.bam").bam 
TEST = BAM[0:3]
print(TEST)

# =================================================================================================
#     Target Rules
# =================================================================================================
rule all:
    input:
        #fastqc = expand("/global/scratch/users/arphillips/qc/fastqc/{sample}_fastqc.zip", sample = SAMPLE),
#        fastp = expand("/global/scratch/users/arphillips/data/trimmed/{sample}.trim.fastq.gz" , sample = SAMPLE),
#        bam = expand("/global/scratch/users/arphillips/data/interm/mark_dups/{sample}.dedup.bam", sample = SAMPLE),
        bamqc = expand("/global/scratch/users/arphillips/reports/bamqc/{sample}_stats/genome_results.txt", sample = SAMPLE)
#        mapdamage = expand("/global/scratch/users/arphillips/reports/mapdamage/{bams}/5pCtoT_freq.txt", bams = TEST)
#        toz19 = expand("/global/scratch/users/arphillips/data/toz19/{bam}.toz19.cov.txt", bam = BAM)

# =================================================================================================
#     Rule Modules
# =================================================================================================
include: "rules/mapping.smk"
include: "rules/ploidy_sex.smk"

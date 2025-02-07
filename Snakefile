configfile: "profiles/config.yaml"
import pandas as pd
from random import randint
import datetime

# Sample names (JGI filenames)
SAMPLE = glob_wildcards("/global/scratch/users/arphillips/raw/jgi_wgs/{sample}.fastq.gz").sample

# BAMs to process
BAM = glob_wildcards("/global/scratch/users/arphillips/data/interm/mark_dups/{bam}.dedup.bam").bam 

# Chromosomes
fai =  pd.read_csv("/global/scratch/projects/fc_moilab/projects/aspen/genome/CAM1604/Populus_tremuloides_var_CAM1604-4_HAP1_V2_release/Populus_tremuloides_var_CAM1604-4/sequences/Populus_tremuloides_var_CAM1604-4_HAP1.mainGenome.fasta.fai", header = None, sep = "\t")
CHR = list(fai[0])
#print(CHROM)

# Date
DATE = datetime.datetime.utcnow().strftime("%Y-%m-%d")

# SNP filters
MAX_DP = ["75"]
MIN_DP = ["10"]

# =================================================================================================
#     Target Rules
# =================================================================================================
rule all:
    input:
      ## Maping
        #fastqc = expand("/global/scratch/users/arphillips/qc/fastqc/{sample}_fastqc.zip", sample = SAMPLE),
#        fastp = expand("/global/scratch/users/arphillips/data/trimmed/{sample}.trim.fastq.gz" , sample = `SAMPLE),
#        bam = expand("/global/scratch/users/arphillips/data/interm/mark_dups/{sample}.dedup.bam", sample = SAMPLE),
#        bamqc = expand("/global/scratch/users/arphillips/reports/bamqc/{sample}_stats/genome_results.txt", sample = SAMPLE),
        mapdamage = expand("/global/scratch/users/arphillips/reports/mapdamage/{bams}/5pCtoT_freq.txt", bams = TEST)
      ## Sex
#        depth = expand("/global/scratch/users/arphillips/data/toz19/{bam}.chr13.cov.txt", bam = BAM),
      ## Calling and filtering
        snp = expand("/global/scratch/users/arphillips/reports/filtering/wgs_aspen.{chr}.table", chr = CHR),
        dp_table = expand("/global/scratch/users/arphillips/reports/filtering/depth/wgs_aspen.{chr}.filtered.nocall.table", chr = CHR),
        vcf = expand("/global/scratch/users/arphillips/data/processed/filtered_snps/wgs_aspen.{chr}.nocall.{min_dp}dp{max_dp}.vcf", chr = CHR, min_dp = MIN_DP, max_dp = MAX_DP)

# =================================================================================================
#     Rule Modules
# =================================================================================================
include: "rules/mapping.smk"
include: "rules/calling.smk"
include: "rules/ploidy_sex.smk"

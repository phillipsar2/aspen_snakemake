configfile: "profiles/config.yaml"
import pandas as pd
from random import randint
import datetime

# Sample names (JGI filenames)
SAMPLE = glob_wildcards("/global/scratch/users/arphillips/raw/jgi_wgs/{sample}.fastq.gz").sample
#print(SAMPLE)

# Other poplar samples (filenames)
OTHER_POP = glob_wildcards("/global/scratch/users/arphillips/raw/other_poplars/{other_pop}_1.fastq.gz").other_pop
#print(OTHER_POP)

# BAMs to process
BAM = glob_wildcards("/global/scratch/projects/fc_moilab/aphillips/aspen_snakemake/data/bams/{bam}.dedup.bam").bam 

# Chromosomes
fai =  pd.read_csv("/global/scratch/projects/fc_moilab/projects/aspen/genome/CAM1604/Populus_tremuloides_var_CAM1604-4_HAP1_V2_release/Populus_tremuloides_var_CAM1604-4/sequences/Populus_tremuloides_var_CAM1604-4_HAP1.mainGenome.fasta.fai", header = None, sep = "\t")
CHR = list(fai[0])
#print(CHROM)

# 1 Mb regions
region_list = pd.read_csv("/global/scratch/projects/fc_moilab/projects/aspen/genome/CAM1604/Populus_tremuloides_var_CAM1604-4_HAP1_V2_release/Populus_tremuloides_var_CAM1604-4/sequences/chr_regions.txt", header = None, sep = "\t")
REGION = list(region_list[0])

# Date
DATE = datetime.datetime.utcnow().strftime("%Y-%m-%d")

# SNP filters
MAX_DP = ["75"]
MIN_DP = ["10"]

# Ploidy range
GENOTYPE_PLOIDY = ["2", "3"]
PLOIDY = ["diploid", "triploid"]

# =================================================================================================
#     Target Rules
# =================================================================================================
rule all:
    input:
      ## Mapping
        #fastqc = expand("/global/scratch/users/arphillips/qc/fastqc/{sample}_fastqc.zip", sample = SAMPLE),
#        fastp = expand("/global/scratch/users/arphillips/data/trimmed/{sample}.trim.fastq.gz" , sample = SAMPLE),
#        bam = expand("/global/scratch/projects/fc_moilab/aphillips/aspen_snakemake/data/bams/{sample}.dedup.bam", sample = SAMPLE),
#        bamqc = expand("/global/scratch/users/arphillips/reports/bamqc/{sample}_stats/genome_results.txt", sample = SAMPLE)
#        mapdamage = expand("/global/scratch/users/arphillips/reports/mapdamage/{bams}/5pCtoT_freq.txt", bams = BAM)
      ## Mapping other poplar
#        mapped = expand("/global/scratch/users/arphillips/data/interm/mapped_bam/{other_pop}.mapped.bam", other_pop = OTHER_POP),
#        bam = expand("/global/scratch/projects/fc_moilab/aphillips/aspen_snakemake/data/bams/{other_pop}.dedup.bam", other_pop = OTHER_POP),
#        bamqc = expand("/global/scratch/users/arphillips/reports/bamqc/{other_pop}_stats/genome_results.txt", other_pop = OTHER_POP)
      ## Sex
#        depth = expand("/global/scratch/users/arphillips/data/toz19/{bam}.chr13.cov.txt", bam = BAM),
      ## Calling and filtering
        raw_vcf = expand("/global/scratch/users/arphillips/data/vcf/wgs_aspen.{bam}.raw.vcf.gz", bam = BAM)
#        snp = expand("/global/scratch/users/arphillips/reports/filtering/wgs_aspen.{chr}.table", chr = REGION),
#        dp_table = expand("/global/scratch/users/arphillips/reports/filtering/depth/wgs_aspen.{chr}.filtered.nocall.table", chr = REGION),
#        vcf = expand("/global/scratch/users/arphillips/data/processed/filtered_snps/wgs_aspen.{chr}.nocall.{min_dp}dp{max_dp}.vcf", chr = REGION, min_dp = MIN_DP, max_dp = MAX_DP)
      ## Plastids
#        map_plas = expand("/global/scratch/users/arphillips/data/interm/mapped_chl/{sample}.mapped_chl.bam", sample = SAMPLE)
#        plas_fastq = expand("/global/scratch/users/arphillips/data/plastid/fastq/{sample}.R1.fastq.gz", sample = SAMPLE)
      ## Genotyping
#        updog = expand("/global/scratch/users/arphillips/data/updog/updog.genomat.{ploidy}.{chr}.{date}.txt", ploidy = "triploid", chr = CHR, date = DATE)

# =================================================================================================
#     Rule Modules
# =================================================================================================
#include: "rules/mapping.smk"
#include: "rules/mapping_otherpoplar.smk"
include: "rules/calling.smk"
#include: "rules/ploidy_sex.smk"
#include: "rules/plastid.smk"
#include: "rules/updog_genotyping.smk"

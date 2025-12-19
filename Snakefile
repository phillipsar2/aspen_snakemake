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
#BAM = glob_wildcards("/global/scratch/projects/fc_moilab/aphillips/aspen_snakemake/data/bams/{bam}.merged.dedup.bam").bam
#print(BAM)

# Chromosomes
fai =  pd.read_csv("/global/scratch/projects/fc_moilab/projects/aspen/genome/CAM1604/Populus_tremuloides_var_CAM1604-4_HAP1_V2_release/Populus_tremuloides_var_CAM1604-4/sequences/Populus_tremuloides_var_CAM1604-4_HAP1.mainGenome.fasta.fai", header = None, sep = "\t")
CHR = list(fai[0])

# 1 Mb regions
region_list = pd.read_csv("/global/scratch/projects/fc_moilab/projects/aspen/genome/CAM1604/Populus_tremuloides_var_CAM1604-4_HAP1_V2_release/Populus_tremuloides_var_CAM1604-4/sequences/chr_regions.txt", header = None, sep = "\t")
REGION = list(region_list[0])
R = glob_wildcards("/global/scratch/users/arphillips/data/updog/updog.genomat.diploid.Ch{r}.vcf.gz").r

# Date
DATE = datetime.datetime.utcnow().strftime("%Y-%m-%d")

# SNP filters
MAX_DP = ["90"]
MIN_DP = ["10"]

# Ploidy range - updog
#GENOTYPE_PLOIDY = ["2", "3"]
#PLOIDY = ["diploid", "triploid"]

# MERGE contains a list of the bams that belong to each genotype (GENO) so they can be merged
file = pd.read_csv("/global/scratch/users/arphillips/reports/filestomerge.08122025.txt", sep = " ", header = 0)
MERGE_A = list(file.Merge_A)
MERGE_B = list(file.Merge_B)
GENO = list(file.Genotype)
#print(GENO)

# Ploidy genotype pairing
pgfile = pd.read_csv("/global/scratch/projects/fc_moilab/aphillips/aspen_snakemake/metadata/ploidy.geno.1127.2025-09-30.csv", sep = ",", header = 0)
GENOTYPE =  list(pgfile["sample"])
GENOTYPE_PLOIDY = list(pgfile.ploidy)

# =================================================================================================
#     Target Rules
# =================================================================================================
rule all:
    input:
      ## Mapping
#        fastp = expand("/global/scratch/users/arphillips/data/trimmed/{sample}.trim.fastq.gz" , sample = SAMPLE),
#        bam = expand("/global/scratch/projects/fc_moilab/aphillips/aspen_snakemake/data/bams/{sample}.dedup.bam", sample = SAMPLE),
#        bamqc = expand("/global/scratch/users/arphillips/reports/bamqc/{sample}_stats/genome_results.txt", sample = SAMPLE)
#        addeam = "/global/scratch/users/arphillips/reports/addeam/plots/damage_report_k3.pdf"
#        merge_bams = expand("/global/scratch/projects/fc_moilab/aphillips/aspen_snakemake/data/bams/{merge_A}_{merge_B}.merged.dedup.bam", zip,  merge_A = MERGE_A, merge_B = MERGE_B)
      ## Mapping other poplar
#        mapped = expand("/global/scratch/users/arphillips/data/interm/mapped_bam/{other_pop}.mapped.bam", other_pop = OTHER_POP),
#        bam = expand("/global/scratch/projects/fc_moilab/aphillips/aspen_snakemake/data/bams/{other_pop}.dedup.bam", other_pop = OTHER_POP),
#        bamqc = expand("/global/scratch/users/arphillips/reports/bamqc/{other_pop}_stats/genome_results.txt", other_pop = OTHER_POP)
      ## Calling and filtering
#        raw_vcf = expand("/global/scratch/users/arphillips/data/vcf/wgs_aspen.{bam}.raw.vcf.gz", bam = BAM),
#        merge_raw = expand("/global/scratch/users/arphillips/data/vcf/wgs_aspen.{region}.raw.merged.vcf.gz", region = REGION),
#        diag = expand("/global/scratch/users/arphillips/reports/filtering/wgs_aspen.{region}.table", region = REGION),
#        snp = expand("/global/scratch/users/arphillips/reports/filtering/wgs_aspen.{chr}.table", chr = REGION),
#        dp_table = expand("/global/scratch/users/arphillips/reports/filtering/depth/wgs_aspen.{chr}.filtered.nocall.table", chr = REGION),
#        filt_vcf = expand("/global/scratch/users/arphillips/data/processed/filtered_snps/wgs_aspen.{region}.nocall.{min_dp}dp{max_dp}.vcf.gz", region = REGION, min_dp = MIN_DP, max_dp = MAX_DP)
#        geno_filt = expand("/global/scratch/users/arphillips/data/processed/filtered_snps/wgs_aspen.{region}.goodg.{min_dp}dp{max_dp}.vcf.gz", region = REGION, min_dp = MIN_DP, max_dp = MAX_DP)
      ## Genotyping
#        merge_filt = expand("/global/scratch/users/arphillips/data/processed/filtered_snps/wgs_aspen.all.nocall.{min_dp}dp{max_dp}.vcf.gz", min_dp = MIN_DP, max_dp = MAX_DP),
#        gbs2ploidy = expand("/global/scratch/users/arphillips/data/gbs2ploidy/{bam}.propOut.csv", bam = BAM)
      ## Plastids
#        map_plas = expand("/global/scratch/users/arphillips/data/interm/mapped_chl/{sample}.mapped_chl.bam", sample = SAMPLE)
#        plas_fastq = expand("/global/scratch/users/arphillips/data/plastid/fastq/{sample}.R1.fastq.gz", sample = SAMPLE)
      ## Genotyping
#         haplotype = expand("/global/scratch/users/arphillips/data/vcf/gatk/called/{geno}.{region}.g.vcf.gz", geno = GENOTYPE, region = CHR)
#        merge_gvcfs = expand("/global/scratch/users/arphillips/data/vcf/gatk/merged/{geno}.g.vcf.gz", geno = GENOTYPE),
        genotyping = expand("/global/scratch/users/arphillips/data/vcf/gatk/called/{geno}_p{geno_ploidy}.g.vcf.gz", zip,geno = GENOTYPE, geno_ploidy = GENOTYPE_PLOIDY)
#        combine)gvcfs = expand("/global/scratch/users/arphillips/data/vcf/gatk/called/wgs_aspen.all.genos.{region}.g.vcf.gz", region = CHR)

# =================================================================================================
#     Rule Modules
# =================================================================================================
#include: "rules/mapping.smk"
#include: "rules/mapping_otherpoplar.smk"
#include: "rules/calling.smk"
#include: "rules/ploidy_sex.smk"
#include: "rules/plastid.smk"
#include: "rules/updog_genotyping.smk"
include: "rules/gatk_genotyping.smk"
#include: "rules/freebayes_genotyping.smk"

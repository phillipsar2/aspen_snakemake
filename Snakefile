configfile: "profiles/config.yaml"
import pandas as pd
from random import randint
import datetime

# Sample names (JGI filenames)
SAMPLE = glob_wildcards("/global/scratch/users/arphillips/data/fastq/{sample}.fastq.gz").sample
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

# Set number of intervals for gatk to 200
INTERVALS = ["{:04d}".format(x) for x in list(range(200))]
wildcard_constraints: 
    intv = r"\d{4}"

# 1 Mb regions
region_list = pd.read_csv("/global/scratch/projects/fc_moilab/projects/aspen/genome/CAM1604/Populus_tremuloides_var_CAM1604-4_HAP1_V2_release/Populus_tremuloides_var_CAM1604-4/sequences/chr_regions.txt", header = None, sep = "\t")
REGION = list(region_list[0])
R = glob_wildcards("/global/scratch/users/arphillips/data/updog/updog.genomat.diploid.Ch{r}.vcf.gz").r

# Date
DATE = datetime.datetime.utcnow().strftime("%Y-%m-%d")

# Ploidy range - updog
#GENOTYPE_PLOIDY = ["2", "3"]
#PLOIDY = ["diploid", "triploid"]

# MERGE contains a list of the bams that belong to each genotype (GENO) so they can be merged
#file = pd.read_csv("reports/filestomerge_or_resolve.08202025.txt", sep = " ", header = 0)
#MERGE_A = list(file.Merge_A)
#MERGE_B = list(file.Merge_B)
#GENO = list(file.Genotype)
#print(GENO)

# Ploidy genotype pairing
pgfile = pd.read_csv("metadata/ploidy.geno.1255.2026-08-07.csv", sep = ",", header = 0)
#pgfile = pd.read_csv("metadata/ploidy.geno.1127.2025-09-30.triploids.csv", sep = ",", header = 0)
#pgfile = pd.read_csv("metadata/california/california.ploidy.geno.01-20-2026.triploids.csv",sep = ",", header = 0)
GENOTYPE =  list(pgfile["sample"])
GENOTYPE_PLOIDY = list(pgfile.ploidy)

# =================================================================================================
#     Target Rules
# =================================================================================================
rule all:
    input:
      ## Mapping
#        bam = expand("/global/scratch/projects/fc_moilab/aphillips/aspen_snakemake/data/bams/{sample}.dedup.bam", sample = SAMPLE),
#        bamqc = expand("reports/bamqc/{sample}_stats/genome_results.txt", sample = SAMPLE),
#        addeam = "/global/scratch/users/arphillips/reports/addeam/plots/damage_report_k3.pdf"
      ## Ploidy calling
#        gbs2ploidy = expand("data/gbs2ploidy/{sample}.propOut.csv", sample = SAMPLE),
      ## Plastids
#        map_plas = expand("/global/scratch/users/arphillips/data/interm/mapped_chl/{sample}.mapped_chl.bam", sample = SAMPLE)
#        plas_fastq = expand("/global/scratch/users/arphillips/data/plastid/fastq/{sample}.R1.fastq.gz", sample = SAMPLE)
      ## Genotyping
#         haplotype = expand("data/vcf/gatk/called/{geno}_p{geno_ploidy}.{region}.haplo.g.vcf.gz", geno = GENOTYPE, region = CHR, geno_ploidy = GENOTYPE_PLOIDY)
        genotyping = expand("data/vcf/gatk/genotyped/{geno}_p{geno_ploidy}.g.vcf.gz", zip,geno = GENOTYPE, geno_ploidy = GENOTYPE_PLOIDY)
#        merge_vcfs = expand("data/vcf/gatk/genotyped/wgs_aspen.all.{chr}.g.vcf.gz", chr = CHR)
      ## Filtering
#        qual = expand("reports/filtering/{geno}_p{geno_ploidy}.table",  zip,geno = GENOTYPE, geno_ploidy = GENOTYPE_PLOIDY)

# =================================================================================================
#     Rule Modules
# =================================================================================================
#include: "rules/mapping.smk"
#include: "rules/bcftools_genotyping.smk"
#include: "rules/plastid.smk"
include: "rules/gatk_genotyping.smk"
#include: "rules/filtering.smk"

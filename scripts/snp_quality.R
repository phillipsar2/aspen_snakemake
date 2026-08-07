## Title: Assess SNP quality per site
## Author: Alyssa Phillips
## Date: 7/1/2026

library(data.table)
library(dplyr)

# load in fai
fai <- read.table("/global/scratch/projects/fc_moilab/projects/aspen/genome/CAM1604/Populus_tremuloides_var_CAM1604-4_HAP1_V2_release/Populus_tremuloides_var_CAM1604-4/sequences/Populus_tremuloides_var_CAM1604-4_HAP1.mainGenome.fasta.fai" )
head(fai)

# load in variant quality tables
path="/global/scratch/projects/fc_moilab/aphillips/aspen_snakemake/reports/filtering/"


## Chr01
chr01_size <- fai$V2[fai$V1 == "Chr01"]
df1 <- fread(paste0(path, "52830.3.464152.AGCTTGAG-AGCTTGAG_p2.table"), 
             sep = "\t", header = T,select = c("CHROM", "POS", "MQ"), 
             nrows = chr01_size, skip = 0)

df2 <- fread(paste0(path, "52830.3.464152.CGGAATAC-CGGAATAC_p2.table"), 
             sep = "\t", header = T,select = c("MQ"), 
             nrows = chr01_size, skip = 0)

df3 <- fread(paste0(path, "52830.4.464202.AAGTCGAG-AAGTCGAG_p2.table"), 
             sep = "\t", header = T,select = c("MQ"), 
             nrows = chr01_size, skip = 0)


# How many SNPs for each individual?
sum(!is.na(df1$MQ)) # 49,222,514 
sum(!is.na(df2$MQ)) # 48,836,420
sum(!is.na(df3$MQ)) # 49,095,311

# How many SNPs overall?
na_count <- cbind(df1$DP, df2$DP, df3$DP) %>%
  is.na() %>%
  rowSums()

sum(na_count < 3) # sites data for at least one individual = 54,773,111


# What is the similarity of MQ scores across samples?



# What is the maximum quality score for each variant?
## When merging vcfs, bcftools selects the maximum QUAL score

### Title: Aspen LFMM
### Author: Alyssa Phillips
### Date: 11/15/25

# modified from https://github.com/moiexpositoalonsolab/grenephase1-paper/blob/main/Genome_Environment_Association/lfmm/lfmm_first_gen/run_lfmm_full_genome_ridge_multiplebiovars.r

###
### LFMM Analysis for Multiple Bioclimatic Variables
###

#
# This script performs LFMM (Latent Factor Mixed Model) analysis for a given
# bioclimatic variable using a specified number of components (K).
#
# Required Input Files:
# 1. /carnegie/nobackup/scratch/tbellagio/gea_grene-net/key_files/delta_p_maf05_mincount05_firstgensamples.csv
#    - Allele frequency data
# 2. /carnegie/nobackup/scratch/tbellagio/gea_grene-net/lfmm_full/env_train_files/env_train_*.csv
#    - Environmental training data for each bioclimatic variable
#
# Output Files:
# 1. /carnegie/nobackup/scratch/tbellagio/gea_grene-net/lfmm_full/lfmm_fullresults_all_k/w_calibration_pvalue_full_genome_*.csv
#    - Calibrated p-values
# 2. /carnegie/nobackup/scratch/tbellagio/gea_grene-net/lfmm_full/lfmm_fullresults_all_k/effect_sizes_simple_full_genome_*.csv
#    - Effect sizes
#

# Required Libraries ----
# remotes::install_github("bcm-uga/lfmm")
library("lfmm")
library("dplyr")
# BiocManager::install("qvalue")
library("qvalue")
library('data.table')
library(vcfR)

# Get Command Line Arguments
# args <- commandArgs(trailingOnly = TRUE)
# 
# # Parse arguments
# samples <- args[1]  # Space-separated list of samples
# bio_name <- args[2]  # Bioclimatic variable name
# k <- as.integer(args[3])  # Number of components



# data("example.data")
# Y <- example.data$genotype

# Clean up samples argument
# samples <- gsub("\\[|\\]|'", "", samples)
# samples_vector <- strsplit(samples, " ")[[1]]
# samples_vector <- trimws(samples_vector)

# Set Up Paths ----
analysis_dir = '/global/scratch/users/arphillips/data/lfmm'
results_dir = '/global/scratch/users/arphillips/data/lfmm/lfmm_results_all_k'
aaf_dir = '/global/scratch/users/arphillips/data/rona/aaf'

# Load Data ----
# Load environmental training data
## needs to be scaled
# env_train_file <- sprintf('%s/env_train_files/env_train_%s.csv', analysis_dir, bio_name)
# env_train <- read.csv(env_train_file, header = TRUE, colClasses = "numeric")
# print(dim(env_train))

clm <- read.csv("/global/scratch/users/arphillips/data/climate/chelsa/chelsav2/GLOBAL/climatologies/1981-2010/bio/CHELSA_climate_data.1197.csv")

metadata <- read.csv("/global/scratch/projects/fc_moilab/aphillips/aspen_snakemake/metadata/megametadata.2025-08-26.csv")

# Load genotype data ----
## n = 818
## SNPs before thinning: 544,839
## SNPs after pruning with PLINK: 77,915

vcf <- read.vcfR("/global/scratch/users/arphillips/data/processed/filtered_snps/wgs_aspen.all.nocall.10dp90.diploids.pruned.vcf.gz")
gt <- extract.gt(vcf, element = "GT", as.numeric = FALSE)
gt <- t(gt) # ind by SNPs

gt[is.na(gt)] <- 0 # replacing all ./. with 0 is wrong but all I've got right now
sum(is.na(gt))

gt[gt == "1/1"] <-2
gt[gt == "0/1"] <- 1

# gt_num <- apply(gt, c(1,2), function(x) {
#   if (x == "1/1") {
#     return(2) # homozygous alt
#   }  else if (x == "0/1"){
#     return(1) # het
#   } else if (x == 0){
#     return(0) # homozygous ref
#   }
# })

dim(gt)
gt[1:5,1:5]

# nacount <- colSums(is.na(gt))
# sum(nacount == 0 )
# hist(nacount)

# gt.imp <- apply(gt, 2, function(x) replace(x, 
#                                            is.na(x), 
#                                            as.numeric( names( which.max(table(x))))))

gt.imp <- apply(gt, 2, as.numeric)

# sum(is.na(gt.imp)) # No NAs
gt.imp[1:10,1:5]
dim(gt.imp)
str(gt.imp)

## Subset env data to vcf
genos <- colnames(vcf@gt)[-1]
length(genos)

metadata$seqname <- gsub(x = metadata$bamlist,pattern = ".dedup.bam", replacement = "") 
metasub <- metadata[metadata$seqname %in% genos,]
dim(metasub)

coor <- metasub %>% # extract and fix coor
  dplyr::select(seqname, Latitude, Longitude) %>%
  unique()
dim(coor)

coor[duplicated(coor$seqname),]
cr_coor <- c(23.524897, 23.524897, -104.684526,  -104.684526)
coor[coor$seqname == "53114.3.593590.AATAGAGATA-ACGCGGCCCT",][2:3] <- cr_coor
coor <- unique(coor)
dim(coor)

# (K) Number of components ----
k <- 4  

# Prep climate data ----
clmsub <- clm %>% # subset clm
  filter(lat %in% coor$Latitude, lon %in% coor$Longitude)
dim(clmsub)

clmsub_ordered <- clmsub[match(coor$Latitude, clmsub$lat), ] # order to match genotypes
clmsub_ordered2 <- clmsub[match(coor$Longitude, clmsub$lon), ]

gt.imp.sub <- gt.imp[-which(is.na(clmsub_ordered2$lat)),]
dim(gt.imp.sub)

clmsub_ordered2.sub <- clmsub_ordered2[complete.cases(clmsub_ordered2),]
dim(clmsub_ordered2.sub)

# Run through bioclim var ----
bioclim <-  paste0("bio", 1:19)

for (i in bioclim){
## Subset climate data ----
biovar <- select(clmsub_ordered2.sub, all_of(i)) 
# %>% scale()

# hist(biovar)
# hist(clmsub_ordered2$bio1)

# Load allele frequency data
# allele_freq <- fread(delta_p_file, select = samples_vector, sep = ',', header = TRUE, check.names = FALSE)
# colnames(allele_freq) <- NULL
# allele_freq <- t(as.matrix(allele_freq))
# print(dim(allele_freq))


## Set Up Output Files ----
output_file_w_calibration_pvalue <- paste0(results_dir, '/w_calibration_pvalue_full_genome_', i, '_k', k, '.csv')
output_file_wo_calibration_pvalue <- paste0(results_dir, '/wo_calibration_pvalue_full_genome_', i, '_k', k, '.csv')
output_file_effect_sizes_simple <- paste0(results_dir, '/effect_sizes_simple_full_genome_', i, '_k', k, '.csv')

## Run LFMM Analysis ----
# Train the model
mod.lfmm <- lfmm_ridge(Y = gt.imp.sub, X = biovar, K = k)
# summary(mod.lfmm)

# Perform association testing
pv <- lfmm_test(Y = gt.imp.sub, X = biovar, lfmm = mod.lfmm, calibrate = "gif")
pv$gif

## Save Results ----
# Save p-values
p_values <- pv$pvalue
# hist(p_values)
write.csv(p_values, file = output_file_wo_calibration_pvalue)

# Save calibrated p-values
calibrated_pvalues <- pv$calibrated.pvalue
# hist(calibrated_pvalues)
write.csv(calibrated_pvalues, file = output_file_w_calibration_pvalue)

# Save effect sizes
beta_i <- pv$B
# hist(beta_i)
write.csv(beta_i, file = output_file_effect_sizes_simple)
}

# Quality plots ----
qqplot(rexp(length(calibrated_pvalues), rate = log(10)),
       -log10(calibrated_pvalues), xlab = "Expected quantile",
       pch = 19, cex = .4)
abline(0,1)

plot(-log10(calibrated_pvalues), 
     pch = 19, 
     cex = .2, 
     xlab = "SNP", ylab = "-Log P",
     col = "grey")
# points(example.data$causal.set, 
#        -log10(calibrated_pvalues)[example.data$causal.set], 
#        type = "h", 
#        col = "blue")


# Alternate Allele Frequency analysis for RONA ----
## Alternate allele frequency
bioclim <-  paste0("bio", 1:19)
# j <- bioclim[1]
for (j in bioclim){
  biovar <- dplyr::select(clmsub_ordered2.sub, all_of(j)) %>%
    scale()
  aaf_df <- as.data.frame(matrix(NA, ncol = 4, nrow = dim(gt.imp.sub)[2]))
  colnames(aaf_df) <- c("site","intercept", "slope", "p")
  for (s in 1:dim(gt.imp.sub)[2] ){
    aaf.mod <- lm(gt.imp.sub[,s] ~ biovar)
    aaf_df$site <- colnames(gt.imp.sub)[s]
    aaf_df[s,2:3] <- aaf.mod$coefficients
    aaf_df$p <- summary(aaf.mod)$coefficients[2,4]
  }
  write.csv(aaf_df, paste0(aaf_dir, "/aaf_model.",j,".csv") )
}

# dim(gt.imp.sub)
# dim(biovar)


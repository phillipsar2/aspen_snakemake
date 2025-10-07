library(tidyverse)

header <- read.table("/global/scratch/projects/fc_moilab/aphillips/aspen_snakemake/samples.in.vcf")
header <- t(header)
dim(header)
head(header)


dir <- "/global/scratch/projects/fc_moilab/aphillips/aspen_snakemake/"
samples_to_drop <- read.csv(paste0(dir, "metadata/samples_to_drop_from_vcf.txt"), header = F)
head(samples_to_drop)

samples_to_drop[,1] %in% header[,1] 

# Exclude those files from header
samples_to_keep <- header[,1][!header[,1] %in% samples_to_drop[,1]] %>%
  as.data.frame()
length(samples_to_keep)
head(samples_to_keep)
dim(samples_to_keep)

# metadata <- read.csv(paste0(dir, "metadata/megametadata.2025-08-26.csv"))
# str(metadata)
# metadata$seq <- gsub(x = metadata$bamlist, pattern = ".dedup.bam", replacement = "")

# sum(samples_to_keep[,1] %in% metadata$seq)
# 
# filter(metadata, seq %in% samples_to_keep[,1]) %>%
#   select(seq, )

# load gbs2ploidy 
ploidy <- read.csv("/global/scratch/users/arphillips/data/gbs2ploidy/flow_cyt_predictions.2025-09-30.csv")
dim(ploidy)

sum(samples_to_keep[,1] %in% ploidy$sample)
sum(!samples_to_keep[,1] %in% ploidy$sample)

ploidy_keep <- filter(ploidy, sample %in% samples_to_keep[,1]) %>%
  select(sample, ploidy_call)


issues <- samples_to_keep[,1][!samples_to_keep[,1] %in% ploidy$sample]
issues <- issues[-c(4:5)] # drop one duplicated sample

manual_ploidy <- c("diploid",
                   "triploid",
                   "triploid",
                   "diploid",
                   "diploid",
                   "diploid"
                   )
fixes <- cbind(issues, manual_ploidy)
colnames(fixes) <- c("sample", "ploidy_call")

ploidy_match <- rbind(ploidy_keep, fixes)
dim(ploidy_match)

ploidy_match$ploidy <- ifelse(ploidy_match$ploidy_call =="diploid", 2, 3)

write.csv(ploidy_match[,c(1,3)], paste0(dir, "metadata/ploidy.geno.1127.", Sys.Date(), ".csv"),row.names = F)

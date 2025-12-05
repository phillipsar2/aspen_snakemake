### Title: Making a mega metadata file
### Date: 2/24/25

library(stringr)
library(forcats)
library(dplyr)

# (1) Load and prep metadata files ----

# Load bamlist of existing bams
dir = "/global/scratch/projects/fc_moilab/aphillips/aspen_snakemake/data/bams/"
bamlist <- list.files(path = dir, pattern = "\\.dedup.bam$", )
bamlist
length(bamlist)

seq_names <- lapply(bamlist, gsub, pattern = "dedup.bam", replacement="fastq.gz") %>%
  unlist() # turn bams into fastq names
length(seq_names)

## Deal with merged bams
merged_files <- seq_names[grep(pattern = "merge", x = seq_names)] 
fixed_seqnames <- paste0(str_split(merged_files, pattern = "_", simplify = T)[,1], ".fastq.gz")
seq_names[grep(pattern = "merge", x = seq_names)] <- fixed_seqnames # replace uncropped names

#  Load meta data files
jgi_meta <- read.csv("/global/scratch/projects/fc_moilab/aphillips/aspen_snakemake/metadata/completed_sequencing_samplereport.d05132025.csv") # JGI metadata file
# jgi_meta$RQC.Seq.Unit.Name

ben_meta <- read.csv("/global/scratch/projects/fc_moilab/aphillips/aspen_snakemake/metadata/all_sites_DO_NOT_SHARE_d02062024.csv")

# Subset metadata
jgi_sub <- jgi_meta[jgi_meta$RQC.Seq.Unit.Name %in% seq_names,]
jgi_sub_ordered <- jgi_sub %>%
  arrange(fct_relevel(RQC.Seq.Unit.Name, seq_names))

dim(jgi_sub_ordered)[1] == length(seq_names)    # correctly subset?
sum(jgi_sub_ordered$RQC.Seq.Unit.Name == seq_names) # correctly ordered?

# (2) Fix bad names ----
bad <- jgi_sub_ordered$Sample.Name[!jgi_sub_ordered$Sample.Name %in% ben_meta$ID] 
length(bad)

issues <- c(            # issues in JGI file
  "POTR-DIXIE-1-178",
  "POTR-DIXIE-7B-231A",
  "POTR-DIXIE-8B-250",
  "POTR-FOCO-3-172A",
  "POTR-FOCO-4-173A",
  "POTR-SJNF-H-B-19",
  "POTR-SJNF-H-D-31",
  "POTR-SJNF-H-E-49",
  "POTR-SJNF-L-E-80",
  "POTR-SJNF-L-G-101B",
  "POTR-SJNF-M-5-288",
  "POTR-SJNF-M-C-216A",
  "POTR-SJNF-M-D-141",
  "POTR-UINTA-8-194A",
  "POTR-UNCO-11-325",
  "POTR-UNCO-3-406",
  "POTR-UNCO-6-200B",
  "POTR-WR-1-203A",
  "POTR-WR-11-237B",
  "POTR-WR-8-210B",
  "POTR-FOCO-1-170A",
  "AB1",
  "AB10",
  "AB3-2",
  "AB6",
  "ESSI-002-1",
  "ESSI-004-1",
  "ESSI-004-3",
  "FORE-004-1",
  "QB7",
  "SP-45"
)

solutions <- c(       # correct names from Ben's file
  "KK-DIXIE-1-178",
  "KK-DIXIE-7B-231A",
  "KK-DIXIE-8B-250",
  "KK-FOCO-3-172B",
  "KK-FOCO-4-173A",
  "KK-SJNF-H-B-19",
  "KK-SJNF-H-D-31",
  "KK-SJNF-H-E-49",
  "KK-SJNF-L-E-80",
  "KK-SJNF-L-G-101B",
  "KK-SJNF-M-5-288",
  "KK-SJNF-M-C-216A",
  "KK-SJNF-M-D-141",
  "KK-UINTA-8-194A",
  "KK-UNCO-11-325",
  "KK-UNCO-3-406",
  "KK-UNCO-6-200B",
  "KK-WR-1-203A",
  "KK-WR-11-237B",
  "KK-WR-8-210A",
  "KK-FOCO-1-170A",
  "AB-1", 
  "AB-10", 
  "AB-3", 
  "AB-6", 
  "ESSI-002",
  "ESSI-004",
  "ESSI-004",
  "FORE-004",
  "QB-7",
  "FP2.5-8"
)

fixes <- cbind(issues, solutions) %>% 
  as.data.frame()
fixes_sort <- fixes[match(bad, fixes$issues), ]

# Replace bad in JGI with good from ben
jgi_sub_ordered$Sample.Name[!jgi_sub_ordered$Sample.Name %in% ben_meta$ID] <- fixes_sort$solutions 

# Zero if fixed
length( jgi_sub_ordered$Sample.Name[!jgi_sub_ordered$Sample.Name %in% ben_meta$ID])

# Check for NAs in Sample.Name - zero if OK
sum(is.na(jgi_sub_ordered$Sample.Name)) 

# Meta sub (duplicates will be represented by one line)
ben_sub <- ben_meta[ben_meta$ID %in% jgi_sub_ordered$Sample.Name,]
dim(ben_sub)

# Merge 
all_meta <- merge(x = jgi_sub_ordered, y = ben_sub, by.x = "Sample.Name", by.y = "ID", all.x = T)
dim(all_meta)
sum(is.na(all_meta$Latitude))
sum(is.na(all_meta$Sample.Name))

# Add on bam names
in_names <- cbind(bamlist, seq_names)
all_meta <- merge(all_meta, in_names, by.x = "RQC.Seq.Unit.Name", by.y = "seq_names")
# all_meta$bams <- lapply(all_meta$RQC.Seq.Unit.Name, gsub, pattern = "fastq.gz", replacement = "dedup.bam") %>%
  # unlist()

# Write meta data
write.csv(all_meta,
          file = paste0("/global/scratch/projects/fc_moilab/aphillips/aspen_snakemake/metadata/megametadata.", Sys.Date(), ".csv"),
          row.names = F)

# Check for duplicate samples (sequenced multiple time)
sum(duplicated(all_meta$Sample.Name))
seq_duplicates <- all_meta$Sample.Name[duplicated(all_meta$Sample.Name)]
duplicated_samples <- all_meta[all_meta$Sample.Name %in% seq_duplicates,]
duplicated_samples <- duplicated_samples %>% distinct(bamlist, .keep_all = TRUE)
# write.table(duplicated_samples, "/global/scratch/users/arphillips/reports/filestomerge_or_resolve.08202025.txt")
# duplicated_samples <- read.table("/global/scratch/users/arphillips/reports/filestomerge_or_resolve.08202025.txt")

# dup_df <- matrix(nrow = length(unique(duplicated_samples$Sample.Name)), ncol = 4)
# dup_df[,1] <- unique(duplicated_samples$Sample.Name)
# for (i in unique(duplicated_samples$Sample.Name)){
#   dup_df[dup_df[,1] == i, 2:4] <- all_meta$RQC.Seq.Unit.Name[all_meta$Sample.Name == i]
# }

dup_list <- lapply(unique(duplicated_samples$Sample.Name), function(x){all_meta$RQC.Seq.Unit.Name[all_meta$Sample.Name == x]})
names(dup_list) <- unique(duplicated_samples$Sample.Name)
dup_list

# colnamdup_listcolnames(dup_df) <- c("Genotype", "Merge_A", "Merge_B")
# dup_df[,2:3] <- gsub(x = dup_df[,2:3], pattern = ".fastq.gz", replacement = "")
# dim(dup_df)

write.table(dup_df, "/global/scratch/users/arphillips/reports/filestomerge.08122025.txt", 
            row.names = F)

View(all_meta_stats[all_meta_stats$Sample.Name %in% seq_duplicates,] )

##################################################  
### Add quality metrics
##################################################

# (1) bamqc metrics
bamqc <- read.csv("/global/scratch/users/arphillips/reports/bamqc/stats.meta.bamqc.2025-06-11.csv")
dim(bamqc)
str(bamqc)

# (2) ploidy calls
gbs2ploidy <- read.csv("/global/scratch/users/arphillips/data/gbs2ploidy/flow_cyt_predictions2025-07-17.csv")
gbs2ploidy$ploidy_call <- as.factor(gbs2ploidy$ploidy_call)
dim(gbs2ploidy)
str(gbs2ploidy)

# (3) Percent missing SNPs per genotype (VCF)
miss <- read.table("/global/scratch/users/arphillips/reports/wgs_aspen.all.10dp90.imiss", header = T)
dim(miss)

# Merge two datasets
tmp <- merge(bamqc[,83:91], gbs2ploidy[,2:3], by.x = "seqname", by.y = "sample")
stats <- merge(tmp, miss[,c(1,5)], by.x = "seqname", by.y = "INDV")
dim(stats)

# Add file ending
stats$seqname <- paste0(stats$seqname, ".fastq.gz")

# Merge with megametadata
all_meta_stats <- merge(all_meta, stats, by.x = "RQC.Seq.Unit.Name", by.y = "seqname")
dim(all_meta_stats)

# Check for duplicate fastq files
duplicates <- all_meta_stats$RQC.Seq.Unit.Name[duplicated(all_meta_stats$RQC.Seq.Unit.Name)]
all_meta_stats <- all_meta_stats[!duplicated(all_meta_stats$RQC.Seq.Unit.Name),]
dim(all_meta_stats)

# Write megametadata_quality file

write.csv(all_meta_stats,
          file = paste0("/global/scratch/projects/fc_moilab/aphillips/aspen_snakemake/metadata/megametadata_quality.", Sys.Date(), ".csv"),
          row.names = F)

##################################################  
### Files to still download
##################################################

to_download_df <- jgi_meta[!jgi_meta$RQC.Seq.Unit.Name %in% seq_names,]
to_download_compl <- to_download_df[to_download_df$Sequencing.Project.Status == "Complete",]
dim(to_download_compl)

seq_list <- to_download_compl %>%
  select(RQC.Seq.Unit.Name, Sample.Name) %>% # subset columns
  filter(RQC.Seq.Unit.Name != "") # some samples don't have a fastq name listed even though complete

dim(seq_list)

write.csv(seq_list, paste0("/global/scratch/projects/fc_moilab/aphillips/aspen_snakemake/metadata/seq_to_download.",Sys.Date(),".csv"), row.names = F)


##################################################  
### CA aspen
##################################################
library(stringr)
library(ggplot2)
all <- read.csv("/global/scratch/projects/fc_moilab/aphillips/aspen_snakemake/metadata/megametadata_quality.2025-08-11.csv")

bamlist <- read.csv("/global/home/users/wandif/metadata/california.bamlist.txt", F)
bamlist$bam <- str_split(bamlist$V1, pattern = "/", simplify = T)[,10]
head(bamlist)

ca <- all[all$bams %in%  bamlist$bam,]
dim(ca)

sum(ca$bams == bamlist$bam)
head(ca$bams)
head(bamlist$bam)


# Assign population IDs
ca$Sample.Name
ca[grep(x = ca$Sample.Name, pattern = "NVO"),c('Latitude', 'Longitude')]

pops <- str_split(ca$Sample.Name, pattern = "[-_]", simplify = T)[,1]
str(pops)
unique(pops) %>% length()

ca$population_ID <- pops
ca$population_ID[ca$population_ID %in% c('CR', 'WR')] <- "CR-WR"

ggplot(ca, aes(x = Longitude, y = Latitude, color = population_ID)) +
  geom_point()

write.csv(ca, "/global/scratch/projects/fc_moilab/aphillips/aspen_snakemake/metadata/california.megametadata.2025-12-03.csv" )


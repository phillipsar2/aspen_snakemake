### Title: Making a mega metadata file
### Date: 2/24/25

library(stringr)
library(forcats)

# Bamlist of existing bams
dir = "/global/scratch/users/arphillips/data/interm/mark_dups/"
bamlist <- list.files(path = dir, pattern = "\\.dedup.bam$", )
bamlist

seq_names <- lapply(bamlist, gsub, pattern = "dedup.bam", replacement="fastq.gz") %>%
  unlist() # turn bams into fastq names
length(seq_names)

# Meta data
jgi_meta <- read.csv("/global/scratch/projects/fc_moilab/aphillips/aspen_snakemake/metadata/completed_sequencing_samplereport_11122024.csv") # JGI metadata file
# jgi_meta$RQC.Seq.Unit.Name

ben_meta <- read.csv("/global/scratch/projects/fc_moilab/aphillips/aspen_snakemake/metadata/all_sites_DO_NOT_SHARE_d02062024.csv")

# Subset metadata
jgi_sub <- jgi_meta[jgi_meta$RQC.Seq.Unit.Name %in% seq_names,]
jgi_sub_ordered <- jgi_sub %>%
  arrange(fct_relevel(RQC.Seq.Unit.Name, seq_names))

dim(jgi_sub_ordered)[1] == length(seq_names)    # correctly subset?
sum(jgi_sub_ordered$RQC.Seq.Unit.Name == seq_names) # correctly ordered?

# Fix bad names
bad <- jgi_sub_ordered$Sample.Name[!jgi_sub_ordered$Sample.Name %in% ben_meta$ID] 
length(bad)

issues <- c(
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
  "POTR-FOCO-1-170A"
)

solutions <- c(
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
  "KK-FOCO-1-170A"
)

fixes <- cbind(issues, solutions) %>% 
  as.data.frame()
fixes_sort <- fixes[match(bad, fixes$issues), ]

# Replace bad in JGI with good from Ben
jgi_sub_ordered$Sample.Name[!jgi_sub_ordered$Sample.Name %in% ben_meta$ID] <- fixes_sort$solutions 

# Zero if fixed
length(jgi_sub_ordered$Sample.Name[!jgi_sub_ordered$Sample.Name %in% ben_meta$ID])

# Meta sub (duplicates will be represented by one line)
ben_sub <- ben_meta[ben_meta$ID %in% jgi_sub_ordered$Sample.Name,]
dim(ben_sub)

# Merge 
all_meta <- merge(x = jgi_sub_ordered, y = ben_sub, by.x = "Sample.Name", by.y = "ID", all.x = T)
# all_meta_ordered <- all_meta %>%
  # arrange(fct_relevel(RQC.Seq.Unit.Name, seq_names))
# dim(all_meta_ordered)
dim(all_meta)

# Add on bam names
all_meta$bams <- lapply(all_meta$RQC.Seq.Unit.Name, gsub, pattern = "fastq.gz", replacement = "dedup.bam") %>%
  unlist()

# Write meta data
write.csv(all_meta,
          file = paste0("/global/scratch/projects/fc_moilab/aphillips/aspen_snakemake/metadata/megametadata.", Sys.Date(), ".csv"),
          row.names = F)

# Title: BAMQC Stats for the WGS data
# Author: Alyssa
# Date: 4/24/25

library(stringr)
library(ggplot2)
library(dplyr)

# Load files and prep
dir = "/global/scratch/users/arphillips/reports/bamqc"
qc_files <- list.files(path = dir, pattern = ".bamqc.txt")
stats <- do.call(cbind, lapply(paste0(dir, "/",qc_files), function(f) read.table(f, header=FALSE)))
# stats <- read.table("~/aspen/qc/bamqc/stats.bamqc.txt", header = F)

colnames(stats) <- gsub(qc_files, pattern = ".bamqc.txt",replacement = "")
# colnames(stats) <- c("bams", "numreads", "numdups", "medianinsertsize", "GC", "meancoverage", "meanMQ", "perreadsmapped")

tmp9 <- str_split(stats$bamlist, pattern = "/", simplify = T)[,9]
tmp10 <- str_split(stats$bamlist, pattern = "/", simplify = T)[,10]
stats$bams <- ifelse(tmp9 == "bams", tmp10, tmp9)
stats$seqname <- str_split(stats$bams, pattern = ".dedup", simplify = T)[,1]

stats[,c("numdups","numreads")] <- lapply(stats[,c("numdups","numreads")], function(x) gsub(",", "", x))
stats$GC <- gsub("%", "", stats$GC)
stats$meancoverage <- gsub("X", "", stats$meancoverage)
stats$perreadsmapped <- gsub("%)", "", stats$perreadsmapped)

stats[,2:8] <- lapply(stats[,2:8], as.numeric)
stats$percentdups <- stats$numdups/stats$numreads * 100

head(stats)
str(stats)
dim(stats)

# Load metadata files
# Megametadata file from scripts/megametadata.R
megameta <- read.csv("/global/scratch/projects/fc_moilab/aphillips/aspen_snakemake/metadata/megametadata.2025-04-24.csv")
str(megameta)
dim(megameta)

meta_sub <- megameta[megameta$bams %in% stats$bams,]
dim(meta_sub)[1] == dim(stats)[1]

# Merge fixed metadata files
stats_meta <- merge(x=meta_sub, y=stats, by.x = "bams", by.y = "bams")
dim(stats_meta)

# Overview of data quality

pdf(paste0(dir, "/bamqc.histograms.",Sys.Date(),".pdf"))

stats %>%
  ggplot(aes(x = numreads/1000000)) +
  geom_density(fill="#69b3a2", color="#e9ecef", alpha=0.9, adjust = 0.8) +
  # xlim(c(0,max(geno_dp, na.rm = T))) +
  theme_bw() +
  labs(x = "Total number of reads (millions)")

stats %>%
  ggplot(aes(x = percentdups)) +
  geom_histogram(fill="#69b3a2", color="#e9ecef", alpha=0.9) +
  # xlim(c(0,max(geno_dp, na.rm = T))) +
  theme_bw() +
  labs(x = "Percent of duplicated reads")

stats %>%
  ggplot(aes(x = medianinsertsize)) +
  geom_histogram(fill="#69b3a2", color="#e9ecef", alpha=0.9) +
  geom_vline(xintercept = 300, color = "black") +
  theme_bw() +
  labs(x = "Median insert size")

stats %>%
  ggplot(aes(x = meancoverage)) +
  geom_histogram(fill="#69b3a2", color="#e9ecef", alpha=0.9) +
  theme_bw() +
  geom_vline(xintercept = mean(stats$meancoverage), color = "black") +
  labs(x = "Average coverage of each genotype")

stats %>%
  ggplot(aes(x = GC)) +
  geom_histogram(fill="#69b3a2", color="#e9ecef", alpha=0.9) +
  theme_bw() +
  labs(x = "Average GC content")
# GC content can be affected by duplicates, bacterial contamination, etc.

stats %>%
  ggplot(aes(x = meanMQ)) +
  geom_histogram(fill="#69b3a2", color="#e9ecef", alpha=0.9) +
  theme_bw() +
  labs(x = "Mean mapping quality (MQ)")

stats %>%
  ggplot(aes(x = perreadsmapped)) +
  geom_histogram(fill="#69b3a2", color="#e9ecef", alpha=0.9) +
  theme_bw() +
  labs(x = "Percent reads mapped")

dev.off()

# Comparing the quality statistics
pairs(stats[,c(2:8,11)], upper.panel = NULL)
cor(stats[,c(2:8,11)])

'Nothing is really standing out as too correlated. Not useful for figuring out what is happening with GC content.'

# Looking at characters of the submitted data
str(stats_meta)
unique(stats_meta$Sample.type)

pdf(paste0(dir,"/bamqc.stats_by_sampletype.",Sys.Date(),".pdf"))

stats_meta %>%
  ggplot(aes(x = Sample.type, y = meancoverage, fill = Sample.type)) +
  geom_violin() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(fill="none") +
  labs(y = "Mean coverage", x = "Sample type")

stats_meta %>%
  ggplot(aes(x = Sample.type, y = GC, fill = Sample.type)) +
  geom_violin() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(fill="none") +
  labs(y = "GC content", x = "Sample type")

# When DNA is highly fragmented, GC content will be higher as they are more stable 
stats_meta %>%
  ggplot(aes(x = Sample.type, y = medianinsertsize, fill = Sample.type)) +
  geom_violin() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(fill="none") +
  labs(y = "Median Insert Size", x = "Sample type")

stats_meta %>%
  ggplot(aes(x = Sample.type, y = percentdups, fill = Sample.type)) +
  geom_violin() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(fill="none") +
  labs(y = "% dups", x = "Sample type")

stats_meta %>%
  ggplot(aes(x = Sample.type, y = perreadsmapped, fill = Sample.type)) +
  geom_violin() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(fill="none") +
  labs(y = "% reads mapped", x = "Sample type")

dev.off()

# Data by plate
pdf(paste0(dir, "/bamqc.stats_by_plate.",Sys.Date(),".pdf"))

stats_meta %>%
  ggplot(aes(x = as.factor(Sample.Plate.Tube.Name) , y = meancoverage)) +
  geom_violin() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(fill="none") +
  labs(y = "Mean coverage", x = "Sequencing Project ID / Plate")

stats_meta %>%
  ggplot(aes(x = as.factor(Sample.Plate.Tube.Name) , y = GC)) +
  geom_violin() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(fill="none") +
  labs(y = "GC content", x = "Plate")

stats_meta %>%
  ggplot(aes(x = Sample.Plate.Tube.Name, y = medianinsertsize)) +
  geom_violin() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(fill="none") +
  labs(y = "Median Insert Size", x = "Plate")

stats_meta %>%
  ggplot(aes(x = Sample.Plate.Tube.Name, y = percentdups)) +
  geom_violin() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(fill="none") +
  labs(y = "% dups", x = "Plate")

stats_meta %>%
  ggplot(aes(x = Sample.Plate.Tube.Name, y = perreadsmapped)) +
  geom_violin() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(fill="none") +
  labs(y = "% reads mapped", x = "Plate")

dev.off()

# Herbarium sample year by date
herb_date <- stats_meta %>%
  filter(Sample.type == "Herbarium sheet") %>%
  ggplot(aes(x = Year, y = medianinsertsize)) +
  geom_point()
herb_date

ggsave(plot = herb_date, filename = paste0(dir, "/herbarium_medianinsertsize_by_date", Sys.Date(), ".jpeg"), units = "in", width = 4, height = 4)

# Samples too low-quality to continue working with
exclude_20 <- stats_meta %>%
  filter(perreadsmapped < 90 | meancoverage < 20) %>%
  select(Sequencing.Project.ID,Sample.type,Sample.Plate.Tube.Name,numreads,medianinsertsize,GC,meancoverage,meanMQ,perreadsmapped,percentdups)

write.csv(exclude_20, file = paste0(dir,"/sample_to_redo.", Sys.Date(), ".MQ20.PR90.csv"))

exclude_10 <- stats_meta %>%
  filter(perreadsmapped < 90 | meancoverage < 10) %>%
  select(Sequencing.Project.ID,Sample.type,Sample.Plate.Tube.Name,numreads,medianinsertsize,GC,meancoverage,meanMQ,perreadsmapped,percentdups)

write.csv(exclude_10, file = paste0(dir, "/sample_to_redo.", Sys.Date(), ".MQ10.PR90.csv"))

# Write full stats file out 
out <- stats_meta %>%
  select(Sequencing.Project.ID,Sample.Name,bams,Sample.type,Sample.Plate.Tube.Name,
         numreads,medianinsertsize,GC,meancoverage,meanMQ,perreadsmapped,percentdups)
write.csv(out, file = paste0(dir,"/stats.meta.bamqc.", Sys.Date(), ".csv"), row.names = F)

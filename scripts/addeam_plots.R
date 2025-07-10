## Title: Parsing AdDeam results
## Author: Alyssa Phillips
## Date: 7/10/25

library(dplyr)
library(stringr)
library(ggplot2)
library(tidyr)

megameta <- read.csv("/global/scratch/projects/fc_moilab/aphillips/aspen_snakemake/metadata/megametadata.2025-06-09.csv")
dim(megameta)
str(megameta)
unique(megameta$Sample.Isolated.From)


meta_sub <- filter(megameta, Sample.Isolated.From %in% c("Herbarium sheet", "Dried leaf + immediate storage in silica, subsequent powdering + freezing"))
dim(meta_sub)
# bamlist <- paste0("/global/scratch/projects/fc_moilab/aphillips/aspen_snakemake/data/bams/", meta_sub$bams)
bamlist <- gsub(meta_sub$bams, pattern = ".bam", replacement = "")

profilesDir = "/global/scratch/users/arphillips/reports/addeam/profiles"

files <- lapply(bamlist, function(x) list.files(path = profilesDir, pattern = x) ) %>% 
  unlist()

files_3p <- files[grep(pattern = "3p", files)]
files_5p <- files[grep(pattern = "5p", files)]


# C > T from '5 end
mat_5p <- matrix(NA, nrow = length(files_5p), ncol = 6) %>%
  as.data.frame()
colnames(mat_5p) <- c("sample", "0", "1", "2", "3", "4")
mat_5p$sample <- files_5p

for (i in 1: length(files_5p)){
  tmp <- read.csv(paste0(profilesDir,"/", files_5p[i]), sep = "\t") %>%
    select("C.T")
  mat_5p[i,2:6] <- t(tmp)
}

mat_5p$sample_type <- meta_sub$Sample.Isolated.From

# mat_5p_long <- melt(setDT(mat_5p), id.vars = c("sample","Sample.Isolated.From"), variable.name = "position")
mat_5p_long <- pivot_longer(mat_5p, cols = `0`:`4`, names_to = "position", values_to = "frequency")

ggplot(mat_5p_long, aes(x = position, y = frequency, color = sample_type, group = sample)) +
  geom_line() +
  guides(color="none") +
  theme_bw() +
  xlab("Position from 5' end") +
  ylab("C>T substitution frequency")

ggsave(paste0("/global/scratch/users/arphillips/reports/addeam/plots/", "CtoT_sub_herbaria_and_freeze_thaw_samples.jpeg"), 
       height = 4, width = 6, unit = "in")

# G > A from 3' end
mat_3p <- matrix(NA, nrow = length(files_3p), ncol = 6) %>%
  as.data.frame()
colnames(mat_3p) <- c("sample", "0", "1", "2", "3", "4")
mat_3p$sample <- files_3p

for (i in 1: length(files_3p)){
  tmp <- read.csv(paste0(profilesDir,"/", files_3p[i]), sep = "\t") %>%
    select("G.A")
  mat_3p[i,2:6] <- t(tmp)
}

mat_3p$sample_type <- meta_sub$Sample.Isolated.From

mat_3p_long <- pivot_longer(mat_3p, cols = `0`:`4`, names_to = "position", values_to = "frequency")

ggplot(mat_3p_long, aes(x = position, y = frequency, color = sample_type, group = sample)) +
  geom_line() +
  guides(color="none") +
  theme_bw() +
  xlab("Position from 3' end") +
  ylab("G>A substitution frequency")

ggsave(paste0("/global/scratch/users/arphillips/reports/addeam/plots/", "GtoA_sub_herbaria_and_freeze_thaw_samples.jpeg"), 
       height = 4, width = 6, unit = "in")


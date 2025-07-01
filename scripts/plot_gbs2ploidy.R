### Title: Plot gbs2ploidy results
### Author: Alyssa Phillips
### Date: 6/12/2025

library(ggplot2)
library(dplyr)
library(tidyr)

# Import propOut values ----
dir <- "/global/scratch/users/arphillips/data/gbs2ploidy/"

# propOut_df <- read.csv("/global/scratch/projects/fc_moilab/aphillips/spectral_aspen/data/gbs2ploidy/flowcyt_predictions.allsnps.propOut.csv")
propOut_files <- Sys.glob("/global/scratch/users/arphillips/data/gbs2ploidy/*.propOut.csv")
propOut_list <- lapply(propOut_files, function(x) read.table(x, sep = ",", header = T))
propOut_df <- do.call(rbind, propOut_list)

head(propOut_df)
str(propOut_df)
dim(propOut_df)
length(unique(propOut_df$sample))

# Plot propOut per sample ----
pdf(paste0(dir, "flowcyt_predictions.plots.pdf"))
par(mfrow=c(4,3))
propOut_df %>%
  dplyr::filter(quantile == '75%') %>%
  ggplot(aes(x = as.factor(allelic_ratio), y = proportion)) +
  geom_point() +
  facet_wrap( ~ sample, nrow = 9) +
  theme_bw()+
  xlab("Allelic ratio") +
  ylab("75th quantile of the posterior distribuion of allelic proportions")
dev.off()


## Inspect allelic proportions for the first nine individuals
# geno_names <- colnames(ad)
# pdf("/global/scratch/projects/fc_moilab/aphillips/spectral_aspen/data/nquack/gbs2ploidy_alleleratios.pdf")
# par(mfrow=c(3,3))
# for(i in 1:9){ # change number from 1 to 9
#   plot(propOut[[i]][,5], ylim=c(0,1), axes=FALSE,
#        xlab="allelic ratios",
#        # ylab="75th quantile of the posterior distribuion of allelic proportions",
#        ylab="proportions",
#        main = geno_names[i])
#   axis(1, at = 1:3,c("1:2","1:1","2:1"))
#   axis(2)
#   box()
#   segments(1:5, propOut[[i]][,1], 1:5,propOut[[i]][,5])
#   # title(main=paste("true ploidy =",dat[[3]][i]))
# }
# dev.off()


# What is the winner ----
winner_prop <- propOut_df %>%
  group_by(sample, quantile) %>%
  filter(proportion == max(proportion)) %>%
  arrange(sample, quantile)

winner_num <- winner_prop %>%
  group_by(sample) %>%
  count(allelic_ratio, name = 'Count')

winner_num$ploidy_call <- ifelse(winner_num$allelic_ratio == '0.5', yes = 'diploid', no = 'triploid')
head(winner_num)

## Check for duplicates
dim(winner_num)[1] == length(unique(propOut_df$sample))

duplicates <- winner_num$sample[duplicated(winner_num$sample)]

filter(winner_num, sample %in% duplicates)
filter(propOut_df, sample %in% duplicates) %>% 
  arrange(quantile)

## Keep samples with over 4 wins
winner_num <- winner_num[winner_num$Count >= 4,]

# Assess quality stats ----

## Read in metadata from bamqc
meta <- read.csv("/global/scratch/users/arphillips/reports/bamqc/stats.meta.bamqc.2025-06-11.csv")
meta$seqnames <- gsub(meta$bams, pattern = ".dedup.bam", replacement = "") %>% sort()

winner_df <- merge(x = winner_num, y = meta, by.x = "sample", by.y = "seqnames")
dim(winner_df)

## Plot winner vs depth
ggplot(winner_df, aes(x = ploidy_call, y = meancoverage)) +
  geom_boxplot()

## PCA & DA to assign ploidies
# library(gbs2ploidy)
# pout <- estploidy(alphas = propOut_list,
#                   het = rep(0.15, 9),
#                   depth = dp[1:9],
#                   nclasses = 2, # number of cytotypes expected
#                   ids = geno_names[1:9], pcs = 1:2)
# 
# plot(pout$pcscrs[,1] ~ pout$pcscrs[,2])

# Plot ploidy on map ----
## Map
library(raster)
library(sf)
library(ggplot2)
library(ggspatial)
library(USAboundaries)
library(purrr)
library(rnaturalearth)
library(rnaturalearthdata)

north_america <- ne_countries(continent = "North America", returnclass = "sf")
# states <- us_states()

events_sf <- winner_df %>% 
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326) 
dim(events_sf)

ploidy_map <- ggplot() + 
  geom_sf(data = north_america, size = 4, color = "black", fill = NA) +
  geom_sf(data = events_sf, size = 3, aes(color = ploidy_call), alpha = 0.4) +
  theme_bw() +
  theme(panel.grid = element_blank(), 
        axis.line = element_blank(),
  ) + 
  ylim(c(10,70)) +
  xlim(c(165, 57)) +
  scale_color_manual(values = c("blue", "red") ) +
  labs(color = "Ploidy")

ggsave(ploidy_map, filename = "/global/scratch/projects/fc_moilab/aphillips/aspen_snakemake/figures/ploidymap.1206.pdf", 
       width = 6, height = 5, unit = "in")

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

# Write class to file
write.csv(winner_num[,c(1,4)], paste0(dir, "flow_cyt_predictions.", Sys.Date(),".csv"))

# winner_num <- read.csv(paste0(dir, "flow_cyt_predictions.csv"))

# Compare calls to sequence data quality to assess bias ----
meta_dir <- "/global/scratch/users/arphillips/reports/filtering/depth"

## Average depth at all sites
dp <- read.csv(paste0(meta_dir, "/meandepthperind.txt.idepth"), sep = "\t")

## Missing data at all sites 
miss <- read.csv(paste0(meta_dir, "/missingdataperind.txt.imiss"), sep = "\t")

## Number of heterozygous sites used
nsites <- read.csv(paste0(meta_dir, "/nsites.hets.csv"))

## Mean reference read depth
refdp <- read.csv(paste0(meta_dir, "/refdepth.hets.csv"))

## Mean alternate read depth
altdp <- read.csv(paste0(meta_dir, "/altdepth.hets.csv"))

## Combine into one dataframe with winner_num
hetstats <- cbind(nsites, refdp[,2], altdp[,2])
colnames(hetstats) <- c("sample", "nsites", "refdp", "altdp")

allstats <- cbind(miss[,c(1,4:5)], dp[,2:3])

stats_df <- merge(x = hetstats, y = allstats, by.x = "sample", by.y = "INDV")
dim(stats_df)

winner_df <- merge(x = winner_num, y = stats_df, by = "sample")
dim(winner_df)
str(winner_df)

winner_df$dp = winner_df$altdp + winner_df$refdp

# Plots ----
# Plot winner vs depth
ggplot(winner_df, aes(x = ploidy_call, y = winner_df$dp )) +
  geom_boxplot()

hist(winner_df$nsites)

plot(winner_df$nsites, winner_df$altdp,
     col = as.factor(winner_df$ploidy_call),
     xlab = "Number of sites",
     ylab = "Mean alternate read depth")

plot(winner_df$nsites, winner_df$dp,
     col = as.factor(winner_df$ploidy_call),
     xlab = "Number of sites",
     ylab = "Mean read depth")

pdf("/global/scratch/projects/fc_moilab/aphillips/aspen_snakemake/figures/gbs2ploidy_qualitymetrics.pdf")
plot(y = winner_df$altdp, x = winner_df$refdp,
     col = as.factor(winner_df$ploidy_call),
     pch = ifelse( winner_df$nsites < 10000, 20, 23),
     xlab = "Mean alternate read depth",
     ylab = "Mean reference read depth")
abline(1,1)
dev.off()

# Ploidy cannot be dtermined for samples with less than 10K sites
winner_df$ploidy_call[winner_df$nsites < 10000] <- "unknown"
write.csv(winner_df[,c(1,3)], paste0(dir, "flow_cyt_predictions", Sys.Date(),".csv"))

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
states <- us_states()

meta <- read.csv("/global/scratch/projects/fc_moilab/aphillips/aspen_snakemake/metadata/megametadata.2025-06-09.csv")
meta$seqname <- gsub(x = meta$RQC.Seq.Unit.Name, pattern = ".fastq.gz", replacement = "")
winner_meta <- merge(winner_df, meta, by.x = "sample", by.y = "seqname")

events_sf <- winner_meta %>% 
  # filter(seqname %in% ca_samples[,1]) %>%
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326) 
dim(events_sf)

ploidy_map <- ggplot() + 
  geom_sf(data = north_america, size = 4, color = "black", fill = NA) +
  geom_sf(data = states, size = 4, fill = NA, color = "black") +
  geom_sf(data = events_sf, size = 2, aes(color = ploidy_call, shape = ploidy_call), alpha = 0.7) +
  theme_bw() +
  theme(panel.grid = element_blank(), 
        axis.line = element_blank(),
  ) + 
  ylim(c(10,70)) +
  xlim(c(165, 57)) +
  # ylim(c(30,50)) +
  # xlim(c(130, 110)) +
  scale_color_manual(values = c("blue", "red","black") ) +
  labs(color = "Ploidy") +
  guides(shape = "none")

ploidy_map

# ggsave(ploidy_map, filename = "/global/scratch/projects/fc_moilab/aphillips/aspen_snakemake/figures/ploidymap.CA.jpg",
       # width = 4, height = 6, unit = "in")
ggsave(ploidy_map, filename = "/global/scratch/projects/fc_moilab/aphillips/aspen_snakemake/figures/ploidymap.1206.pdf", 
       width = 6, height = 5, unit = "in")

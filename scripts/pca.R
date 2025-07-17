# Title: WGS aspen PCA
# Author: Alyssa Phillips
# Date: 1/28/25

# install.packages('ldsep')
library(ldsep)
library(vcfR)
library(stringr)
library(ggplot2)
library(dplyr)
library(viridis)

# (1) Load genotype matrices from updog ----
dir = "/global/scratch/users/arphillips/data/updog/"
files <- list.files(path = dir, pattern = "updog.genomat.diploid*")

# dir = "/global/scratch/users/arphillips/data/updog/"
# files_dip <- list.files(path = dir, pattern = "updog.genomat.diploid*")
# files_trip <- list.files(path = dir, pattern = "updog.genomat.triploid*")

genos <-  lapply(paste0(dir, files), read.table)
genos <- dplyr::bind_rows(genos)
dim(genos)
head(genos)

paste0("Starting SNPs: ",dim(genos)[1])

# genos_dip <-  lapply(paste0(dir, files_dip), read.table)
# genos_dip <- dplyr::bind_rows(genos_dip)
# dim(genos_dip)
# 
# genos_dip_sub <- genos_dip[,1:87]

# genos_trip <-  lapply(paste0(dir, files_trip), read.table)
# genos_trip <- dplyr::bind_rows(genos_trip)
# dim(genos_trip)
# 
# genos_trip_sub <- genos_trip[,88:175]

# (1b) read in genotypes from bcftools filtered vcfs ----
# vcf_files <- Sys.glob("/global/scratch/users/arphillips/data/processed/filtered_snps/wgs_aspen.*.nocall.10dp90.vcf")
# gt_list <- lapply(vcf_files, function(x){
#   vcf <- read.vcfR(x, verbose = F)
#   extract.gt(vcf)
# } )
# gt_mat <- do.call(rbind, gt_list)
# dim(gt_mat)

vcf <- read.vcfR("/global/scratch/users/arphillips/data/processed/filtered_snps/wgs_aspen.all.nocall.10dp90.vcf.gz")
gt <- extract.gt(vcf)
rm(vcf)

# (2) Load in metadata ----
meta <- read.csv("/global/scratch/projects/fc_moilab/aphillips/aspen_snakemake/metadata/megametadata.2025-06-09.csv")
str(meta)
dim(meta)

# cov <- read.table("/global/scratch/users/arphillips/reports/bamqc/stats.bamqc.txt", header = F) %>%
#   select(V1, V7)
# names(cov) <- c("bams", "cov")
# # cov$cov <- gsub(cov$cov, pattern = "X", replacement = "") %>% as.numeric
# cov$bams <- str_split(cov$bams, pattern = "/", simplify=T)[,9] 
# head(cov)
# 
# meta_cov <- merge(meta, cov, by = "bams", all = T )
# head(meta_cov)
# dim(meta_cov)

# (3) Prep the data ----
## Merge dataframes to shared sites
# sum(rownames(genos_dip_sub) %in% rownames(genos_trip_sub)) # number of sites
# dip_sub <- genos_dip_sub[rownames(genos_dip_sub) %in% rownames(genos_trip_sub),]
# trip_sub <- genos_trip_sub[rownames(genos_trip_sub) %in% rownames(dip_sub),]
# 
# genos <- cbind(dip_sub, trip_sub) %>% as.matrix()
# dim(genos)
# 
# paste0("Starting SNPs: ", dim(genos)[1])

## Missing data across sites
is.na(genos) %>% sum

# hist(as.matrix(dip_sub))
# hist(as.matrix(trip_sub))
hist(as.matrix(genos))

# (4) LD thin ----
## Thin SNPs by one every 500 bp
# genos <- genos[grep("Chr02", rownames(genos)), ] # subset to one grom

pos <- str_split(rownames(genos), pattern = "_", simplify = T) %>%
  as.data.frame()
colnames(pos) <- c("chr", "pos")

pos$int <- cut(as.numeric(pos[,2]),
    seq( from = min(as.numeric(pos[,2])),
         to = max(as.numeric(pos[,2])),
         by = 500 )
    )

keep_pos <- pos %>%
  group_by(chr, int) %>% # group by the chr and interval, 4743 intervals
  slice_sample(n = 1) %>% # randomly select one
  mutate(chr_pos = paste0(chr, "_", pos)) # remerge names

genos_thin <- filter(as.data.frame(genos), row.names(genos) %in% keep_pos$chr_pos)
paste0("Number of SNPs after thinning: ", dim(genos_thin)[1])

# (5) MAF filter ----
## convert genotypes if from vcf not updog
# genos_thin_cont <- apply(genos_thin, c(1,2), function(x) {
#   if (is.na(x)){      # is.na needs to be the first test
#     "<NA>"               # assign to homozygous ref and remove the first alleles
#   } else if (x == "0/0"){
#     "0" 
#   } else if (x == "0/1"){
#     "1"
#   } else if (x == "1/1"){
#     "2"}})

## not right for mixed-ploidy pops
hist(as.matrix(genos_thin)[,2])
q2 <- rowSums(genos_thin == 0, na.rm = T) / dim(genos_thin)[2] # q^2, minor allele
af <- sqrt(q2)
hist(as.matrix(af))

genos_maf <- genos_thin[af > 0.05 & af < 0.95,]
dim(genos_maf)
paste0("Number of SNPs after MAF & invariant site filter: ", dim(genos_maf)[1])

# (5) Run PCA ----
pca_out <- prcomp(genos_maf, scale. = T, center = T)
pve <- summary(pca_out)$importance[2,1:10] *100

# (6) Plot PCA ----
## Make sure meta is ordered correctly
geno_names <- gsub(x = colnames(genos_maf), pattern = "X", replacement = "" ) 
meta_names <- gsub(meta_cov$bams, pattern = ".dedup.bam", replacement = "") %>% sort()

meta_cov <- dplyr::arrange(meta_cov, bams)
# meta_cov$ploidy <- c(rep("diploid", 87), rep("triploid", 88))

tail(geno_names)
tail(meta_names)

pca_df <- data.frame(meta_cov,
                     pca_out$rotation[,1:10])
head(pca_df)
# write.csv(pca_df, "~/test.csv", col.names = T, row.names = F)
# read.csv("~/test.csv")

# ggsave("/global/scratch/projects/fc_moilab/aphillips/aspen_snakemake/data/pca/pca_by_latitude.jpeg",
       # width = 6, height = 4, unit = "in")

pdf(paste0("/global/scratch/users/arphillips/data/pca/updog.diploidgenotypes.pca.", Sys.Date(), ".pdf"))

ggplot(pca_df, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = Latitude)) +
  theme_bw() +
  # guides(color="none") +
  xlab(paste0("PC1 (", round(pve[1],2), "%)")) +
  ylab(paste0("PC2 (", round(pve[2],2), "%)")) +
  scale_color_viridis()

ggplot(pca_df, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = Longitude)) +
  theme_bw() +
  # guides(color="none") +
  xlab(paste0("PC1 (", round(pve[1],2), "%)")) +
  ylab(paste0("PC2 (", round(pve[2],2), "%)"))+
  scale_color_viridis(option = "B")

ggplot(pca_df, aes(x = PC1, y = PC2, color = cov)) +
  geom_point() +
  theme_bw() +
  xlab(paste0("PC1 (", round(pve[1],2), "%)")) +
  ylab(paste0("PC2 (", round(pve[2],2), "%)")) +
  scale_color_continuous(name = "Mean coverage") +
  scale_color_viridis()

ggplot(pca_df, aes(x = PC1, y = PC2, color = Cytotype)) +
  geom_point() +
  theme_bw() +
  xlab(paste0("PC1 (", round(pve[1],2), "%)")) +
  ylab(paste0("PC2 (", round(pve[2],2), "%)")) +
  scale_color_viridis(discrete = TRUE)

ggplot(pca_df, aes(x = PC1, y = PC2, color = Year)) +
  geom_point(alpha = 0.7) +
  theme_bw() +
  xlab(paste0("PC1 (", round(pve[1],2), "%)")) +
  ylab(paste0("PC2 (", round(pve[2],2), "%)")) +
  scale_color_viridis()

ggplot(pca_df, aes(x = PC1, y = PC2, color = Investigators)) +
  geom_point(alpha = 0.7) +
  theme_bw() +
  xlab(paste0("PC1 (", round(pve[1],2), "%)")) +
  ylab(paste0("PC2 (", round(pve[2],2), "%)"))

ggplot(pca_df, aes(x = PC1, y = PC2, color = Sample.Isolated.From)) +
  geom_point(alpha = 0.7) +
  theme_bw() +
  xlab(paste0("PC1 (", round(pve[1],2), "%)")) +
  ylab(paste0("PC2 (", round(pve[2],2), "%)")) + 
  theme(legend.title = element_text(size = 3), 
        legend.text = element_text(size = 3))

dev.off()

## Plot on map
ggplot(pca_df, aes(y = Latitude, x = Longitude, color = PC1)) +
  geom_point(alpha = 0.7) +
  theme_bw() +
  scale_color_viridis()

pca_df$Longitude[pca_df$PC1 < 0.07]

ggplot(pca_df, aes(y = Latitude, x = Longitude, color = PC2)) +
  geom_point(alpha = 0.7) +
  theme_bw() +
  scale_color_viridis(option = "C")


####### PLINK PCA ----
library(tidyverse)
library(ggpubr)
library(viridis)

# pca <- read_table("/global/scratch/users/arphillips/data/pca/wgs_aspen.all.10dp90.thin.eigenvec", col_names = FALSE)
# eigenval <- scan("/global/scratch/users/arphillips/data/pca/wgs_aspen.all.10dp90.thin.eigenval")

## CA PCA
pca <- read_table("/global/scratch/users/arphillips/data/pca/wgs_aspen.CA.10dp90.thin.eigenvec", col_names = FALSE)
eigenval <- scan("/global/scratch/users/arphillips/data/pca/wgs_aspen.CA.10dp90.thin.eigenval")

# sort out the pca data
# remove nuisance column
pca <- pca[,-1]
# set names
names(pca)[1] <- "ind"
names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))

# Sort metadata
# (2) Load in metadata ----
meta <- read.csv("/global/scratch/projects/fc_moilab/aphillips/aspen_snakemake/metadata/megametadata.2025-06-09.csv")
# meta <- read.csv("/global/scratch/users/arphillips/reports/bamqc/stats.meta.bamqc.2025-06-11.csv")
str(meta)
dim(meta)

meta$seqnames <- gsub(meta$bams, pattern = ".dedup.bam", replacement = "")

pca_df <- merge(x = pca, y = meta, by.x = "ind", by.y = "seqnames")
dim(pca_df)

# first convert to percentage variance explained
pve <- data.frame(PC = 1:20, pve = eigenval/sum(eigenval)*100)

# make plot
a <- ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity")
a + ylab("Percentage variance explained") + theme_light()

# calculate the cumulative sum of the percentage variance explained
cumsum(pve$pve)

# remake data.frame
# pca <- as.tibble(data.frame(pca, spp, loc, spp_loc))

# plot pca
pca_p <- ggplot(pca_df, aes(PC1, PC2, col = Latitude)) + 
  geom_point(size = 3) +
  theme_bw() + 
  xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + 
  ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)")) +
  scale_color_viridis()

ggsave(pca_p, filename = "/global/scratch/projects/fc_moilab/aphillips/aspen_snakemake/figures/pca.CA.plink.pdf", 
       width = 5, height = 4, unit = "in")

ggplot(pca_df, aes(PC1, PC2, col = Longitude)) + 
  geom_point(size = 3) +
  theme_bw() + 
  xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + 
  ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)")) +
  scale_color_viridis(option = "C")

ggplot(pca_df, aes(PC3, PC2, col = PC3)) + 
  geom_point(size = 3, alpha = 0.7) +
  theme_bw() + 
  ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)")) + 
  xlab(paste0("PC3 (", signif(pve$pve[3], 3), "%)")) +
  scale_color_viridis(option = "F")

ggsave(filename = "/global/scratch/projects/fc_moilab/aphillips/aspen_snakemake/figures/pca.CA.PC2_PC3.plink.pdf", 
       width = 5, height = 4, unit = "in")

# ggplot(pca_df, aes(y = Latitude, x = Longitude, color = PC3)) +
#   geom_point(alpha = 0.7) +
#   theme_bw() +
#   scale_color_viridis()

# ggplot(pca_df, aes(y = Latitude, x = Longitude, color = PC2)) +
#   geom_point(alpha = 0.7) +
#   theme_bw() +
#   scale_color_viridis(option = "C")
# 
# ggplot(pca_df, aes(y = Latitude, x = Longitude, color = PC3)) +
#   geom_point(alpha = 0.7) +
#   theme_bw() +
#   scale_color_viridis(option = "F")

## Map
library(raster)
library(sf)
library(ggplot2)
library(ggspatial)
library(USAboundaries)
library(purrr)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggpubr)

north_america <- ne_countries(continent = "North America", returnclass = "sf")
states <- us_states()

events_sf <- pca_df %>% 
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326) 
dim(events_sf)

pc1_p <- ggplot() + 
  geom_sf(data = north_america, size = 4, color = "black", fill = NA) +
  geom_sf(data = states, size = 4, fill = NA, color = "black") +
  geom_sf(data = events_sf, size = 3, aes(color = PC1), alpha = 0.5) +
  theme_bw() +
  theme(panel.grid = element_blank(), 
        axis.line = element_blank(),
  ) + 
  # ylim(c(10,70)) +
  # xlim(c(165, 57)) +
  ylim(c(30,50)) +
  xlim(c(130, 110)) +
  scale_color_viridis()

pc2_p <- ggplot() + 
  geom_sf(data = north_america, size = 4, color = "black", fill = NA) +
  geom_sf(data = events_sf, size = 3, aes(color = PC2), alpha = 0.5) +
  theme_bw() +
  theme(panel.grid = element_blank(), 
        y.axis.text = element_blank(),
        axis.line = element_blank(),
  ) + 
  # ylim(c(10,70)) +
  # xlim(c(165, 57)) +
  ylim(c(30,50)) +
  xlim(c(130, 110)) +
  scale_color_viridis(option = "C")

pc3_p <- ggplot() + 
  geom_sf(data = north_america, size = 4, color = "black", fill = NA) +
  geom_sf(data = events_sf, size = 3, aes(color = PC3), alpha = 0.5) +
  theme_bw() +
  theme(panel.grid = element_blank(), 
        # axis.text = element_blank(),
        axis.line = element_blank(),
        # axis.ticks = element_blank(),
        # legend.title = element_blank(),
        # legend.text = element_blank()
  ) + 
  # ylim(c(10,70)) +
  # xlim(c(165, 57)) +
  ylim(c(30,50)) +
  xlim(c(130, 110)) +
  scale_color_viridis(option = "F")

ggsave(pc3_p, filename = "/global/scratch/projects/fc_moilab/aphillips/aspen_snakemake/figures/pcamap.1206.pc3.plink.pdf", 
       width = 6, height = 5, unit = "in")

pcs_p <- ggarrange(pc1_p, pc2_p, ncol = 2, nrow = 1)

# ggsave(pcs_p, filename = "/global/scratch/projects/fc_moilab/aphillips/aspen_snakemake/figures/pcamap.1206.plink.pdf", 
#        width = 10, height = 5, unit = "in")
ggsave(pcs_p, filename = "/global/scratch/projects/fc_moilab/aphillips/aspen_snakemake/figures/pcamap.CA.plink.pdf",
       width = 10, height = 5, unit = "in")

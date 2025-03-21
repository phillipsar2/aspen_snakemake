### Title: Plot single read PCA analysis
### Author: Alyssa Phillips
### Date: 02/07/2025

library(stringr)
library(ggplot2)
library(dplyr)
library(forcats)
library(viridis)

# (1) Load data ----
C <- as.matrix(read.table("/global/scratch/users/arphillips/data/angsd/pca/singlepca.3dp70.chr1.covMat"))

# (2) Load metadata ----
bamlist <- read.csv("/global/scratch/projects/fc_moilab/aphillips/aspen_snakemake/data/bams/bamlist.02242025.txt", header = F, col.names = c("bam")) # bamlist ANGSD input

seq_names <- lapply(bamlist$bam, gsub, pattern = "dedup.bam", replacement="fastq.gz") %>%
  unlist()

meta <- read.csv("/global/scratch/projects/fc_moilab/aphillips/aspen_snakemake/metadata/megametadata.2025-02-24.csv")
str(meta)
dim(meta)

meta <- arrange(meta, bams)

tail(seq_names)
tail(meta$bams)

# (4)  Compute eigenvalues ----
# Compute eigenvalues and corresponding eigenvectors of S
# these are the principal components of S
e <- eigen(C)

# Print proportion of total var explained by the components
pc <- (e$values / sum(e$values)) *100
pc[1:10]

# (5) Plot ----
df <- as.data.frame(cbind(e$vectors[,1:10], meta))
names(df) <- c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", colnames(meta))
dim(df)

# df <- filter(df, PC1 < 0.25, PC2 < 0.25, PC2 > -0.025)

pdf(paste0("/global/scratch/users/arphillips/data/angsd/pca/PCA_plot.", Sys.Date(), ".", chr ,".pdf"),
    width = 7, height = 5)

df %>%
  ggplot(aes(x = PC1, y = PC2)) + 
  geom_point(aes(color =  Sample.type),
             size = 2) +
  xlab(paste0("PC1 (", round(pc[1],1), "%)")) +
  ylab(paste0("PC2 (", round(pc[2], 1), "%)")) +
  theme_bw(base_size = 12) +
  theme(
    # legend.position="none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())  +
  ggtitle("JGI Sample.Plate.Tube.Name")
# guides(color = 'none') 

df %>%
  ggplot(aes(x = PC1, y = PC2)) + 
  geom_point(aes(color =  Investigators),
             size = 2) +
  xlab(paste0("PC1 (", round(pc[1],1), "%)")) +
  ylab(paste0("PC2 (", round(pc[2], 1), "%)")) +
  theme_bw(base_size = 12) +
  theme(
    # legend.position="none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    text = element_text(size = 6))  +
  ggtitle("Investigator")
# guides(color = 'none') 

df %>%
  ggplot(aes(x = PC1, y = PC2)) + 
  geom_point(aes(color =  Cytotype),
             size = 2) +
  xlab(paste0("PC1 (", round(pc[1],1), "%)")) +
  ylab(paste0("PC2 (", round(pc[2], 1), "%)")) +
  theme_bw(base_size = 12) +
  theme(
    # legend.position="none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())  +
  ggtitle("Cytotype")
# guides(color = 'none') 

df %>%
  ggplot(aes(x = PC1, y = PC2)) + 
  geom_point(aes(color =  Longitude),
             size = 2) +
  xlab(paste0("PC1 (", round(pc[1],1), "%)")) +
  ylab(paste0("PC2 (", round(pc[2], 1), "%)")) +
  theme_bw(base_size = 12) +
  theme(
    # legend.position="none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())  +
  ggtitle("Longitude") +
  scale_color_viridis_c(option = "magma", direction = 1)
# guides(color = 'none') 

df %>%
  ggplot(aes(x = PC1, y = PC2)) + 
  geom_point(aes(color =  Latitude),
             size = 2, alpha = 0.7) +
  xlab(paste0("PC1 (", round(pc[1],1), "%)")) +
  ylab(paste0("PC2 (", round(pc[2], 1), "%)")) +
  theme_bw(base_size = 12) +
  theme(
    # legend.position="none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())  +
  ggtitle("Latitude") +
  scale_color_viridis()


df %>%
  ggplot(aes(x = PC1, y = PC2)) + 
  geom_point(aes(color =  depth),
             size = 2) +
  xlab(paste0("PC1 (", round(pc[1],1), "%)")) +
  ylab(paste0("PC2 (", round(pc[2], 1), "%)")) +
  theme_bw(base_size = 12) +
  theme(
    # legend.position="none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())  +
  scale_color_viridis_c(option = "magma", direction = -1) +
  ggtitle("Depth")

dev.off()


## Plot on map
df %>%
  ggplot(aes(y = Latitude, x = Longitude)) + 
  geom_point(aes(color = PC2),
             size = 2, alpha = 1) +
  ylab("Latitude") +
  xlab("Longitude") +
  theme_bw(base_size = 12) +
  scale_color_viridis_c(option = "magma", direction = -1) +
  theme(
    # legend.position="none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()) 
# guides(color = guide_legend(title = "PC1"))
# shape = guide_legend(title="Ploidy"))


## Plot map ----

library(ggmap)
# library(readxl)
library("rnaturalearth")
library("rnaturalearthdata")


## Get map & Specify bounding box for samples
na_bbox <- make_bbox(lat = Latitude, lon = Longitude, data = df)
# na_big <- get_map(location = na_bbox, maptype = "terrain", source = "google")
na <- ne_countries(scale = "medium", 
                   returnclass = "sf", 
                   continent = "north america")

# ggmap(na_big)+ 
#   geom_point(data = done_seqs, 
#              mapping = aes(y = Latitude, x = Longitude))

ggplot(na) +
  geom_sf() +
  coord_sf(xlim = c(na_bbox[1], na_bbox[3]), 
           ylim = c(na_bbox[2], na_bbox[4]), expand = FALSE) + 
  geom_point(data = df,
             mapping = aes(y = Latitude, 
                           x = Longitude,
                           color = PC1), size = 2) +
  scale_color_viridis_c(option = "magma", direction = -1) +
  theme_bw()

ggsave(paste0("/global/scratch/users/arphillips/data/angsd/pca/PCA_map.PC1.", Sys.Date(), ".", chr ,".pdf"),
       width = 7, height = 5, unit = "in") 

ggplot(na) +
  geom_sf() +
  coord_sf(xlim = c(na_bbox[1], na_bbox[3]), 
           ylim = c(na_bbox[2], na_bbox[4]), expand = FALSE) + 
  geom_point(data = df,
             mapping = aes(y = Latitude, 
                           x = Longitude,
                           color = PC2)) +
  scale_color_viridis_c(option = "magma", direction = -1) +
  theme_bw()

ggsave(paste0("/global/scratch/users/arphillips/data/angsd/pca/PCA_map.PC2.", Sys.Date(), ".", chr ,".pdf"),
       width = 7, height = 5, unit = "in") 



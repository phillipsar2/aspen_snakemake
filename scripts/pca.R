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

# (2) Load in metadata ----
meta <- read.csv("/global/scratch/projects/fc_moilab/aphillips/aspen_snakemake/metadata/megametadata.2025-02-24.csv")
str(meta)
dim(meta)

cov <- read.table("/global/scratch/users/arphillips/reports/bamqc/stats.bamqc.txt", header = F) %>%
  select(V1, V7)
names(cov) <- c("bams", "cov")
# cov$cov <- gsub(cov$cov, pattern = "X", replacement = "") %>% as.numeric
cov$bams <- str_split(cov$bams, pattern = "/", simplify=T)[,9] 
head(cov)

meta_cov <- merge(meta, cov, by = "bams", all = T )
head(meta_cov)
dim(meta_cov)

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
## The naive LD estimates (using just sample moments) with the posterior mean
# mom_ld <- mldest(geno = as.matrix(genos)[1:1000,], K = 2, type = "comp")
# plot(mom_ld)
# mom_ld$chri <- str_split(mom_ld$snpi, pattern = "_", simplify = T)[,1]
# mom_ld$chrj <- str_split(mom_ld$snpj, pattern = "_", simplify = T)[,1]
# mom_ld$posi <- str_split(mom_ld$snpi, pattern = "_", simplify = T)[,2]
# mom_ld$posj <- str_split(mom_ld$snpj, pattern = "_", simplify = T)[,2]
# 
# mom_ld %>%
#   dplyr::filter(chri == "Chr01") %>% 
#   dplyr::filter(chrj == "Chr01") %>%
#   ggplot(aes(x = posi, y = posj, color = r2)) +
#   geom_tile() +
#   scale_color_viridis()
  

# dist <- abs(as.numeric(mom_ld$posi[mom_ld$r2 > 0.9]) - as.numeric(mom_ld$posj[mom_ld$r2 > 0.9]))/1000
# hist(dist, breaks = 50, xlab = "Kbp")
# 50*1000

## Alternative option is to thin by position
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
## not right for mixed-ploidy pops
hist(as.matrix(genos_thin)[,2])
q2 <- rowSums(genos_thin == 0) / dim(genos_thin)[2] # q^2, minor allele
af <- sqrt(q2)
hist(as.matrix(af))

genos_maf <- genos_thin[af > 0.05 & af < 1,]
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

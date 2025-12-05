### Title: Aspen RONA
### Author: Alyssa Phillips
### Date: 11/17/2025

library(dplyr)
library(vcfR)

# Load data ----
bioclim = "bio1"
k = 4

##LFMM results ----
indir <- '/global/scratch/users/arphillips/data/lfmm/lfmm_results_all_k'
pvalues <- read.csv(paste0(indir, '/w_calibration_pvalue_full_genome_', bioclim, '_k', k, '.csv'))
betas <- read.csv(paste0(indir, '/effect_sizes_simple_full_genome_', bioclim, '_k', k, '.csv'))



# Drop p-values below threshold
p_threshold = c(0.1,0.01, 0.001)

# is.na(pvalues$bio1) %>% sum()
q_threshold <- quantile(pvalues$bio1, 0.001, na.rm = T)
keep_pvalues <- pvalues[pvalues$bio1 < q_threshold, ]
dim(keep_pvalues)

plot(-log10(pvalues$bio1), 
     pch = 19, 
     cex = .2, 
     xlab = "SNP", ylab = "-Log P",
     col = "grey")
abline(h = -log10(q_threshold), col = "blue")


## VCF ----
vcf <- read.vcfR("/global/scratch/users/arphillips/data/processed/filtered_snps/wgs_aspen.all.nocall.10dp90.diploids.pruned.vcf.gz")
gt <- vcfR::extract.gt(vcf, element = "GT", as.numeric = F)
gt[1:5,1:5]
gt <- t(gt) # ind by SNPs
dim(gt)

gt[is.na(gt)] <- 0 # replacing all ./. with 0 is wrong but all I've got right now
sum(is.na(gt))

gt[gt == "1/1"] <-2
gt[gt == "0/1"] <- 1

# gt.imp <- apply(gt, 2, function(x) replace(x, 
#                                            is.na(x), 
#                                            as.numeric( names( which.max(table(x))))))

gt.imp <- apply(gt, 2, as.numeric)
dim(gt.imp)
gt.imp[1:5,1:5]

sum(colnames(gt.imp) == betas$X)

pc <- prcomp(gt.imp)
plot(pc$sdev[1:20]^2, xlab = 'PC', ylab = "Variance explained")

## Past climate data ----
clm <- read.csv("/global/scratch/users/arphillips/data/climate/chelsa/chelsav2/GLOBAL/climatologies/1981-2010/bio/CHELSA_climate_data.1197.csv")
genos <- colnames(vcf@gt)[-1]
length(genos)

metadata <- read.csv("/global/scratch/projects/fc_moilab/aphillips/aspen_snakemake/metadata/megametadata.2025-08-26.csv")
metadata$seqname <- gsub(x = metadata$bamlist,pattern = ".dedup.bam", replacement = "") 
metasub <- metadata[metadata$seqname %in% genos,]
dim(metasub)

coor <- metasub %>% # extract and fix coor
  dplyr::select(seqname, Latitude, Longitude) %>%
  unique()
dim(coor)

coor[duplicated(coor$seqname),]
cr_coor <- c(23.524897, 23.524897, -104.684526,  -104.684526)
coor[coor$seqname == "53114.3.593590.AATAGAGATA-ACGCGGCCCT",][2:3] <- cr_coor
coor <- unique(coor)
dim(coor)

clmsub <- clm %>% # subset clm
  filter(lat %in% coor$Latitude, lon %in% coor$Longitude)
dim(clmsub)

clmsub_ordered <- clmsub[match(coor$Latitude, clmsub$lat), ] # order to match genotypes
clmsub_ordered2 <- clmsub[match(coor$Longitude, clmsub$lon), ]

gt.imp.sub <- gt.imp[-which(is.na(clmsub_ordered2$lat)),] # subset genotype data to match
dim(gt.imp.sub)

clmsub_ordered2.sub <- clmsub_ordered2[complete.cases(clmsub_ordered2),]
dim(clmsub_ordered2.sub)

clmsub_ordered2.sub[,4:22] <- apply(clmsub_ordered2.sub[,4:22], 2, scale)

# library(lfmm)
# mod.lfmm <- lfmm_ridge(Y = gt.imp.sub, X = clmsub_ordered2.sub$bio1, K = k)
# summary(mod.lfmm)
# mod.lfmm$U

## Future climate data ----
future <- read.csv("/global/scratch/users/arphillips/data/climate/chelsa/chelsav2/GLOBAL/climatologies/2071-2100/IPSL-CM6A-LR/ssp370/CHELSA_ipsl-cm6a-lr_ssp370_climate_data.1192.2025-11-19.csv")
head(future)

### Subset future climate data and sort to match the genotypes
futsub <- future %>% # subset clm
  filter(lat %in% clmsub_ordered2.sub$lat, lon %in% clmsub_ordered2.sub$lon)
head(futsub)

futsub_ordered <- futsub[match(clmsub_ordered2.sub$lat, futsub$lat), ] # order to match genotypes
futsub_ordered2 <- futsub[match(clmsub_ordered2.sub$lon, futsub$lon), ]
head(futsub_ordered2)
dim(clmsub_ordered2.sub)

futsub_ordered2[,5:23] <- apply(futsub_ordered2[,5:23], 2, scale)

# Calculate RONA ----
## idk what the intercept is

## one snp
dim(betas)
dim(future)
hist(gt.imp.sub[,2])

plot(x = gt.imp.sub[,260], y = clmsub_ordered2.sub$bio1)
abline(a = 0, b = betas$bio1[260] )

plot(x = gt.imp.sub[,119], y = clmsub_ordered2.sub$bio1)
abline(a = 0, b = betas$bio1[119] )


# rona <- (betas$bio1[1] * futsub_ordered2$bio01) - gt.imp.sub[,1]

rona_mat <- matrix(nrow = length(betas$bio1), ncol = 813 )
for (b in 1:length(betas$bio1)){
  rona_mat[b,] <- (betas$bio1[b] * futsub_ordered2$bio01) - gt.imp.sub[,b]
}
rona_mat[1:15,1:15]

rona_avg <- colMeans(abs(rona_mat), na.rm = T)
hist(rona_avg)
length(rona_avg)

library(ggmap)
# library(readxl)
library("rnaturalearth")
library("rnaturalearthdata")

df <- cbind(rona_avg, futsub_ordered2)
dim(df)

## Get map & Specify bounding box for samples
na_bbox <- make_bbox(lat = lat, lon = lon, data = df)
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
             mapping = aes(y = lat, 
                           x = lon,
                           color = rona), size = 2) +
  scale_color_viridis_c(option = "magma", direction = -1) +
  theme_bw()

# ggsave(paste0("/global/scratch/users/arphillips/data/rona_map.bio1.", Sys.Date(), ".", chr ,".pdf"),
       # width = 7, height = 5, unit = "in") 

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


# Make a map ----
#http://membres-timc.imag.fr/Olivier.Francois/LEA/files/Spatial_Prediction_Genomic_Offset.html
library(terra)
library(geodata)
library(fields)
library(maps)

coordinates <- clmsub_ordered2.sub[,c(3,2)]
plot(coordinates, cex = .4, col = "darkblue", ## verify points
     xlab = "Longitude", ylab = "Latitude",
     main = "Sample coordinates", las = 1)
maps::map(add = TRUE, interior = FALSE, col = "grey40")

s <- rast("/global/scratch/users/arphillips/data/climate/chelsa/chelsav2/GLOBAL/climatologies/1981-2010/bio/CHELSA_bio1_1981-2010_V.2.1.tif")
f <- rast("/global/scratch/users/arphillips/data/climate/chelsa/chelsav2/GLOBAL/climatologies/2071-2100/IPSL-CM6A-LR/ssp370/CHELSA_ipsl-cm6a-lr_ssp370_bio01_2071-2100_V.2.1.tif")
X.env = terra::extract(x = s, ## load historical climate data
                       y = data.frame(coordinates), 
                       cells = FALSE)
X.env = X.env[,-1] # remove IDs

my.colors = colorRampPalette(c("lightblue3", "orange2", "red3"))(100)

# future data
nc = 200 ## nc = resolution, higher is better but slower 

# range of longitude for NA (deg E)
long.mat <- seq(-160, -60, length = nc)

# range of latitude for NA (deg N)
lat.mat <- seq(15, 70, length = nc)

# matrix of cells for NA (nc times nc)
coord.mat <- NULL
for (x in long.mat){
  for (y in lat.mat) {coord.mat <- rbind(coord.mat, c(x,y))}
} 
  

## bins extreme values above .1 - see histogram
# go2  = go
# go2[go2 > .1] = .1

# Extract historical climate
env.new = terra::extract(x = s, 
                         y = data.frame(coord.mat), 
                         cells = FALSE)
env.new = env.new[,-1] 

# Extract future climate
env.pred = terra::extract(x = f, 
                          y = data.frame(coord.mat), 
                          cells=FALSE)
env.pred = env.pred[,-1]

# scale
m.x <- mean(X.env)
sd.x <- sd(X.env)
# m.x <- apply(X.env, 2, FUN = function(x) mean(x, na.rm = TRUE))
# sd.x <- apply(X.env, 2, function(x) sd(x, na.rm = TRUE))

env.new <- t(t(env.new) - m.x) %*% 1/sd.x
env.pred <- t(t(env.pred) - m.x) %*% 1/sd.x

dim(env.new)
mean(((env.new[1,] - env.pred[1,])  %*% betas$bio1[betas$X %in% keep_pvalues$X])^2, na.rm = TRUE)


# (betas$bio1[1] * futsub_ordered2$bio01) - gt.imp.sub[,1]
gg = NULL
for (i in 1:nrow(env.new)){
  gg[i] = mean(((env.new[i,] - env.pred[i,])  %*% betas$bio1[betas$X %in% keep_pvalues$X])^2, na.rm = TRUE)
} # could also run with just keep_pvalues
go = t(matrix(gg, byrow =  FALSE, ncol = nc))
hist(as.numeric(go), 
     main = "Histogram of GO values",
     xlab = "Geometric GO")

my.colors = colorRampPalette(c("lightblue3", "orange2", "red3"))(100)
fields::image.plot(long.mat, lat.mat, go, 
                   col = my.colors,
                   las = 1,
                   xlab = "Longitude",
                   ylab = "Latitude")

maps::map(add = TRUE, interior = FALSE, col = "grey40")
points(coordinates, col = "grey40", cex = .3)

zoom_xlim <- c(-125, -103)
zoom_ylim <- c(30, 45)

image.plot(long.mat, lat.mat, go, 
           col = my.colors,
           las = 1,
           xlim = zoom_xlim, 
           ylim = zoom_ylim)
maps::map(add = TRUE, interior = FALSE, col = "grey40")
points(coordinates, col = "grey40", cex = .3)


# v <- terra::vect(coordinates, crs="+proj=longlat")
# p <- terra::project(v, crs(r),crs="" )
library("rnaturalearth")
library("rnaturalearthdata")
library(ggmap)


## Get map & Specify bounding box for samples
na_bbox <- make_bbox(lat = lat, lon = lon, data = clmsub_ordered2.sub)
# na_big <- get_map(location = na_bbox, maptype = "terrain", source = "google")
na <- ne_countries(scale = "medium", 
                   returnclass = "sf", 
                   continent = "north america")

library(raster)
r <- rast(t(go)[nc:1,])
ext(r) <- c(-160, -60, 15, 70 )

terra::plot(r, 
             col = my.colors, 
             las = 1,
             xlab = "xcells",
             ylab = "ycells")
maps::map(add = TRUE, interior = FALSE, col = "grey40")
points(coordinates, col = "grey40", cex = .3)

new_extent <- maps::map(add = TRUE, interior = FALSE, col = "grey40")

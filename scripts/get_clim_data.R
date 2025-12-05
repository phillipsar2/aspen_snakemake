### Title: Download climate data
### Author: Alyssa Phillips
### Date: 08/14/2025

library(dplyr)
library("climenv")
library("fs")
library("raster")
library("dplyr")
library("sf")
library("terra")
library(sp)

# Make spatial points object
metadata <- read.csv("/global/scratch/projects/fc_moilab/aphillips/aspen_snakemake/metadata/megametadata.2025-08-26.csv")
head(metadata)
dim(metadata)

unique(metadata$Sample.Name) %>% length()
sum(is.na(metadata$Sample.Name)) # should be zero

n = unique(metadata$Sample.Name) %>% length() # number of samples-ish
n

## Subset to complete dataset
sub_meta <- metadata[!is.na(metadata[,c('Latitude')]),]
dim(sub_meta)

## Spatial points object
elev_df<- read.csv("/global/scratch/users/arphillips/data/climate/gee/gee_elevation_data.1192.fixed.09122025.csv")[,-1]
# pts <- st_as_sf(sub_meta, coords = c( 'Longitude', 'Latitude') ) # lat long coord
pts <- st_as_sf(elev_df, coords = c( 'Longitude', 'Latitude') ) # lat long coord
pts <- st_set_crs(pts, 4326) 

#################################
# present - Extract CHELSA bioclim data ----
#################################
# chelsa v2 climatology 1981-2010 bioclim variables

## Set up pts
chelsa_path <- "/global/scratch/users/arphillips/data/climate/chelsa/chelsav2/GLOBAL/climatologies/1981-2010/bio/"

# r <- rast(paste0(chelsa_path, "CHELSA_bio1_1981-2010_V.2.1.tif")) # WGS84

pts_transform <- st_transform(pts, "WGS84") # transform points to projection of raster

# plot(r)
# plot(pts_transform, add = T)

## Sample point data from rasters
bioclim <- paste0("bio", seq(1, 19, 1))

clim_mat <- matrix(NA, nrow = dim(pts_transform)[1], ncol = 22)
colnames(clim_mat) <- c("Sample.Name", "lat", "lon", bioclim)

# clim_mat[,1] <- sub_meta$Latitude
# clim_mat[,2] <- sub_meta$Longitude
clim_mat[,1] <- elev_df$Sample.Name
clim_mat[,2] <- elev_df$Latitude
clim_mat[,3] <- elev_df$Longitude
clim_df <- as.data.frame(clim_mat)
# head(clim_df)

for (i in bioclim[1:19]){
  r <- rast(paste0(chelsa_path, "CHELSA_", i, "_1981-2010_V.2.1.tif"))
  clim_df[,i] <- extract(r, pts_transform)[,2]
  print(paste0(i, " DONE"))
}

# lapply(clim_df[,3:19], function(x){hist(x)})

head(clim_df)
tail(clim_df)

write.csv(clim_df, paste0(chelsa_path, "CHELSA_climate_data.",n,".", Sys.Date(),".csv"))

# offsets <- c(-273.15, 0,0,0,-273.15,-273.15,0,-273.15,-273.15,-273.15,-273.15,0,0,0,0,0,0,0,0)
## offsets from https://chelsa-climate.org/wp-admin/download-page/CHELSA_tech_specification_V2.pdf
# scale = 0.1

#################################
# future -  Extract CHELSA bioclim data ----
#################################

## Set up pts
chelsa_path <- "/global/scratch/users/arphillips/data/climate/chelsa/chelsav2/GLOBAL/climatologies/2071-2100/IPSL-CM6A-LR/ssp370/"

pts_transform <- st_transform(pts, "WGS84") # transform points to projection of raster

## Sample point data from rasters
bioclim <- paste0("bio", c(paste0("0",seq(1, 9, 1)), seq(10, 19, 1)))

clim_mat <- matrix(NA, nrow = dim(pts_transform)[1], ncol = 22)
colnames(clim_mat) <- c("Sample.Name", "lat", "lon", bioclim)

# clim_mat[,1] <- sub_meta$Latitude
# clim_mat[,2] <- sub_meta$Longitude
clim_mat[,1] <- elev_df$Sample.Name
clim_mat[,2] <- elev_df$Latitude
clim_mat[,3] <- elev_df$Longitude
clim_df <- as.data.frame(clim_mat)
# head(clim_df)

for (i in bioclim[1:19]){
  r <- rast(paste0(chelsa_path, "CHELSA_ipsl-cm6a-lr_ssp370_", i, "_2071-2100_V.2.1.tif"))
  clim_df[,i] <- extract(r, pts_transform)[,2]
  print(paste0(i, " DONE"))
}

head(clim_df)
tail(clim_df)

write.csv(clim_df, paste0(chelsa_path, "CHELSA_ipsl-cm6a-lr_ssp370_climate_data.",n,".", Sys.Date(),".csv"))

#########################################
# Get elevation, aspect, and slope ----
#########################################
## Export pts as shapefile
st_write(pts[,c('X', 'Y')], # just the coordinate columns
         paste0("/global/scratch/users/arphillips/data/climate/coordinates.", n ,".shp"),
         append=FALSE) # overwrite

## Run GEE script `get_slope_aspect_elev`

## Import elevation calls
gee <- read.csv("/global/scratch/users/arphillips/data/climate/gee/elv_to_pts.csv")
gee_sub <- gee %>% dplyr::select('X', 'Y', elevation, Smpl_Nm)
dim(gee_sub)
head(gee_sub)

elev_df <- merge(metadata, gee, by.x = "Sample.Name", by.y = "Smpl_Nm", all = T, no.dups = T) %>%
  dplyr::select(Sample.Name, Latitude, Longitude, elevation) %>%
  distinct(.keep_all = T)
dim(elev_df)
head(elev_df)

## Extract samples that failed getting their data
# fail <- pts[!pts$Sample.Name %in% gee$Smpl_Nm,] %>%
#   unique()
# 
# elev_df[is.na(elev_df$X.x),]
# elev_df[is.na(elev_df$elevation),]
# 
# fail_df <- metadata[metadata$Sample.Name %in% fail$Sample.Name,] 
# fail_min <- dplyr::select(fail_df, Sample.Name, Latitude, Longitude, Investigators ) %>%
#   unique()
# dim(fail_min)

# write.table(fail_min, "/global/scratch/users/arphillips/data/climate/gee/coordinate_issues.08292025.txt", sep = "\t", row.names = F)

## Export the datapoints and elev
write.csv(elev_df, paste0("/global/scratch/users/arphillips/data/climate/gee/gee_elevation_data.",n,".", Sys.Date(),".csv"))

## Manually fix some datapoints with resolved_coordinate_issues sheet
elev_df <- read.csv("/global/scratch/users/arphillips/data/climate/gee/gee_elevation_data.1192.fixed.09122025.csv")[,-1]
head(elev_df)

#################################
# Extract ClimateNA data ----
#################################
library(data.table) # V 1.14.6
# install.packages('/global/scratch/users/arphillips/toolz/', repos=NULL, type='source')
library(ClimateNAr) # V 3.1.0
library(dplyr)

clm_dir = "/global/scratch/users/arphillips/data/climate/climatena/"

## Generate input file
# format is coordinates and elevations of sites
input_file <- '/global/scratch/users/arphillips/data/climate/gee/gee_elevation_data.1192.fixed.09122025.csv'

# Get data from ClimateNA 'NA'
period <- c(
  # "Normal_1961_1990.nrm", 
            "Normal_1971_2000.nrm", 
            "Normal_1981_2010.nrm"
            )
for (p in period){
  clm <- ClimateNAr::ClimateNA_API2(ClimateBC_NA='NA', 
                                    inputFile = input_file, 
                                    period = p, 
                                    MSY= 'Y')
  write.csv(clm, 
            file = paste0(clm_dir, p, ".Y.",Sys.Date(),".n",n,".csv"),
            row.names = F, col.names = T)
}

head(clm)

#################################
# Compare climate data ----
#################################
clm_dir = "/global/scratch/users/arphillips/data/climate/climatena/"
clm_files <- list.files(clm_dir)

n1961 <- read.csv(paste0(clm_dir, clm_files[1]))
n1971 <- read.csv(paste0(clm_dir, clm_files[2]))
n1981 <- read.csv(paste0(clm_dir, clm_files[3]))

chelsa_path <- "/global/scratch/users/arphillips/data/climate/chelsa/chelsav2/GLOBAL/climatologies/1981-2010/bio/"
list.files(chelsa_path, pattern = "CHELSA_climate_data.*")

chelsa <- read.csv(paste0(chelsa_path,"CHELSA_climate_data.1192.2025-09-15.csv"))[,-1]
head(chelsa)

# bind the files together

# plot the data comparisions
## bio12 to MAP
## bio1 to MAT
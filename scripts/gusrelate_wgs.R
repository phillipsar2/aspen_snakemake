### Title: GUSrelate - WGS
### Author: Alyssa Phillips
### Date: 1/27/2025

# library("devtools")
# devtools::install_github("tpbilton/GUSbase")
# devtools::install_github("tpbilton/GUSrelate")
library("GUSrelate")
library(vcfR)
library(stringr)
# install.packages("pheatmap")
library(pheatmap)
# install.packages("egg")
library(egg)
library(dplyr)
library(pheatmap)

# (1) Load data ----
dataset <- "CA"
file <- as.character("/global/scratch/users/arphillips/data/processed/filtered_snps/wgs_aspen.CA.nocall.10dp90.vcf.gz", verbose = FALSE )

vcf <- read.vcfR(file, verbose = FALSE )
vcf

# (2) Filter VCF ----
## LD Thinning
distance_thin <- function(vcfR,
                          min.distance=NULL){
  
  #if vcfR is not class vcfR, fail gracefully
  if (!inherits(vcfR, what = "vcfR")){
    stop("specified vcfR object must be of class 'vcfR'")
  }
  
  #function is useless if no distance is provided
  if (is.null(min.distance)){
    stop("filtering distance must be provided")
  }
  
  #logical test specifying that the minimum distance between SNPs for filtering must be at least 1 base pair or the logic doesn't work
  if (is.numeric(min.distance) != TRUE){
    stop("specified filtering distance must be numeric")
  }
  
  #logical test specifying that the minimum distance between SNPs for filtering must be at least 1 base pair or the logic doesn't work
  if (min.distance < 1){
    stop("filtering distance must be >= 1 base pair")
  }
  
  #logical test to ensure formatting of input vcf will allow accurate analysis of position in genome
  if (colnames(vcfR@fix)[1] != "CHROM"){
    stop("vcfR incorrectly formatted. vcfR@fix column 1 must be 'CHROM'")
  }
  
  #logical test to ensure formatting of input vcf will allow accurate analysis of position in genome
  if (colnames(vcfR@fix)[2] != "POS"){
    stop("vcfR incorrectly formatted. vcfR@fix column 2 must be 'POS'")
  }
  
  #set min distance specified by user
  j=min.distance
  
  #generate dataframe containing information for chromosome and bp locality of each SNP
  df<-as.data.frame(vcfR@fix[,1:2])
  
  #write test to identify and remove duplicated SNPs in input vcf
  if (length(unique(paste(df$CHROM,df$POS))) < nrow(df)){
    #remove duplicated SNPs
    vcfR<-vcfR[!duplicated(paste(df$CHROM,df$POS)),]
    #report to user
    message(nrow(df) - length(unique(paste(df$CHROM,df$POS)))," duplicated SNPs removed from input vcf")
    #regenerate df without duplicate inputs
    df<-as.data.frame(vcfR@fix[,1:2])
  }
  
  #generate list of all unique chromosomes in alphabetical order
  chroms<-levels(as.factor(df$CHROM))
  
  #intialize empty df to hold filtering
  keep.df<-data.frame()
  
  #make progress bar
  pb <- utils::txtProgressBar(min = 0, max = length(chroms), style = 3)
  
  #begin tracker
  pbtrack<-1
  
  #loop over each chromosome
  #for t in vector containing the name of each chromosome
  for (t in chroms){
    
    #isolate the SNP positions on the given chromosome
    fix.sub<-as.numeric(df$POS[df$CHROM == t])
    
    #order the positions numerically
    fix.sub<-fix.sub[order(fix.sub)]
    
    #set the first position
    prev<-fix.sub[1]
    
    #initialize empty vector
    k<-c()
    
    #always keep first SNP on the chromosome
    k[1]<-TRUE
    
    #loop to decide whether to keep each following SNP
    if (length(fix.sub) < 2){
      #if chrom only has 1 SNP, do nothing
    } else{
      
      #else, use a for loop to determine which SNPs to keep that satisfy our distance criteria
      for (i in 2:length(fix.sub)){
        
        #store logical indicating whether this SNP is greater than j base pairs from the previous SNP
        k[i]<- fix.sub[i] > prev+j
        
        #if it is, then we keep this SNP, making it the new 'previous' for assessing the next point.
        #If we don't keep the SNP, we don't update the closest point
        if (fix.sub[i] > prev+j){
          prev<-fix.sub[i]
        }
        
        #close for loop
      }
      #close else statement
    }
    
    #make a dataframe with the precise info for each SNP for this chromosome
    chrom.df<-data.frame(CHROM=rep(t, times=length(fix.sub)), POS=fix.sub, keep=k)
    
    #now we rbind in the information for this chromosome to the overall df
    keep.df<-rbind(keep.df,chrom.df)
    
    #empty df for this chrom to prepare for the next one
    chrom.df<-NULL
    
    #update progress bar
    utils::setTxtProgressBar(pb, pbtrack)
    
    #update tracker
    pbtrack<-pbtrack+1
    
  } #close for loop, start over on next chromosome
  
  #close progress bar
  close(pb)
  
  #order the dataframe to match the order of the input vcf file
  #remove scientific notation
  keep.df$POS<-format(keep.df$POS,scientific = FALSE)
  df$POS<-format(df$POS,scientific = FALSE)
  #make sure class matches between the columns you're trying to merge
  keep.df$POS<-as.numeric(as.character(keep.df$POS))
  df$POS<-as.numeric(as.character(df$POS))
  #make sure class matches between the columns you're trying to merge
  keep.df$CHROM<-as.character(keep.df$CHROM)
  df$CHROM<-as.character(df$CHROM)
  #add tracking column
  df$id<-c(1:nrow(df))
  #merge
  order.df<-merge(keep.df,df)
  #order based on tracking column
  order.df<-order.df[order(order.df$id),]
  
  #order.df<-keep.df[match(paste(df$CHROM,format(df$POS,scientific = FALSE)), paste(keep.df$CHROM,format(keep.df$POS,scientific = FALSE))),]
  #note: the position vector must be stripped of scientific notation otherwise identical numbers will not be recognized as matches, giving NA values
  #note: this old approach using match() has been deprecated because there were too many formatting issues causing
  #match to not be able to correctly match the order between 'df' and 'keep.df'. Now we use merge() which seems to be more robust
  
  #write a test to catch if this internal dataset is not able to merge correctly
  if (sum(is.na(order.df)) > .5){
    stop("internal error with the merge function. Please email a copy of your input vcf to devonderaad@gmail.com for a bug fix")
  }
  
  #write a test to catch if this internal dataset is not able to merge correctly
  if (sum(order.df$id != c(1:nrow(df))) > .5){
    stop("internal error with the merge function. Please email a copy of your input vcf to devonderaad@gmail.com for a bug fix")
  }
  
  #subset vcfR locus info based on the logical column from our dataframe
  #vcfR@fix<-vcfR@fix[order.df$keep,]
  #subset genotypes based on logical
  #vcfR@gt<-vcfR@gt[order.df$keep,]
  
  #realized there is no need to do this subsetting separately
  vcfR<-vcfR[order.df$keep,]
  
  #calculate number of total SNPs input
  z<-nrow(keep.df)
  
  #calculate total SNPs retained
  p<-nrow(vcfR@fix)
  
  #print info to screen
  message(p," out of ",z," input SNPs were not located within ",j," base-pairs of another SNP and were retained despite filtering")
  
  #return vcfR
  return(vcfR)
}
  
# vcf.thin <- distance_thin(vcf, min.distance = 500)
# vcf.thin

keep_sites <- read.csv("/global/scratch/users/arphillips/data/pca/wgs_aspen.CA.10dp90.thin.prune.in", header = F)
# sites <- str_split(keep_sites[,1], pattern = ":", simplify = T)
# colnames(sites) <- c("CHROM", "POS")
head(sites)

vcf_sites <- as.data.frame(vcf@fix) %>%
  select(CHROM, POS)
vcf_sites$chrpos <- paste0(vcf_sites$CHROM, ":", vcf_sites$POS)

vcf_site_index <- which(vcf_sites$chrpos %in% keep_sites[,1])
length(vcf_site_index)

vcf@fix <- vcf@fix[vcf_site_index,]
vcf@gt <- vcf@gt[vcf_site_index,]
vcf
dim(vcf@gt)

## Save filtered vcf
filt_thin_file <- paste0("/global/scratch/users/arphillips/data/vcf/", dataset, ".thinned.vcf.gz")
write.vcf(vcf, file = filt_thin_file)

# filt_file <- as.character("/global/scratch/projects/fc_moilab/aphillips/obv_aspen/data/vcf/obv_aspen.filt.LD.MAC.vcf")
# write.vcf(vcf.maf, file = filt_file)

# (2) Load in ploidy information ----

meta <- read.csv("/global/scratch/users/arphillips/data/gbs2ploidy/flow_cyt_predictions.csv")
str(meta)
dim(meta)

# (3) Create required genotype metadata file ----
## Assign numeric ploidy ID
meta$ploidy <- ifelse(meta$ploidy_call == "diploid", yes = 2, no = 3)

## Remove samples without ploidy info
sum(colnames(vcf@gt) %in% meta$sample)
samples_keep <- colnames(vcf@gt) %in% meta$sample
samples_keep[1] <- TRUE # keep FORMAT column
vcf.sub <- vcf[,samples_keep]
vcf.sub

meta$sample[!meta$sample %in% colnames(vcf@gt)]

## Write filtered vcf to file for GUSrelate
sub_file <- paste0("/global/scratch/users/arphillips/data/vcf/",dataset,".genowithploidy.vcf.gz")
write.vcf(vcf.sub, file = sub_file)

## Make csv files for gusrelate
geno_info_sub <- dplyr::select(meta, sample, ploidy) %>%
  filter(sample %in% colnames(vcf.sub@gt))
colnames(geno_info_sub) <- c("ID", "Ploidy")
dim(geno_info_sub)

geno_sub_file <- as.character(paste0("/global/scratch/users/arphillips/data/gusrelate/",dataset,".subset_ploidy.meta.csv"))
write.csv(geno_info_sub, geno_sub_file, row.names = F, quote = F) # ID, Ploidy, Group


# (3) Create GRM object ----
ra <- VCFtoRA(sub_file) # convert VCF to RA format
radata <- readRA(ra) # process RA file to RA object
grm <- makeGRM(RAobj = radata, 
               samfile = geno_sub_file,
               filter=list(MAF=0.05)
               )

# (4) Compute genotype relatedness matrix (GRM) ----
# seq_error = 0.01
grm$computeGRM(name = "GRM_VR", method="VanRaden")

# (5) PCA of GRM ----
grm$PCA(name = "GRM_VR", colour="Ploidy", shape=NULL, )
ggsave(filename = paste0("/global/scratch/users/arphillips/data/gusrelate/",dataset,".grm.PCA.ploidy.",Sys.Date(),".jpeg"),
       width = 4, height = 5, unit = "in")

# (6) Plot GRM ----
grm_mat <- grm$extractGRM(name = "GRM_VR")

# grm_mat <- read.table("/global/scratch/projects/fc_moilab/aphillips/obv_aspen/data/gusrelate/grm.filt.LD.2025-01-27.csv", sep = ",", header =T, row.names = 1)

metadata <- read.csv("/global/scratch/projects/fc_moilab/aphillips/aspen_snakemake/metadata/megametadata.2025-06-09.csv")
metadata$seqname <- str_split(metadata$bams, ".dedup.bam", simplify = T)[,1]
metadata_sub <- metadata[metadata$seqname %in% geno_info_sub$ID,]
metadata_sub <- arrange(metadata_sub, seqname)
dim(metadata_sub)

sum( metadata_sub$seqname == colnames(grm_mat) )

colnames(grm_mat) <- metadata_sub$Sample.Name
rownames(grm_mat) <- metadata_sub$Sample.Name

dim(grm_mat)
heat <- pheatmap(grm_mat, fontsize = 4)
ggsave(heat, filename = paste0("/global/scratch/users/arphillips/data/gusrelate/",dataset,".grm.heatmap.",Sys.Date(),".jpeg"), 
       height = 10, width = 10, units = "in")

## Plot with no diagonal
diag(grm_mat) <- NA
heat <- pheatmap(grm_mat, fontsize = 4)
ggsave(heat, filename = paste0("/global/scratch/users/arphillips/data/gusrelate/",dataset,".grm.heatmap.nodiag.",Sys.Date(),".jpeg"), 
       height = 20, width = 20, units = "in")

## Plot with less than zero as zero
grm_mat[grm_mat < 0] <- 0
heat <- pheatmap(grm_mat, fontsize = 10)
ggsave(heat, filename = paste0("/global/scratch/users/arphillips/data/gusrelate/",dataset,".grm.heatmap.nonegs.withdiag.",Sys.Date(),".jpeg"), 
       height = 10, width = 10, units = "in")

# Look at inbreeding ----
grm_mat <- grm$extractGRM(name = "GRM_VR")
colnames(grm_mat) <- metadata_sub$Sample.Name
rownames(grm_mat) <- metadata_sub$Sample.Name

d <- diag(grm_mat)
hist(d,
     main = NULL, xlab = "Inbreeding coefficient",
     breaks = 15,
     xlim = c(0,1))

# (7) Pull out clones ----
kin <- which(grm_mat > 0.8, arr.ind = TRUE) %>% 
  as.data.frame() # rownames are colnames

kin$colnames <- colnames(grm_mat)[kin$col] 
kin$rownames <- rownames(grm_mat)[kin$row]

# Expected self-relatedness ~ 0.9-1.25
# Parent-offspring or full-sibiling or grandparent-offspring ~ 0.4-0.8
# Half siblings or half cousins pibling-half nibling ~ 0.2-0.35
# Unrelated = 0 - 0.2

# > Plot by geographic position ----
metadata <- read.csv("/global/scratch/projects/fc_moilab/aphillips/aspen_snakemake/metadata/megametadata.2025-06-09.csv")
metadata$seqname <- str_split(metadata$bams, ".dedup.bam", simplify = T)[,1]

meatadata_kin <- metadata[metadata$seqname %in% kin$colnames,]


kin_nodups$colSite <- str_split(kin_nodups$colnames, pattern = "_", simplify = T)[,1]
kin_nodups$rowSite <- str_split(kin_nodups$rownames, pattern = "_", simplify = T)[,1]

merge_col <- merge(x = kin_nodups, y = meta_RMBL, by.x = "colSite", by.y = "Site_Code") %>%
  select(colSite, Latitude, Longitude, relationship, values)
merge_row <- merge(x = kin_nodups, y = meta_RMBL, by.x = "rowSite", by.y = "Site_Code")  %>%
  select(rowSite, Latitude, Longitude)

segments_df <- cbind(merge_col, merge_row)
colnames(segments_df) <- c("colSite", "colLatitude", "colLongitude", "relationship", "values", "rowSite", "rowLatitude", "rowLongitude")
head(segments_df)


ggplot(meta_RMBL, aes(x = Longitude, y = Latitude)) + 
  geom_point(size = 2)  +
  # geom_point(data = start, aes(x=X, y=Y, color=as.character(group)), size=2) +
  geom_segment(data = segments_df[segments_df$relationship == "clone",], 
               aes(x = colLongitude, y = colLatitude, 
                   xend = rowLongitude, yend = rowLatitude,
                   color = relationship 
                   ) 
               )


## NEED to somehow put the start and stop coordinates for each segment into one df
# https://stackoverflow.com/questions/56862791/ggplot-line-segments-from-one-point-to-many-from-different-dataframes


plot(heat$tree_row) # cut heatmap into clusters by branhc
abline(h=2.5, col="red", lty=2, lwd=2)

heat_groups <- cutree(heat$tree_row, h=2.5) %>% as.data.frame() # designate groups
colnames(heat_groups) <- "group"
unique(heat_groups$group) %>% length()

heat_groups$Site_Code <- str_split(rownames(heat_groups), pattern = "_", simplify = T)[,1]
heat_meta <- merge(meta_RMBL, heat_groups, by = "Site_Code")
heat_meta$group <- as.factor(heat_meta$group)

ggplot(heat_meta, aes(x = Longitude, y = Latitude, shape = group, col = group, group = group)) +
  geom_point() +
  # geom_line() +
  theme_bw() +
  scale_shape_manual(values = 1:nlevels(heat_meta$group))
  # theme(legend.position="none")

group4 <- heat_meta$Site_Code[heat_meta$group == "4"]
group4_tmp <- kin_nodups[kin_nodups$colSite %in% group4,]
group4_df <- kin_nodups[kin_nodups$rowSite %in% group4,]
dim(group4_df)



## Write relationships to file
colnames(clones) <- NULL
rownames(clones) <- NULL

write.csv(clones, paste0("/global/scratch/projects/fc_moilab/aphillips/spectral_aspen/data/gusrelate/", dataset,".gusrelate.clonelist.",Sys.Date(),".csv"),
          row.names = F)  

# (8) Export GRM ----
grm_file <- paste0("/global/scratch/projects/fc_moilab/aphillips/spectral_aspen/data/gusrelate/", dataset,".grm.filt.",Sys.Date(),".csv")
grm$writeGRM(name = "GRM_VR", filename = grm_file, IDvar=NULL)

# grm_mat <- read.table("/global/scratch/projects/fc_moilab/aphillips/spectral_aspen/data/gusrelate/spectral_aspen.grm.filt.2025-04-21.csv",
# header = T, sep = ",", row.names = 1)
# grm_mat <- read.table("/global/scratch/projects/fc_moilab/aphillips/spectral_aspen/data/gusrelate/RMBLgrm.filt.2025-05-12.csv",
# header = T, sep = ",", row.names = 1)

# (9) Jaccard Similarity Index ----
# jaccard <- function(a, b) {
#   intersection = length(intersect(a, b))
#   union = length(a) + length(b) - intersection
#   return (intersection/union)
# }

## Extract genotypes
gt <- extract.gt(vcf.sub)
gt[1:5,1:5]

## Swap genotypes to 0,1,2
gt_cont <- NULL
gt_cont <- apply(gt, c(1,2), function(x) {
  if (is.na(x)){      # is.na needs to be the first test
    NA               # assign to homozygous ref and remove the first alleles
  } else if (x == "0/0"){
    "0" 
  } else if (x == "0/1"){
    "1"
  } else if (x == "1/1"){
    "2"}})

gt_num <- apply(gt_cont, 1, as.numeric)
gt_num[1:30,1:30]
dim(gt_num)

## Remove invariant sites
gt_var <- gt_num[,!colSums(gt_num, na.rm = T) == 0]
dim(gt_var)

## Estimate JSI
# library(devtools)
# install_github("ramhiser/clusteval")
library(clusteval)
# cluster_similarity(gt_var[,1], gt_var[,2], similarity="jaccard")
# cluster_similarity(gt_var[,100], gt_var[,2], similarity="jaccard")

n_geno <- dim(gt_var)[1]
j_mat <- matrix(nrow = n_geno, ncol = n_geno)

for (i in 1:n_geno){
  for (j in 1:n_geno){
    j_mat[i,j] <- cluster_similarity(gt_var[i,], gt_var[j,], similarity="jaccard")
  }
}

colnames(j_mat) <- colnames(gt_cont)
rownames(j_mat) <- colnames(gt_cont)

write.csv(j_mat, paste0("/global/scratch/projects/fc_moilab/aphillips/spectral_aspen/data/gusrelate/", dataset, ".jaccard_mat.", Sys.Date(), ".csv"))

# (10) Plot Jaccard Similarity Matirx ----
j_mat <- read.table(paste0("/global/scratch/projects/fc_moilab/aphillips/spectral_aspen/data/gusrelate/WUSG.jaccard_mat.05132025.csv"), 
                    sep = ",", header = T, row.names = 1)
j_mat[1:5,1:5]

jheat <- pheatmap(j_mat, fontsize = 4)
ggsave(jheat, filename = paste0("/global/scratch/projects/fc_moilab/aphillips/spectral_aspen/data/gusrelate/",dataset,".jsi.heatmap.",Sys.Date(),".jpeg"), 
       height = 20, width = 20, units = "in")

# Histograms of the two matrices
# pdf("/global/scratch/projects/fc_moilab/aphillips/spectral_aspen/data/gusrelate/WUSG.grm_and_jsi_hist.05132025.pdf")
par(mfrow=c(1,2))
hist(j_mat[lower.tri(j_mat)], xlab = "Jaccard Similarity Index", main = "WUSG")
hist(grm_mat[lower.tri(grm_mat)], xlab = "Relatedness", main = NULL)
dev.off()

# Compare values from two matrices
j_df <- t(combn(colnames(j_mat), 2))
jsi_df <- data.frame(j_df, jsi=j_mat[j_df])
jsi_df$combo <- paste0(jsi_df$X1, "-", jsi_df$X2)

g_df <- t(combn(colnames(grm_mat), 2))
grm_df <- data.frame(g_df, grm=grm_mat[g_df])
grm_df$combo <- paste0(grm_df$X1, "-", grm_df$X2)

grm_jsi_df <- merge(grm_df, jsi_df, by = "combo", all.x = T)

plot(grm_jsi_df$grm, grm_jsi_df$jsi, xlab = "Genetic relatedness", ylab = "Jaccard Similarity Index")
cor(grm_jsi_df$grm, grm_jsi_df$jsi)

# Title: SNP Quality Assessment - WGS
# Author: Alyssa Phillips
# Date: 1/21/24

library(dplyr)
library(ggplot2)

##############
### Raw snps ----
##############

# Load quality table for one chromosome
chr <- "Chr01"
# qual <- read.table(paste0("/global/scratch/users/arphillips/reports/filtering/wgs_aspen.", chr ,".table") , header = T)
# head(qual)
# dim(qual)

chr_files <- Sys.glob(paste0("/global/scratch/users/arphillips/reports/filtering/wgs_aspen.", chr ,"*.table"))
chr_list <- lapply(chr_files, function(x) read.table(x, header = T))
qual <- do.call(rbind, chr_list)

genos = "1207"

pdf(paste0("/global/scratch/users/arphillips/reports/filtering/alignment_quality.", Sys.Date(), ".", chr ,".pdf"))

qual %>%
  ggplot(aes(x=QUAL)) +
  geom_density(fill="#69b3a2", color="#e9ecef", alpha=0.9, adjust = 0.8) +
  xlim(c(0,150)) +
  theme_bw() +
  geom_vline(xintercept = 20, color = "black")+
  geom_vline(xintercept = 30, color = "black") +
  xlab("Base quality (QUAL)") +
  ggtitle(paste0(chr))

qual %>%
  ggplot(aes(x=MQ)) +
  geom_histogram(fill="#69b3a2", color="#e9ecef", alpha=0.9) +
  xlim(c(0,65)) + # bwa-mem doesn't go above MQ of 60
  theme_bw() +
  geom_vline(xintercept = 20, color = "black")+
  geom_vline(xintercept = 30, color = "black") +
  xlab("Mapping quality (MQ)") +
  ggtitle(paste0(chr, ", n = ", genos ))

qual %>%
  ggplot(aes(x=DP)) +
  geom_histogram(fill="#69b3a2", color="#e9ecef") +
  xlim(c(0,10000)) +
  theme_bw() +
  xlab("Total depth at each site (DP)") +
  ggtitle(paste0(chr, ", n = ", genos ))

dev.off()

# Playing around with filter thresholds
dim(qual)

dplyr::filter(qual, QUAL > 20, MQ > 20) %>%
  dim()

dplyr::filter(qual, QUAL > 30, MQ > 30) %>%
  dim()

##############
### Filtering for Depth ----
##############
# Load quality table
chr <- "Chr02"
# qual <- read.table(paste0("/global/scratch/users/arphillips/reports/filtering/depth/wgs_aspen.", chr ,".filtered.nocall.table") , header = T)

chr_files <- Sys.glob(paste0("/global/scratch/users/arphillips/reports/filtering/depth/wgs_aspen.", chr ,"*.table"))
chr_list <- lapply(chr_files[1], function(x) read.table(x, header = T))
qual <- do.call(rbind, chr_list)

# head(qual)
# str(qual)
dim(qual)

qual[1:10,1:10]

genotypes <- dim(qual)[2] - 2

# Prep table
# qual <- qual[qual$CHROM != "CHROM",] # shouldn't be necessary if files are properly joined

qual[,3:dim(qual)[2]] <- lapply(qual[,3:dim(qual)[2]], as.numeric)
str(qual)
dim(qual)

# Replace 0s with NA
qual[qual == 0] <- "NA"
dim(qual)

# Estimate mean genotype depth across sites
geno_dp <- colMeans(qual[3:dim(qual)[2]], na.rm = T)

## Plot mean genotype depth
p_gdp <- reshape2::melt(geno_dp) %>%
  ggplot(aes(x=value)) +
  geom_density(fill="#69b3a2", color="#e9ecef", alpha=0.9, adjust = 0.8) +
  xlim(c(0,max(geno_dp, na.rm = T))) +
  theme_bw() +
  geom_vline(xintercept = mean(geno_dp, na.rm = T), color = "black") +
  labs(x = "Average genotype depth", 
       title = paste0("Average genotype depth for ", genotypes, " genotypes"))

ggsave(plot = p_gdp, filename = paste0("/global/scratch/users/arphillips/reports/filtering/depth/genotype_depth.", Sys.Date(), ".", chr ,".jpeg"),
       width = 5, height = 4, units = "in")

# Plot genotype depth distributions
qlist <- matrix(nrow = genotypes, ncol = 3) # define number of samples (10 samples here)
qlist <- data.frame(qlist, row.names=colnames(qual)[-c(1:2)])

pdf(paste0("/global/scratch/users/arphillips/reports/filtering/depth/genotype_depth_distributions.",Sys.Date(),".pdf"))
par(mfrow=c(4,3)) # mfrow sets max number of plots (rows by columns), automatically makes new pages if PDF

for (i in 3:dim(qual)[2]) {
  qlist[i-2,] <- quantile(qual[,i], c(.05, .1, .99), na.rm=T)
  d <- density(qual[,i], from=0, to=100, bw=1, na.rm = T)
  plot(d, xlim = c(0,50), main=rownames(qlist)[i-2], col="blue", lwd=2, xlab = NULL, cex.main=0.5)
  abline(v=qlist[i-2,c(1,3)], col='red', lwd=3)
}

dev.off()

# # Test out some depth + missingness filters
max = 90
min = 10
miss = 0.1

qual[qual < min] <- "NA"
qual[qual > max] <- "NA"

# genos_with_data <- rowSums(is.na(qual[,3:dim(qual)[2]])) # number of genotypes to be removed at each site
hist(rowSums(is.na(qual[,3:dim(qual)[2]])))

# Depth per site
# site_dp <- rowSums(qual[3:dim(qual)[2]], na.rm = T)
# hist(site_dp, xlim = c(0, 40000))

# Genotypes with data (not zero)
genos_with_data <- rowSums(is.na(qual[,3:dim(qual)[2]]))
hist(genos_with_data, 
     xlab = "Number of genotypes sequenced", 
     main = "Number of genotypes sequenced per site", 
     ylab = "Number of sites")


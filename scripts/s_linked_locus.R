# Title: Locating the sex-linked locus
# Author: Alyssa Phillips
# Date: 1/28/25

library(data.table)
library(stringr)
library(dplyr)
library(ggplot2)

# Get file names
dir <- as.character("/global/scratch/users/arphillips/data/toz19")
files <- list.files(path = dir, pattern = ".cov.txt")[1:175]

# test <- read.table(paste0(dir, "/", files[1])) %>%
#   data.table()
# dim(test)
# 
# # Create intervals
# int_length <- 10000
# chr13_length <- 18426903
# interval_start <- seq(from = 1, to = chr13_length, by = int_length) # mark intervals
# pos_group <- rep(interval_start, 10000) %>% sort() # assign group to each position
# test$pos_group <- pos_group[1:chr13_length] # truncate last window
# 
# # Estimate mean and SD in these intervals
# chr13_mean = test[,
#                   mean( V1 ), # function
#                   by = pos_group # group
# ]
# 
# chr13_sd = test[,
#                   mean( V1 ), # function
#                   by = pos_group # group
# ]
# 
# dim(chr13_mean)
# head(chr13_mean)
#
## Plot
# plot(x = chr13_mean$pos_group, y = chr13_mean$V1, type ="l" )

###
### Write a for loop to do this for all the individuals
###

# Create intervals
int_length <- 10000
chr13_length <- 18426903
interval_start <- seq(from = 1, to = chr13_length, by = int_length) # mark intervals
pos_group <- rep(interval_start, each = 10000) %>% sort() # assign group to each position
pos_group_crop <- pos_group[1:chr13_length] # truncate last window

# Create empty matrix
chr13_means <- matrix(nrow = length(interval_start), ncol = length(files))
colnames(chr13_means) <- str_split(files, pattern = ".chr13.", simplify = T)[,1]

# Estimate means
for (f in 1:length(files)){
  cov <- read.table(paste0(dir, "/", files[f])) %>%
    data.table()
  
  cov$pos <- pos_group_crop 
  
  # Estimate mean and SD in these intervals
  tmp = cov[,
                    mean( V1 ), # function
                    by = pos # group
                    ]
  chr13_means[,f] <- tmp$V1
}

####
# Depth from vcf
####
chr <- "Chr13"
qual <- read.table(paste0("/global/scratch/users/arphillips/reports/filtering/depth/wgs_aspen.", chr ,".filtered.nocall.table") , header = T)
qual[1:5,1:5]
genos <- dim(qual)[2]-2

# Create intervals
int_length <- 10000
chr13_length <- 19195830 # from fai
interval_start <- seq(from = 1, to = chr13_length, by = int_length) # mark intervals

interval_stop <- interval_start - 1
interval_stop[1] <- chr13_length
interval_stop <- sort(interval_stop)

intervals <- cbind(interval_start, interval_stop)
intervals
head(intervals)

# pos_group <- rep(interval_start, each = 10000) %>% sort() # assign group to each position
# pos_group_crop <- pos_group[1:chr13_length] # truncate last window

# Estimate genotype mean depth in 10k windows
chr13_means <- matrix(nrow = length(interval_start), ncol = genos)
names <- str_split(colnames(qual), pattern = ".DP", simplify = T)[,1]
colnames(chr13_means) <- names[-c(1:2)]

for (i in 1: length(interval_start)){
  chr13_means[i,] <- qual %>%
    filter(POS >= interval_start[i], POS <= interval_stop[i]) %>%
    select(-POS, -CHROM) %>%
    colMeans()
}

# Estimate mean depth of each genotype across sites
geno_dp <- colMeans(qual[3:dim(qual)[2]], na.rm = T)
hist(geno_dp)

# Standardize the window coverage by the mean genome-wide coverage
chr13_stand <- mapply(`/`, data.frame(chr13_means), geno_dp)

int_df <- data.frame(intervals, chr13_stand)
str(int_df)

# Plot on individual
plot(x = int_df$interval_start, y=int_df$X52832.2.464561.GTGTCAAC.GTATCGTG, "l")
abline(a=1, b=0, lty = 2, col = "red")

# Plot all
## Melt wide to long
int_df_long <- melt(setDT(int_df), id.vars = c("interval_start","interval_stop"), variable.name = "genotype")
head(int_df_long)

## ggplot
ggplot(int_df_long, aes(x = interval_start, y = value, color = genotype)) +
  geom_line(alpha = 0.5) +
  theme_bw() +
  guides(color="none") +
  geom_hline(yintercept = 1, linetype="dashed") +
  xlab("Chr13 positon") +
  ylab("Coverage relative to the mean") +
  xlim(c(12000000,12030001)) +
  ylim(c(-1,2)) +
  geom_segment(aes(x = 12003994, y = 0, xend = 12013044, yend = 0), color = "purple") +
  geom_segment(aes(x = 12003994, y = 0, xend = 12009181, yend = 0), color = "blue")

# Plot ind with known sex
known <- read.csv("/global/scratch/projects/fc_moilab/aphillips/aspen_snakemake/sex_data/test_sex_data.csv")
head(known)
known$ID <- str_split(known$RQC.Seq.Unit.Name, ".fastq.gz", simplify = T)[,1]
known$IDX <- paste0("X", known$ID)

bamlist <- read.csv("")

sum(known$ID %in% colnames(int_df))

colnames(int_df)

# Title: Updog genotyping
# Author: Alyssa Phillips
# Date: 1/15/2025

#devtools::install_github("dcgerard/updog")
library("updog") # v2.1.5
#devtools::install_github(repo="knausb/vcfR")
library(vcfR)
library(dplyr)
library(stringr)
library("argparser")
library(ldsep) # v2.1.5
library(corrplot)

# Argument name
ap <- arg_parser("Updog genotyping")

# add mandatory positional arguments (filename)
ap <- add_argument(ap, "vcf", help = "Input VCF")

# add additional arguments
ap <- add_argument(ap, "--meta", help = "Metadata file")
ap <- add_argument(ap, "--ploidy", help = "choose 'diploid' or 'triploid'")
ap <- add_argument(ap, "--cores", help = "Number of available cores")
ap <- add_argument(ap, "--outdir", help = "Output directory")
ap <- add_argument(ap, "--chr", help = "Region")

# parse arguments
argv <- parse_args(ap)

ploidy_level <- as.character(argv$ploidy)
chr <- as.character(argv$chr)
cores <- as.numeric(argv$cores)
outdir <- as.character(argv$outdir)

## test arguments
# ploidy_level <- as.character("diploid")
# chr <- as.character("Chr02")
# cores <- as.numeric(4)
# outdir <- as.character("/global/scratch/users/arphillips/data/updog")
# vcf <- read.vcfR("/global/scratch/users/arphillips/data/processed/filtered_snps/wgs_aspen.all.nocall.10dp90.vcf.gz", verbose = FALSE, nrows = 1000 )
# meta <- read.csv("/global/scratch/users/arphillips/data/gbs2ploidy/flow_cyt_predictions.csv")

# (1) Load and process VCF ----
vcf <- read.vcfR(as.character(argv$vcf), verbose = FALSE)
vcf

# (2) Load and process metadata ----
meta <- read.csv(as.character(argv$meta))
head(meta)

# (3) Subset vcf to only ploidy members ----
keep_mat <- colnames(vcf@gt) %in% meta$sample[meta$ploidy_call == ploidy_level]
keep_mat[1] <- TRUE
vcf@gt <- vcf@gt[,keep_mat]

paste0("Number of ", ploidy_level, " genotypes: ", sum(keep_mat)-1)

# sites <- str_split(colnames(refmat), pattern = "_", simplify = T)[,1]

# key <- cbind(colnames(refmat), sites) %>% as.data.frame()
# colnames(key) <- c("tree", "site")

# key_ploidy <- merge(x = key, y = meta, by.x = "site", by.y = "Site_Code") %>%
# dplyr::select(site, tree, Ploidy_level)

# (4) Extract matrices ----
# Extract allele depth
ad <- extract.gt(vcf, element = "AD", as.numeric = F)

# Create input matrix of depth of reference read
refmat <- masplit(ad, record = 1, sort = 0)

# Create input matrix of total read depth
alt <- masplit(ad, record = 2, sort = 0)

sizemat <- refmat + alt

# (Alt) Specify ploidy for all genotypes
#ploidy_level <- as.character('diploid')

# (5) Genotype ----
if (ploidy_level == 'diploid'){
  mout <- multidog(refmat = refmat,
                   sizemat = sizemat,
                   ploidy = 2,
                   model = "norm",
                   nc = cores)

  # write.table(genomat, file=paste0(outdir, "/updog.genomat.diploid.", Sys.Date(),".txt"), quote = F)

} else if (ploidy_level == 'triploid'){

  mout <- multidog(refmat = refmat, 
                   sizemat = sizemat, 
                   ploidy = 3, 
                   model = "norm",
                   nc = cores)
  
} else {
  print("ERROR: --ploidy can only be diploid or triploid")
}

saveRDS(mout, paste0(outdir,"/mout.",ploidy_level,".raw.rds"))

# (6) Filter SNPs based on updog recommendations ----
mout_cleaned <- filter_snp(mout, prop_mis < 0.05 & bias > 0.5 & bias < 2)

# Extract genotype matrix
genomat <- format_multidog(mout_cleaned, varname = "geno")

# Save filtered genotype matrix
write.table(genomat, 
            file = paste0(outdir, "/updog.genomat.", ploidy_level, ".", chr, ".",Sys.Date(),".txt"), 
            quote = F)

# Examine genotype and SNP Quality
pdf(paste0(outdir,"/genotype_depth_distributions.",ploidy_level,".", chr, ".",Sys.Date(),".pdf"))
## The (posterior) proportion of individuals mis-genotyped at each site
hist(mout$snpdf$prop_mis, main = "Proportion of ind mis-genotyped at each SNP")

## Overdispersion of each snp - simulations suggest dropping > 0.05
hist(mout$snpdf$od, main = "Overdispersion of each SNP")

## Bias - simulations suggest filtering 0.5 < x > 2
hist(mout$snpdf$bias, main = "Bias at each SNP")

dev.off()

# (7) Look at linkage between sites ----
varnames <- names(mout_cleaned$inddf)[grepl(x = names(mout_cleaned$inddf), pattern = "Pr_")]

# larray <- format_multidog(x = mout_cleaned, varname = varnames)
gp <- format_multidog(x = mout_cleaned, varname = varnames)

# like_ld <- mldest(geno = larray, 
#                   K = ploidy_level, 
#                   type = "comp", 
#                   nc = cores)
ldout <- ldfast(gp = gp, type = "r2")

pdf(paste0(outdir,"/LD_corrplot.",ploidy_level,".", chr, ".",Sys.Date(),".pdf"))
corrplot(corr = ldout$ldmat, 
         method = "color", 
         type = "upper",
         diag = FALSE,
         tl.pos = "n")
dev.off()

write.table(ldout$ldmat, 
            file = paste0(outdir, "/ldsep.r2.", ploidy_level, ".", chr, ".",Sys.Date(),".txt"), 
            quote = F)

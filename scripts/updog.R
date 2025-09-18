# Title: Updog genotyping - WGS
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
library(BIGr)

# Argument name
ap <- arg_parser("Updog genotyping")

# add mandatory positional arguments (filename)
ap <- add_argument(ap, "vcf", help = "Input VCF")

# add additional arguments
ap <- add_argument(ap, "--meta", help = "Metadata file")
ap <- add_argument(ap, "--ploidy", help = "choose 'diploid' or 'triploid'")
ap <- add_argument(ap, "--cores", help = "Number of available cores")
# ap <- add_argument(ap, "--outdir", help = "Output directory")
# ap <- add_argument(ap, "--chr", help = "Region")
ap <- add_argument(ap, "--prefix", help = "Prefix for output files")

# parse arguments
argv <- parse_args(ap)

ploidy_level <- as.character(argv$ploidy)
#chr <- as.character(argv$chr)
cores <- as.numeric(argv$cores)
#outdir <- as.character(argv$outdir)
prefix <- as.character(argv$prefix)

## test arguments
# ploidy_level <- as.character("diploid")
# # chr <- as.character("Chr02")
# cores <- as.numeric(1)
# # outdir <- as.character("/global/scratch/users/arphillips/data/updog")
# prefix <- as.character("/global/scratch/users/arphillips/data/updog/test")
# vcf <- read.vcfR("/global/scratch/users/arphillips/data/processed/filtered_snps/wgs_aspen.Chr09:3000000-4000000.goodg.10dp90.vcf.gz",
#         verbose = FALSE, nrows = 1000 )
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

# saveRDS(mout, paste0(outdir,"/mout.",ploidy_level,".raw.rds"))
saveRDS(mout, paste0(prefix,".raw.rds"))

# (6) Filter SNPs based on updog recommendations ----
mout_cleaned <- filter_snp(mout, prop_mis < 0.05 & bias > 0.5 & bias < 2)

# Extract genotype matrix
genomat <- format_multidog(mout_cleaned, varname = "geno")

# Save filtered genotype matrix
# write.table(genomat, 
#             file = paste0(outdir, "/updog.genomat.", ploidy_level, ".", chr, ".",Sys.Date(),".txt"), 
#             quote = F)
write.table(genomat, 
            file = paste0(prefix,".txt"), # what the snakemake rule asks for
            quote = F)

# Examine genotype and SNP Quality
# pdf(paste0(outdir,"/genotype_depth_distributions.",ploidy_level,".", chr, ".",Sys.Date(),".pdf"))
pdf(paste0(prefix,".genotype_depth_distributions",".pdf"))
## The (posterior) proportion of individuals mis-genotyped at each site
hist(mout$snpdf$prop_mis, main = "Proportion of ind mis-genotyped at each SNP")

## Overdispersion of each snp - simulations suggest dropping > 0.05
hist(mout$snpdf$od, main = "Overdispersion of each SNP")

## Bias - simulations suggest filtering 0.5 < x > 2
hist(mout$snpdf$bias, main = "Bias at each SNP")

dev.off()

# (6b) Export as VCF ----
updog2vcf <- function (multidog.object, output.file, updog_version = NULL, 
                       compress = TRUE) {
  mout <- multidog.object
  ploidy <- as.numeric(unique(multidog.object$snpdf$ploidy))
  if (!grepl(".vcf", output.file)) 
    output.file <- paste0(output.file, ".vcf")
  model_selected <- unique(multidog.object$snpdf$model)
  updog_meta <- paste0("##UpdogCommandLine.multidog=<ID=Multidog,Version=\"", 
                       updog_version, "\",CommandLine=\"> multidog(refmat = matrices$ref_matrix, sizemat = matrices$size_matrix, ploidy = ", 
                       ploidy, ", model = ", model_selected, ")\">")
  bigr_meta <- paste0("##BIGrCommandLine.updog2vcf=<ID=updog2vcf,Version=\"", 
                      packageVersion("BIGr"), "\",Data=\"", Sys.time(), "\", CommandLine=\"> updog2vcf(", 
                      deparse(substitute(multidog.object)), ",", output.file, 
                      ",", updog_version, ")\">")
  vcf_header <- c("##fileformat=VCFv4.3", "##reference=NA", 
                  "##contig=<ID=NA,length=NA>", "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">", 
                  "##INFO=<ID=ADS,Number=R,Type=Integer,Description=\"Depths for the ref and each alt allele in the order listed\">", 
                  "##INFO=<ID=BIAS,Number=1,Type=Float,Description=\"The estimated allele bias of the SNP from updog\">", 
                  "##INFO=<ID=OD,Number=1,Type=Float,Description=\"The estimated overdispersion parameter of the SNP from updog\">", 
                  "##INFO=<ID=PMC,Number=1,Type=Float,Description=\"The estimated proportion of individuals misclassified in the SNP from updog\">", 
                  "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype, where 1 is the count of alternate alleles\">", 
                  "##FORMAT=<ID=UD,Number=1,Type=Integer,Description=\"Dosage count of reference alleles from updog, where 0 = homozygous alternate\">", 
                  "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read depth\">", 
                  "##FORMAT=<ID=RA,Number=1,Type=Integer,Description=\"Reference allele read depth\">", 
                  "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed\">", 
                  "##FORMAT=<ID=MPP,Number=1,Type=Float,Description=\"Maximum posterior probability for that dosage call from updog\">", 
                  updog_meta, bigr_meta)
  if (any(model_selected %in% c("f1", "f1pp", "s1", "s1pp"))) {
    mout$snpdf$maxpostprob <- apply(mout$snpdf[, grep("Pr_", 
                                                      colnames(mout$snpdf))], 1, max)
    if (model_selected == "f1" | model_selected == "f1pp") {
      sele_col <- c("snp", "p1ref", "p1size", "p2ref", 
                    "p2size", "p1geno", "p2geno", "maxpostprob")
      parents <- pivot_longer(mout$snpdf[, sele_col], cols = c("p1geno", 
                                                               "p2geno"), names_to = "ind", values_to = "geno")
      parents.depth <- apply(parents, 1, function(x) {
        if (grepl("p1", x[7])) 
          return(data.frame(ref = x[2], size = x[3], 
                            check.names = FALSE))
        else if (grepl("p2", x[7])) 
          return(data.frame(ref = x[4], size = x[5], 
                            check.names = FALSE))
        else return(data.frame(ref = NA, size = NA, check.names = FALSE))
      })
      parents.depth <- do.call(rbind, parents.depth)
      parents <- cbind(parents[, c(1, 6:8)], parents.depth)
      parents$ind <- gsub("p1geno", "parent1", parents$ind)
      parents$ind <- gsub("p2geno", "parent2", parents$ind)
    }
    else {
      sele_col <- c("snp", "maxpostprob", "pgeno", "p1ref", 
                    "p1size")
      parents <- data.frame(mout$snpdf[, sele_col[1:2]], 
                            ind = "parent", mout$snpdf[, sele_col[3:5]], 
                            check.names = FALSE)
      colnames(parents)[4:6] <- c("geno", "ref", "size")
    }
    inddf <- mout$inddf[, c(1, 7, 2, 3, 4, 5)]
    inddf <- rbind(parents, inddf)
    inddf$ref <- as.numeric(inddf$ref)
    inddf$size <- as.numeric(inddf$size)
  }
  else {
    inddf <- mout$inddf[, c(1, 7, 2, 3, 4, 5)]
  }
  depth_df <- inddf %>% group_by(snp) %>% summarize(total_ref = sum(ref), 
                                                    total_size = sum(size), total_alt = sum(size) - sum(ref))
  depth_df <- depth_df %>% arrange(match(snp, mout$snpdf$snp))
  new_df <- mout$snpdf %>% tidyr::separate(snp, into = c("CHROM", 
                                                         "POS", "other"), sep = "_") %>% select(CHROM, POS)
  new_df$POS <- sub("^0+", "", new_df$POS)
  vcf_df <- data.frame(CHROM = new_df$CHROM, POS = new_df$POS, 
                       # ID = mout$snpdf$snp,
                       ID = ".",
                       ### NEED TO CHANGE TO REAL  REF ALT ----
                       REF = "A", ALT = "B", QUAL = ".", 
                       FILTER = ".", INFO = NA, FORMAT = NA, check.names = FALSE)
  vcf_df$INFO <- paste0("DP=", depth_df$total_size, ";", "ADS=", 
                        depth_df$total_ref, ",", depth_df$total_alt, ";", "BIAS=", 
                        mout$snpdf$bias, ";", "OD=", mout$snpdf$od, ";", "PMC=", 
                        mout$snpdf$prop_mis)
  vcf_df$FORMAT <- paste("GT", "UD", "DP", "RA", "AD", "MPP", 
                         sep = ":")
  convert_dosage_to_genotype <- function(dosage, ploidy) {
    if (is.na(dosage)) {
      return("./.")
    }
    ref_count <- dosage
    alt_count <- ploidy - dosage
    genotype <- paste(c(rep("0", ref_count), rep("1", alt_count)), 
                      collapse = "/")
    return(genotype)
  }
  geno_df <- inddf[, c("snp", "geno")] %>% mutate(genotype = sapply(geno, 
                                                                    convert_dosage_to_genotype, ploidy = as.numeric(ploidy)))
  format_df <- data.frame(snp = inddf$snp, ind = inddf$ind, 
                          format = paste0(geno_df$genotype, ":", inddf$geno, ":", 
                                          inddf$size, ":", inddf$ref, ":", inddf$ref, ",", 
                                          (inddf$size - inddf$ref), ":", inddf$maxpostprob), 
                          check.names = FALSE)
  format_wide <- format_df %>% pivot_wider(names_from = ind, 
                                           values_from = format)
  vcf_df <- cbind(vcf_df, format_wide[, -1])
  colnames(vcf_df)[1] <- "#CHROM"
  if (!compress) {
    writeLines(vcf_header, con = output.file)
    suppressWarnings(write.table(vcf_df, file = output.file, 
                                 sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE, 
                                 append = TRUE))
  }
  else {
    temp_loc <- tempfile()
    writeLines(vcf_header, con = temp_loc)
    suppressWarnings(write.table(vcf_df, file = temp_loc, 
                                 sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE, 
                                 append = TRUE))
    suppressMessages(Rsamtools::bgzip(temp_loc, dest = paste0(output.file, 
                                                              ".gz"), overwrite = FALSE))
  }
}
updog2vcf(mout_cleaned, 
          paste0(prefix, ".vcf"),
          updog_version = '2.1.5', compress = F)

# (7) Look at linkage between sites ----
varnames <- names(mout_cleaned$inddf)[grepl(x = names(mout_cleaned$inddf), pattern = "Pr_")]

# larray <- format_multidog(x = mout_cleaned, varname = varnames)
gp <- format_multidog(x = mout_cleaned, varname = varnames)

# like_ld <- mldest(geno = larray, 
#                   K = ploidy_level, 
#                   type = "comp", 
#                   nc = cores)
ldout <- ldfast(gp = gp, type = "r2")

# pdf(paste0(outdir,"/LD_corrplot.",ploidy_level,".", chr, ".",Sys.Date(),".pdf"))
pdf(paste0(prefix,"LD_corrplot.",".pdf"))
corrplot(corr = ldout$ldmat, 
         method = "color", 
         type = "upper",
         diag = FALSE,
         tl.pos = "n")
dev.off()

write.table(ldout$ldmat, 
            file = paste0(prefix, ".ldsep.r2",".txt"),
            # file = paste0(outdir, "/ldsep.r2.", ploidy_level, ".", chr, ".",Sys.Date(),".txt"), 
            quote = F)

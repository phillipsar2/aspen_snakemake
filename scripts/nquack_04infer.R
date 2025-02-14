#devtools::install_github("mgaynor1/nQuack")
#devtools::install("/global/scratch/users/arphillips/toolz/nQuack")
library(nQuack)
library("argparser")

getwd()

# Argument name
ap <- arg_parser("nQuack")

# add mandatory positional arguments (filename)
ap <- add_argument(ap, "samp", help = "CSV file from step 3")

# parse arguments
argv <- parse_args(ap)

samp <- as.character(argv$samp)
print(samp)

inpath <- "/global/scratch/users/arphillips/data/nquack/processed/"
outpath <- "/global/scratch/users/arphillips/data/nquack/model_inference/"


### One at a time
temp <- as.matrix(read.csv(paste0(inpath, samp, ".csv")))
out1 <- quackNormal(xm = temp, samplename = samp, cores = 10, parallel = FALSE)
out2 <- quackBeta(xm = temp, samplename = samp, cores = 10, parallel = FALSE)
out3 <- quackBetaBinom(xm = temp, samplename = samp, cores = 10, parallel = FALSE)
allout <- rbind(out1, out2, out3)
write.csv(allout, 
            file = paste0(outpath, samp, ".csv"),
            row.names = FALSE)

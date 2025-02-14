### Title: MapDamage
### Author: Alyssa Phillips
### Date: 2/10/24

library(dplyr)

dir = "/global/scratch/users/arphillips/reports/mapdamage/"
filelist <- list.files(dir)

# Occurrences for mutation types at each position 
mut_occr <- read.table(paste0(dir,"52832.2.464561.GTGTCAAC-GTATCGTG/misincorporation.txt"),
                       header = T)
head(mut_occr)
dim(mut_occr)

###
### Frequencies of C -> T from 5' end ----
###

plot(NA, ylim = c(0, 0.05), xlim = c(0, 25),
     xlab = "Position",
     ylab = "Frequency",
     main = "Frequencies of C -> T from 5' end")

for (i in filelist){
  ctot <- read.table(paste0(dir, i, "/5pCtoT_freq.txt"),
                    header = T)
  
  lines(x = ctot$pos, y = ctot$X5pC.T)
}

# At first base
ctot_1 <- c()
for (i in filelist){
  ctot <- read.table(paste0(dir, i, "/5pCtoT_freq.txt"),
                     header = T)
  
  ctot_1 <- c(ctot_1, ctot$X5pC.T[ctot$pos == 1])
}

hist(ctot_1)

###
### Frequencies of A & G from 3' end ----
###

plot(NA, ylim = c(0, 0.05), xlim = c(0, 25),
     xlab = "Position",
     ylab = "Frequency")

for (i in filelist){
  pur <- read.table(paste0(dir, i, "/3pGtoA_freq.txt"),
                      header = T)

  lines(x = pur$pos, y = pur$X3pG.A)
}


###
### Fragment size distribution ----
###
plot(NA, ylim = c(0, 3000), xlim = c(30, 180),
     xlab = "Length",
     ylab = "Occurences")

for (i in filelist[1:10]){
  fsize <- read.table(paste0(dir, i, "/lgdistribution.txt"),
                      header = T) %>% arrange(Length)
  
  fsize_minus <- filter(fsize, Std == "-")
  fsize_plus <- filter(fsize, Std == "+")
  
  lines(x = fsize_minus$Length, y = fsize_minus$Occurences, col = "red")
  # lines(x = fsize_plus$Length, y = fsize_plus$Occurences, col = "blue")
}


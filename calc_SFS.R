# Calculate SFS

library(vcfR)
library(reshape2)
library(pinfsc50)
library(reshape2)
library(ggplot2)
library(ape)
library(stringr)
library(fuzzyjoin)
library(dplyr)
library(inbreedR)

# load vcf_file
vcf_file <- "data/inbreeding/nes_filtered.recode.vcf"
# read vcf
nes_vcf <- read.vcfR(vcf_file, verbose = FALSE )
nes_vcf
gt <- as.data.frame(extract.gt(nes_vcf), stringsAsFactors = FALSE)

ind_names <- names(gt)
gt <- t(gt)
row.names(gt) <- ind_names
# NA handling
#gt[gt == "./."] <- NA
# split columns
snp_geno <- do.call(cbind, apply(gt, 2, function(x) colsplit(x, "/", c("a","b"))))




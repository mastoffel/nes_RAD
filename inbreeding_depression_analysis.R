# inbreeding analysis

# Quantifying inbreeding
# Filter in a variety of ways and visualise inbreeding ~ missingness
# sMLH, IBCS, relatedness (IBD)

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
options(scipen=999)
library(reshape2)
library(inbreedR)
library(vcfR)
library(purrr)
source("martin.R")
# system("mkdir data/inbreeding")


# Count number of SNPs in vcf file
system("grep -v '#' data/nes_rad.vcf | wc -l")

#~~ Filter raw vcf file
#~~~~~~~~~~~~~~~~~~~~~~#
# 5% geno  # maf 0.01 # minDP 10
#~~~~~~~~~~~~~~~~~~~~~~# 
# 
# system("/home/martin/bin/vcftools --vcf data/nes_rad.vcf --out data/inbreeding/nes_filtered --min-meanDP 10 --max-meanDP 20 --minDP 10 --maxDP 30 --max-missing 0.05 --remove-indels --maf 0.01 --min-alleles 2 --max-alleles 2 --recode --recode-INFO-all")
system("grep -v '#' data/inbreeding/nes_filtered.recode.vcf | wc -l")

# 
# load vcf_file
vcf_file <- "data/inbreeding/nes_filtered.recode.vcf"
# read vcf
nes_vcf <- read.vcfR(vcf_file, verbose = FALSE )

# genotypes
nes_gt <- vcfR2loci(nes_vcf)
# convert to character
nes_gt[] <- map(nes_gt, as.character)
# individuals are row names
row.names(nes_gt)

calc_het <- function(x) {
  case_when(
    x == "0/0" | x == "1/0" ~ 0,
    x == "0/1" | x == "1/0" ~ 1
    #is.na(gt_GT) ~ NA_real_
  )
}
nes_gt_het <- nes_gt %>% 
              mutate_all(funs(calc_het))

g2_nes <- g2_snps(nes_gt_het, nperm = 1000)

missing_vals <- rowSums(is.na(nes_gt_het))
nes_het <- sMLH(nes_gt_het)
# missing data vs. heterozygosity check
plot(missing_vals, nes_het)

# create plink files from vcf
system("/home/martin/bin/vcftools --vcf data/inbreeding/nes_filtered.recode.vcf --plink --out data/inbreeding/nes_filtered")

# create plink raw file
system("/home/martin/bin/plink --file data/inbreeding/nes_filtered --make-bed --recodeAD --out data/inbreeding/nes_filtered")


get_sMLH_from_plinkraw <- function(file) {
  
  x <- fread(file, colClasses = "character")
  ids <- x$IID
  x <- dplyr::select(x, contains("HET"))
  row.names(x) <- ids
  NAs <- apply(x, 1, function(x) sum(is.na(x)))
  
  sMLH <- as.data.frame(sMLH(x))
  sMLH$ANIMAL <- ids
  sMLH$NAS <- NAs
  colnames(sMLH) <- c("sMLH", "ANIMAL", "NAs")
  sMLH 
  
}

raw_files <- paste("data/inbreeding/", list.files(path = "data/inbreeding", pattern="*.raw"), sep = "")
# sMLH dataframe
nes_sMLH <- get_sMLH_from_plinkraw(raw_files)



#~~ g2 (boot over loci)

get_g2_from_plinkraw <- function(file) {
  
  x <- fread(file, colClasses = "numeric")
  
  ids <- x$IID
  x <- dplyr::select(x, contains("HET"))
  row.names(x) <- ids

  g2 <- g2_snps(x, nperm = 0, nboot = 1000, CI = 0.95)
  g2
  
}

nes_g2 <- get_g2_from_plinkraw(raw_files)
plot(nes_g2, col = "grey")
g2
# var(filter(sMLH, Run == "GQ_DP_maf05_G_ldi_ME_SSB_hwe_LD")$sMLH)



#~~ Get Fhats

# recode bim files for GCTA

recode_bim_1chr <- function(file){
  file <- fread(file) %>%
    mutate(V1 = 1) %>%
    fwrite(file, quote = F, row.names = F,
           col.names = F, sep = " ")
}

bim_files <- paste("data/inbreeding/", list.files(path = "data/inbreeding", pattern="*.bim"), sep = "")
recode_bim_1chr(bim_files)


# get fhats using gcta
plink_files <- paste("data/inbreeding/", list.files(path = "data/inbreeding", pattern="*.ped"), sep = "")
plink_files <- lapply(plink_files, function(x) gsub(".ped", "", x))

for (i in 1:length(plink_files)){
  system(paste0("/home/martin/bin/gcta64 --bfile ", plink_files[i]," --autosome --ibc --out ", plink_files[i]," --thread-num 10"))
}

# load fhats

load_fhats <- function(file) {
  
  fhats <- fread(file, header = T)
  
}

ibc_file <- paste("data/inbreeding/", list.files(path = "data/inbreeding", pattern="*.ibc"), sep = "")
fhats <- load_fhats(ibc_file)
names(fhats)[1] <- "ANIMAL"

ibcs <- fhats %>% left_join(nes_sMLH, by = c("ANIMAL"))

# ~~ Plot g2

library(ggthemr)

ggthemr(palette = "pale", layout = "clean",
        line_weight = 0.7, text_size = 20, type = "outer")
swatch()
to_swap <- swatch()[3:4]

nes_g2$g2 / var(ibcs$Fhat1)
g2$g2 / var(filter(ibcs, Run == "GQ_DP_maf05_G_ldi_ME_SSB_hwe_LD")$sMLH)
g2$g2 / var(filter(ibcs, Run == "GQ_DP_maf05_G_ldi_ME_SSB_hwe_LD")$Fhat1)
g2$g2 / var(filter(ibcs, Run == "GQ_DP_maf05_G_ldi_ME_SSB_hwe_LD")$Fhat2)
g2$g2 / var(filter(ibcs, Run == "GQ_DP_maf05_G_ldi_ME_SSB_hwe_LD")$Fhat3)

#View(filter(ibcs, Run == "GQ_DP_maf05_G_ldi_ME_SSB_LD"))

g2_plot <- data.frame(nes_g2$g2_boot)
lcl <- nes_g2$CI_boot[1]
ucl <- nes_g2$CI_boot[2]
g2_boot_summary <- data.frame(lcl, ucl)

ibcs_vars_summary <- data.frame(c(nes_g2$g2,
                                  var(ibcs$sMLH),
                                  var(ibcs$Fhat1),
                                  var(ibcs$Fhat2),
                                  var(ibcs$Fhat3)),
                                c("g2", "sMLH", "Fhat1", "Fhat2", "Fhat3"))

colnames(ibcs_vars_summary) <- c("val", "var")

plot(nes_g2)


# g2 bootstrapping distribution showing empirical g2 with CIs

require(gridExtra)
library(sitools)
cbPalette <- c( "#1B9E77", "#66A61E", "#E6AB02", "black", "#7570B3", "#D95F02", "#E7298A")

#png("figs/g2_boot.png", units = "in", res = 300, width = 8, height = 7)
# remove F2 and F1
ibcs_vars_summary <- ibcs_vars_summary %>% 
                        filter(!(var %in% c("Fhat1", "Fhat2")))
cbPalette <- c( "#1B9E77", "#66A61E", "#E6AB02", "black", "#7570B3")
g2_CI_plot <-
  ggplot(g2_plot, aes(nes_g2$g2_boot)) +
  geom_histogram(colour = "grey45", fill = "grey45") +
  geom_errorbarh(aes(xmin = g2_boot_summary$lcl , xmax = g2_boot_summary$ucl , y = 90),
                 size = 0.8, color = "black", linetype = "solid", height = 0) +
  geom_vline(data = ibcs_vars_summary, aes(xintercept = val, colour = var), size = 0.8,
             linetype = c("dashed", "solid", "solid"), show.legend = T) +
  scale_colour_manual(values = cbPalette, name = "",
                      breaks = c("g2","Fhat3", "sMLH"),
                      labels = c(expression(italic(g[2])),
                                 expression("var"(italic(hat(F)["III"]))),
                                 expression("var"("sMLH")))) +
  labs(y = "Counts", x = expression(italic(g[2]))) +
  theme_martin() +
  theme(axis.text.x = element_text(face = "plain"),
        axis.text.y = element_text(face = "plain"),
        axis.title.x = element_text(face = "plain"),
        axis.title.y = element_text(face = "plain"),
        plot.title=element_text(hjust=0, size = 18, face = "plain")) 


g2_CI_plot
#dev.off()





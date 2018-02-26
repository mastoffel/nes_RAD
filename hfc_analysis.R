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

nes_mingeno8 <- scan("nesF_min8")
nes_mingeno10 <- scan("nesF_min10")
nes_mingeno12 <- scan("nesF_min12")


nesF86_filt <- nesF86[!(nesF86 > 0.01)]
het_filt <- het[!(nesF86 > 0.01)]
plot(nesF86_filt, het_filt)

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
# convert
snp_genotypes <- inbreedR::convert_raw(snp_geno)
# check data
inbreedR::check_data(snp_genotypes)
nes <- snp_genotypes

missing_vals <- rowSums(is.na(nes))
het <- sMLH(nes)

plot(het, het_old)

g2_nes <- g2_snps(snp_genotypes, nboot = 100, nperm = 100)
g2_nes

var(het)
var(fhats$Fhat3)

plot(missing_vals, het)
plot(missing_vals, nes1)
hist(missing_vals, breaks = 50)

plot(het, het_1)

# vcf to plink
system("vcftools --vcf data/inbreeding/nes_filtered.recode.vcf --plink --out data/inbreeding/relatedness")
system("plink --file data/inbreeding/relatedness --make-bed --out data/inbreeding/relatedness")

#~~ Get Fhats

# recode bim files for GCTA
recode_bim_1chr <- function(file){
  file <- fread(file) %>%
    mutate(V1 = 1) %>%
    fwrite(file, quote = F, row.names = F,
           col.names = F, sep = " ")
}

bim_files <- paste("data/inbreeding/", list.files(path = "data/inbreeding", pattern="*.bim"), sep = "")
lapply(bim_files, recode_bim_1chr)

# get fhats using gcta

plink_files <- paste("data/inbreeding/", list.files(path = "data/inbreeding", pattern="*.ped"), sep = "")

plink_files <- lapply(plink_files, function(x) gsub(".ped", "", x))

system(paste0("/home/martin/bin/gcta64 --bfile ", plink_files[[1]]," --autosome --ibc --out ", plink_files[i]," --thread-num 10"))

# load fhats
load_fhats <- function(file) {
  fhats <- fread(file, header = T)
}
ibc_files <- paste("data/inbreeding/", list.files(path = "data/inbreeding", pattern="*.ibc"), sep = "")
fhats <- load_fhats(ibc_files)

plot(fhats$Fhat3, nesF86)
nesF86_2_filt <- nesF86_2[!(nesF86_2>0.010)]
fhat3_filt <- fhats$Fhat3[!(nesF86_2>0.010)]
plot(nesF86_2_filt, fhat3_filt)

ibcs <- rbindlist(fhats, idcol = "Run") %>%
  left_join(sMLH, by = c("IID", "Run")) %>%
  left_join(ms_sMLH, by  = "IID")





# get individual data
library(inbreedR)
library(readxl)
library(stringr)
sample_ids <- read_excel("data/sample_ids.xlsx")
# concatenate
sample_ids$full_id <- str_c(sample_ids$library, sample_ids$barcode, sep = "_")
sample_ids <- sample_ids %>% mutate(full_id = str_c(library, barcode, sep = "_")) 
ind_names <- str_replace_all(ind_names, "merged_sample_", "")

# find the correct sequence to match genotypes to data
new_sequence <- match(ind_names, sample_ids$full_id)

# resort data to match real names
nes <- nes[new_sequence, ]
# rename rows with real ids
rownames(nes) <- sample_ids$sample_id


# check in how many individuals a snp was called
typed_snp <- colSums(!is.na(nes))
hist(typed_snp)

seals <- data.frame("id" = rownames(nes), "mlh" = het)
rownames(seals) <- 1:nrow(seals)

# read in raw data
seal_data <- read_excel("data/nes_rad_data.xlsx")

seal_data <- seal_data[seal_data$`Animals ID` %in% seals$id, ]

seals <- seals %>% dplyr::mutate(groups = seal_data$group) %>% 
                   dplyr::mutate(groups_sum = if_else(groups %in% c("control1", "control2", "control3", "control4"), "control",
                                        if_else(groups == "bacteria", "bacteria", if_else(groups == "worms", "worms", "")))) %>% 
                   dplyr::mutate(blubber = seal_data$`Relative Blubber depth`) %>% 
                   dplyr::mutate(weight = as.numeric(seal_data$`death weight`)) %>% 
                   dplyr::mutate(msat_het = as.numeric(seal_data$HL)) 
        
ggplot(seals, aes(x = groups_sum, y=mlh)) + 
  geom_boxplot() +
  stat_summary(fun.y=mean, geom="point", shape=5, size=4) +
  geom_jitter()


ggplot(seals, aes(x = mlh, y=msat_het)) + geom_point() + geom_smooth(method = "lm")


plot(seal_data$`st het_Obs`, seals$mlh)
mod <- lm(seal_data$`st het_Obs` ~ seals$mlh)
summary(mod)
plot(mod)
abline()
summary(lm(mlh~groups, data = seals))


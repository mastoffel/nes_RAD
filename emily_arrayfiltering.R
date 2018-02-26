# Filtering RAD SNPs for SNP array
# May 2017

library(data.table)
library(dplyr)
library(seqinr)
library(tidyr)
library(plyr)
source("scripts/getsnpflanking.R")
source("scripts/arraydesign_seq.R")
options(scipen=999)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 1. Load raw vcf file                                   #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

getwd()

# Count number of SNPs in vcf file

system("grep -v '#' /external/data/Emily/Seal_RADseq/cebitec_2016/AFS/GATK/ArcGaz_genotype_gvcf.vcf | wc -l")

# Update vcf to include only biallelic snps

system("bcftools view -m2 -M2 -v snps /external/data/Emily/Seal_RADseq/cebitec_2016/AFS/GATK/ArcGaz_genotype_gvcf.vcf > data/ArcGaz_biallelic.vcf")
system("grep -v '#' data/ArcGaz_biallelic.vcf | wc -l") # 797768

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 2. Recode sample names                                 #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# current sample names

system("bcftools query -l data/ArcGaz_biallelic.vcf")

# use list of new names in correct order to recode

system("bcftools reheader -s data/raw/sample_names.txt -o data/ArcGaz_biallelic_rename.vcf data/ArcGaz_biallelic.vcf")

# new sample names

system("bcftools query -l data/ArcGaz_biallelic_rename.vcf")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 3. VCFtools filtering                                  #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# filter for individuals in correct triads

system("vcftools --vcf data/ArcGaz_biallelic_rename.vcf --out data/ArcGaz_biallelic_sub --keep data/raw/individuals_to_keep.txt --recode --recode-INFO-all")
system("bcftools query -l data/ArcGaz_biallelic_sub.recode.vcf")

# get depth stats
system("vcftools --vcf data/ArcGaz_biallelic_sub.recode.vcf --geno-depth --out data/ArcGaz_biallelic_sub")

depth <- fread("data/ArcGaz_biallelic_sub.gdepth", colClasses = "numeric") %>%
  unite(SNP, CHROM, POS)

depth <- dplyr::select(depth, -SNP)

depth <- depth %>%
  dplyr::mutate(mean = rowMeans(.), 
                sum = rowSums(.), 
                logmeandepth = log10(mean))

CI <- 0.90
CI_meanDP <- stats::quantile(depth$mean, c((1-CI)/2,1-(1-CI)/2), na.rm=TRUE)
CI_sumDP <- stats::quantile(depth$sum, c((1-CI)/2,1-(1-CI)/2), na.rm=TRUE)
CI_meanDP
CI_sumDP

hist(depth$logmeandepth)
hist(log10(depth$sum))
hist(depth$sum)

good_depth <- filter(depth, sum > 400)

# for the purpose of axiom scores sub_qual_maf_miss:
# min mean depth 5 (because of inbreeding analysis)
# max mean depth 18
# max missing data 60%
# maf 0.05

# SNP chip
system("vcftools --vcf data/ArcGaz_biallelic_sub.recode.vcf --out data/ArcGaz_biallelic_sub_qual_maf_miss --min-meanDP 5 --max-meanDP 18 --maf 0.05 --max-missing 0.6 --recode --recode-INFO-all")
system("vcftools --vcf data/ArcGaz_biallelic_sub_qual_maf_miss.recode.vcf --plink --out data/ArcGaz_biallelic_sub_qual_maf_miss")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 4. Get Mendel error SNPs        #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# recode map file to include contig in chr column

system("scripts/recode_full_map.sh data/ArcGaz_biallelic_sub_qual_maf_miss.map data/ArcGaz_biallelic_sub_qual_maf_miss.map")

# Make bed

system("plink --file data/ArcGaz_biallelic_sub_qual_maf_miss --make-bed --out data/ArcGaz_biallelic_sub_qual_maf_miss --allow-extra-chr --debug")

# Determine Working Pedigree for PLINK
# Pedigree code from Susan Johnston: https://goo.gl/iRw3Ro 

pedigree <- read.table("data/raw/new_pedigree.txt", header = T, 
                       colClasses = c("character", "character", "character"))
pedigree$MOTHER <- as.character(pedigree$MOTHER)
pedigree$FATHER <- as.character(pedigree$FATHER)

famped <- NULL

for(i in pedigree$ANIMAL){
  ped1 <- pedigree[which(pedigree$ANIMAL == i),]
  {
    ped2 <- ped1[,1:3]
    ped2 <- rbind(data.frame(ANIMAL = c(unlist(ped2[1, 2:3])), FATHER = 0, MOTHER = 0), ped2)
    famped <- rbind(famped, ped2)
    famped <- famped[which(famped[,1] != 0),]
    rm(ped2)
  }
}

fatherID <- famped$FATHER
fatherID <- fatherID[fatherID != 0]
fatherID <- na.omit(fatherID)
fatherID <- rep(fatherID, each = 3)
popID <- as.character(famped[c(73:98),1])

famIDs <- c(fatherID, popID)
famped$Family <- famIDs

famped <- famped %>%
  dplyr::mutate(Phenotype = 0,
                Sex = c(rep(c(2,1,"unknown"),24), rep("unknown", 26))) %>%
  dplyr::mutate(FATHER = ifelse(is.na(FATHER),0,FATHER),
                MOTHER = ifelse(is.na(MOTHER),0,MOTHER))

famped <- famped %>%
  dplyr::distinct(ANIMAL, .keep_all = T)

# Write over fam file
# Arrange rows to match plink
# COLS: FAMID, IID, FATHER, MOTHER, SEX, PHENOTYPE

fam <- fread("data/ArcGaz_biallelic_sub_qual_maf_miss.fam")


new_famped <- fam %>%
  left_join(famped, by = c("V2" = "ANIMAL"))


write.table(new_famped[c(9,1,7,8,11,10)], "data/ArcGaz_biallelic_sub_qual_maf_miss.fam", col.names = F,
            row.names = F, quote = F)

# write fam IDs for other analyses

write.table(famped[c(4,1)], "data/processed/ArcGaz_fam.txt", col.names = F,
            row.names = F, quote = F)

# write full pedigree file

write.table(famped[c(1,2,3)], "data/processed/full_pedigree.txt", col.names = F,
            row.names = F, quote = F)

# Run PLINK --mendel

system("mkdir data/mendel")
system("plink --bfile data/ArcGaz_biallelic_sub_qual_maf_miss --mendel --out data/mendel/ArcGaz_biallelic_sub_qual_maf_miss --allow-extra-chr --debug")

# Process Mendelian inconsistencies

# Load mendel output and NA errors due to missing parental genotypes

menderr_plink <- fread("data/mendel/ArcGaz_biallelic_sub_qual_maf_miss.mendel") %>%
  mutate(V6 = gsub("\\*", NA, V6), V8 = gsub("\\*", NA, V8)) %>%
  mutate(V6 = gsub("NA\\/NA", NA, V6), V8 = gsub("NA\\/NA", NA, V8)) %>%
  na.omit()

# Determine per SNP error rate

error_rate <- menderr_plink %>%
  group_by(V4) %>%
  dplyr::summarise(length(V4)) %>%
  `colnames<-`(c("V4", "N")) %>%
  filter(N > 0) %>%
  mutate(rate = (N/24)*100) # 24 N of triads

hist(error_rate$rate)

# write out list of Mendel Error SNPs

# strict
mendel_error_strict <- error_rate %>%
  filter(N > 1) %>% # wrong in one or more triads (4%)
  dplyr::select(V4) %>%
  separate(V4, c("Contig", "Position")) 

write.table(mendel_error_strict, "data/mendel/MendelErrorStrict_SNPs.txt", col.names = F,
            row.names = F, quote = F, sep = ":")
write.table(mendel_error_strict, "data/mendel/MendelErrorStrict_SNPs_tab.txt", col.names = F,
            row.names = F, quote = F, sep = "\t")

mendel_error_mod <- error_rate %>%
  filter(N > 4) %>% # wrong in more than 4 (16%)
  dplyr::select(V4) %>%
  separate(V4, c("Contig", "Position"))

write.table(mendel_error_mod, "data/mendel/MendelErrorMod_SNPs.txt", col.names = F,
            row.names = F, quote = F, sep = ":")
write.table(mendel_error_mod, "data/mendel/MendelErrorMod_SNPs_tab.txt", col.names = F,
            row.names = F, quote = F, sep = "/t")

system("wc -l data/mendel/MendelErrorMod_SNPs.txt") # 225
system("wc -l data/mendel/MendelErrorStrict_SNPs.txt") # 12942

#~~ Filter plink files for mendel errors

#system("plink --bfile data/ArcGaz_biallelic_sub_qual_maf_miss --exclude data/mendel/MendelError_SNPs.txt --nonfounders --recode --make-bed --out data/mendel/ArcGaz_biallelic_sub_qual_maf_miss_me --allow-extra-chr")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#    Select SNPs for chip      #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 1. Remove SNPs close to start or end of scaffold       #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# load scaffold length data
lengths <- fread("data/raw/PBJelly_lengths.txt", header = F) %>%
  `colnames<-`(c("Contig", "Length"))

# get SNP positions:

snps <- fread("data/ArcGaz_biallelic_sub_qual_maf_miss.bim") %>%
  select(V1,V4,V5,V6) %>%
  `colnames<-`(c("Contig", "Position", "Ref", "Alt"))


get_axiomflanks <- function(x, y){
  snp.pos <- x
  snp.pos$one <- snp.pos$Position - 36   
  snp.pos$two <- snp.pos$Position + 35
  snp.pos <- left_join(snp.pos, y, by = "Contig", sort = F)
  
  for (i in 1:length(snp.pos$two)){
    if (snp.pos$two[i] > snp.pos$Length[i]){
      snp.pos$two[i] <- snp.pos$Length[i]}
    if (snp.pos$one[i] < 1){
      snp.pos$one[i] <- 1
    }
  }
  
  snp.pos$SNPenclose <- (snp.pos$Position - snp.pos$one)
  snp.pos <- snp.pos # [-5]
}

axiom <- get_axiomflanks(snps, lengths) %>%
  filter(SNPenclose == 36 & Position < Length & two - Position == 35)
length(axiom[,1])


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 2. Extract flanking sequences from genome            #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

system("mkdir data/bedtools")

# write bed file with flanking sequence coordinates

write.table(axiom[c(1,5,6)], "data/bedtools/axiomflanks.bed", 
            quote = F, col.names = F, row.names = F, sep = "\t")

# run bed tools

bedCmd <- paste("bedtools getfasta -fi ~/assemblies/ArcGaz_PBJelly/final.assembly.ArcGaz002_PBJelly.fasta -bed data/bedtools/axiomflanks.bed -fo data/bedtools/axiomflanks.fasta")
system(bedCmd)

# get list of flanks
seqs <- read.fasta("data/bedtools/axiomflanks.fasta", as.string = T, forceDNAtolower = F) %>%
  lapply(function(x) x[[1]]) %>%
  as.character(unlist(.)) 

axiom <- axiom %>%
  mutate(seq = seqs)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 4. Make file for Axiom p-convert scores       #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

system("mkdir data/affy")
x <- axiom
x$seq <- toupper(x$seq) # make upper case for string matching

# get snp according to sequence
axiom_snps <- data.frame(Refdf = substr(x$seq, x$SNPenclose, x$SNPenclose))

# get alt alleles according to allele in seq
for (i in 1:length(x$seq)){
  if (axiom_snps$Refdf[i] == x$Ref[i]){
    axiom_snps$Altdf[i] <- as.character(x$Alt[i])}
  else {
    axiom_snps$Altdf[i] <- as.character(x$Ref[i])
  }
}

# alphabetical
axiom_snps <- as.data.frame(t(apply(axiom_snps, 1, sort)))
colnames(axiom_snps) <- c("Refdf", "Altdf")
x <- data.frame(x, axiom_snps)

# enclose SNP
for (i in 1:length(x$seq)){
  substr(x$seq[i], x$SNPenclose[i], x$SNPenclose[i]) <- "]"
}

to <- paste("[", x$Refdf, "/", x$Altdf, "]", sep = "")

for(i in 1:length(x$seq))
  x$seq[i]<-gsub("]", to[i], x$seq[i])

x <- mutate(x, SNPid = paste(Contig, Position, sep = '_'))

axiom_df <- data.frame(Organism = "Arctocephalus_gazella",
                       Locus_Name = x$SNPid,
                       SEQ = x$seq,
                       SNP_PRIORITY = 1,
                       CHR = "unknown",
                       CHR_TYPE = "autosomal",
                       SNP_VAL = 0)

write.table(axiom_df, "data/affy/axiomdesign_Seq_FurSeal_RAD.txt",
            quote = F, col.names = T, row.names = F, sep = "\t")


# cat together with transcriptome SNP axiom file from transcriptome paper

trans_snps <- fread("data/raw/axiomdesign_Seq_FurSeal.txt")
all_snps_axiom <- rbind(axiom_df, trans_snps)

write.table(all_snps_axiom, "data/affy/axiomdesign_Seq_FurSeal_RAD_transcriptome.txt",
            quote = F, col.names = T, row.names = F, sep = "\t")




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#  3. BLAST SNP flanking sequences to genome      #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# write file of query seqs

system("mkdir data/blast")

get.fasta <- function(x){
  names <- do.call(paste, c(x[c(1,2)], sep = "_"))
  names <- data.frame(paste(">", names, sep = ""))
  fasta <- cbind(names, x$seq)
  fasta <- as.vector(t(fasta))
  fasta <- as.data.frame(fasta)
  write.table(fasta, paste("data/blast/", substitute(x), ".fasta", sep = ""), quote = F, col.names = F, row.names = F)
}

get.fasta(axiom)

# BLAST query seqs to fur seal genome
# fur seal blastdb in data/blast/

infilenames <- list.files(path = "data/blast", pattern = "*.fasta")
outnames <- paste(unlist(sapply(infilenames, strsplit, split = "*.fasta")),
                  "SNPs2genome", sep = "")

# function to create the commands
cmdCreate <- function(infile, outfile){
  paste("nohup ~/programs/blastn -db data/blast/ArcGaz_PBJelly -outfmt 6 -num_threads 32 -evalue 1e-12 -query data/blast/",infile, " -out data/blast/",
        outfile," &", sep = "")
}

# create the commands
cmds <- mapply(FUN = cmdCreate, infile = infilenames, outfile = outnames)

# run the blasts. this will obviously take a while when using the full genome

sapply(cmds, system)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#      Transcriptome Data BLAST        #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# blast transcriptome seqs to new genome
# extracting seqs from axiom file: 31590

trans_fasta <- trans_snps %>%
  select(Locus_Name, SEQ) %>%
  mutate(SEQ = gsub("\\[", "", SEQ)) %>%
  mutate(SEQ = gsub("\\/[A-Z]\\]", "", SEQ)) %>%
  mutate(Locus_Name = paste0(">", Locus_Name)) 

trans_fasta <- as.vector(t(trans_fasta))
write.table(trans_fasta, "data/transcriptome/flanks_for_blast.fasta", quote = F, col.names = F, row.names = F)


infilenames <- list.files(path = "data/transcriptome", pattern = "flanks_for_blast.fasta")
outnames <- paste(unlist(sapply(infilenames, strsplit, split = "*.fasta")),
                  "SNPs2genome", sep = "")

# function to create the commands
cmdCreate <- function(infile, outfile){
  paste("nohup ~/programs/blastn -db data/blast/ArcGaz_PBJelly -outfmt 6 -num_threads 32 -evalue 1e-12 -query data/transcriptome/",infile, " -out data/transcriptome/",
        outfile," &", sep = "")
}

# create the commands
cmds <- mapply(FUN = cmdCreate, infile = infilenames, outfile = outnames)

# run the blasts. this will obviously take a while when using the full genome

sapply(cmds, system)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 4. Prepare summary blast df            #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Read in blast output and prepare for genomic context filtering

blasts <-  fread("data/blast/axiomSNPs2genome") %>%
  `colnames<-`(c("SNPid","Subject","PercentIdentity", "AlignmentLength", 
                 "Mismatches", "Gap_Opening", "QueryStart", "QueryEnd", 
                 "SubjectStart", "SubjectEnd", "E.value", "BitScore")) 

# Count number of hits

hits <- blasts %>%
  mutate(Count = 1) %>%
  dplyr::group_by(SNPid) %>%
  dplyr::summarise(hits = sum(Count))

blasts <- dplyr::inner_join(hits, blasts, by = "SNPid")

# Extract top blast hit

snps <- blasts[!duplicated(blasts$SNPid),]

# Identify SNPs mapping uniquely and completely

snps <- snps %>%
  mutate(good_mapping = ifelse(hits == 1 & AlignmentLength == 71, 1, 0)) %>%
  separate(SNPid, c("Contig", "Position"))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#  Prepare summary Transcriptome blast df          #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

Tblasts <-  fread("data/transcriptome/flanks_for_blastSNPs2genome") %>%
  `colnames<-`(c("SNPid","Subject","PercentIdentity", "AlignmentLength", 
                 "Mismatches", "Gap_Opening", "QueryStart", "QueryEnd", 
                 "SubjectStart", "SubjectEnd", "E.value", "BitScore")) 

# Count number of hits

Thits <- Tblasts %>%
  mutate(Count = 1) %>%
  dplyr::group_by(SNPid) %>%
  dplyr::summarise(hits = sum(Count))

Tblasts <- dplyr::inner_join(Thits, Tblasts, by = "SNPid")

# Extract top blast hit

Tsnps <- Tblasts[!duplicated(Tblasts$SNPid),] 

# Filter for SNPs mapping uniquely and completely

Tsnps <-Tsnps %>%
  mutate(good_mapping = ifelse(hits == 1 & AlignmentLength == 71, 1, 0)) %>%
  mutate(SNPid = gsub("_v1.1", ".v1.1", SNPid)) %>%
  separate(SNPid, c("Contig", "Position"), sep = "_") %>%
  mutate(Contig = gsub(".v1.1", "_v1.1", Contig))

# Summary blast dataframe for RAD and transcriptome seqs
# Including unique mapping score (good_mapping) :

all_snps <- rbind(snps, Tsnps) %>%
  mutate(Position = as.numeric(Position))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#  2. Write list of SNPs with SNPs in flanks   #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

all_snps <- read.table("data/affy/axiomdesign_Seq_FurSeal_RAD_transcriptome.txt", header = T) %>%
  mutate(Locus_Name = gsub("_v1.1", ".v1.1", Locus_Name)) %>%
  separate(Locus_Name, c("Contig", "Position"), sep = "_") %>%
  mutate(Contig = gsub(".v1.1", "_v1.1", Contig)) %>%
  mutate(Position = as.numeric(Position)) %>%
  left_join(all_snps, by =c("Contig", "Position")) %>%
  select(-c(SEQ, CHR, CHR_TYPE, SNP_VAL))

# Thin snps to those without variants within 20 bp (Axiom recommendations) 
# If one flank is clean then SNP can still be used

# split df into list of Contigs
spacing <- split(all_snps , f = all_snps$Contig)

# arrange by SNP order within each contig
spacing <- lapply(spacing, function(x) 
  arrange(x, Position)) 

thinned <- lapply(spacing, function(x) 
  mutate(x, distleft = (x$Position - lag(x$Position))))

thinned <- lapply(thinned, function(x) 
  mutate(x, distright = (lead(x$Position) - x$Position)))

# lag to get distance of first SNP on contig
thinned <- lapply(thinned, function(x)
  mutate(x, distleft = ifelse(is.na(distleft),distright,distleft),
         distright = ifelse(is.na(distright),distleft,distright)))

library(dplyr)
thinned <- ldply(thinned)


# combine distance data with full dataframe
# add columns indicating presence of snps in left & right flanks\

all_snps <- thinned %>%
  mutate(good_leftflank = ifelse(distleft > 20 | is.na(distleft), 1, 0),
         good_rightflank = ifelse(distright > 20| is.na(distright), 1, 0)) %>%
  #unite(SNPid, Contig, Position) %>%
  left_join(all_snps, by = c("Contig", "Position"))


#~~ Write new plink files

#filter <- thinned %>%
#  mutate(id = paste(Contig, Position, sep = ":")) %>%
#  select(id)

#write.table(filter, "data/plink/thinned.txt",
#            quote = F, col.names = F, row.names = F, sep = "\t")

#system("plink --bfile data/plink/ArcGaz_biallelic_rename_sub_fixfilt_mendel_maf --extract data/plink/thinned.txt --recode --make-bed --out data/plink/ArcGaz_biallelic_rename_sub_fixfilt_mendel_maf_thinned --allow-extra-chr")



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#      Identify SNPs with Mendel Errors        #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

me_mod <- fread("data/mendel/MendelErrorMod_SNPs.txt") %>%
  #unite(SNPid, V1, V2) %>%
  `colnames<-`(c("Contig", "Position")) %>%
  mutate(me_modPass = 0)

me_strict <- fread("data/mendel/MendelErrorStrict_SNPs.txt") %>%
  #unite(SNPid, V1, V2) %>%
  `colnames<-`(c("Contig", "Position")) %>%
  mutate(me_strictPass = 0)

# add column to summary df stating whether SNP is ME

snp_info_df <- all_snps %>%
  left_join(me_mod, by = c("Contig", "Position")) %>%
  left_join(me_strict, by = c("Contig", "Position")) %>%
  mutate(me_modPass = ifelse(is.na(me_modPass), 1 ,me_modPass)) %>%
  mutate(me_strictPass = ifelse(is.na(me_strictPass), 1 ,me_strictPass)) %>%
  select(Contig, Position, good_leftflank, good_rightflank, good_mapping.x, me_modPass, me_strictPass)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#        Identify A/T and G/C SNPs             #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~ A/T & C/G SNPs take up twice as much room on array
#~~ Ideally we avoid them

extract <- "\\[[A-Z]\\/[A-Z]\\]"

#~~ get snp alleles and prepare dataframe

library(stringr)
snp_types <- read.table("data/affy/axiomdesign_Seq_FurSeal_RAD_transcriptome.txt", header = T) %>%
  select(Locus_Name, SEQ) %>%
  mutate(SEQ = str_extract(SEQ, extract)) %>%
  mutate(SEQ = gsub("\\[", "", SEQ)) %>%
  mutate(SEQ = gsub("\\]", "", SEQ)) %>%
  separate(SEQ, c("Ref", "Alt")) %>%
  mutate(Locus_Name = gsub("_v1.1", ".v1.1", Locus_Name)) %>%
  separate(Locus_Name, c("Contig", "Position"), sep = "_") %>%
  mutate(Contig = gsub(".v1.1", "_v1.1", Contig)) %>%
  mutate(Position = as.numeric(Position))

snp_info_df <- snp_info_df %>%
  left_join(snp_types, by = c("Contig", "Position")) 

snp_info_df <- snp_info_df %>%
  mutate(AT = ifelse(Ref == "A" & Alt == "T", 1, 0),
         CG = ifelse(Ref == "C" & Alt == "G", 1, 0))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#          Incorporate Affy Scores             #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# read 182,652 SNP flanking seqs scored by Affy

pcon <- read.table("data/affy/axiomdesign_Seq_FurSeal_RAD_transcriptome_scored.txt", header = T) %>%
  select(Snpid, forwardPconvert, forwardRecommendation, reversePconvert, reverseRecommendation) %>%
  mutate(Snpid = gsub("_v1.1", ".v1.1", Snpid)) %>%
  separate(Snpid, c("Contig", "Position"), sep = "_") %>%
  mutate(Contig = gsub(".v1.1", "_v1.1", Contig)) %>%
  mutate(Position = as.numeric(Position))

# combine with summary df
# note: flanking seqs that didn't map to the genome are included but have NAs for genome mapping

snp_info_df <-  pcon %>%
  left_join(snp_info_df, by = c("Contig", "Position"))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#          Incorporate Transcriptome Annotations        #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

annotations <- read.table("data/processed/Transcript_annotations.txt", header =T)

snp_info_df <- snp_info_df %>%
  left_join(annotations, by = c("Contig" = "Contig_Name")) %>%
  mutate(immune = ifelse(is.na(immune),0,immune),
         growth = ifelse(is.na(growth),0,growth),
         metabolism = ifelse(is.na(metabolism),0,metabolism))


# save file
system("mkdir data/processed")
write.table(snp_info_df, "data/processed/SNP_info_df.txt", 
            row.names = F, col.names = T, quote = F)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#              Identify SNPs for SNP chip               #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

nrow(filter(snp_info_df, me_modPass == 1 & good_mapping.x == 1 & (good_leftflank == 1 | good_rightflank ==1) & AT == 0 & CG == 0))
nrow(filter(snp_info_df, good_mapping.x == 1 & (good_leftflank == 1 | good_rightflank ==1) & AT == 0 & CG == 0))

nrow(filter(snp_info_df, good_mapping.x == 1 & (good_leftflank == 1 | good_rightflank ==1) & AT == 0 & CG == 0 & (forwardRecommendation == "recommended" | reverseRecommendation == "recommended")))
nrow(filter(snp_info_df, (forwardRecommendation == "recommended" | reverseRecommendation == "recommended")))
nrow(filter(snp_info_df, forwardRecommendation == "recommended" & reverseRecommendation == "recommended"))


# investigating

recommended <- filter(snp_info_df, (forwardRecommendation == "recommended" | reverseRecommendation == "recommended"))
not_recommended <- filter(snp_info_df, (forwardRecommendation != "recommended" & reverseRecommendation != "recommended"))

sum(recommended$good_mapping.x, na.rm = T) / length(recommended$good_mapping.x)
sum(not_recommended$good_mapping.x, na.rm = T) / length(not_recommended$good_mapping.x)

table(not_recommended$good_mapping.x)
table(recommended$good_mapping.x)






#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#      Prepare RAD SNP seqs for validation         #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Get 100 bp flanks

validation <- snp_info_df %>%
  #filter(grepl("^Contig", Contig)) %>%
  mutate(Position = as.numeric(Position)) %>%
  left_join(axiom[c(1:2)], by = c("Contig", "Position"))

get_validation_flanks <- function(x, y){
  snp.pos <- x
  snp.pos$one <- snp.pos$Position - 101   
  snp.pos$two <- snp.pos$Position + 100
  snp.pos <- left_join(snp.pos, y, by = "Contig", sort = F)
  
  for (i in 1:length(snp.pos$two)){
    if (snp.pos$two[i] > snp.pos$Length[i]){
      snp.pos$two[i] <- snp.pos$Length[i]}
    if (snp.pos$one[i] < 1){
      snp.pos$one[i] <- 1
    }
  }
  
  snp.pos$SNPenclose <- (snp.pos$Position - snp.pos$one)
  snp.pos <- snp.pos # [-5]
}


validation <- get_validation_flanks(validation, lengths) %>%
  filter(SNPenclose == 101 & Position < Length & two - Position == 100)

length(validation[,1])

system("mkdir data/validation")

# write bed file with flanking sequence coordinates

write.table(validation[c(1,11,12)], "data/validation/validationflanks.bed", 
            quote = F, col.names = F, row.names = F, sep = "\t")

# run bed tools

bedCmd <- paste("bedtools getfasta -fi ~/assemblies/ArcGaz_PBJelly/final.assembly.ArcGaz002_PBJelly.fasta -bed data/validation/validationflanks.bed -fo data/validation/validationflanks.fasta")
system(bedCmd)

# get list of flanks
val_seqs <- read.fasta("data/validation/validationflanks.fasta", as.string = T, forceDNAtolower = F) %>%
  lapply(function(x) x[[1]]) %>%
  as.character(unlist(.)) 

validation <- validation %>%
  mutate(seq = val_seqs)

# enclose SNP for primer design

for (i in 1:length(validation$seq)){
  substr(validation$seq[i], validation$SNPenclose[i], validation$SNPenclose[i]) <- "]"
}

id <- paste("[", validation$Ref, "]", sep = "")

for(i in 1:length(validation$seq))
  validation$seq[i]<-gsub("]", id[i], validation$seq[i])


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#         Select SNPs for validation           #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Good SNPs
# Unique
# no AT or CG
# Pass mendel moderate
# No secondary SNPs in left of right flank

good_snps_ME <- filter(validation, me_modPass == 1 & good_mapping == 1 & (good_leftflank == 1 | good_rightflank ==1) & AT == 0 & CG == 0) 
length(validation$me_modPass == 1)

# select random number

set.seed(100)
sample <- sample(1:nrow(good_snps_ME), 50)
for_validation <- good_snps_ME[sample,] %>%
  write.table("data/validation/RAD_validation_SNPs.txt", quote = F, col.names = F, row.names = F)

good_snps_ME[sample,] %>%
  unite(ID, Contig, Position, sep = ":") %>%
  .[1] %>%
  write.table("data/validation/validationIDs.txt", quote = F, col.names = F, row.names = F)

# extract selected SNPs from plink files in het format

system("plink --file data/ArcGaz_biallelic_sub_qual_maf_miss --extract data/validation/validationIDs.txt --out data/validation/validation --nonfounders --recodeA --allow-extra-chr --debug")

# get genotypes 0/1



# one SNP flank had too high GC content. Replace.

write.table(good_snps_ME[1,], "data/validation/RAD_validation_SNPs_extra.txt", quote = F, col.names = F, row.names = F)

system("plink --file data/ArcGaz_biallelic_sub_qual_maf_miss --extract data/validation/validationIDs_extra.txt --out data/validation/validation_extra --recodeA --allow-extra-chr --debug")










#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#       Prepare tSNPs for validation           #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Extract longer flanks for df of transcriptome SNPs

trans_val <- select(Tsnps, SNPid) %>% 
  separate(SNPid, c("Contig", "Position"), sep = ":") %>%
  separate(Position, c("start", "end"), remove = F) %>%
  mutate(start = as.numeric(start),
         end = as.numeric(end)) %>%
  mutate(start = start - 25,
         end = end + 25) %>%
  unite(SNPid, Contig, Position, sep = ":") %>%
  left_join(Tsnps) %>%
  filter(start > 0) %>%
  separate(SNPid, c("Contig", "Position"), sep = ":")

write.table(trans_val[c(1,3,4)], "data/transcriptome/trans_val.bed", 
            quote = F, row.names = F, col.names = F, sep = "\t")


# Get longer flanking sequence from transcriptome.
# run bed tools

bedCmd <- paste("bedtools getfasta -fi ~/assemblies/ArcGaz_transcriptome/joined_transcriptome.fasta -bed data/transcriptome/trans_val.bed -fo data/transcriptome/trans_val.fasta")
system(bedCmd)


# Some flanks are longer than transcript lengths
# Get list of flanks

trans_val_seqs <- read.fasta("data/transcriptome/trans_val.fasta", as.string = T, forceDNAtolower = F) %>%
  lapply(function(x) x[[1]])

library(reshape2)
trans_val_seqs <- melt(trans_val_seqs)

trans_val <- trans_val %>%
  unite(newflank, start, end, sep = "-") %>%
  unite(L1, Contig, newflank, sep = ":", remove = F) %>%
  left_join(trans_val_seqs)

# enclose SNP

trans_val$SNPenclose <- 61

trans_refs <- paste0("[",str_sub(trans_val$value, 61, 61),"]")
trans_val$value <- as.character(trans_val$value)


for (i in 1:length(trans_val$value)){
  substr(trans_val$value[i], trans_val$SNPenclose[i], trans_val$SNPenclose[i]) <- "]"
}

for(i in 1:length(trans_val$value))
  trans_val$value[i]<-gsub("]", trans_refs[i], trans_val$value[i])


# Randomly select 5 good 5 bad

head(trans_val)






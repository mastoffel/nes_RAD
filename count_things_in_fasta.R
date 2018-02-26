# count stuff in fasta file

library(Biostrings)
nes_ref <- readDNAStringSet("../angsd_analysis/nes_catalogue.fa")
head(nes_ref)

nes_ranges <- nes_ref@ranges
mean_length <- mean(nes_ranges@width) #489.9654
median_length <- median(nes_ranges@width) # 571
quantile(nes_ranges@width, probs = seq(0, 1, 0.25), na.rm = FALSE,
         names = TRUE)

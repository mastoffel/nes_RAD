# stacks parameter ranges 80% rule

library(readr)
library(stringr)
library(dplyr)
run <- "M1"
extract_df_from_poplog <- function(run){
  # read file
  pop_file <- readLines(paste0("/external/data2/martin/NES/stacks_runs/stacks_", run, "/populations/populations.log"))
  # extract part
  start_extract <- which(pop_file == "# Distribution of the number of SNPs per catalog locus after filtering.")
  # formatting
  pop_file <- pop_file[start_extract:(length(pop_file)-2)][-c(1,2)]
  pop_file_split <- str_split(pop_file, "\t")
  df <- do.call(rbind, lapply(pop_file_split, function(x) data.frame("number_snps" = x[[1]], "number_loci" = x[[2]])))
  df <- data.frame(apply(df, 2, function(x) as.numeric(as.character(x))))
  df <- data.frame("run" = run, df)
}

all_runs <- c("M1", "M2", "M3", "M4", "M5", "M6", "M7", "M8")

# gather data from all runs
all_runs_df <- do.call(rbind, lapply(all_runs, extract_df_from_poplog))

# plot all kept sites

all_runs_sum <- all_runs_df %>% dplyr::group_by(run) %>% 
  filter(number_snps != 0) %>% 
  dplyr::summarise(number_loci_sum = sum(number_loci))

ggplot(all_runs_sum, aes(run, number_loci_sum)) + geom_point()










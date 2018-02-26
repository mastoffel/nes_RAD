# collect vars from bootstrap replicates

library(readr)
library(purrr)
library(stringr)
library(parallel)
options(scipen=999)
# how many bootstraps
nboot <- 100
num_sims <- 10
setwd("/home/martin/nes/nes_RAD")


# folders to collect from
# bootstrap folders
folders_boot <- list.files("fastsimcoal_analyses/analyse_bootstrap_sfs/", pattern = "boot")
# simulations per bootstrap replicate folders
# folders_sim <- list.files("fastsimcoal_analyses/analyse_bootstrap_sfs/boot_1/", pattern = "fscrun[0-9]")
folders_sim <- paste0("fscrun", 1:num_sims)
# combine
folders_df <- data.frame("boot" = rep(folders_boot, each = length(folders_sim)), "sims" = rep(folders_sim, times = length(folders_boot)))

# collect likelihoods 
collect_vars <- function(folder_boot, folder_sim){
  fsc_pars <- read_delim(paste0("fastsimcoal_analyses/analyse_bootstrap_sfs/", folder_boot, "/", folder_sim, "/nes/nes.bestlhoods"), delim = "\t")
}
# collect everything (usually just take maximum likelihood per bootstrap)
all_vars <- map2(folders_df$boot, folders_df$sims, possibly(collect_vars, NA_real_)) 
all_vars_working <- do.call(rbind, all_vars)

quantile(all_vars_working$NCUR)

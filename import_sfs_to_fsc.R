# imports an sfs from angsd to the fastsimcoal_files folder
sfs <- scan("../angsd_analysis/SFS37/nes37.sfs")
# create names
sfs_names <- sapply(1:length(sfs), function(x) paste0("d0_", x))

sink("../nes_RAD/fastsimcoal_files/nes_MAFpop0.obs")
cat("1 observations")
cat("\n")
cat(sfs_names)
cat("\n")
cat(sfs)
sink()



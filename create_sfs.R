# create obs sfs for fastsimcoal
sfs <- scan("../../angsd_analysis/SFS25/nes25.sfs")

# create names
sfs_names <- sapply(1:length(sfs), function(x) paste0("d0_", x))

sink("/home/martin/nes/nes_RAD/fsc_run/nes_MAFpop0.obs")
cat("1 observations")
cat("\n")
cat(sfs_names)
cat("\n")
cat(sfs)
sink()

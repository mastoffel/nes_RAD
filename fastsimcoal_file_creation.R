# test run with fastsimcoal

# read site frequency spectrum 
sfs <- scan("smallFolded.sfs")[97:1]

# format sfs as obs file for fastsimcoal2 ----
sfs_df <- as.data.frame(matrix(sfs, ncol = length(sfs)))
names(sfs_df) <- paste0("d0_", c(1:97))
# non-scientific notation
options(scipen = 999)
# create obs file
write.table("1 observation", file = "nes_MAFpop0.obs", quote = FALSE, 
            col.names = FALSE, row.names = FALSE)
# add sfs
write.table(sfs_df, file = "nes_MAFpop0.obs", append = TRUE,quote = FALSE, 
            col.names = TRUE, row.names = FALSE)

library(stringr)

sink(file = "nes.tpl")
cat(c("//Number of population samples (demes)",
      "1",
      "//Population effective sizes (number of genes)", 
      "NCUR", 
      "//Sample sizes", 
      "96",
      "//Growth rates : negative growth implies population expansion", 
      "0", 
      "//Number of migration matrices : 0 implies no migration between demes",  
      "0", 
      "//historical event: time, source, sink, migrants, new size, new growth rate, migr. matrix", 
      "2 historical event",  
      "TBOT 0 0 0 RESBOT 0 0",    
      "TENDBOT 0 0 0 RESENDBOT 0 0",
      "//Number of independent loci [chromosome]",
      "1 0",
      "//Per chromosome: Number of linkage blocks", 
      "1", 
      "//per Block: data type, num loci, rec. rate and mut rate + optional parameters", "\n",
      "FREQ 1 0 2.5e-8"),
    sep = "\n"
)
sink()


# est file

sink(file = "nes.est")
cat(c("// Priors and rules file",
      "// *********************",
      "[PARAMETERS]", 
      "1 NCUR unif 10 100000 output", 
      "1 NANC unif 10 100000 output", 
      "1 NBOT unif 10 100000 output",
      "1 TBOT unif 10 1000 output", 
      "[RULES]", 
      "[COMPLEX PARAMETERS]",  
      "0 RESBOT = NBOT/NCUR hide", 
      "0 RESENDBOT = NANC/NBOT hide",
      "1 TENDBOT = TBOT+100 hide"),
    sep = "\n"
)
sink()


      
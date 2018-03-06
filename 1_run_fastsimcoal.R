# run several iterations of fastsimcoal

library(parallel)

setwd("/home/martin/nes/nes_RAD/fastsimcoal_analyses/fsc_run")
if (length(list.files("./")) != 0) { 
    system("rm -r *")
  }
# system("~/bin/fsc26")
#folders <- list.files("fsc_run/", pattern = "fsc_run[0-9]")

run_fsc <- function(run_num){
  # create directory for fsc run
  system(paste0("mkdir fscrun", run_num))
  
  # paste all relevant files into directory
  system(paste0("cp ../fastsimcoal_files/nes* fscrun", run_num))
  
  # change to directory
  setwd(paste0("/home/martin/nes/nes_RAD/fastsimcoal_analyses/fsc_run/fscrun", run_num))
  
  #if (run_num < 51) {
    # write fsc_run file 
    writeLines(c(paste0("/home/martin/nes/nes_RAD/fastsimcoal_analyses/fsc_run/fscrun", run_num),
                 "-t nes.tpl -n 100000 -m -e nes.est -M -L 40 -q -w 0.01 --foldedSFS -x -C 10 --nosingleton"),
               paste0("/home/martin/nes/nes_RAD/fastsimcoal_analyses/fsc_run/fscrun", run_num, "/fsc_run.txt"))
    
  #}

  # run fsc
  system("~/bin/fsc26")
  
  # change back
  setwd(paste0("/home/martin/nes/nes_RAD/fastsimcoal_analyses/fsc_run/"))
}

# run all
cl <- makeCluster(getOption("cl.cores", 30))
parLapply(cl, 1:50, run_fsc)
stopCluster(cl)

setwd("/home/martin/nes/nes_RAD/")

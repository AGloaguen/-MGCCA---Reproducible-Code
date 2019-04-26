rm(list = ls())

######################
#####  LIBRARIES #####
######################

library(pracma)
library(R.matlab)

script_location = pwd()
katri_rao       = TRUE
output_dir      = paste0(script_location, paste0("/simu_", gsub(pattern = "  | |:", replacement = "_", x = date())))

dir.create(file.path(output_dir), showWarnings = FALSE)

#######################
#####  PARAMETERS #####
#######################
n_simu             = 10
SNR                = c(0.1, 0.5, 1)
n_random_start     = 2
nstart_next_comps  = 2
tol                = 10^(-8)
init_nfac          = 5
SNR_to_plot        = "SNR_0.5"
methods_to_compare = c("PARAFAC_nbcomp_5", "CMTF_nbcomp_5", "RGCCA_", "MGCCA_dim_2")



#Parameters MGCCA
ncomp_MGCCA          = c(2, 2)
ncomp_MGCCA          = paste(ncomp_MGCCA, collapse = ",")
deflation_mode_MGCCA = c("NULL", "dim_2")
#Parameters RGCCA
ncomp_RGCCA          = c(2, 2)
ncomp_RGCCA          = paste(ncomp_RGCCA, collapse = ",")
#Parameters CMTF
ncomp_CMTF           = 5
mc.cores_CMTF        = 2
#Parameters PARAFAC
ncomp_PARAFAC        = 5
mc.cores_PARAFAC     = 1


##########################
#####  GENERATE DATA #####
##########################
fileConn<-file(paste0(output_dir, "/", "commands.txt"), "w")
writeLines(paste0("Rscript ", script_location,
                  "/generate_simulations.R ", SNR, " ", init_nfac, " ",
                  output_dir, " ", katri_rao, " ", n_simu, " "), fileConn)

##################
#####  MGCCA #####
##################
for (deflation_mode in deflation_mode_MGCCA){
  writeLines(paste0(paste("Rscript", paste0(script_location, "/mapper_mgcca.R"),
                     SNR, n_random_start, ncomp_MGCCA, deflation_mode, output_dir, nstart_next_comps, 
                     sep = " "), " "), fileConn)
  writeLines(paste0(paste("Rscript", paste0(script_location, "/reducer_mgcca.R"),
                     deflation_mode, output_dir, sep = " "), " "), fileConn)
}

##################
#####  RGCCA #####
##################
writeLines(paste0(paste("Rscript", paste0(script_location, "/mapper_rgcca.R"),
                   SNR, n_random_start, ncomp_RGCCA,  output_dir, nstart_next_comps,sep = " "), " "), 
      fileConn)
writeLines(paste0(paste("Rscript", paste0(script_location, "/reducer_rgcca.R"),
                   output_dir, sep = " "), " "), fileConn)

##################
#####  CMTF  #####
##################
writeMat(paste(output_dir, "/mapper_cmtf_parameter_file.mat", sep = ""),
         nstart   = n_random_start*nstart_next_comps, 
         path     = output_dir,
         nb_cores = mc.cores_CMTF,
         tol      = tol)

for (snr in SNR){
  for (R in ncomp_CMTF){
  cmd = paste("\"clear all; close all; load(\'",
              paste(output_dir, "/mapper_cmtf_parameter_file.mat", sep = ""),"\');",
              "SNR=", snr, ";",
              "R=", R, ";",
              "run(\'", paste0(script_location, "/mapper_cmtf.m\');exit;\""), " ",
              sep = "")
  writeLines(paste("matlab", "-nodisplay", "-nosplash", "-nodesktop", "-r", cmd), fileConn)
  writeLines(paste0(paste("Rscript", paste0(script_location, "/reducer_cmtf.R"),
                   snr,  R, output_dir, sep = " "), " "), fileConn)
  }
}

####################
#####  PARAFAC #####
####################
for (ncomp in ncomp_PARAFAC){
  writeLines(paste0(paste("Rscript", paste0(script_location, "/mapper_parafac.R"),
                     SNR, ncomp, n_random_start*nstart_next_comps, mc.cores_PARAFAC,  output_dir, sep = " "), " "), fileConn)
  writeLines(paste0(paste("Rscript", paste0(script_location, "/reducer_parafac.R"),
                     ncomp,  output_dir, sep = " "), " "), fileConn)
}

##########################
#####  FINAL REDUCER #####
##########################
methods_to_compare = paste(methods_to_compare, collapse = ",")
writeLines(paste0("Rscript ", script_location,
           "/final_reducer.R ", output_dir, " ", n_simu, " ", SNR_to_plot, " ", methods_to_compare, " "), fileConn)

close(fileConn)
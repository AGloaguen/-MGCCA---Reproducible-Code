rm(list = ls())

args = commandArgs(trailingOnly = TRUE)

library(stringr)

if(length(args) != 5){
  stop("SNR/Number of Components/Number of random starts/Number of cores/current path have to be specified")
}else{
  SNR          = as.numeric(args[1])
  nfac         = as.numeric(args[2])
  nstart       = as.numeric(args[3])
  mc.cores     = as.numeric(args[4])
  current_path = args[5]
}

library(multiway)
library(pbmcapply)

norm2 <-function(vec){
  a <- sqrt(sum(vec^2))
  if(a==0) a <- .05
  return(a)
}

Simu_Folds        = list.files(path = paste0(current_path, "/SNR_", SNR, "/"), pattern = "Simu_Fold_", full.names = T)
Truth             = list.files(path = paste0(current_path, "/SNR_", SNR, "/"), pattern = "W_", full.names = T)
Already_computed  = list.files(path = paste0(current_path, "/SNR_", SNR, "/", "Parafac/")
                               , pattern = "parafac.results_Fold_", full.names = F)

Already_computed = sapply(Already_computed, function(x) as.numeric(str_match(x, "parafac.results_Fold_(.*?).Rdata")[, 2]))
Already_computed = unname(Already_computed)

idx_Simu_Folds = sapply(Simu_Folds, function(x) as.numeric(str_match(x, "Simu_Fold_(.*?).Rdata")[, 2]))
idx_Simu_Folds = unname(idx_Simu_Folds)

Simu_Folds_to_rm = match(Already_computed, idx_Simu_Folds)

if (length(Simu_Folds_to_rm) != 0){
  Simu_Folds = Simu_Folds[-Simu_Folds_to_rm]
}

idx_Simu_Folds = sapply(Simu_Folds, function(x) as.numeric(str_match(x, "Simu_Fold_(.*?).Rdata")[, 2]))
idx_Simu_Folds = unname(idx_Simu_Folds)

for (i in 1:length(Truth)){
  load(file = Truth[i])
}

parafac.results = list()

compute_parafac = function(i){
  load(file = Simu_Folds[i])
  
  #Parafac sur X_1
  PF_1      = parafac(X = A[[1]], nfac = nfac, nstart = 1, 
                      Bstart = INIT$init_B1[, 1:nfac], Cstart = INIT$init_C1[, 1:nfac], ctol = 1e-8)
  PF_1_rand = parafac(X = A[[1]], nfac = nfac, nstart = nstart, ctol = 1e-8)

  if (PF_1_rand$SSE < PF_1$SSE){
    PF_1        = PF_1_rand
    init_bloc_1 = "random"
  }else{
    init_bloc_1 = "svd"
  }
  parafac.results[[1]] = PF_1
  
  #Parafac sur X_2
  PF_2      = parafac(X = A[[2]], nfac = nfac, nstart = 1, 
                      Bstart = INIT$init_B2[, 1:nfac], Cstart = INIT$init_C2[, 1:nfac], ctol = 1e-8)
  PF_2_rand = parafac(X = A[[2]], nfac = nfac, nstart = nstart, ctol = 1e-8)
  
  if (PF_2_rand$SSE < PF_2$SSE){
    PF_2        = PF_2_rand
    init_bloc_2 = "random"
  }else{
    init_bloc_2 = "svd"
  }

  parafac.results[[2]] = PF_2
  parafac.results[[3]] = init_bloc_1
  parafac.results[[4]] = init_bloc_2

  names(parafac.results) = c("raw_results_1", "raw_results_2", "init_bloc_1", "init_bloc_2")

  save(file = paste0(current_path, "/SNR_", SNR, "/", "Parafac/parafac.results_nbcomp_", nfac, "_Fold_", idx_Simu_Folds[i], ".Rdata"), parafac.results)
}

res = mclapply(1:length(Simu_Folds), compute_parafac, mc.cores = mc.cores)


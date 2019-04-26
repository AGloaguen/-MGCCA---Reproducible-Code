rm(list = ls())

args = commandArgs(trailingOnly = TRUE)

if(length(args) != 6){
  stop("SNR/Number of random starts/Number of comp/deflation mode/current path/Number of random starts for next components
       have to be specified")
}else{
  SNR                = as.numeric(args[1])
  nstart             = as.numeric(args[2])
  ncomp              = as.numeric(strsplit(args[3], ",")[[1]])
  deflation_mode_str = args[4]
  current_path       = args[5]
  nstart_next_comps  = args[6]
}

if (deflation_mode_str == "NULL"){
  deflation_mode = NULL
}else{
  deflation_mode = deflation_mode_str
}

if (require("devtools")){
  load_all(paste0(current_path, "/../../RGCCA/."))
}else{
  file.sources = list.files(path = paste0(current_path, "/../../RGCCA/R", sep = ""), pattern = "*.R", full.names = T)
  sapply(file.sources, source, .GlobalEnv)
}


tau              = c(1, 1)
bias             = T
scale            = F
C                = matrix(1, nrow = 2, ncol = 2) - diag(2)
scheme           = "factorial"
init             = "svd"
center           = F
nstart_at_comp_2 = T
M_regularisation = rep("kronecker_Identity_RGCCA", 2)

mgcca.results = list()

Simu_Folds = list.files(path = paste0(current_path, "/SNR_", SNR, "/"), pattern = "Simu_Fold_", full.names = T)
Truth      = list.files(path = paste0(current_path, "/SNR_", SNR, "/"), pattern = "W_", full.names = T)

for (i in 1:length(Truth)){
  load(file = Truth[i])
}


for (i in 1:length(Simu_Folds)){
  load(file = Simu_Folds[i])
  
  init         = list()
  init[[1]]    = mapply(function(x, y) x %x% y, as.data.frame(INIT$init_C1), 
                        as.data.frame(INIT$init_B1))
  init[[2]]    = mapply(function(x, y) x %x% y, as.data.frame(INIT$init_C2), 
                        as.data.frame(INIT$init_B2))
  mgcca.result = mgcca_array(A = A, C = C, tau = tau, ncomp = ncomp, scheme = scheme, 
                             center = center, scale = scale, bias = bias, init = init,
                             verbose = F, deflation_mode = deflation_mode, nstart = 1,
                             M_regularisation = M_regularisation)
  
  mgcca.results[[1]]      = list()
  mgcca.results[[1]][[1]] = mgcca.result
  
  for (j in 1:nstart){
    mgcca.result = mgcca_array(A = A, C = C, tau = tau, ncomp = ncomp, scheme = scheme, 
                               center = center, scale = scale, bias = bias, init = "random", 
                               verbose = F, deflation_mode = deflation_mode, nstart = nstart_next_comps, 
                               M_regularisation = M_regularisation, nstart_at_comp_2 = nstart_at_comp_2)
    
    mgcca.results[[1]][[j+1]] = mgcca.result
  }
  
  best_start         = which.max(unlist(lapply(mgcca.results[[1]], function(x){
    if(is.list(x$crit)){
      return(Reduce("sum", lapply(x$crit, max)))
    }else{
      return(max(x$crit))
    }
  } )))
  mgcca.result       = mgcca.results[[1]][[best_start]]
  mgcca.results[[1]] = mgcca.result
  
  if(best_start == 1){
    mgcca.results[[2]] = "svd"
  }else{
    mgcca.results[[2]] = "random"
  }
  
  names(mgcca.results) = c("raw_results", "best_start")
  
  save(file = paste0(current_path, "/SNR_", SNR, "/", "MGCCA/mgcca.results_", deflation_mode_str, "_Fold_", i, ".Rdata"), mgcca.results)
}


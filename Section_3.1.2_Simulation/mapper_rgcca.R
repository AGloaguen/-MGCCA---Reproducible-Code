rm(list = ls())

args = commandArgs(trailingOnly = TRUE)

if(length(args) != 5){
  stop("SNR/Number of random starts/Number of components/current path/Number of random starts for next components
       have to be specified")
}else{
  SNR               = as.numeric(args[1])
  nstart            = as.numeric(args[2])
  ncomp             = as.numeric(strsplit(args[3], ",")[[1]])
  current_path      = args[4]
  nstart_next_comps = args[5]
  
}

if (require("devtools")){
  load_all(paste0(current_path, "/../../RGCCA/."))
}else{
  file.sources = list.files(path = paste0(current_path, "/../../RGCCA/R", sep = ""), pattern = "*.R", full.names = T)
  sapply(file.sources, source, .GlobalEnv)
}


tau             = c(1, 1)
bias            = T
scale           = F
C               = matrix(1, nrow = 2, ncol = 2) - diag(2)
scheme          = "factorial"
scale_size_bloc = F

rgcca.results = list()

Simu_Folds = list.files(path = paste0(current_path, "/SNR_", SNR, "/"), pattern = "Simu_Fold_", full.names = T)
Truth      = list.files(path = paste0(current_path, "/SNR_", SNR, "/"), pattern = "W_", full.names = T)

for (i in 1:length(Truth)){
  load(file = Truth[i])
}


for (i in 1:length(Simu_Folds)){
  load(file = Simu_Folds[i])
  
  A_2          = lapply(A, function(x) t(apply(x, 1, c)))
  init         = list()
  init[[1]]    = init[[2]] = INIT[[1]]
  rgcca.result = rgcca_nstart(A = A_2, C = C, tau = tau, ncomp = ncomp,
                              scheme = scheme, scale = scale, bias = bias, 
                              init = init, verbose = F, nstart_next_comps = 1)
  
  rgcca.results[[1]]      = list()
  rgcca.results[[1]][[1]] = rgcca.result
  
  for (j in 1:nstart){
    rgcca.result = rgcca_nstart(A = A_2, C = C, tau = tau, ncomp = ncomp, 
                                scheme = scheme, scale = scale, bias = bias, 
                                init = "random", verbose = F, nstart_next_comps = nstart_next_comps)
    
    rgcca.results[[1]][[j+1]] = rgcca.result
  }
  print(i)
  best_start         = which.max(unlist(lapply(rgcca.results[[1]], function(x){
    if(is.list(x$crit)){
      return(Reduce("sum", lapply(x$crit, max)))
    }else{
      return(max(x$crit))
    }
  } )))
  rgcca.result       = rgcca.results[[1]][[best_start]]
  rgcca.results[[1]] = rgcca.result
  
  if(best_start == 1){
    rgcca.results[[2]] = "svd"
  }else{
    rgcca.results[[2]] = "random"
  }
  
  names(rgcca.results) = c("raw_results", "best_start")
  
  save(file = paste0(current_path, "/SNR_", SNR, "/", "RGCCA/rgcca.results_Fold_", i, ".Rdata"), rgcca.results)
}


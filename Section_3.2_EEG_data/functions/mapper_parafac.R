args = commandArgs(trailingOnly = TRUE)

library(multiway)
library(pbmcapply)

#Get parametrs from command line
bootstrap_sample_parameters_file = args[1]
output_file                      = args[2]
bootstrap_sample_number          = strtoi(args[3])
pathtoscript                     = args[4]

if (require("devtools")){
  load_all(paste0(pathtoscript, "../../RGCCA/."))
}else{
  file.sources = list.files(path = paste0(pathtoscript, "../../RGCCA/R", sep = ""), pattern = "*.R", full.names = T)
  sapply(file.sources, source, .GlobalEnv)
}

#Load parameters for bootstrap
load(bootstrap_sample_parameters_file)

#Load data file
load(data_path)
N = nrow(X[[1]])
L = length(X) - 1 #parafac works with only one bloc

full_data_run = full_data_run$parafac

check_for_sign = function(x, y){
  if (length(x) == 1){
    return(list(drop(sign(x*y))))
  }else{
    return(list(drop(sign(diag(cor(x, y))))))
  }
}

redefine_sign = function(x, y){
  if(length(y) == 1){
    return(list(x*y))
  }else{
    return(list(x %*% diag(y)))
  }
}

par_parafac = function(full_data_run, X, idx, nfac, nstart, init, tol){
  #Create boostrap sample
  X.boot        = scale3_array(X[[1]][idx, , ], center = T, scale = F, bias = F)
  
  if (init == "svd"){
    Bstart = apply(X.boot, 2, c)
    Bstart = svd(Bstart,nu=0,nv=nfac)$v
    Cstart = apply(X.boot, 3, c)
    Cstart = svd(Cstart,nu=0,nv=nfac)$v
  }else{
    Bstart = Cstart = NULL
  }
  
  #Compute RGCCA and save results
  result.parafac = parafac(X      = X.boot, 
                           nfac   = nfac, 
                           nstart = nstart, 
                           Bstart = Bstart, 
                           Cstart = Cstart, 
                           ctol   = tol)
  
  Y_signes     = mapply(function(x, y) check_for_sign(x, y), list(result.parafac$A), list(full_data_run$A))
  B_signes     = mapply(function(x, y) check_for_sign(x, y), list(result.parafac$B), list(full_data_run$B))
  C_signes     = mapply(function(x, y) check_for_sign(x, y), list(result.parafac$C), list(full_data_run$C))
  
  Ystar = mapply(function(x, y) redefine_sign(x, y), list(result.parafac$A), Y_signes)
  Bstar = mapply(function(x, y) redefine_sign(x, y), list(result.parafac$B), B_signes)
  Cstar = mapply(function(x, y) redefine_sign(x, y), list(result.parafac$C), C_signes)
  
  return(list(Ystar = Ystar, Bstar = Bstar, Cstar = Cstar))
}

nb_boot = dim(idx_boot)[1]
res     = mclapply(1:nb_boot, function(i) par_parafac(full_data_run = full_data_run, 
                                                    X             = X, 
                                                    idx           = idx_boot[i, ], 
                                                    nfac          = nfac, 
                                                    nstart        = nstart,
                                                    init          = init,
                                                    tol           = tol), 
                   mc.cores = nb_cores)

Ystar = lapply(res[[1]]$Ystar, function(x) array(0, dim = c(dim(x), nb_boot)))
Bstar = lapply(res[[1]]$Bstar, function(x) array(0, dim = c(dim(x), nb_boot)))
Cstar = lapply(res[[1]]$Cstar, function(x) array(0, dim = c(dim(x), nb_boot)))
for(i in 1:nb_boot){
  for(j in 1:L){
    Ystar[[j]][, , i] = res[[i]]$Ystar[[j]]
    if (j %in% B_3D){
      Bstar[[j]][, , i] = res[[i]]$Bstar[[j]]
      Cstar[[j]][, , i] = res[[i]]$Cstar[[j]]
    }
  }
}

save(Ystar = Ystar, Bstar = Bstar, Cstar = Cstar, 
     bootstrap_sample_number = bootstrap_sample_number, file = output_file)


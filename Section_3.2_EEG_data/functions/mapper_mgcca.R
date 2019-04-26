
args = commandArgs(trailingOnly = TRUE)

library(pbmcapply)

#Get parametrs from command line
bootstrap_sample_parameters_file = args[1]
output_file                      = args[2]
bootstrap_sample_number          = strtoi(args[3])
pathtoscript                     = args[4]
design                           = args[5]
scheme                           = args[6]

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
L = length(X)

#Generate design matrix C
if (design == "hierarchical") { 
  C       = matrix(0, ncol = L, nrow = L) 
  C[L, ]  = C[, L] = 1
  C[L, L] = 0
} else if (design == "complete") { 
  C = matrix(1, ncol = L, nrow = L) - diag(L)
}

full_data_run = full_data_run$mgcca

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

par_mgcca = function(full_data_run, X, idx, B_3D, tau, ncomp, C, sheme, init, L, tol){
  
  #Create boostrap sample
  X.boot        = rep(list(0), L)
  X.boot[B_3D]  = lapply(X[B_3D], function(x) scale3_array(x[idx, , ], center = T, scale = F, bias = F))
  X.boot[-B_3D] = lapply(X[-B_3D], function(x) scale(as.matrix(as.matrix(x)[idx, ]), center = T, scale = F))
  names(X.boot) = names(X)
  
  #Compute RGCCA and save results
  result.mgcca = mgcca_array(A       = X.boot, 
                             tau     = tau, 
                             ncomp   = ncomp, 
                             C       = C, 
                             scheme  = scheme,     
                             verbose = F, 
                             init    = init, 
                             scale   = F, 
                             center  = F, 
                             bias    = F,
                             tol     = tol,
                             M_regularisation = M_regularisation)
  
  Y_signes     = mapply(function(x, y) check_for_sign(x, y), result.mgcca$Y, full_data_run$Y)
  Astar_signes = mapply(function(x, y) check_for_sign(x, y), result.mgcca$astar, full_data_run$astar)
  B_signes     = mapply(function(x, y) check_for_sign(x, y), result.mgcca$b, full_data_run$b)
  C_signes     = mapply(function(x, y) check_for_sign(x, y), result.mgcca$c, full_data_run$c)
  
  Ystar = mapply(function(x, y) redefine_sign(x, y), result.mgcca$Y, Y_signes)
  Astar = mapply(function(x, y) redefine_sign(x, y), result.mgcca$astar, Astar_signes)
  Bstar = mapply(function(x, y) redefine_sign(x, y), result.mgcca$b, B_signes)
  Cstar = mapply(function(x, y) redefine_sign(x, y), result.mgcca$c, C_signes)

  return(list(Ystar = Ystar, Astar = Astar, Bstar = Bstar, Cstar = Cstar))
}

nb_boot = dim(idx_boot)[1]
res     = mclapply(1:nb_boot, function(i) par_mgcca(full_data_run = full_data_run, 
                                                    X             = X, 
                                                    idx           = idx_boot[i, ], 
                                                    B_3D          = B_3D, 
                                                    tau           = tau, 
                                                    ncomp         = ncomp, 
                                                    C             = C, 
                                                    sheme         = sheme, 
                                                    init          = init, 
                                                    L             = L, 
                                                    tol           = tol), 
                   mc.cores = nb_cores)

Ystar = lapply(res[[1]]$Ystar, function(x) array(0, dim = c(dim(x), nb_boot)))
Astar = lapply(res[[1]]$Astar, function(x) array(0, dim = c(dim(x), nb_boot)))
Bstar = lapply(res[[1]]$Bstar, function(x) array(0, dim = c(dim(x), nb_boot)))
Cstar = lapply(res[[1]]$Cstar, function(x) array(0, dim = c(dim(x), nb_boot)))
for(i in 1:nb_boot){
  for(j in 1:L){
    Ystar[[j]][, , i] = res[[i]]$Ystar[[j]]
    Astar[[j]][, , i] = res[[i]]$Astar[[j]]
    if (j %in% B_3D){
      Bstar[[j]][, , i] = res[[i]]$Bstar[[j]]
      Cstar[[j]][, , i] = res[[i]]$Cstar[[j]]
    }
  }
}

save(Ystar = Ystar, Astar = Astar, Bstar = Bstar, Cstar = Cstar, 
     bootstrap_sample_number = bootstrap_sample_number, file = output_file)


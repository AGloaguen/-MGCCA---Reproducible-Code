full_data_run = list()

if (require("devtools")){
  load_all(paste0(pathtoscript, "../../RGCCA/."))
}else{
  file.sources = list.files(path = paste0(current_path, "../../RGCCA/R", sep = ""), pattern = "*.R", full.names = T)
  sapply(file.sources, source, .GlobalEnv)
}

X[B_3D]  = lapply(X[B_3D], function(x) scale3_array(x, center = T, scale = F, bias = F))
X[-B_3D] = lapply(X[-B_3D], function(x) scale(as.matrix(as.matrix(x)), center = T, scale = F))

j = 1
full_data_run_names = c()
for(i in 1:length(mappers)){
  mapper = mappers[i]
  if (mapper == "mgcca"){    
    design = method.cfg$option.mapper[[mapper]]$design
    scheme = method.cfg$option.mapper[[mapper]]$scheme
    
    #Generate design matrix C
    if (design == "hierarchical") { 
      C       = matrix(0, ncol = L, nrow = L) 
      C[L, ]  = C[, L] = 1
      C[L, L] = 0
    } else if (design == "complete") { 
      C = matrix(1, ncol = L, nrow = L) - diag(L)
    }
    
    #Compute RGCCA and save results
    reswr = mgcca_array(A       = X, 
                        tau     = tau, 
                        ncomp   = ncomp, 
                        C       = C, 
                        scheme  = scheme,     
                        verbose = F, 
                        init    = init, 
                        scale   = F, 
                        center  = F, 
                        bias    = F,
                        M_regularisation = M_regularisation)
    full_data_run[[j]]  = reswr
    j                   = j + 1
    full_data_run_names = c(full_data_run_names, mapper)
  }
  
  if (mapper == "parafac"){
    library(multiway)
    if (init == "svd"){
      Bstart = apply(X[[1]], 2, c)
      Bstart = svd(Bstart,nu=0,nv=nfac)$v
      Cstart = apply(X[[1]], 3, c)
      Cstart = svd(Cstart,nu=0,nv=nfac)$v
    }else{
      Bstart = Cstart = NULL
    }
    reswr = parafac(X      = X[[1]], 
                    nfac   = nfac, 
                    nstart = nstart, 
                    Bstart = Bstart, 
                    Cstart = Cstart, 
                    ctol   = tol)
    full_data_run[[j]]  = reswr
    j                   = j + 1
    full_data_run_names = c(full_data_run_names, mapper)
  }
  
  if (mapper == "cmtf"){
    full_data_parameter_file = paste(pathtowritefile, method.cfg$option.mapper[[mapper]]$full_data_parameter_file, sep = "")
    full_data_file           = paste(pathtowritefile, method.cfg$option.mapper[[mapper]]$full_data_file, sep = "")
    writeMat(full_data_parameter_file,
             nfac            = nfac,
             parfor_bool     = F,
             init            = init, 
             tol             = tol,
             data_path       = paste(file_path_sans_ext(data_path), "mat", sep = "."))

    cmd = paste("\"clear all; close all; parameters_file = ",
                paste0("\'", full_data_parameter_file, "\'"), "; output_file = ",
                paste0("\'", full_data_file, "\'"), "; file_number = 10;",
                paste0("run(\'", pathtoscript, method.cfg$model.mapper[[mapper]], "\');exit;\""), " ",
                sep = "")
  system(paste("matlab", "-nodisplay", "-nosplash", "-nodesktop", "-r", cmd))
  }
}

names(full_data_run) = full_data_run_names




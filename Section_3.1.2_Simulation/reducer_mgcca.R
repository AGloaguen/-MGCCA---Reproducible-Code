rm(list = ls())

args = commandArgs(trailingOnly = TRUE)

library(abind)

if(length(args) != 2){
  stop("Deflation mode/current path have to be specified")
}else{
  deflation_mode = args[1]
  current_path   = args[2]
}

directories = list.dirs(current_path, recursive = F, full.names = F)
SNR_dirs    = directories[ grepl("SNR_", directories) ]

idx_SNR = sapply(SNR_dirs, function(x) as.numeric(gsub("SNR_","",x)))
idx_SNR = unname(idx_SNR)

for (j in 1:length(idx_SNR)){
  
  current_SNR_path = paste0(current_path, "/SNR_", idx_SNR[j])
  
  load(file = paste0(current_SNR_path, "/W_B_1.Rdata"))
  load(file = paste0(current_SNR_path, "/W_B_2.Rdata"))
  load(file = paste0(current_SNR_path, "/W_C_1.Rdata"))
  load(file = paste0(current_SNR_path, "/W_C_2.Rdata"))
  
  KRON_W_1 = sapply(1:dim(W_C_1)[2], function(x) W_C_1[, x] %x% W_B_1[, x])
  KRON_W_2 = sapply(1:dim(W_C_2)[2], function(x) W_C_2[, x] %x% W_B_2[, x])
  
  MGCCA_files = list.files(path = paste0(current_SNR_path, "/MGCCA"), pattern = paste0("mgcca.results_", deflation_mode, "_Fold_"), full.names = T)
  
  MGCCA_cor_1 = c()
  MGCCA_cor_2 = c()
  MGCCA       = list()
  n_simu      = length(MGCCA_files)
  
  for (i in 1:n_simu){
    load(MGCCA_files[i])
    
    #Gestion du signe des corrélations
    MGCCA_sign_1        = drop(sign(diag(cor(mgcca.results$raw_results$b[[1]], W_B_1))))
    MGCCA_sign_2        = drop(sign(diag(cor(mgcca.results$raw_results$b[[2]], W_B_2))))
    MGCCA_sign_3        = drop(sign(diag(cor(mgcca.results$raw_results$c[[1]], W_C_1))))
    MGCCA_sign_4        = drop(sign(diag(cor(mgcca.results$raw_results$c[[2]], W_C_2))))
    
    mgcca.result        = mgcca.results$raw_results
    mgcca.result$b[[1]] = mgcca.results$raw_results$b[[1]] %*% diag(MGCCA_sign_1)
    mgcca.result$b[[2]] = mgcca.results$raw_results$b[[2]] %*% diag(MGCCA_sign_2)
    mgcca.result$c[[1]] = mgcca.results$raw_results$c[[1]] %*% diag(MGCCA_sign_3)
    mgcca.result$c[[2]] = mgcca.results$raw_results$c[[2]] %*% diag(MGCCA_sign_4)
    
    mgcca.results[[3]] = mgcca.result
    
    #Evaluation de la correlation avec la vérité
    cor_1   = cor(c(mgcca.result$b[[1]], mgcca.result$b[[2]], mgcca.result$c[[1]], mgcca.result$c[[2]]),
                  c(W_B_1, W_B_2, W_C_1, W_C_2))
    cor_2 = cor(c(sapply(1:dim(mgcca.result$c[[1]])[2], function(x) mgcca.result$c[[1]][, x] %x% mgcca.result$b[[1]][, x]),
                  sapply(1:dim(mgcca.result$c[[2]])[2], function(x) mgcca.result$c[[2]][, x] %x% mgcca.result$b[[2]][, x])), 
                c(KRON_W_1, KRON_W_2))
    
    mgcca.results[[4]] = cor_1
    mgcca.results[[5]] = cor_2
    
    cor_1_by_comp = diag(cor(rbind(mgcca.result$b[[1]], mgcca.result$b[[2]], mgcca.result$c[[1]], mgcca.result$c[[2]]),
                             rbind(W_B_1, W_B_2, W_C_1, W_C_2)))
    cor_2_by_comp = diag(cor(rbind(sapply(1:dim(mgcca.result$c[[1]])[2], function(x) mgcca.result$c[[1]][, x] %x% mgcca.result$b[[1]][, x]),
                                   sapply(1:dim(mgcca.result$c[[2]])[2], function(x) mgcca.result$c[[2]][, x] %x% mgcca.result$b[[2]][, x])), 
                             rbind(KRON_W_1, KRON_W_2)))
    
    mgcca.results[[6]] = cor_1_by_comp
    mgcca.results[[7]] = cor_2_by_comp
    
    MGCCA_cor_1 = c(MGCCA_cor_1, cor_1)
    MGCCA_cor_2 = c(MGCCA_cor_2, cor_2)
    if(i == 1){
      dim_b_1             = dim(mgcca.result$b[[1]])
      dim_c_1             = dim(mgcca.result$c[[1]])
      dim_b_2             = dim(mgcca.result$b[[2]])
      dim_c_2             = dim(mgcca.result$c[[2]])
      MGCCA$b[[1]]        = array(data = mgcca.result$b[[1]], dim = c(dim_b_1, 1))
      MGCCA$c[[1]]        = array(data = mgcca.result$c[[1]], dim = c(dim_c_1, 1))
      MGCCA$b[[2]]        = array(data = mgcca.result$b[[2]], dim = c(dim_b_2, 1))
      MGCCA$c[[2]]        = array(data = mgcca.result$c[[2]], dim = c(dim_c_2, 1))
      crit_list           = mgcca.results$raw_results$crit
      crit                = ifelse(is.list(crit_list),
                                   Reduce("sum", lapply(crit_list, max)), 
                                   max(crit_list))
      MGCCA$crit          = crit
      MGCCA_cor_1_by_comp = cor_1_by_comp
      MGCCA_cor_2_by_comp = cor_2_by_comp
    }else{
      MGCCA$b[[1]]        = abind(MGCCA$b[[1]], array(data = mgcca.result$b[[1]], dim = c(dim_b_1, 1)), along = 3)
      MGCCA$c[[1]]        = abind(MGCCA$c[[1]], array(data = mgcca.result$c[[1]], dim = c(dim_c_1, 1)), along = 3)
      MGCCA$b[[2]]        = abind(MGCCA$b[[2]], array(data = mgcca.result$b[[2]], dim = c(dim_b_2, 1)), along = 3)
      MGCCA$c[[2]]        = abind(MGCCA$c[[2]], array(data = mgcca.result$c[[2]], dim = c(dim_c_2, 1)), along = 3)
      crit_list           = mgcca.results$raw_results$crit
      crit                = ifelse(is.list(crit_list), 
                                   Reduce("sum", lapply(crit_list, max)), 
                                   max(crit_list))
      MGCCA$crit          = c(MGCCA$crit, crit)
      MGCCA_cor_1_by_comp = rbind(MGCCA_cor_1_by_comp, cor_1_by_comp)
      MGCCA_cor_2_by_comp = rbind(MGCCA_cor_2_by_comp, cor_2_by_comp)
    }
    
    names(mgcca.results) = c("raw_results", "best_start", "sign_changed_results", 
                             "cor_1", "cor_2", "cor_1_by_comp", "cor_2_by_comp")
    
    save(file = MGCCA_files[i], mgcca.results)
  }
  MGCCA$cor_1         = MGCCA_cor_1
  MGCCA$cor_2         = MGCCA_cor_2
  MGCCA$cor_1_by_comp = MGCCA_cor_1_by_comp
  MGCCA$cor_2_by_comp = MGCCA_cor_2_by_comp
  
  save(file = paste0(current_SNR_path, "/MGCCA_reduced_", deflation_mode, "_SNR_", idx_SNR[j], ".Rdata"), MGCCA)
}

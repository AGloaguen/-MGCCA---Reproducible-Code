rm(list = ls())

args = commandArgs(trailingOnly = TRUE)

if(length(args) != 1){
  stop("Current path have to be specified")
}else{
  current_path   = args[1]
}

library(abind)

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
  
  RGCCA_files          = list.files(path = paste0(current_SNR_path, "/RGCCA"), pattern = paste0("rgcca.results_Fold_"), full.names = T)
  
  RGCCA_cor_2 = c()
  RGCCA       = list()
  n_simu      = length(RGCCA_files)
  
  for (i in 1:n_simu){
    load(RGCCA_files[i])
    
    #Gestion du signe des corrélations
    RGCCA_sign_1        = drop(sign(diag(cor(rgcca.results$raw_results$a[[1]], KRON_W_1))))
    RGCCA_sign_2        = drop(sign(diag(cor(rgcca.results$raw_results$a[[2]], KRON_W_2))))
    
    rgcca.result        = rgcca.results$raw_results
    rgcca.result$a[[1]] = rgcca.results$raw_results$a[[1]] %*% diag(RGCCA_sign_1)
    rgcca.result$a[[2]] = rgcca.results$raw_results$a[[2]] %*% diag(RGCCA_sign_2)
    
    rgcca.results[[3]] = rgcca.result
    
    #Evaluation de la correlation avec la vérité
    cor_2 = cor(c(rgcca.result$a[[1]], rgcca.result$a[[2]]), c(KRON_W_1, KRON_W_2))
    
    rgcca.results[[4]] = cor_2
    
    cor_2_by_comp = diag(cor(rbind(rgcca.result$a[[1]], rgcca.result$a[[2]]), rbind(KRON_W_1, KRON_W_2)))
    
    rgcca.results[[5]] = cor_2_by_comp
    
    RGCCA_cor_2 = c(RGCCA_cor_2, cor_2)
    if(i == 1){
      dim_a_1             = dim(rgcca.result$a[[1]])
      dim_a_2             = dim(rgcca.result$a[[2]])
      RGCCA$a[[1]]        = array(data = rgcca.result$a[[1]], dim = c(dim_a_1, 1))
      RGCCA$a[[2]]        = array(data = rgcca.result$a[[2]], dim = c(dim_a_2, 1))
      crit_list           = rgcca.results$raw_results$crit
      crit                = ifelse(is.list(crit_list), 
                                   Reduce("sum", lapply(crit_list, max)),
                                   max(crit_list))
      RGCCA$crit          = crit
      RGCCA_cor_2_by_comp = cor_2_by_comp
    }else{
      RGCCA$a[[1]]        = abind(RGCCA$a[[1]], array(data = rgcca.result$a[[1]], dim = c(dim_a_1, 1)), along = 3)
      RGCCA$a[[2]]        = abind(RGCCA$a[[2]], array(data = rgcca.result$a[[2]], dim = c(dim_a_2, 1)), along = 3)
      crit_list           = rgcca.results$raw_results$crit
      crit                = ifelse(is.list(crit_list), 
                                   Reduce("sum", lapply(crit_list, max)), 
                                   max(crit_list))
      RGCCA$crit          = c(RGCCA$crit, crit)
      RGCCA_cor_2_by_comp = rbind(RGCCA_cor_2_by_comp, cor_2_by_comp)
    }
    
    names(rgcca.results) = c("raw_results", "best_start", "sign_changed_results", "cor_2", "cor_2_by_comp")
    
    save(file = RGCCA_files[i], rgcca.results)
  }
  RGCCA$cor_2         = RGCCA_cor_2
  RGCCA$cor_2_by_comp = RGCCA_cor_2_by_comp
  
  save(file = paste0(current_SNR_path, "/RGCCA_reduced_", "_SNR_", idx_SNR[j], ".Rdata"), RGCCA)
}

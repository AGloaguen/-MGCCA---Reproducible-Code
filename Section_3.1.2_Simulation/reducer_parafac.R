rm(list = ls())

args = commandArgs(trailingOnly = TRUE)

library(abind)

if(length(args) != 2){
  stop("Nummber components/current path have to be specified")
}else{
  nfac           = as.numeric(args[1])
  current_path   = args[2]
}

norm2 <-function(vec){
  a <- sqrt(sum(vec^2))
  if(a==0) a <- .05
  return(a)
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
  
  nb_comp_truth_1 = dim(W_C_1)[2]
  nb_comp_truth_2 = dim(W_C_2)[2]
  
  KRON_W_1 = sapply(1:nb_comp_truth_1, function(x) W_C_1[, x] %x% W_B_1[, x])
  KRON_W_2 = sapply(1:nb_comp_truth_2, function(x) W_C_2[, x] %x% W_B_2[, x])
  
  PARAFAC_files = list.files(path = paste0(current_SNR_path, "/Parafac"), pattern = paste0("parafac.results_nbcomp_", nfac, "_Fold_"), full.names = T)
  
  comb_bloc_1 = t(combn(x = 1:nfac, m = nb_comp_truth_1))
  comb_bloc_1 = rbind(comb_bloc_1, comb_bloc_1[ , nb_comp_truth_1:1])
  comb_bloc_2 = t(combn(x = 1:nfac, m = nb_comp_truth_2))
  comb_bloc_2 = rbind(comb_bloc_2, comb_bloc_2[ , nb_comp_truth_2:1])
  
  cross_comb_idx = expand.grid(1:dim(comb_bloc_1)[1], 1:dim(comb_bloc_2)[1])
  cross_comb     = cbind(comb_bloc_1[cross_comb_idx[, 1], ], comb_bloc_2[cross_comb_idx[, 2], ])
  
  PARAFAC_cor_1 = c()
  PARAFAC_cor_2 = c()
  PARAFAC       = list()
  n_simu        = length(PARAFAC_files)
  
  for (i in 1:n_simu){
    load(PARAFAC_files[i])
    
    PF_1  = parafac.results$raw_results_1
    PF_2  = parafac.results$raw_results_2
    
    PF_1$B = apply(PF_1$B, 2, function(x) x/norm2(x))
    PF_1$C = apply(PF_1$C, 2, function(x) x/norm2(x))
    PF_2$B = apply(PF_2$B, 2, function(x) x/norm2(x))
    PF_2$C = apply(PF_2$C, 2, function(x) x/norm2(x))
    
    cor_comb = apply(cross_comb, 1, function(x){
      PF_1_B  = PF_1$B[, x[1:2]]
      PF_2_B  = PF_2$B[, x[3:4]]
      PF_1_C  = PF_1$C[, x[1:2]]
      PF_2_C  = PF_2$C[, x[3:4]]
      
      PARAFAC_sign_1_B = drop(sign(diag(cor(PF_1_B, W_B_1))))
      PARAFAC_sign_2_B = drop(sign(diag(cor(PF_2_B, W_B_2))))
      PARAFAC_sign_1_C = drop(sign(diag(cor(PF_1_C, W_C_1))))
      PARAFAC_sign_2_C = drop(sign(diag(cor(PF_2_C, W_C_2))))
      
      PF_1_B = PF_1_B %*% diag(PARAFAC_sign_1_B)
      PF_2_B = PF_2_B %*% diag(PARAFAC_sign_2_B)
      PF_1_C = PF_1_C %*% diag(PARAFAC_sign_1_C)
      PF_2_C = PF_2_C %*% diag(PARAFAC_sign_2_C)
      
      KRON_1 = mapply(function(y, z) y %x% z, as.data.frame(PF_1_C), as.data.frame(PF_1_B))
      KRON_2 = mapply(function(y, z) y %x% z, as.data.frame(PF_2_C), as.data.frame(PF_2_B))
      
      cor_1 = cor(c(W_B_1, W_B_2, W_C_1, W_C_2), c(PF_1_B, PF_2_B, PF_1_C, PF_2_C))
      cor_2 = cor(c(KRON_1, KRON_2), c(KRON_W_1, KRON_W_2))
      
      return(c(cor_1, cor_2))
    })
    
    cor_comb  = t(cor_comb)
    idx_cor_1 = cross_comb[which.max(cor_comb[, 1]), ]
    idx_cor_2 = cross_comb[which.max(cor_comb[, 2]), ]
    
    PF_1_B  = PF_1$B[, idx_cor_2[1:2]]
    PF_2_B  = PF_2$B[, idx_cor_2[3:4]]
    PF_1_C  = PF_1$C[, idx_cor_2[1:2]]
    PF_2_C  = PF_2$C[, idx_cor_2[3:4]]
    
    PARAFAC_sign_1_B = drop(sign(diag(cor(PF_1_B, W_B_1))))
    PARAFAC_sign_2_B = drop(sign(diag(cor(PF_2_B, W_B_2))))
    PARAFAC_sign_1_C = drop(sign(diag(cor(PF_1_C, W_C_1))))
    PARAFAC_sign_2_C = drop(sign(diag(cor(PF_2_C, W_C_2))))
    
    PF_1_B = PF_1_B %*% diag(PARAFAC_sign_1_B)
    PF_2_B = PF_2_B %*% diag(PARAFAC_sign_2_B)
    PF_1_C = PF_1_C %*% diag(PARAFAC_sign_1_C)
    PF_2_C = PF_2_C %*% diag(PARAFAC_sign_2_C)
    
    parafac.results[[5]]  = idx_cor_1
    parafac.results[[6]]  = idx_cor_2
    cor_1                 = max(cor_comb[, 1])
    cor_2                 = max(cor_comb[, 2])
    parafac.results[[7]]  = cor_1
    parafac.results[[8]]  = cor_2
    KRON_1                = mapply(function(y, z) y %x% z, as.data.frame(PF_1_C), as.data.frame(PF_1_B))
    KRON_2                = mapply(function(y, z) y %x% z, as.data.frame(PF_2_C), as.data.frame(PF_2_B))
    cor_1_by_comp         = diag(cor(rbind(W_B_1, W_B_2, W_C_1, W_C_2), rbind(PF_1_B, PF_2_B, PF_1_C, PF_2_C)))
    cor_2_by_comp         = diag(cor(rbind(KRON_1, KRON_2), rbind(KRON_W_1, KRON_W_2)))
    parafac.results[[9]]  = cor_1_by_comp
    parafac.results[[10]] = cor_2_by_comp
    
    PARAFAC_cor_1 = c(PARAFAC_cor_1, cor_1)
    PARAFAC_cor_2 = c(PARAFAC_cor_2, cor_2)
    if(i == 1){
      dim_b_1                   = dim(PF_1_B)
      dim_c_1                   = dim(PF_1_C)
      dim_b_2                   = dim(PF_2_B)
      dim_c_2                   = dim(PF_2_C)
      PARAFAC$b[[1]]            = array(data = PF_1_B, dim = c(dim_b_1, 1))
      PARAFAC$c[[1]]            = array(data = PF_1_C, dim = c(dim_c_1, 1))
      PARAFAC$b[[2]]            = array(data = PF_2_B, dim = c(dim_b_2, 1))
      PARAFAC$c[[2]]            = array(data = PF_2_C, dim = c(dim_c_2, 1))
      PARAFAC$best_comp_NO_mway = idx_cor_2
      PARAFAC_cor_1_by_comp     = cor_1_by_comp
      PARAFAC_cor_2_by_comp     = cor_2_by_comp
    }else{
      PARAFAC$b[[1]]            = abind(PARAFAC$b[[1]], array(data = PF_1_B, dim = c(dim_b_1, 1)), along = 3)
      PARAFAC$c[[1]]            = abind(PARAFAC$c[[1]], array(data = PF_1_C, dim = c(dim_c_1, 1)), along = 3)
      PARAFAC$b[[2]]            = abind(PARAFAC$b[[2]], array(data = PF_2_B, dim = c(dim_b_2, 1)), along = 3)
      PARAFAC$c[[2]]            = abind(PARAFAC$c[[2]], array(data = PF_2_C, dim = c(dim_c_2, 1)), along = 3)
      PARAFAC$best_comp_NO_mway = rbind(PARAFAC$best_comp_NO_mway, idx_cor_2)
      PARAFAC_cor_1_by_comp     = rbind(PARAFAC_cor_1_by_comp, cor_1_by_comp)
      PARAFAC_cor_2_by_comp     = rbind(PARAFAC_cor_2_by_comp, cor_2_by_comp)
    }
    print(i)
    names(parafac.results) = c("raw_results_1", "raw_results_2", "init_bloc_1", "init_bloc_2",
                               "best_comps_cor_1", "best_comps_cor_2", 
                               "cor_1", "cor_2", "cor_1_by_comp", "cor_2_by_comp")
    
    save(file = PARAFAC_files[i], parafac.results)
  }
  PARAFAC$cor_1         = PARAFAC_cor_1
  PARAFAC$cor_2         = PARAFAC_cor_2
  PARAFAC$cor_1_by_comp = PARAFAC_cor_1_by_comp
  PARAFAC$cor_2_by_comp = PARAFAC_cor_2_by_comp
  
  save(file = paste0(current_SNR_path, "/PARAFAC_reduced_nbcomp_", nfac, "_SNR_", idx_SNR[j], ".Rdata"), PARAFAC)
}

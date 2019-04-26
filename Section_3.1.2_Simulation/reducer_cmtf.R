rm(list = ls())

args = commandArgs(trailingOnly = TRUE)

library(R.matlab)
library(abind)

if(length(args) != 3){
  stop("SNR/Number of components/Current path have to be specified")
}else{
  SNR          = as.numeric(args[1])
  nfac         = as.numeric(args[2])
  current_path = args[3]
}

norm2 <-function(vec){
  a <- sqrt(sum(vec^2))
  if(a==0) a <- .05
  return(a)
}

cmtf.results = list()

res_CMTF_Folds = list.files(path = paste0(current_path, "/SNR_", SNR, "/CMTF"), pattern = paste0("cmtf.results_nbcomp_", nfac, "_Fold_"), full.names = T)
Truth          = list.files(path = paste0(current_path, "/SNR_", SNR, "/"), pattern = "W_", full.names = T)

for (i in 1:length(Truth)){
  load(file = Truth[i])
}

nb_comp_truth_1 = dim(W_C_1)[2]
nb_comp_truth_2 = dim(W_C_2)[2]

KRON_W_1 = sapply(1:nb_comp_truth_1, function(x) W_C_1[, x] %x% W_B_1[, x])
KRON_W_2 = sapply(1:nb_comp_truth_2, function(x) W_C_2[, x] %x% W_B_2[, x])

cross_comb     = t(combn(x = 1:nfac, m = nb_comp_truth_1))
cross_comb     = rbind(cross_comb, cross_comb[ , nb_comp_truth_1:1])

CMTF_cor_1 = c()
CMTF_cor_2 = c()
CMTF       = list()
for (i in 1:length(res_CMTF_Folds)){
  cur_CMTF = readMat(res_CMTF_Folds[i])
  cur_CMTF = cur_CMTF$x
  
  cur_CMTF$u[[2]][[1]] = apply(cur_CMTF$u[[2]][[1]], 2, function(x) x/norm2(x))
  cur_CMTF$u[[3]][[1]] = apply(cur_CMTF$u[[3]][[1]], 2, function(x) x/norm2(x))
  cur_CMTF$u[[4]][[1]] = apply(cur_CMTF$u[[4]][[1]], 2, function(x) x/norm2(x))
  cur_CMTF$u[[5]][[1]] = apply(cur_CMTF$u[[5]][[1]], 2, function(x) x/norm2(x))
  
  cor_comb = apply(cross_comb, 1, function(x){
    CMTF_1_B  = cur_CMTF$u[[2]][[1]][, x]
    CMTF_2_B  = cur_CMTF$u[[4]][[1]][, x]
    CMTF_1_C  = cur_CMTF$u[[3]][[1]][, x]
    CMTF_2_C  = cur_CMTF$u[[5]][[1]][, x]
    
    CMTF_sign_1_B = drop(sign(diag(cor(CMTF_1_B, W_B_1))))
    CMTF_sign_2_B = drop(sign(diag(cor(CMTF_2_B, W_B_2))))
    CMTF_sign_1_C = drop(sign(diag(cor(CMTF_1_C, W_C_1))))
    CMTF_sign_2_C = drop(sign(diag(cor(CMTF_2_C, W_C_2))))
    
    CMTF_1_B = CMTF_1_B %*% diag(CMTF_sign_1_B)
    CMTF_2_B = CMTF_2_B %*% diag(CMTF_sign_2_B)
    CMTF_1_C = CMTF_1_C %*% diag(CMTF_sign_1_C)
    CMTF_2_C = CMTF_2_C %*% diag(CMTF_sign_2_C)
    
    KRON_1 = mapply(function(y, z) y %x% z, as.data.frame(CMTF_1_C), as.data.frame(CMTF_1_B))
    KRON_2 = mapply(function(y, z) y %x% z, as.data.frame(CMTF_2_C), as.data.frame(CMTF_2_B))
    
    cor_1 = cor(c(W_B_1, W_B_2, W_C_1, W_C_2), c(CMTF_1_B, CMTF_2_B, CMTF_1_C, CMTF_2_C))
    cor_2 = cor(c(KRON_1, KRON_2), c(KRON_W_1, KRON_W_2))
    
    return(c(cor_1, cor_2))
  })
  
  cor_comb  = t(cor_comb)
  idx_cor_1 = cross_comb[which.max(cor_comb[, 1]), ]
  idx_cor_2 = cross_comb[which.max(cor_comb[, 2]), ]
  cor_1     = max(cor_comb[, 1])
  cor_2     = max(cor_comb[, 2])
  
  CMTF_1_B  = cur_CMTF$u[[2]][[1]][, idx_cor_2]
  CMTF_2_B  = cur_CMTF$u[[4]][[1]][, idx_cor_2]
  CMTF_1_C  = cur_CMTF$u[[3]][[1]][, idx_cor_2]
  CMTF_2_C  = cur_CMTF$u[[5]][[1]][, idx_cor_2]
  
  CMTF_sign_1_B = drop(sign(diag(cor(CMTF_1_B, W_B_1))))
  CMTF_sign_2_B = drop(sign(diag(cor(CMTF_2_B, W_B_2))))
  CMTF_sign_1_C = drop(sign(diag(cor(CMTF_1_C, W_C_1))))
  CMTF_sign_2_C = drop(sign(diag(cor(CMTF_2_C, W_C_2))))
  
  CMTF_1_B = CMTF_1_B %*% diag(CMTF_sign_1_B)
  CMTF_2_B = CMTF_2_B %*% diag(CMTF_sign_2_B)
  CMTF_1_C = CMTF_1_C %*% diag(CMTF_sign_1_C)
  CMTF_2_C = CMTF_2_C %*% diag(CMTF_sign_2_C)
  
  cmtf.results[[1]] = cur_CMTF
  cmtf.results[[2]] = list(CMTF_1_B, CMTF_1_C, CMTF_2_B, CMTF_2_C)
  cmtf.results[[3]] = idx_cor_1
  cmtf.results[[4]] = idx_cor_2
  cmtf.results[[5]] = cor_1
  cmtf.results[[6]] = cor_2
  
  KRON_1            = mapply(function(y, z) y %x% z, as.data.frame(CMTF_1_C), as.data.frame(CMTF_1_B))
  KRON_2            = mapply(function(y, z) y %x% z, as.data.frame(CMTF_2_C), as.data.frame(CMTF_2_B))
  cor_1_by_comp     = diag(cor(rbind(W_B_1, W_B_2, W_C_1, W_C_2), rbind(CMTF_1_B, CMTF_2_B, CMTF_1_C, CMTF_2_C)))
  cor_2_by_comp     = diag(cor(rbind(KRON_1, KRON_2), rbind(KRON_W_1, KRON_W_2)))
  cmtf.results[[7]] = cor_1_by_comp
  cmtf.results[[8]] = cor_2_by_comp
  
  names(cmtf.results) = c("raw_results", "used_components", "best_comp_cor_1", "best_comp_cor_2",
                          "cor_1", "cor_2", "cor_1_by_comp", "cor_2_by_comp")
  
  names(cmtf.results[[2]]) = c("B_1", "C_1", "B_2", "C_2")
  
  save(file = paste0(current_path, "/SNR_", SNR, "/CMTF/cmtf.results_Fold_", i, ".Rdata"), cmtf.results)
  
  CMTF_cor_1 = c(CMTF_cor_1, cor_1)
  CMTF_cor_2 = c(CMTF_cor_2, cor_2)
  if(i == 1){
    dim_b_1                = dim(CMTF_1_B)
    dim_c_1                = dim(CMTF_1_C)
    dim_b_2                = dim(CMTF_2_B)
    dim_c_2                = dim(CMTF_2_C)
    CMTF$b[[1]]            = array(data = CMTF_1_B, dim = c(dim_b_1, 1))
    CMTF$c[[1]]            = array(data = CMTF_1_C, dim = c(dim_c_1, 1))
    CMTF$b[[2]]            = array(data = CMTF_2_B, dim = c(dim_b_2, 1))
    CMTF$c[[2]]            = array(data = CMTF_2_C, dim = c(dim_c_2, 1))
    CMTF$best_comp_NO_mway = idx_cor_2
    CMTF_cor_1_by_comp     = cor_1_by_comp
    CMTF_cor_2_by_comp     = cor_2_by_comp
  }else{
    CMTF$b[[1]]            = abind(CMTF$b[[1]], array(data = CMTF_1_B, dim = c(dim_b_1, 1)), along = 3)
    CMTF$c[[1]]            = abind(CMTF$c[[1]], array(data = CMTF_1_C, dim = c(dim_c_1, 1)), along = 3)
    CMTF$b[[2]]            = abind(CMTF$b[[2]], array(data = CMTF_2_B, dim = c(dim_b_2, 1)), along = 3)
    CMTF$c[[2]]            = abind(CMTF$c[[2]], array(data = CMTF_2_C, dim = c(dim_c_2, 1)), along = 3)
    CMTF$best_comp_NO_mway = rbind(CMTF$best_comp_NO_mway, idx_cor_2)
    CMTF_cor_1_by_comp     = rbind(CMTF_cor_1_by_comp, cor_1_by_comp)
    CMTF_cor_2_by_comp     = rbind(CMTF_cor_2_by_comp, cor_2_by_comp)
  }
}

CMTF$cor_1         = CMTF_cor_1
CMTF$cor_2         = CMTF_cor_2
CMTF$cor_1_by_comp = CMTF_cor_1_by_comp
CMTF$cor_2_by_comp = CMTF_cor_2_by_comp

save(file = paste0(current_path, "/SNR_", SNR, "/CMTF_reduced_nbcomp_", nfac, "_SNR_", SNR, ".Rdata"), CMTF)

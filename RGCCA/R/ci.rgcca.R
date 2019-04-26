ci.rgcca = function(object, A, B, alpha = 0.05, ndim = 1, verbose = FALSE, plot = FALSE, nb_cores = 4) {
  n = NROW(A[[1]])
  p = sapply(A, ncol)
  J = length(A)
  
  ##############################################
  # Initialization of the outer weight vectors #
  ##############################################
  if(any(object$ncomp!=1)){
    Yinit <- sapply(object$Y, function(x) x[, ndim])
  }else{
    Yinit <- sapply(object$Y, cbind)
  }  
  W = lapply(object$a, function(x) x[, ndim])
  
  ########################################
  # Construction of the Boostrap samples #
  ########################################
  if (verbose){
    boot_b  = pbmcapply::pbmclapply(1:B, function(z) bootstrap(n = n, J = J, A = A, object = object, W = W, ndim = ndim), mc.cores = nb_cores)
  }else{
    boot_b  = parallel::mclapply(1:B, function(z) bootstrap(n = n, J = J, A = A, object = object, W = W, ndim = ndim), mc.cores = nb_cores)
  }
  
  ########################################
  # Construction of the Coefficeints'CI  #
  ########################################
  Astar = NULL
  for (i in 1:J){
    Astar[[i]] = sapply(1:B, FUN = function(x) boot_b[[x]][[2]][[i]])
  }
  
  M1 = lapply(Astar, function(w) apply(w, 1,  function(x) c(mean(x), sd(x))))
  
  mat  = list()
  tail = qnorm(1-alpha/(2))
  for (j in 1:J){
    mat[[j]]           <- cbind(W[[j]], M1[[j]][1, ]-tail*M1[[j]][2, ], M1[[j]][1, ]+tail*M1[[j]][2, ])
    rownames(mat[[j]]) <- colnames(A[[j]])
    colnames(mat[[j]]) <- c("Initial weights","Lower Bound","Upper Bound")
  }
  
  ########################################
  # Construction of the Correlations'CI  #
  ########################################
  MAT_COR = t(sapply(1:B, FUN = function(x) boot_b[[x]][[1]]))
  M2      = apply(MAT_COR, 2,  function(x) c(mean(x), sd(x)))
  
  inner_relation           = matrix(0, 3, (J*(J-1))/2)
  connection               = cbind(rep(paste0("X", 1:J), J), rep(paste0("X", 1:J), each = J))
  colnames(inner_relation) = matrix(paste0(connection[, 1], "-", connection[, 2]), J, J)[upper.tri(diag(J))]
  rownames(inner_relation) = c("Initial Coorelation","Lower Bound","Upper Bound")
  inner_relation[1, ]      = (cor(Yinit)*object$C)[upper.tri(diag(J))]
  inner_relation[2, ]      = (matrix(M2[1, ]-tail*M2[2, ], J, J)*object$C)[upper.tri(diag(J))]
  inner_relation[3, ]      = (matrix(M2[1, ]+tail*M2[2, ], J, J)*object$C)[upper.tri(diag(J))]
  inner_relation           = t(inner_relation)
  
  
  ########################################
  # Construction of a vizualisation      #
  ########################################
  if (plot){
    mat = lapply(mat, function(x) x[order(x[, 1], decreasing = TRUE), ])
    
    par(cex = .8)
    for (j in 1:J){
      color = rep("red", nrow(mat[[j]])) ; color[which(mat[[j]][, 2]/mat[[j]][, 3]>0)] = "green3"
      r     = barplot(mat[[j]][, 1], col = color, ylim = c(min(0, min(mat[[j]][, 2])), max(0, mat[[j]][, 3])), las = 2,  omi = c(50, 4, 4, 4))
      segments(r, mat[[j]][, 2], r, mat[[j]][, 3])
      segments(r-0.1, mat[[j]][, 2], r+0.1, mat[[j]][, 2])
      segments(r-0.1, mat[[j]][, 3], r+0.1, mat[[j]][, 3])
    }
  }
  return(list(a_boot = Astar, CI = mat, inner_relation = inner_relation))
}

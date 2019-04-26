generate_simu_data_KKT = function(N, J, K, L, mean_y, mean_noise, sd_noise, Sigma){
  
  A = list()
  Y = mvrnorm(n = N, mu = mean_y, Sigma = Sigma, empirical = T)
  
  for (i in 1:L){
    mydim  = c(N,J[i],K[i])
    nf     = 1
    Amat   = Y[, i]
    Bmat   = matrix(runif(mydim[2]*nf),mydim[2],nf)
    Cmat   = matrix(runif(mydim[3]*nf),mydim[3],nf)
    Xmat   = array(tcrossprod(Amat,krprod(Cmat,Bmat)),dim=mydim)
    Emat   = array(sapply(sd_noise[[i]], function(x) rnorm(prod(mydim[1:2]), mean = mean_noise, sd = x)),dim=mydim)
    X      = Xmat + Emat
    A[[i]] = X
  }
  return(A)
}
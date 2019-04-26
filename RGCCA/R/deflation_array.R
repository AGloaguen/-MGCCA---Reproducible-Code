deflation_array <- function(X, y, b, c){
  # Computation of the residual matrix R
  X_hat <- array(c(as.matrix(c %x% y %x% b)), dim(X))
  R     <- X - X_hat
  return(list(R=R))
}

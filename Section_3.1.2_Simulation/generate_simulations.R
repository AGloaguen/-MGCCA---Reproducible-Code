rm(list = ls())

args = commandArgs(trailingOnly = TRUE)

if(length(args) != 5){
  stop("SNR/Number of components/Current path/Kathri Rao have to be specified")
}else{
  SNR          = as.numeric(args[1])
  nfac         = as.numeric(args[2])
  current_path = args[3]
  katri_rao    = eval(parse(text=args[4]))
  n_simu       = as.numeric(args[5])
}

library(mvtnorm)
library(multiway)
library(abind)
library(pracma)
library(R.matlab)

dir.create(file.path(current_path, paste0("SNR_", SNR)), showWarnings = FALSE)
dir.create(file.path(current_path, paste0("SNR_", SNR, "/Parafac")), showWarnings = FALSE)
dir.create(file.path(current_path, paste0("SNR_", SNR, "/MGCCA")), showWarnings = FALSE)
dir.create(file.path(current_path, paste0("SNR_", SNR, "/RGCCA")), showWarnings = FALSE)
dir.create(file.path(current_path, paste0("SNR_", SNR, "/CMTF")), showWarnings = FALSE)
dir.create(file.path(current_path, paste0("SNR_", SNR, "/matfiles")), showWarnings = FALSE)


if (require("devtools")){
  load_all(paste0(current_path, "/../../RGCCA/."))
}else{
  file.sources = list.files(path = paste0(current_path, "/../../RGCCA/R", sep = ""), pattern = "*.R", full.names = T)
  sapply(file.sources, source, .GlobalEnv)
}

N      = 50
J      = c(30, 30)
K      = c(30, 30)
L      = length(K)

#Initialisation de la matrice de variance-covariance
Sigma = matrix(0, ncol = sum(J*K), nrow = sum(J*K))

#Zones d'intéraction entre X_1 et X_2
Interaction_comp1 = list(list(6:15, 11:15), list(c(16:25), c(6:8, 13:15)))
Interaction_comp2 = list(list(c(1:5, 16:20), 26:30), list(c(1:10), c(21:25)))
alpha             = 0.7

random_DP = function(p){
  random_MAT    = matrix(rnorm(n = p*p), nrow = p, ncol = p)
  random_SVD    = svd(x = random_MAT)
  eigen_values  = runif(n = p, min = 10^-10, max = 1)
  DP            = random_SVD$u %*% diag(eigen_values) %*% t(random_SVD$u)
  return(list(matrix = DP, eigen_vectors = random_SVD$u, eigen_values = eigen_values))
}

Sigma_11 = Sigma_22 = matrix(0, ncol = sum(J*K), nrow = sum(J*K))
if(katri_rao){
  random_DP_K_1 = random_DP(K[1])
  random_DP_J_1 = random_DP(J[1])
  Sigma_11_tmp  = mapply(function(x, y, z) list(z * x %x% y %*% t(x %x% y)), 
                         as.data.frame(random_DP_K_1$eigen_vectors), 
                         as.data.frame(random_DP_J_1$eigen_vectors),
                         as.data.frame(t(random_DP_K_1$eigen_values)))
  
  random_DP_K_2 = random_DP(K[2])
  random_DP_J_2 = random_DP(J[2])
  Sigma_22_tmp  = mapply(function(x, y, z) list(z * x %x% y %*% t(x %x% y)), 
                         as.data.frame(random_DP_K_2$eigen_vectors), 
                         as.data.frame(random_DP_J_2$eigen_vectors),
                         as.data.frame(t(random_DP_K_2$eigen_values)))
  
  
  Sigma_11[1:900, 1:900]       = Reduce("+", Sigma_11_tmp)
  Sigma_22[901:1800, 901:1800] = Reduce("+", Sigma_22_tmp)
}else{
  Sigma_11[1:900, 1:900]       = random_DP(J[1]*K[1])$matrix
  Sigma_22[901:1800, 901:1800] = random_DP(J[2]*K[2])$matrix
}

#Ajout de la variance de l'intéraction. Les vecterus associées à l'intéraction sont gardées
# afin d'être confrontés plus tard aux résultats.
W = list()

ind_1             = Interaction_comp1[[1]]
xJ1               = matrix(0, ncol = 1, nrow = J[1])
xK1               = matrix(0, ncol = 1, nrow = K[1])
xJ1[ind_1[[1]], ] = 1
xK1[ind_1[[2]], ] = 1

ind_2             = Interaction_comp1[[2]]
xJ2               = matrix(0, ncol = 1, nrow = J[2])
xK2               = matrix(0, ncol = 1, nrow = K[2])
xJ2[ind_2[[1]], ] = 1
xK2[ind_2[[2]], ] = 1

W[[1]]         = list(xJ1, xK1)
W[[2]]         = list(xJ2, xK2)

ind_1                   = Interaction_comp2[[1]]
xJ1                     = matrix(0, ncol = 1, nrow = J[1])
xK1                     = matrix(0, ncol = 1, nrow = K[1])
xJ1[ind_1[[1]][1:5], ]  = 1
xJ1[ind_1[[1]][6:10], ] = -1
xK1[ind_1[[2]], ]       = 1

ind_2                   = Interaction_comp2[[2]]
xJ2                     = matrix(0, ncol = 1, nrow = J[2])
xK2                     = matrix(0, ncol = 1, nrow = K[2])
xJ2[ind_2[[1]][1:5], ]  = 1
xJ2[ind_2[[1]][6:10], ] = -1
xK2[ind_2[[2]], ]       = 1

W[[1]][[1]] = cbind(W[[1]][[1]], xJ1)
W[[1]][[2]] = cbind(W[[1]][[2]], xK1) 
W[[2]][[1]] = cbind(W[[2]][[1]], xJ2)
W[[2]][[2]] = cbind(W[[2]][[2]], xK2) 

W[[1]][[1]] = apply(W[[1]][[1]], 2, function(x) x/norm2(x))
W[[1]][[2]] = apply(W[[1]][[2]], 2, function(x) x/norm2(x))
W[[2]][[1]] = apply(W[[2]][[1]], 2, function(x) x/norm2(x))
W[[2]][[2]] = apply(W[[2]][[2]], 2, function(x) x/norm2(x))

x = c(W[[1]][[2]][, 1] %x% W[[1]][[1]][, 1], W[[2]][[2]][, 1] %x% W[[2]][[1]][, 1])
x = cbind(x, c(W[[1]][[2]][, 2] %x% W[[1]][[1]][, 2], W[[2]][[2]][, 2] %x% W[[2]][[1]][, 2])) 

SIGNAL = x %*% diag(c(alpha, 1-alpha)) %*% t(x)
NOISE  = Sigma_11 + Sigma_22

SNR    = SNR
noise  = drop(sqrt(crossprod(c(NOISE))))
signal = drop(sqrt(crossprod(c(SIGNAL))))

Sigma = (signal/noise/SNR)*NOISE + SIGNAL

#Normalisation de poids
W_B_1 = W[[1]][[1]]
W_B_2 = W[[2]][[1]]
W_C_1 = W[[1]][[2]]
W_C_2 = W[[2]][[2]]

save(file = paste0(current_path, "/SNR_", SNR, "/W_B_1.Rdata"), W_B_1)
save(file = paste0(current_path, "/SNR_", SNR, "/W_B_2.Rdata"), W_B_2)
save(file = paste0(current_path, "/SNR_", SNR, "/W_C_1.Rdata"), W_C_1)
save(file = paste0(current_path, "/SNR_", SNR, "/W_C_2.Rdata"), W_C_2)

writeMat(paste0(current_path, "/SNR_", SNR, "/matfiles/W_B_1.mat"), W_B_1 = W_B_1)
writeMat(paste0(current_path, "/SNR_", SNR, "/matfiles/W_B_2.mat"), W_B_2 = W_B_2)
writeMat(paste0(current_path, "/SNR_", SNR, "/matfiles/W_C_1.mat"), W_C_1 = W_C_1)
writeMat(paste0(current_path, "/SNR_", SNR, "/matfiles/W_C_2.mat"), W_C_2 = W_C_2)

#Generation des données
set.seed(1000)
X = mvrnorm(n = n_simu*N, mu = rep(0, sum(J*K)), Sigma = Sigma, empirical = F)

for (i in 1:n_simu){
  ind_simu = 1:N
  ind_simu = ind_simu + (i - 1)*N

  #Séparation en 2 blocs
  Xm_1 = X[ind_simu, 1:(J*K)[2]]
  Xm_2 = X[ind_simu, ((J*K)[2] + 1):sum(J*K)]
  #Création des tenseurs
  Xu_1 = array(c(Xm_1), dim = c(N, J[1], K[1]))
  Xu_2 = array(c(Xm_2), dim = c(N, J[2], K[2]))

  A       = list()
  A[[1]]  = Xu_1
  A[[2]]  = Xu_2

  A = lapply(1:L, function(x) scale3_array(A[[x]], bias = T, scale = F, center = T))
  A = lapply(A, function(x) x/sd(c(x)))
  
  INIT        = list()
  INIT[[1]]   = svd(cbind(t(apply(A[[1]], 1, c)), t(apply(A[[2]], 1, c)) ), nu = nfac, nv = 0)$u
  INIT[[2]]   = svd(t(apply(A[[1]], 2, c)), nu = nfac, nv = 0)$u
  INIT[[3]]   = svd(t(apply(A[[1]], 3, c)), nu = nfac, nv = 0)$u
  INIT[[4]]   = svd(t(apply(A[[2]], 2, c)), nu = nfac, nv = 0)$u
  INIT[[5]]   = svd(t(apply(A[[2]], 3, c)), nu = nfac, nv = 0)$u
  names(INIT) = c("init_A", "init_B1", "init_C1", "init_B2", "init_C2")

  save(file = paste0(current_path, "/SNR_", SNR, "/Simu_Fold_", i, ".Rdata"), A, INIT)

  writeMat(paste0(current_path, "/SNR_", SNR, "/matfiles/Simu_Fold_", i, "_T1.mat"), T1 = A[[1]])
  writeMat(paste0(current_path, "/SNR_", SNR, "/matfiles/Simu_Fold_", i, "_T2.mat"), T2 = A[[2]])
  writeMat(paste0(current_path, "/SNR_", SNR, "/matfiles/INIT_simu_Fold_", i, ".mat"), INIT = INIT)

  print(i)
}






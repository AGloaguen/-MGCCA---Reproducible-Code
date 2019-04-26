rm(list = ls())

###########################
#####    LIBRARIES    #####
###########################
library(Matrix)
library(multiway)

if (require("devtools")){
  load_all("../RGCCA/.")
}else{
  file.sources = list.files(path = paste0(current_path, "../RGCCA/R", sep = ""), pattern = "*.R", full.names = T)
  sapply(file.sources, source, .GlobalEnv)
}
source("generate_simu_data_KKT.R")

###########################
#####    PARAMETERS   #####
###########################
#To simulate data set
J          = c(200, 500, 1000)
K          = c(5, 10, 10)
N          = 90
L          = length(J)
mean_y     = rep(0, L)
mean_noise = rep(0, L)
mean_noise = rep(2, L)
Sigma      = matrix(c(1, 0, 0.7, 0, 1, 0.7, 0.7, 0.7, 1), nrow = 3, byrow = T)
#To apply MGCCA
C       = matrix(c(0, 0, 1, 0, 0, 1, 1, 1, 0), nrow = 3, byrow = T)
tau     = rep(1, L)
ncomp   = rep(1, L)
scheme  = "factorial"
center  = T
scale   = T 
init    = "random"
bias    = T
verbose = F
M_regularisation = rep("kronecker_Identity_RGCCA", L)
#For the simulation
nb_iteration_of_mgcca_algo   = 39
nb_application_of_mgcca_algo = 200
KKT                          = list()
tresholds_KKT                = c(1e-1, 1e-3, 1e-5, 1e-7, 1e-9, 1e-11, 1e-13, 1e-14)

###########################
##### DATA GENERATION #####
###########################
X = generate_simu_data_KKT(N          = N, 
                           J          = J, 
                           K          = K, 
                           L          = L,
                           mean_y     = mean_y, 
                           mean_noise = mean_noise, 
                           sd_noise   = mean_noise,
                           Sigma      = Sigma)
###########################
#####   APPLY MGCCA   #####
###########################
set.seed(53)
for (i in 1:nb_application_of_mgcca_algo){
  res = mgcca_array_KKT(A       = X, 
                        C       = C, 
                        tau     = tau, 
                        ncomp   = ncomp, 
                        scheme  = scheme, 
                        center  = center, 
                        scale   = scale, 
                        init    = init, 
                        bias    = bias, 
                        verbose = verbose,
                        M_regularisation = M_regularisation, 
                        nb_iter          = nb_iteration_of_mgcca_algo)
  KKT[[i]] = res$kkt
  print(i)
}

KKT                   = Reduce("rbind", KKT)
results_KKT           = sapply(tresholds_KKT, function(x) min(which(apply(KKT, 2, function(y) all(y<x)))))
results_KKT           = rbind(tresholds_KKT, results_KKT)
rownames(results_KKT) = c("KKT", "s")
colnames(results_KKT) = rep("", length(tresholds_KKT))
results_KKT

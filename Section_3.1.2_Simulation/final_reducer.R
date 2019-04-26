rm(list = ls())

args = commandArgs(trailingOnly = TRUE)

if(length(args) != 4){
  stop("Current path, number of simulations, SNR to plot and methds to compare have to be sprecified")
}else{
  current_path       = args[1]
  n_simu             = as.numeric(args[2])
  SNR_to_plot        = args[3]
  methods_to_compare = strsplit(args[4], ",")[[1]]
}

library(ggplot2)
library(grid)
library(gridExtra)

print(methods_to_compare)

width  = 60
height = 42

directories = list.dirs(current_path, recursive = F, full.names = F)
SNR_dirs    = directories[ grepl("SNR_", directories) ]

idx_SNR = sapply(SNR_dirs, function(x) as.numeric(gsub("SNR_","",x)))
idx_SNR = unname(idx_SNR)

for (j in 1:length(SNR_dirs)){
  
  current_SNR_path = paste0(current_path, "/", SNR_dirs[j])
  patter_to_rm     = paste0("_reduced|_", SNR_dirs[j], "|.Rdata")
    
  Parafac_files = list.files(path = current_SNR_path, pattern = "PARAFAC_reduced_", full.names = T)
  for (i in 1:length(Parafac_files)){
    load(Parafac_files[i])
    name = gsub(pattern = patter_to_rm, replacement = "", x = basename(Parafac_files[i]))
    if (i == 1 & j == 1){
      nbcomp = dim(PARAFAC$b[[1]])[2]
      df_4   = data.frame(y1 = c(PARAFAC$cor_2_by_comp), 
                          y2 = rep(name, nbcomp*n_simu),
                          y3 = rep(idx_SNR[j], nbcomp*n_simu),
                          y4 = rep("Parafac", nbcomp*n_simu),
                          y5 = rep(paste0("Component n°", 1:nbcomp), each = n_simu))
    }else{
      df_4   = rbind(df_4, 
                     data.frame(y1 = c(PARAFAC$cor_2_by_comp), 
                                y2 = rep(name, nbcomp*n_simu),
                                y3 = rep(idx_SNR[j], nbcomp*n_simu),
                                y4 = rep("Parafac", nbcomp*n_simu),
                                y5 = rep(paste0("Component n°", 1:nbcomp), each = n_simu)))
    }
  }
  
  CMTF_files    = list.files(path = current_SNR_path, pattern = "CMTF_reduced_", full.names = T)
  for (i in 1:length(Parafac_files)){
    load(CMTF_files[i])
    name   = gsub(pattern = patter_to_rm, replacement = "", x = basename(CMTF_files[i]))
    df_4   = rbind(df_4, 
                   data.frame(y1 = c(CMTF$cor_2_by_comp), 
                              y2 = rep(name, nbcomp*n_simu),
                              y3 = rep(idx_SNR[j], nbcomp*n_simu),
                              y4 = rep("CMTF", nbcomp*n_simu),
                              y5 = rep(paste0("Component n°", 1:nbcomp), each = n_simu)))
  }
  
  RGCCA_files   = list.files(path = current_SNR_path, pattern = "RGCCA_reduced_", full.names = T)
  for (i in 1:length(RGCCA_files)){
    load(RGCCA_files[i])
    name   = gsub(pattern = patter_to_rm, replacement = "", x = basename(RGCCA_files[i]))
    df_4   = rbind(df_4, 
                   data.frame(y1 = c(RGCCA$cor_2_by_comp), 
                              y2 = rep(name, nbcomp*n_simu),
                              y3 = rep(idx_SNR[j], nbcomp*n_simu),
                              y4 = rep("RGCCA", nbcomp*n_simu),
                              y5 = rep(paste0("Component n°", 1:nbcomp), each = n_simu)))
  }

  MGCCA_files   = list.files(path = current_SNR_path, pattern = "MGCCA_reduced_", full.names = T)
  for (i in 1:length(MGCCA_files)){
    load(MGCCA_files[i])
    name   = gsub(pattern = patter_to_rm, replacement = "", x = basename(MGCCA_files[i]))
    df_4   = rbind(df_4, 
                   data.frame(y1 = c(MGCCA$cor_2_by_comp), 
                              y2 = rep(name, nbcomp*n_simu),
                              y3 = rep(idx_SNR[j], nbcomp*n_simu),
                              y4 = rep("MGCCA", nbcomp*n_simu),
                              y5 = rep(paste0("Component n°", 1:nbcomp), each = n_simu)))
  }
}


###################################################################################
#######################           Display Results           #######################
###################################################################################
fig_to_plot        = list("Parafac", "CMTF", "RGCCA", "MGCCA")
fig_dir            = paste0(current_path, "/figures")
comp_to_plot       = 2

dir.create(file.path(fig_dir), showWarnings = FALSE)

df_4$y2 = as.factor(df_4$y2)
df_4$y4 = as.factor(df_4$y4)
df_4$y5 = as.factor(df_4$y5)

fig_final_4 = ggplot(df_4[df_4$y2 == methods_to_compare, ], aes(x = 20*log10(y3), y = y1, group = interaction(factor(y3), y2), fill = y4)) +
  geom_boxplot() + theme_bw() + scale_fill_grey(start = 0.5, end = 1)+
  ggtitle("") + theme(plot.title = element_text(hjust = 0.5)) +
  xlab("")+ylab("")+ theme(plot.title = element_text(size = 40, hjust = 0.5),
                              axis.text.x =element_text(size=40, angle = 0),
                              axis.text.y =element_text(size=40),
                              axis.title.x =element_text(size=60),
                              axis.title.y =element_text(size=60), 
                              legend.key.size = unit(3,"line"),
                              legend.text = element_text(size = 40), 
                              legend.title = element_text(size = 0),
                              legend.position = c(0.35, 0.2), 
                              strip.text.x = element_text(size = 40),
                              legend.background = element_rect(colour = NA, fill = NA)) +  facet_wrap(~ y5)


ggsave(filename = paste0(fig_dir, "/All_Methods_comp_separated_NoAxisTitle_DB.pdf"),
       plot = fig_final_4, width = width, height = height, units ="cm")


source(paste0(current_path, "/../display_simu_envelop_V3.R"))
j       = which(SNR_dirs %in% SNR_to_plot)
figName = paste0(fig_dir, "/fig_Weights_OneTensor_OneMode_SNR_", idx_SNR[j], ".pdf")

current_SNR_path = paste0(current_path, "/SNR_", idx_SNR[j])

load(file = paste0(current_SNR_path, "/W_B_1.Rdata"))
load(file = paste0(current_SNR_path, "/W_B_2.Rdata"))
load(file = paste0(current_SNR_path, "/W_C_1.Rdata"))
load(file = paste0(current_SNR_path, "/W_C_2.Rdata"))

TRUTH        = list()
TRUTH$b[[1]] = W_B_1
TRUTH$c[[1]] = W_C_1
TRUTH$b[[2]] = W_B_2
TRUTH$c[[2]] = W_C_2

Parafac_files = list.files(path = current_SNR_path, pattern = "PARAFAC_reduced_", full.names = T)
CMTF_files    = list.files(path = current_SNR_path, pattern = "CMTF_reduced_", full.names = T)
MGCCA_files   = list.files(path = current_SNR_path, pattern = "MGCCA_reduced_", full.names = T)
RGCCA_files   = list.files(path = current_SNR_path, pattern = "RGCCA_reduced_", full.names = T)

patter_to_rm     = paste0("_reduced|_", SNR_dirs[j], "|.Rdata")
idx_Parafac_file = which(sapply(Parafac_files, function(x) gsub(pattern = patter_to_rm, replacement = "", x = basename(x))) 
                         %in% methods_to_compare)
idx_CMTF_file    = which(sapply(CMTF_files, function(x) gsub(pattern = patter_to_rm, replacement = "", x = basename(x))) 
                         %in% methods_to_compare)
idx_MGCCA_file   = which(sapply(MGCCA_files, function(x) gsub(pattern = patter_to_rm, replacement = "", x = basename(x))) 
                         %in% methods_to_compare)

load(Parafac_files[idx_Parafac_file])
load(CMTF_files[idx_CMTF_file])
load(MGCCA_files[idx_MGCCA_file])

idx_MGCCA   = which(MGCCA$cor_2_by_comp[, comp_to_plot] > quantile(MGCCA$cor_2_by_comp[, comp_to_plot], probs = 0.1))
idx_PARAFAC = which(PARAFAC$cor_2_by_comp[, comp_to_plot] > quantile(PARAFAC$cor_2_by_comp[, comp_to_plot], probs = 0.1))
idx_CMTF    = which(CMTF$cor_2_by_comp[, comp_to_plot] > quantile(CMTF$cor_2_by_comp[, comp_to_plot], probs = 0.1))

width  = 60
height = 20

J = unlist(lapply(PARAFAC$b, function(x) dim(x)[1]))
K = unlist(lapply(PARAFAC$c, function(x) dim(x)[1]))
display_simu_envelop_V3(J = J, K = K, 
                        TRUTH   = lapply(TRUTH,   function(x) lapply(x, function(y) y[, comp_to_plot])), 
                        PARAFAC = lapply(PARAFAC[c("b", "c")], function(x) lapply(x, function(y) 
                          apply(y[, comp_to_plot, idx_PARAFAC], 1, function(z) c(mean(z), median(z), min(z), max(z))))), 
                        MGCCA   = lapply(MGCCA[c("b", "c")],   function(x) lapply(x, function(y) 
                          apply(y[, comp_to_plot, idx_MGCCA], 1, function(z) c(mean(z), median(z), min(z), max(z))))), 
                        CMTF    = lapply(CMTF[c("b", "c")],    function(x) lapply(x, function(y) 
                          apply(y[, comp_to_plot, idx_CMTF], 1, function(z) c(mean(z), median(z), min(z), max(z))))), 
                        ft_size = 30, width = width, height = height, figName = figName)
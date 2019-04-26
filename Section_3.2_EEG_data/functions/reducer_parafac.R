args = commandArgs(trailingOnly = TRUE)

library(abind)

pathtoReadFile         = args[1]
map_res_file_pattern   = args[2]
bootstrap_samples      = strtoi(args[3])
output_file            = args[4]
parameter_file         = args[5]


GlobalMapRes_files = list.files(pathtoReadFile, pattern = map_res_file_pattern, full.names = T)

load(GlobalMapRes_files[1])
Ystar_tmp = Ystar
Bstar_tmp = Bstar
Cstar_tmp = Cstar
for(i in 2:length(GlobalMapRes_files)){
  load(GlobalMapRes_files[i])
  Ystar_tmp = mapply(function(x, y) list(abind(x, y, along = 3)), Ystar_tmp, Ystar)
  Bstar_tmp = mapply(function(x, y) list(abind(x, y, along = 3)), Bstar_tmp, Bstar)
  Cstar_tmp = mapply(function(x, y) list(abind(x, y, along = 3)), Cstar_tmp, Cstar)
}

load(parameter_file)
load(data_path)

compute_CI = function(y, alpha){
  mean = mean(y)
  sd   = 1.4826*(quantile(y, probs = 0.75) - quantile(y, probs = 0.25))/2
  Q    = qnorm(p    = c(alpha/2, 1-alpha/2),
               mean = mean, 
               sd   = sd)
  return(c(mean, Q))
}

Ystar = lapply(Ystar_tmp, function(x) apply(x, c(1, 2), function(y) compute_CI(y, alpha = alpha)))
Bstar = lapply(Bstar_tmp, function(x) apply(x, c(1, 2), function(y) compute_CI(y, alpha = alpha)))
Cstar = lapply(Cstar_tmp, function(x) apply(x, c(1, 2), function(y) compute_CI(y, alpha = alpha)))

save(Ystar_tmp = Ystar_tmp,
     Bstar_tmp = Bstar_tmp,
     Cstar_tmp = Cstar_tmp, 
     Ystar     = Ystar,
     Bstar     = Bstar,
     Cstar     = Cstar, 
     file      = output_file)

output_channel_figure_file    = paste(pathtofigurefile, "channel_figure_eeglab_parafac", sep = "")
output_time_figure_file       = paste(pathtofigurefile, "time_figure_ggplot_parafac", sep = "")
channel_figure_parameter_file = paste(pathtofigurefile, "channel_figure_parameter_file_parafac", sep = "")

for(s in to_source_in_reducer) { source(paste(pathtoscript, s, sep="")) }


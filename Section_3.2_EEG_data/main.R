rm(list=ls())

library(rjson)
library(CMA)
library(tools)
library(pracma)
library(R.matlab)

DATE = gsub(pattern     = "  | |:", 
            replacement = "_", 
            x           = date())

right_file_name = function(filename, DATE, name_of_computation, pathtowritefile){
  filename = paste(paste(file_path_sans_ext(filename), 
                         name_of_computation, 
                         DATE, 
                         sep = "_"), 
                   file_ext(filename), 
                   sep = ".")
  filename = paste(pathtowritefile,
                   filename,
                   sep = "")
  return(filename)
}

path_of_main = paste0(pwd(), "/")

#########################################################################
###
###  reading JSON configuration file, parameters setting and data loading
###
#########################################################################
args   = commandArgs(trailingOnly = TRUE)

# Set working directory.
data_json_file   = args[1]
method_json_file = args[2]
data.cfg         = fromJSON(paste(readLines(data_json_file), collapse=""))
method.cfg       = fromJSON(paste(readLines(method_json_file), collapse=""))


#Paths
name_of_computation          = data.cfg$name_of_computation
pathtowritefile_basename     = data.cfg$pathtowritefile
pathtowritefile              = paste0(path_of_main, data.cfg$pathtowritefile)
dir.create(pathtowritefile, showWarnings = FALSE)
pathtowritefile              = paste0(pathtowritefile, 
                                      name_of_computation, 
                                      "_", 
                                      DATE, 
                                      "/")
pathtoscript                 = paste0(path_of_main, method.cfg$pathtoscript)
data_path                    = paste0(path_of_main, data.cfg$data_path)
sources                      = method.cfg$to_source_in_main
command_file                 = data.cfg$command_file
command_file                 = right_file_name(filename            = command_file, 
                                               DATE                = DATE, 
                                               name_of_computation = name_of_computation, 
                                               pathtowritefile     = pathtowritefile)
channels_file                = paste0(path_of_main, data.cfg$channels_file)
pathtofigurefile             = data.cfg$pathtofigurefile
channel_figure_pathtoscript  = paste0(path_of_main, data.cfg$channel_figure_pathtoscript)

# Create pathtowritefile
dir.create(pathtowritefile, showWarnings = FALSE)
dir.create(paste0(pathtowritefile, "bootstrap_parameters_files"), showWarnings = FALSE)
pathtofigurefile = paste(pathtowritefile, pathtofigurefile)
dir.create(pathtofigurefile, showWarnings = FALSE)

# CV parameters
scale                     = method.cfg$scale
bootsrap_samples          = method.cfg$bootsrap_samples
nb_of_bootstrap_group     = method.cfg$nb_of_bootstrap_group
nb_cores                  = method.cfg$nb_cores
scale_size_bloc           = method.cfg$scale_size_bloc
bias                      = method.cfg$bias
tau                       = method.cfg$tau
ncomp                     = method.cfg$ncomp
init                      = method.cfg$init
cv_method                 = method.cfg$cv.method
nstart                    = method.cfg$nstart
tol                       = method.cfg$tol
nfac                      = method.cfg$nfac
parfor_bool               = method.cfg$parfor_bool

if ((bootsrap_samples/nb_of_bootstrap_group) %% nb_cores != 0) stop("nb_of_bootstrap_group is not divisible by nb_cores")

#Parameters reducers
alpha                = method.cfg$alpha
to_source_in_reducer = method.cfg$to_source_in_reducer

#Mappers
mappers      = names(method.cfg$model.mapper)

load(data_path)
N    = nrow(X[[1]])
L    = length(X)
DIM  = lapply(X, dim)
LEN  = unlist(lapply(DIM, length))
B_3D = which(LEN == 3)

M_regularisation = rep(method.cfg$M_regularisation, L)
# source required scripts
for(s in sources) { source(paste(pathtoscript, s, sep="")) }

#########################################################################
###
###  prediction performances assessment
###
#########################################################################
### set parameter for CV
### Generation of the training and testing sets with a method from package CMA.
if(cv_method == "paired_bootstrap"){
  trainmat.outer = CMA::GenerateLearningsets(n      = N/2, 
                                             method = "bootstrap", 
                                             niter  = bootsrap_samples)@learnmatrix
  trainmat.outer = cbind(trainmat.outer, trainmat.outer + N/2)
}else{
  trainmat.outer = CMA::GenerateLearningsets(n      = N, 
                                             method = cv_method, 
                                             niter  = bootsrap_samples)@learnmatrix
}



trainmat.outer_list = tapply(1:bootsrap_samples, 
                             cut(1:bootsrap_samples, 
                                 breaks = nb_of_bootstrap_group), 
                             function(x) list(trainmat.outer[x, ]))


### Cross validation loop
invisible(file.create(command_file))
for (i in 1:nb_of_bootstrap_group){

  idx_boot = trainmat.outer_list[[i]]
  # save data for the fold
  bootstrap_sample_parameters_file = paste(pathtowritefile, "bootstrap_parameters_files/bootstrap_sample_parameters_file_", i, ".Rdata", sep="")
  save(full_data_run, 
       idx_boot, 
       B_3D,
       nb_cores,
       bias, 
       tau, 
       cv_method, 
       scale, 
       scale_size_bloc, 
       ncomp, 
       nfac,
       parfor_bool,
       init, 
       nstart,
       tol,
       M_regularisation,
       data_path, 
       file = bootstrap_sample_parameters_file)
  
  writeMat(paste(file_path_sans_ext(bootstrap_sample_parameters_file), "mat", sep = "."),
           idx_boot        = idx_boot, 
           nb_cores        = nb_cores,
           nfac            = nfac,
           parfor_bool     = parfor_bool,
           init            = init, 
           nstart          = nstart,
           tol             = tol,
           data_path       = paste(file_path_sans_ext(data_path), "mat", sep = "."),
           full_data_file  = full_data_file)
           
  
  ### tau tune for the fold i 
  write(sapply(mappers, 
                   function(x) write_mapper_command(model.mapper_name    = x, 
                                                    model.mapper_script  = method.cfg$model.mapper[[x]],
                                                    model.mapper_options = method.cfg$option.mapper[[x]],
                                                    i                    = i, 
                                                    pathtowritefile      = pathtowritefile,
                                                    pathtoscript         = pathtoscript, 
                                                    parameter_file       = bootstrap_sample_parameters_file)),
        file = command_file, append = T)
}


reducers_parameter_file = paste(pathtowritefile, "reducers_parameter_file", ".Rdata", sep="")
save(full_data_run, 
     B_3D,
     to_source_in_reducer, 
     alpha, 
     pathtoscript,
     data_path, 
     channels_file,
     pathtofigurefile,
     channel_figure_pathtoscript,
     file = reducers_parameter_file)

## Reducer call
write(sapply(mappers, 
                 function(x) write_reducer_command(model.mapper_name    = x, 
                                                   model.mapper_options = method.cfg$option.mapper[[x]],
                                                   model.reducer_script = method.cfg$model.reducer[[x]],
                                                   nf                   = bootsrap_samples, 
                                                   pathtowritefile      = pathtowritefile,
                                                   pathtoscript         = pathtoscript, 
                                                   parameter_file       = reducers_parameter_file)), 
      file = command_file, append = T)

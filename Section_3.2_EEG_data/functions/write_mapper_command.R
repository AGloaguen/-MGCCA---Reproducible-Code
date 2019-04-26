

write_mapper_command = function(model.mapper_name, model.mapper_script, model.mapper_options, i, pathtowritefile, pathtoscript, parameter_file){

  switch(model.mapper_name,
         
         ###############################
         #             MGCCA           #
         ###############################
         "mgcca" = 
         {
           design = model.mapper_options$design
           scheme = model.mapper_options$scheme
           for (j in 1:length(design)){
             output_file = paste(pathtowritefile, 
                                 "map_res_mgcca_"       , design[j], 
                                 "_scheme_"             , scheme,
                                 "_fold_"                , i, 
                                 ".Rdata", sep="")
             return(paste("Rscript", 
                       paste(pathtoscript, model.mapper_script, sep=""), 
                       parameter_file, 
                       output_file, 
                       i, 
                       pathtoscript, 
                       design[j], 
                       scheme, 
                       " ",
                       sep = " "))
           }},
         
         ###############################
         #             RGCCA           #
         ###############################
         "rgcca" = 
         {
           design    = model.mapper_options$design
           HP        = model.mapper_options$hyperparameter_parallelized
           scheme    = model.mapper_options$scheme
           unfolding = model.mapper_options$unfolding
           for (j in 1:length(design)){
             for (k in 1:length(unfolding)){
               if(HP){
                 for (hyperparameter in 1:nb_hyperparameters){
                   output_file = paste(pathtowritefile, "map_res_rgcca_", design[j], "_scheme_", scheme,
                                       "_unfolding_", unfolding[k],"_hyperparameter_", hyperparameter, 
                                       "_fold_", i, ".Rdata", sep="")
                   cat(paste("Rscript", paste(pathtoscript, "mapper.rgcca.R", sep=""), parameter_file, output_file, 
                             hyperparameter, i, pathtoscript, design[j], scheme, unfolding[k]), "\n", sep = " ")
                 }
               }else{
                 output_file = paste(pathtowritefile, "map_res_rgcca_", design[j], "_scheme_", scheme,
                                     "_unfolding_", unfolding[k],"_ALL_hyperparameters_", "fold_", i, ".Rdata", sep="")
                 cat(paste("Rscript", paste(pathtoscript, "mapper.rgcca.R", sep=""), parameter_file, output_file, "ALL", 
                           i, pathtoscript, design[j], scheme, unfolding[k]), "\n", sep = " ")
               } 
             }
           }},
         
         ###############################
         #   RGCCA    Boostrap         #
         ###############################   
         "rgcca_boot" = 
         {
           design    = model.mapper_options$design
           scheme    = model.mapper_options$scheme
           for (j in 1:length(design)){
               output_file = paste(pathtowritefile, "map_res_rgcca_boot_", design[j], "_scheme_", scheme,
                                   "_bootstrapSample_", i, ".Rdata", sep="")
               return(paste(paste("Rscript", paste(pathtoscript, "mapper.rgcca.bootstrap.R", sep=""), parameter_file, output_file, 
                         i, pathtoscript, design[j], scheme), "\n", sep = " "))
           }},
         
         ###############################
         #             SGCCA           #
         ###############################
         "sgcca" = 
         {
           design    = model.mapper_options$design
           HP        = model.mapper_options$hyperparameter_parallelized
           scheme    = model.mapper_options$scheme
           for (j in 1:length(design)){
             if(HP){
               for (hyperparameter in 1:nb_hyperparameters){
                 output_file = paste(pathtowritefile, "map_res_sgcca_", design[j], "_scheme_", scheme,
                                     "_hyperparameter_", hyperparameter, "_fold_", i, ".Rdata", sep="")
                 cat(paste("Rscript", paste(pathtoscript, "mapper.sgcca.R", sep=""), parameter_file, output_file, 
                           hyperparameter, i, pathtoscript, design[j], scheme), "\n", sep = " ")
               }
             }else{
               output_file = paste(pathtowritefile, "map_res_sgcca_", design[j], "_scheme_", scheme,
                                   "_ALL_hyperparameters_", "fold_", i, ".Rdata", sep="")
               cat(paste("Rscript", paste(pathtoscript, "mapper.sgcca.R", sep=""), parameter_file, output_file, "ALL", 
                         i, pathtoscript, design[j], scheme), "\n", sep = " ")
             } 
           }},
         
         ###############################
         #           Parafac           #
         ###############################
         "parafac" = 
         {
           output_file = paste(pathtowritefile, 
                               "map_res_parafac_", 
                               "fold_"          , i, 
                               ".Rdata", sep="")
           return(paste("Rscript", 
                        paste(pathtoscript, model.mapper_script, sep=""), 
                        parameter_file, 
                        output_file, 
                        i, 
                        pathtoscript,
                        " ",
                        sep = " "))
         },
         
         ###############################
         #            CMTF             #
         ###############################
         "cmtf" = 
         {
           output_file = paste(pathtowritefile, 
                               "map_res_cmtf_", 
                               "fold_" , i, 
                               ".mat", sep="")
           cmd = paste("\"clear all; close all; parameters_file = ",
                       paste0("\'", file_path_sans_ext(parameter_file), ".mat", "\'"), "; output_file = ",
                       paste0("\'", output_file, "\'"), "; file_number = ", i, ";",
                       paste0("run(\'", pathtoscript, model.mapper_script, "\');exit;\""), " ",
                       sep = "")
          return(paste("matlab", "-nodisplay", "-nosplash", "-nodesktop", "-r", cmd))
         },
         
         ###############################
         #           NPLS              #
         ###############################
         "npls" = 
         {
           HP          = model.mapper_options$hyperparameter_parallelized
           if(HP){
             for (hyperparameter in 1:nb_hyperparameters){
               output_file = paste(pathtowritefile, "map_res_npls_", "hyperparameter_", hyperparameter, 
                                   "_fold_", i, ".Rdata", sep="")
               cat(paste("Rscript", paste(pathtoscript, "mapper.npls.R", sep=""), parameter_file, output_file, 
                         hyperparameter, i, pathtoscript), "\n", sep = " ")
             }
           }else{
             output_file = paste(pathtowritefile, "map_res_npls_", "ALL_hyperparameters_", "fold_", i, ".Rdata", sep="")
             cat(paste("Rscript", paste(pathtoscript, "mapper.npls.R", sep=""), parameter_file, output_file, "ALL", 
                       i, pathtoscript), "\n", sep = " ")
           }
         }
  ) 
}

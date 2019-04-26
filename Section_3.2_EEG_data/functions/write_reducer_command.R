

write_reducer_command = function(model.mapper_name, model.reducer_script, model.mapper_options, 
                                 nf, pathtowritefile, pathtoscript, parameter_file){
  
  switch(model.mapper_name,
         
         ###############################
         #             MGCCA           #
         ###############################
         "mgcca" = 
         {
           design = model.mapper_options$design
           scheme = model.mapper_options$scheme
           for (j in 1:length(design)){
             map_res_pattern = paste("map_res_mgcca", design[j], 
                                     "scheme"       , scheme, 
                                     sep = "_")
             red_res_file    = paste0(pathtowritefile, 
                                      "red_res_mgcca_", design[j], 
                                      "_scheme_"      , scheme, 
                                      ".Rdata") 
             return(paste("Rscript", 
                       paste(pathtoscript, model.reducer_script, sep=""), 
                       pathtowritefile, 
                       map_res_pattern, 
                       nf, 
                       red_res_file, 
                       design[j],
                       parameter_file,
                       " ",
                       sep=" "))
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
               map_res_pattern = paste("map_res_rgcca", design[j], "scheme", scheme, "unfolding", unfolding[k], sep = "_")
               red_res_file    = paste0(pathtowritefile, "red_res_rgcca_", design[j], "_scheme_", scheme, "_unfolding_", unfolding[k], ".Rdata") 
               cat(paste("Rscript", paste(pathtoscript, "reducer.rgcca.R", sep=""), pathtowritefile, map_res_pattern, nb_hyperparameters, nf, 
                         red_res_file, PLOT, design[j], unfolding[k], HP, "\n",sep=" "))
             }
           }},
         
         ###############################
         #     RGCCA    bootstrap      #
         ###############################
         "rgcca_boot" = 
         {
           design    = model.mapper_options$design
           scheme    = model.mapper_options$scheme
           for (j in 1:length(design)){
             map_res_pattern = paste("map_res_rgcca_boot", design[j], "scheme", scheme, sep = "_")
             red_res_file    = paste0(pathtowritefile, "red_res_rgcca_boot_", design[j], "_scheme_", scheme, ".Rdata") 
             return(paste("Rscript", paste(pathtoscript, "reducer.rgcca.bootstrap.R", sep=""), pathtowritefile, map_res_pattern, 
                       nf, red_res_file, design[j], "\n",sep=" "))
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
             map_res_pattern = paste("map_res_sgcca", design[j], "scheme", scheme, sep = "_")
             red_res_file    = paste0(pathtowritefile, "red_res_sgcca_", design[j], "_scheme_", scheme, ".Rdata") 
             cat(paste("Rscript", paste(pathtoscript, "reducer.sgcca.R", sep=""), pathtowritefile, map_res_pattern, nb_hyperparameters, nf, 
                       red_res_file, design[j], HP, "\n",sep=" "))
           }},
         
         ###############################
         #           Parafac           #
         ###############################
         "parafac" = 
         {
           map_res_pattern = "map_res_parafac_"
           red_res_file    = paste0(pathtowritefile, 
                                    "red_res_parafac.Rdata") 
           return(paste("Rscript", 
                        paste(pathtoscript, model.reducer_script, sep=""), 
                        pathtowritefile, 
                        map_res_pattern, 
                        nf, 
                        red_res_file, 
                        parameter_file,
                        " ",
                        sep=" "))
         },
         
         ###############################
         #             CMTF            #
         ###############################
         "cmtf" = 
         {
           map_res_pattern = "map_res_cmtf_"
           red_res_file    = paste0(pathtowritefile, 
                                    "red_res_cmtf.Rdata") 
           return(paste("Rscript", 
                        paste(pathtoscript, model.reducer_script, sep=""), 
                        pathtowritefile, 
                        map_res_pattern, 
                        nf, 
                        red_res_file, 
                        parameter_file,
                        " ",
                        sep=" "))
         },
         
         ###############################
         #           NPLS              #
         ###############################
         "npls" = 
         {
           map_res_pattern = "map_res_npls_";
           red_res_file    = paste(pathtowritefile, "red_res_npls.Rdata" , sep="") 
           cat(paste("Rscript", paste(pathtoscript, "reducer.npls.R", sep=""), pathtowritefile, map_res_pattern, nb_hyperparameters, nf, 
                     red_res_file, PLOT, "\n",sep=" "))
         }
  ) 
}

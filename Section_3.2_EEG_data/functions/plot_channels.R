library(R.matlab)


channels_file_new       = paste(dirname(channels_file), "channel_file_tmp.xyz", sep = "/")
plotrad                 = 0.5
intrad                  = 0.5
show_electrodes         = "on"
numcontour              = 1
color_map               = "gray"
cluster_channels_symbol = "o"
cluster_channels_color  = "k"
cluster_channels_size   = 4

nfac = dim(Cstar[[1]])[3]

for (fac in 1:nfac){
  elec_kept                     = which(apply(Cstar[[1]], c(2, 3), function(x) prod(x[2:3])) > 0)
  coef_channels                 = Cstar[[1]][1, , fac]
  cluster_channels              = elec_kept

  writeMat(paste(channel_figure_parameter_file, "_comp_", fac, ".mat", sep = ""),
           coef_channels              = coef_channels,
           channels_file              = channels_file,
           channels_file_new          = channels_file_new,
           plotrad                    = plotrad,
           intrad                     = intrad,
           show_electrodes            = show_electrodes,
           numcontour                 = numcontour,
           color_map                  = color_map,
           cluster_channels           = cluster_channels,
           cluster_channels_symbol    = cluster_channels_symbol,
           cluster_channels_color     = cluster_channels_color,
           cluster_channels_size      = cluster_channels_size,
           channels_to_remove         = c(17, 125, 126, 127, 128),
           output_channel_figure_file = paste(output_channel_figure_file, "_comp_", fac, ".jpg", sep = ""))
  
  cmd = paste("\"load(\'",
              paste(channel_figure_parameter_file, "_comp_", fac, ".mat", sep = ""),
              "\');run(\'",
              channel_figure_pathtoscript,
              "\');exit;\"",
              sep = "")
  system(paste("matlab", "-nodisplay", "-nosplash", "-nodesktop", "-r", cmd))
}

library(ggplot2)

nfac = dim(Bstar[[1]])[3]

for (fac in 1:nfac){
  ft_size             = 20
  Mean                = Bstar[[1]][1, , fac]
  Q25                 = Bstar[[1]][2, , fac]
  Q975                = Bstar[[1]][3, , fac]
  print(length(Mean))
  print(length(Q25))
  print(length(Q975))
  print(length(time))
  df                  = data.frame(x = time, Mean = Mean, Q25 = Q25, Q975 = Q975)
  
  B_Q_sig      = abs(diff(c(-1, sign(df$Q25*df$Q975))))/2
  B_Q_sig      = which(B_Q_sig != 0)
  nb_breaks_qt = length(B_Q_sig)
  if (nb_breaks_qt !=0){
    if((nb_breaks_qt %% 2)){B_Q_sig = c(B_Q_sig, length(time))}
    B_Q_sig = matrix(time[B_Q_sig], byrow = T, ncol = 2)
    B_Q_sig = data.frame(x_min = B_Q_sig[, 1], x_max = B_Q_sig[, 2])
    
    fig_1 = ggplot() +
      geom_rect(data = B_Q_sig, aes(xmin = x_min, xmax = x_max, ymin = -Inf, ymax = Inf, colour = "gray"), 
                alpha = 0.4, col = "gray") +
      geom_line(data=df, aes(x=time, y=Mean)) +
      geom_line(data=df, aes(x=time, y=Q25), lty = "dashed") +
      geom_line(data=df, aes(x=time, y=Q975), lty = "dashed") +
      theme_bw() +
      xlab("Time (s)") +
      theme(axis.title.x=element_text(size=25),
            axis.title.y=element_text(size=25),
            axis.text.x = element_text(size = 25),
            axis.text.y = element_text(size = 25)) + 
      geom_hline(yintercept = 0) + geom_vline(xintercept = c(0, -1.8, -1.2, -0.6)) +
      scale_x_continuous(breaks = round(seq(min(df$x), max(df$x), by = 0.5),1))
      print(B_Q_sig)
  }else{
    fig_1 = ggplot(data=df, aes(x=time, y=Mean, ymin=Q25, ymax=Q975)) +
      geom_line() +
      geom_ribbon(alpha=0.3)+
      xlab("Time (s)") +
      theme_bw() +
      theme(axis.title.x=element_text(size=25),
            axis.title.y=element_text(size=25),
            axis.text.x = element_text(size = 25),
            axis.text.y = element_text(size = 25)) + 
      geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
      scale_x_continuous(breaks = round(seq(min(df$x), max(df$x), by = 0.5),1))
  }
  
  ggsave(filename = paste(output_time_figure_file, "_comp_", fac, ".jpg", sep = ""), 
         plot     = fig_1, 
         device   = "jpeg", 
         height   = 7, 
         width    = 21)
}



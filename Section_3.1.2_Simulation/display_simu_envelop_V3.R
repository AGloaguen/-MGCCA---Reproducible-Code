display_simu_envelop_V3 = function(J, K, TRUTH, PARAFAC, MGCCA, CMTF, ft_size, width = 29.7, height = 7, figName){
  
  color_legends = c("Truth"="black", "Parafac"="black", "MGCCA"="black", "CMTF" = "black")
  
  #############################################################
  ###################   MODE 2 : B coeffs   ###################
  #############################################################
  
  
  #########################
  ####     TENSOR 1    ####
  #########################
  
  x  <- 1:J[1]
  y1 = TRUTH$b[[1]]
  Y2 = PARAFAC$b[[1]]
  Y3 = MGCCA$b[[1]]
  Y4 = CMTF$b[[1]]
  df <- data.frame(x, y1 = y1, y2 = Y2[1, ], y3 = Y3[1, ], y4 = Y4[1, ], 
                   y2_min = Y2[3, ], y3_min = Y3[3, ], y4_min = Y4[3, ],
                   y2_max = Y2[4, ], y3_max = Y3[4, ], y4_max = Y4[4, ])
  
  B_T1_PF <-  ggplot(df, aes(x)) + theme_bw() +
    geom_line(aes(y=y1), colour="black", size = 1) +
    geom_point(aes(y=y2), colour="black", size = 2, shape = 5) +
    geom_ribbon(data=df, aes(x=x, ymin=y2_min, ymax=y2_max),fill="black", alpha=0.2) +
    xlab("") + ylab("") +
    scale_colour_manual(name="Methods",values=color_legends ) + theme(legend.position="bottom") +
    scale_y_continuous(limits = c(-0.6, 0.6), breaks = c(-0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6)) +
    theme(axis.text=element_text(size=ft_size),
          axis.title.x=element_text(size=ft_size), 
          title = element_text(size=ft_size), 
          plot.title = element_text(hjust = 0.5),
          axis.title.y = element_text(size=ft_size, angle = 0, vjust = 0.5))
  
  B_T1_CMTF <-  ggplot(df, aes(x)) + theme_bw() +
    geom_line(aes(y=y1), colour="black", size = 1) +
    geom_point(aes(y=y4), colour="black", size = 2, shape = 5) +
    geom_ribbon(data=df, aes(x=x, ymin=y4_min, ymax=y4_max),fill="black", alpha=0.2) +
    xlab("") + ylab("") +
    scale_colour_manual(name="Methods",values=color_legends ) + theme(legend.position="bottom") +
    scale_y_continuous(limits = c(-0.6, 0.6), breaks = c(-0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6)) +
    theme(axis.text=element_text(size=ft_size),
          title = element_text(size=ft_size),
          axis.title=element_text(size=ft_size),
          plot.title = element_text(hjust = 0.5))
  
  B_T1_MGCCA <-  ggplot(df, aes(x)) + theme_bw() +
    geom_line(aes(y=y1), colour="black", size = 1) +
    geom_point(aes(y=y3), colour="black", size = 2, shape = 5) +
    geom_ribbon(data=df, aes(x=x, ymin=y3_min, ymax=y3_max),fill="black", alpha=0.2) +
    xlab("") + ylab("") +
    scale_colour_manual(name="Methods",values=color_legends ) + theme(legend.position="bottom") +
    scale_y_continuous(limits = c(-0.6, 0.6), breaks = c(-0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6)) +
    theme(axis.text=element_text(size=ft_size),
          title = element_text(size=ft_size),
          axis.title=element_text(size=ft_size),
          plot.title = element_text(hjust = 0.5))
  
  
  #############################################################
  ###################         Truth         ###################
  #############################################################
  plot <- grid.arrange(arrangeGrob(B_T1_PF),  
                       arrangeGrob(B_T1_CMTF), 
                       arrangeGrob(B_T1_MGCCA),
                       nrow=1,heights = 10, widths = c(10, 10, 10))
  
  ggsave(filename = figName, plot = plot, width = width, height = height, units ="cm")
}

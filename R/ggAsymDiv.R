ggAsymDiv <- function(output){
  cbPalette <- rev(c("#999999", "#E69F00", "#56B4E9", "#009E73", 
                     "#330066", "#CC79A7", "#0072B2", "#D55E00"))
  ggplot(output, aes(x=order, y=qD, colour=Site, lty = method)) +
    geom_line(size=1.2) +
    scale_colour_manual(values = cbPalette) +
    geom_ribbon(data = output %>% filter(method=="Estimated"),
                aes(ymin=qD.LCL, ymax=qD.UCL, fill=Site), alpha=0.2, linetype=0) +
    scale_fill_manual(values = cbPalette) +
    labs(x="Order q", y="Species diversity") +
    theme_bw(base_size = 18) +
    theme(text=element_text(size=18)) +
    theme(legend.position="bottom", legend.box = "vertical", 
          legend.key.width = unit(1.2,"cm"), 
          # plot.margin = unit(c(1.5,0.3,1.2,0.3), "lines"),
          legend.title=element_blank(),
          legend.margin=margin(0,0,0,0),
          legend.box.margin = margin(-10,-10,-5,-10)) +
    scale_linetype_manual(values = c(2,1), breaks=c("Estimated", "Empirical"),
                          labels=c("Asymptotic", "Empirical"))
}

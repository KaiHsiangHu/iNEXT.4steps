Evenness <- function(x, q=seq(0,2,0.2), datatype = "abundance") {
  name = colnames(x)
  est = estimateD(x, q=seq(0,2,0.1), datatype, base="coverage", level=NULL, nboot=0)
  maxC = min(unique(est[,"SC"]))
  even = cbind(est[,c("site","order")],
               Evenness=as.numeric(sapply(name, function(i) {tmp=(est %>% filter(Community==i))$qD; tmp/tmp[1]}))
  )
  return(even)
}

ggEven <- function(output) {
  cbPalette <- rev(c("#999999", "#E69F00", "#56B4E9", "#009E73", 
                     "#330066", "#CC79A7", "#0072B2", "#D55E00"))
  ggplot(output, aes(x=order, y=Evenness, colour=site)) +
    geom_line(size=1.2) +
    scale_colour_manual(values = cbPalette) +
    labs(x="Order q", y="Evenness") +
    theme_bw(base_size = 18) +
    theme(text=element_text(size=18)) +
    theme(legend.position="bottom", legend.box = "vertical", 
          legend.key.width = unit(1.2,"cm"), 
          # plot.margin = unit(c(1.5,0.3,1.2,0.3), "lines"),
          legend.title=element_blank())
}

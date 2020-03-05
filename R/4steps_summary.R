summary.deal <- function(table, step, Pielou=NULL) {
  if (step==1) {
    tmp = (table %>% filter(order %in% c(0,1,2)))[,c("order","qD","Site")]
    out = acast(tmp, Site~order, value.var="qD")
    colnames(out) = paste("q=", c(0,1,2), sep="")
  }
  if (step==2){
    tmp = (table %>% filter(order %in% c(0,1,2)))[,c("order","qD","Site","method")]
    out = sapply(unique(tmp$Site), function(k) {
      tmp = tmp %>% filter(Site==k)
      tmp2 = acast(tmp,  method~order, value.var="qD")
      colnames(tmp2) = paste("q=", c(0,1,2), sep="")
      tmp2 = rbind(tmp2, "Undetected"=tmp2[2,]-tmp2[1,])
    }, simplify = "array")
    dimnames(out)[[3]] = unique(tmp$Site)
  }
  if (step==3){
    tmp = (table %>% filter(order %in% c(0,1,2)))[,c("order","qD","site")]
    Cmax = round(min(table$SC), 3)
    out = dcast(tmp, site~order, value.var="qD")
    colnames(out) = c(paste("'maxC=", Cmax, "'", sep=""),  
                      paste("q=", c(0,1,2), sep=""))
  }
  if (step==4){
    tmp = table %>% filter(order %in% c(0,1,2))
    out = acast(tmp, site~order, value.var="Evenness")
    
    logD = (Pielou %>% filter(order == 1))[,c("site","qD")]
    logS = (Pielou %>% filter(order == 0))[,c("site","qD")]
    out[,1] = log(logD$qD)/log(logS$qD)
    colnames(out) = c("Pielou J'", paste("q=", c(1,2), sep=""))
  }
  
  return(out)
}

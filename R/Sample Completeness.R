SC <- function (x, q = seq(0, 2, 0.2), datatype = "abundance", nboot = 50, 
                conf = 0.95)
{
  TYPE <- c("abundance", "incidence", "incidence_freq", "incidence_raw")
  if (is.na(pmatch(datatype, TYPE))) 
    stop("invalid datatype")
  if (pmatch(datatype, TYPE) == -1) 
    stop("ambiguous datatype")
  datatype <- match.arg(datatype, TYPE)
  class_x <- class(x)[1]
  if (datatype == "incidence") {
    stop("datatype=\"incidence\" was no longer supported after v2.0.8, \n         please try datatype=\"incidence_freq\".")
  }
  if (datatype == "incidence_raw") {
    if (class_x == "list") {
      x <- lapply(x, as.incfreq)
    }
    else {
      x <- as.incfreq(x)
    }
    datatype <- "incidence"
  }
  if (datatype == "incidence_freq") 
    datatype <- "incidence"
  if (class(x) == "numeric" | class(x) == "integer") {
    x <- list(data = x)
  }
  if (class(x) == "data.frame" | class(x) == "matrix") {
    datalist <- lapply(1:ncol(x), function(i) x[, i])
    if (is.null(colnames(x))) 
      names(datalist) <- paste0("data", 1:ncol(x))
    else names(datalist) <- colnames(x)
    x <- datalist
  }
  if (datatype == "abundance") {
    out <- lapply(1:length(x), function(i) {
      dq <- sample_completeness(x[[i]], q, "abundance")
      if (nboot > 1) {
        Prob.hat <- iNEXT:::EstiBootComm.Ind(x[[i]])
        Abun.Mat <- rmultinom(nboot, sum(x[[i]]), Prob.hat)
        error <- qnorm(1-(1-conf)/2)*
          apply(apply(Abun.Mat,  2, function(xb) sample_completeness(xb, q, "abundance")), 
                1, sd, na.rm = TRUE)
      }
      else {
        error = 0
      }
      out <- data.frame(order = q, qD = dq, qD.LCL = dq-error, qD.UCL = dq+error, Site = names(x)[i], 
                        method = rep("Estimated", length(q)))
      out$qD.LCL[out$qD.LCL < 0] <- 0
      out$qD.UCL[out$qD.UCL > 1] <- 1
      out
    })
    out <- do.call(rbind, out)
  }
  else if (datatype == "incidence") {
    out <- lapply(1:length(x), function(i) {
      dq <- sample_completeness(x[[i]], q, "incidence_freq")
      if (nboot > 1) {
        nT <- x[[i]][1]
        Prob.hat <- iNEXT:::EstiBootComm.Sam(x[[i]])
        Incid.Mat <- t(sapply(Prob.hat, function(p) rbinom(nboot, nT, p)))
        Incid.Mat <- matrix(c(rbind(nT, Incid.Mat)), ncol = nboot)
        tmp <- which(colSums(Incid.Mat) == nT)
        if (length(tmp) > 0) 
          Incid.Mat <- Incid.Mat[, -tmp]
        if (ncol(Incid.Mat) == 0) {
          error = 0
          warning("Insufficient data to compute bootstrap s.e.")
        }
        else {
          error <- qnorm(1-(1-conf)/2)*
            apply( apply(Incid.Mat, 2, function(yb) sample_completeness(yb, q, "incidence_freq")), 
                   1, sd, na.rm = TRUE)
        }
      }
      else {
        error = 0
      }
      out <- data.frame(order = q, qD = dq, qD.LCL = dq-error, qD.UCL = dq+error, Site = names(x)[i], 
                        method = rep("Estimated", length(q)))
      out$qD.LCL[out$qD.LCL < 0] <- 0
      out$qD.UCL[out$qD.UCL > 1] <- 1
      out
    })
    out <- do.call(rbind, out)
  }
}

ggSC <- function(output) {
  cbPalette <- rev(c("#999999", "#E69F00", "#56B4E9", "#009E73", 
                     "#330066", "#CC79A7", "#0072B2", "#D55E00"))
  ggplot(output, aes(x=order, y=qD, colour=Site))+
    geom_line(size=1.2) +
    scale_colour_manual(values = cbPalette) + 
    geom_ribbon(aes(ymin=qD.LCL, ymax=qD.UCL, fill=Site, colour=NULL), alpha=0.2) + 
    scale_fill_manual(values = cbPalette) +
    labs(x="Order q", y="Sample completeness") + 
    theme_bw(base_size = 18) +
    theme(text=element_text(size=18)) +
    theme(legend.position="bottom", legend.box = "vertical", 
          legend.key.width = unit(1.2,"cm"), 
          # plot.margin = unit(c(1.5,0.3,1.2,0.3), "lines"),
          legend.title=element_blank())
}

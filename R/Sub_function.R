####
# Generate the summary table from the deatils of each function.
#
# \code{summary.deal} Generate the summary table of four lists from the deatils of each function.
#
# @param table a table
# @param step a numerical value of which step summary do you need
# @param Pielou a logical variable whether you want to generate step 4 summary
# @return a matrix or a list according different steps

summary.deal <- function(table, step, Pielou=NULL) {
  if (step==1) {
    tmp = (table %>% filter(order %in% c(0,1,2)))[,c("order","SC","Site")]
    out = acast(tmp, Site~order, value.var="SC")
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
    tmp = table[,c("order","qD","site")]
    Cmax = round(min(table$SC), 3)
    out = dcast(tmp, site~order, value.var="qD")
    colnames(out) = c(paste("'maxC=", Cmax, "'", sep=""),
                      paste("q=", c(0,1,2), sep=""))
  }
  if (step==4){
    tmp = (table[[1]] %>%
             filter( (order %in% c(0,1,2)) & method %in% "Estimated") )[,c("order","Evenness","Site")]
    out = acast(tmp, Site~order, value.var="Evenness")

    D = (Pielou %>% filter(order == 1))[,c("site","qD")]
    S = (Pielou %>% filter(order == 0))[,c("site","qD")]
    out[,1] = sapply(rownames(out), function(x) log(D[D$site==x,"qD"])/log(S[S$site==x,"qD"]))
    colnames(out) = c("Pielou J'", paste("q=", c(1,2), sep=""))
  }

  return(out)
}

#
####
# Sample Completeness main function
#
# \code{SC} Estimation of Sample Completeness with order q
#
# @param x data.frame or matrix of data
# @param q a integer vector of the order of Hill number
# @param datatype a binary choose with 'abundance' or 'incidence_freq'
# @param nboot the number of bootstrap resampling times, default is 50
# @param conf a integer value between 0 to 1 for confidence interval
# @return a matrix of estimated sample completeness with order q

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
      out <- data.frame(order = q, SC = dq, SC.LCL = dq-error, SC.UCL = dq+error, Site = names(x)[i],
                        method = rep("Estimated", length(q)))
      out$SC.LCL[out$SC.LCL < 0] <- 0
      out$SC.UCL[out$SC.UCL > 1] <- 1
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
      out <- data.frame(order = q, SC = dq, SC.LCL = dq-error, SC.UCL = dq+error, Site = names(x)[i],
                        method = rep("Estimated", length(q)))
      out$SC.LCL[out$SC.LCL < 0] <- 0
      out$SC.UCL[out$SC.UCL > 1] <- 1
      out
    })
    out <- do.call(rbind, out)
  }
}

#
####
# Sample Completeness
#
# \code{sample_completeness} Estimation of Sample Completeness with order q
#
# @param x a vector of data
# @param q a integer vector of the order of Hill number
# @param datatype a binary choose with 'abundance' or 'incidence_freq'
# @return a vector of estimated sample completeness with order q

sample_completeness = function(x, q, datatype = c("abundance","incidence_freq")){
  if(datatype=="abundance"){
    x = x[x>0]
    n = sum(x)
    f1 = sum(x==1); f2 = sum(x==2)
    A = ifelse(f2>0,2*f2/((n-1)*f1+2*f2),ifelse(f1>0,2/((n-1)*(f1-1)+2),1))

    sc.abund = function(q){
      if (q==0){
        S_obs = length(x)
        f0_hat = ifelse(f2==0, ((n-1)/n)*(f1*(f1-1)/2), ((n-1)/n)*((f1^2)/(2*f2)))
        c_hat = S_obs/(S_obs+ceiling(f0_hat))
        return(c_hat)
      } else if (q==1){
        c_hat = 1-(f1/n)*(1-A)
        return(c_hat)
      } else if (q==2){
        x = x[x>=2]
        c_hat = 1-(f1/n)*((A*(1-A))/sum(x*(x-1)/(n*(n-1))))
        return(c_hat)
      } else {
        r <- 0:(n-1)
        sort.data = sort(unique(x))
        tab = table(x)
        term = sapply(sort.data,function(z){
          k=0:(n-z)
          sum(choose(k-q,k)*exp(lchoose(n-k-1,z-1)-lchoose(n,z)))
        })
        lambda_hat = sum(tab*term) + ifelse(f1==0|A==1, 0, f1/n*(1-A)^(1-n)*(A^(q-1)-sum(choose(q-1,r)*(A-1)^r)))
        c_hat = 1-((f1/n)*(A^(q-1))*(1-A)/lambda_hat)
        return(c_hat)
      }
    }
    est=sapply(q, sc.abund)
  } else {
    t = x[1]
    x = x[-1]; x = x[x>0]
    u = sum(x)
    Q1 = sum(x==1); Q2 = sum(x==2)
    B = ifelse(Q2>0,2*Q2/((t-1)*Q1+2*Q2), ifelse(Q1>0, 2/((t-1)*(Q1-1)+2), 1))

    sc.incid = function(q){
      if (q==0){
        S_obs = length(x)
        f0_hat = ifelse(Q2==0, ((t-1)/t)*(Q1*(Q1-1)/2), ((t-1)/t)*((Q1^2)/(2*Q2)))
        c_hat = S_obs/(S_obs+f0_hat)
        return(c_hat)
      } else if (q==1){
        c_hat = 1-(Q1/u)*(1-B)
        return(c_hat)
      } else if (q==2){
        x = x[x>=2]
        c_hat = 1-(t-1)*Q1*((B*(1-B))/sum(x*(x-1)))
        return(c_hat)
      } else {
        r <- 0:(t-1)
        sort.data = sort(unique(x))
        tab = table(x)
        term = sapply(sort.data,function(z){
          k=0:(t-z)
          sum(choose(k-q,k)*exp(lchoose(t-k-1,z-1)-lchoose(t,z)))
        })
        phi_hat = sum(tab*term) + ifelse(Q1==0|B==1, 0, Q1/t*(1-B)^(1-t)*(B^(q-1)-sum(choose(q-1,r)*(B-1)^r)))
        c_hat = 1-((Q1/t)*(B^(q-1))*(1-B)/phi_hat)
        return(c_hat)
      }
    }
    est=sapply(q, sc.incid)
  }
  return(est)
}

#
####
# ggplot for Sample Completeness
#
# \code{SC} The figure for estimation of Sample Completeness with order q
#
# @param output a table generated from SC function
# @return a figure of estimated sample completeness with order q

ggSC <- function(output) {
  cbPalette <- rev(c("#999999", "#E69F00", "#56B4E9", "#009E73",
                     "#330066", "#CC79A7", "#0072B2", "#D55E00"))
  ggplot(output, aes(x=order, y=SC, colour=Site))+
    geom_line(size=1.2) +
    scale_colour_manual(values = cbPalette) +
    geom_ribbon(aes(ymin=SC.LCL, ymax=SC.UCL, fill=Site, colour=NULL), alpha=0.2) +
    scale_fill_manual(values = cbPalette) +
    labs(x="Order q", y="Sample completeness") +
    theme_bw(base_size = 18) +
    theme(text=element_text(size=18)) +
    theme(legend.position="bottom", legend.box = "vertical",
          legend.key.width = unit(1.2,"cm"),
          # plot.margin = unit(c(1.5,0.3,1.2,0.3), "lines"),
          legend.title=element_blank())
}

#
####
# ggplot for Asymptotic diversity
#
# \code{ggAsymDiv} The figure for estimation of Asymptotic diversity with order q
#
# @param output a table generated from AsymDiv function
# @return a figure of estimated sample completeness with order q

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

#
####
# Calculate six classes for Evenness
#
# @param q a integer vector for the order of Hill number
# @param qD a vector for diversity
# @param S a integer value for species
# @param E a integer value between 1 to 6
# @return a vector for evenness value

even = function(q, qD, S, E) {
  tmp = c()
  if (E == 1)
    tmp = ifelse(q!=1, (1-qD^(1-q))/(1-S^(1-q)), log(qD)/log(S))
  if (E == 2)
    tmp = ifelse(q!=1, (1-qD^(q-1))/(1-S^(q-1)), log(qD)/log(S))
  if (E == 3)
    tmp = (qD-1)/(S-1)
  if (E == 4)
    tmp = (1-1/qD)/(1-1/S)
  if (E == 5)
    tmp = log(qD)/log(S)
  if (E == 6) {
    p = x/sum(x)
    tmp = sapply(q, function(k) {
      if(k==0){
        nu <- abs(p - (1/S))
        nu <- nu[nu > 0]
        sub1 <- (sum(log(abs(nu)))/sum(nu>0)-(log(1-1/S)+(1-S)*log(S))/S)
        1-exp(sub1)
      }else{
        1-(sum(abs(p-1/S)^k)/((1-1/S)^k+(S-1)*S^(-k)))^(1/k)
      }
    })
    tmp
  }
  return(tmp)
}

#
####
# Evenness profile
#
# \code{Evenness.profile} Estimation or Empirical of Evenness with order q
#
# @param x vector for data.
# @param q a integer vector for the order of Hill number.
# @param datatype a binary choose with 'abundance' or 'incidence_freq'
# @param method the number of bootstrap resampling times, default is 50
# @param E.type a integer vector between 1 to 6
# @return a list of estimated(empirical) evenness with order q, each list is combined with a matrix

Evenness.profile <- function(x, q, datatype=c("abundance","incidence_freq"),
                             method=c("Estimated", "Empirical"), E.type) {
  x = x[x>0]
  if (method == "Estimated") {
    estqD = estimateD(x, q, datatype, base="coverage", level=NULL, nboot=0)
    estS = estimateD(x, 0, datatype, base="coverage", level=NULL, nboot=0)

    out = lapply(E.type, function(i) even(q, estqD$qD, estS$qD, i))
    out
  } else if (method == "Empirical") {
    if (datatype == "abundance") {
      empqD = iNEXT:::Diversity_profile_MLE(x, q)
      empS = iNEXT:::Diversity_profile_MLE(x, 0)
    } else if (datatype == "incidence_freq") {
      empqD = iNEXT:::Diversity_profile_MLE.inc(x, q)
      empS = iNEXT:::Diversity_profile_MLE.inc(x, 0)
    }

    out = lapply(E.type, function(i) even(q, empqD, empS, i))
    out
  }
  return(out)
}

#
####
# Evenness main function
#
# \code{SC} Estimation (Empirical) of Evenness with order q
#
## R scipts "Evenness" for Chao and Ricotta (2019) Ecology paper.
## This R code is for computing Figures 2, 3 and 4 of Chao and Ricotta (2019) paper.
# installed and loaded before running the scripts.

# @param x data.frame or matrix of data
# @param q a integer vector of the order of Hill number
# @param datatype a binary choose with 'abundance' or 'incidence_freq'
# @param nboot the number of bootstrap resampling times, default is 50
# @param conf a integer value between 0 to 1 for confidence interval
# @param E.type a integer vector between 1 to 6
# @return A list of estimated(empirical) evenness with order q.
#         Different lists represents different classes of Evenness.
#         Each list is combined with order.q and sites.

Evenness <- function (x, q = seq(0, 2, 0.2), datatype = "abundance", nboot = 50,
                      conf = 0.95, E.type = c(1:5))
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

  out2 <- lapply(E.type, function(j) {
    if (datatype == "abundance") {
      out1 <- lapply(1:length(x), function(i) {
        dq <- c(unlist(Evenness.profile(x[[i]], q, "abundance", "Estimated", j)),
                unlist(Evenness.profile(x[[i]], q, "abundance", "Empirical", j)))
        if (nboot > 1) {
          Prob.hat <- iNEXT:::EstiBootComm.Ind(x[[i]])
          Abun.Mat <- rmultinom(nboot, sum(x[[i]]), Prob.hat)
          error <- qnorm(1-(1-conf)/2)*apply(
            apply(Abun.Mat, 2, function(xb) c(unlist(Evenness.profile(xb, q, "abundance", "Estimated", j)),
                                              unlist(Evenness.profile(xb, q, "abundance", "Empirical", j))))
            , 1, sd, na.rm = TRUE)
        } else {
          error = 0
        }
        out <- data.frame(order = rep(q, 2), Evenness = dq, Even.LCL = dq - error,
                          Even.UCL = dq + error, Site = names(x)[i],
                          method = rep(c("Estimated", "Empirical"), each = length(q)))
        out$Even.LCL[out$Even.LCL < 0] <- 0
        out
      })
    } else if (datatype == "incidence") {
      out1 <- lapply(1:length(x), function(i) {
        dq <- c(unlist(Evenness.profile(x[[i]], q, "incidence_freq", "Estimated", j)),
                unlist(Evenness.profile(x[[i]], q, "incidence_freq", "Empirical", j)))
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
            error <- qnorm(1-(1-conf)/2)*apply(
              apply(Incid.Mat, 2, function(yb) c(unlist(Evenness.profile(yb, q, "incidence_freq", "Estimated", j)),
                                                 unlist(Evenness.profile(yb, q, "incidence_freq", "Empirical", j))))
              , 1, sd, na.rm = TRUE)
          }
        } else {
          error = 0
        }
        out <- data.frame(order = rep(q, 2), Evenness = dq, Even.LCL = dq - error,
                          Even.UCL = dq + error, Site = names(x)[i],
                          method = rep(c("Estimated", "Empirical"), each = length(q)))
        out$Even.LCL[out$Even.LCL < 0] <- 0
        out
      })
    }
    out1 <- do.call(rbind, out1)
  })

  names(out2) <- paste("E", E.type, sep="")
  return(out2)
}


#
####
# ggplot for Evenness
#
# \code{ggAsymDiv} The figure for estimation of Evenness with order q
#
# @param output a table generated from Evenness function
# @return a figure of estimated sample completeness with order q

ggEven <- function(output) {
  cbPalette <- rev(c("#999999", "#E69F00", "#56B4E9", "#009E73",
                     "#330066", "#CC79A7", "#0072B2", "#D55E00"))
  fig <- list()
  for (i in 1:length(output)) {
    fig[[i]] = ggplot(output[[i]], aes(x=order, y=Evenness, colour=Site, lty = method)) +
      geom_line(size=1.2) +
      scale_colour_manual(values = cbPalette) +
      # geom_ribbon(data = output[[1]] %>% filter(method=="Estimated"),
      #             aes(ymin=Even.LCL, ymax=Even.UCL, fill=Site),
      #             alpha=0.2, linetype=0) +
      scale_fill_manual(values = cbPalette) +
      labs(x="Order q", y="Evenness", title=names(output[i])) +
      theme_bw(base_size = 18) +
      theme(text=element_text(size=18)) +
      theme(legend.position="bottom", legend.box = "vertical",
            legend.key.width = unit(1.2,"cm"),
            # plot.margin = unit(c(1.5,0.3,1.2,0.3), "lines"),
            legend.title=element_blank(),
            legend.margin=margin(0,0,0,0),
            legend.box.margin = margin(-10,-10,-5,-10),
            text=element_text(size=12),
            plot.margin = unit(c(5.5,5.5,5.5,5.5), "pt"),
            plot.title = element_text(size=20, colour='purple', face="bold.italic",hjust = 0.5))
  }
  fig[[(length(output)+1)]] = ggarrange(plotlist=fig)
  return(fig)
}


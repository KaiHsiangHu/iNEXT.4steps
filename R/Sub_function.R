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
    tmp = (table %>% filter(Order.q %in% c(0,1,2)))[,c("Order.q","Estimate.SC","Assemblage")]
    out = dcast(tmp, Assemblage~Order.q, value.var = "Estimate.SC")
    colnames(out)[-1] = paste("q=", c(0,1,2), sep="")
  }
  if (step==2){
    tmp1 = (table %>% filter((Order.q %in% c(0,1,2)) & (method == "Empirical")))[,c("Assemblage", "Order.q", "qD")]
    names(tmp1) = c("Assemblage", "Diversity", "Observed diversity")
    tmp2 = (table %>% filter((Order.q %in% c(0,1,2)) & (method == "Estimated")))[,c("qD", "s.e.", "qD.LCL", "qD.UCL")]
    names(tmp2) = c("Asymptotic diversity estimate", "s.e.", "LCL", "UCL")
    out = cbind(tmp1, tmp2)
    out$Diversity[out$Diversity == c(0,1,2)] = c("Species richness", "Shannon diversity", "Simpson diversity")
  }
  if (step==3){
    tmp = table[,c("Assemblage","Order.q","qD")]
    C = round(min(table$SC), 3)
    out = dcast(tmp, Assemblage~Order.q, value.var="qD")
    colnames(out) = c(paste("maxC=", C, sep=""),
                      paste("q=", c(0,1,2), sep=""))
  }
  if (step==4){
    if (names(table[1]) == "C")  table = table[-1]
    tmp = (table[[1]] %>%
             filter(Order.q %in% c(0,1,2)))[,c("Order.q","Evenness","Assemblage")]
    out = acast(tmp, Assemblage~Order.q, value.var="Evenness")

    D = (Pielou %>% filter(Order.q == 1))[,c("Assemblage","qD")]
    S = (Pielou %>% filter(Order.q == 0))[,c("Assemblage","qD")]
    out[,1] = sapply(rownames(out), function(x) log(D[D$Assemblage==x,"qD"])/log(S[S$Assemblage==x,"qD"]))
    colnames(out) = c("Pielou J'", paste("q=", c(1,2), sep=""))
  }

  return(out)
}

#
####
#' Sample Completeness main function
#'
#' \code{SC} Estimation of Sample Completeness with order q
#'
#' @param x a matrix/data.frame/list/vector of abundances-based/incidences-based species data.\cr
#' @param q a integer vector for the order of Hill number\cr
#' @param datatype data type of input data: individual-based abundance data (\code{datatype = "abundance"}),
#' sampling-unit-based incidence frequencies data (\code{datatype = "incidence_freq"}) or species by sampling-units incidence matrix (\code{datatype = "incidence_raw"}).\cr
#' @param nboot an integer specifying the number of bootstrap replications, default is 30.\cr
#' @param conf  positive number < 1 specifying the level of confidence interval, default is 0.95.\cr\cr
#' @return a matrix of estimated sample completeness with order q: \cr\cr
#'
#' @examples
#' ## Type (1) example for abundance based data (data.frame)
#' ## Ex.1
#' data(Spider)
#' out1 <- SC(x = Spider, datatype = "abundance")
#' out1
#'
#' ## Type (2) example for incidence based data (list of data.frame)
#' ## Ex.2
#' data(woody_incid)
#' out2 <- SC(x = woody_incid[,c(1,4)], datatype = "incidence_freq")
#' out2
#'
#' @references
#' Chao,A.,Y.Kubota,D.Zelený,C.-H.Chiu.
#' Quantifying sample completeness and comparing diversities among assemblages.
#' @export

SC <- function (x, q = seq(0, 2, 0.2), datatype = "abundance", nboot = 30,
                conf = 0.95)
{
  TYPE <- c("abundance", "incidence", "incidence_freq", "incidence_raw")
  if (is.na(pmatch(datatype, TYPE)))
    stop("invalid datatype")
  if (pmatch(datatype, TYPE) == -1)
    stop("ambiguous datatype")
  datatype <- match.arg(datatype, TYPE)
  class_x <- class(x)[1]
  if (datatype == "incidence_raw") {
    if (class_x == "list") {
      x <- lapply(x, as.incfreq)
    }
    else {
      x <- as.incfreq(x)
    }
    datatype <- "incidence"
  }
  if (datatype %in% c("incidence_freq", "incidence_raw"))
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
        se <- apply( matrix(apply(Abun.Mat,  2, function(xb) sample_completeness(xb, q, "abundance")), nrow=length(q)),
                1, sd, na.rm = TRUE)
      }
      else {
        se = NA
      }
      out <- data.frame(Order.q = q, Estimate.SC = dq,
                        SC.LCL = dq-qnorm(1-(1-conf)/2)*se, SC.UCL = dq+qnorm(1-(1-conf)/2)*se,
                        Assemblage = names(x)[i], method = rep("Estimated", length(q)))
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
          se = 0
          warning("Insufficient data to compute bootstrap s.e.")
        }
        else {
          se <- apply( matrix(apply(Incid.Mat, 2, function(yb) sample_completeness(yb, q, "incidence_freq")), nrow=length(q)),
                   1, sd, na.rm = TRUE)
        }
      }
      else {
        se = NA
      }
      out <- data.frame(Order.q = q, Estimate.SC = dq,
                        SC.LCL = dq-qnorm(1-(1-conf)/2)*se, SC.UCL = dq+qnorm(1-(1-conf)/2)*se,
                        Assemblage = names(x)[i], method = rep("Estimated", length(q)))
      out$SC.LCL[out$SC.LCL < 0] <- 0
      out$SC.UCL[out$SC.UCL > 1] <- 1
      out
    })
    out <- do.call(rbind, out)
  }
  return(out)
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
  } else if (datatype == "incidence_freq") {
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
#' ggplot for Sample Completeness
#'
#' \code{SC} The figure for estimation of Sample Completeness with order q
#'
#' @param output a table generated from SC function
#' @return a figure of estimated sample completeness with order q
#'
#' @examples
#' ## Type (1) example for abundance based data (data.frame)
#' ## Ex.1
#' data(Spider)
#' out1 <- SC(x = Spider, datatype = "abundance")
#' ggSC(out1)
#'
#' ## Type (2) example for incidence based data (list of data.frame)
#' ## Ex.2
#' data(woody_incid)
#' out2 <- SC(x = woody_incid[,c(1,4)], datatype = "incidence_freq")
#' ggSC(out2)
#'
#' @references
#' Chao,A.,Y.Kubota,D.Zelený,C.-H.Chiu.
#' Quantifying sample completeness and comparing diversities among assemblages.
#' @export

ggSC <- function(output) {
  cbPalette <- rev(c("#999999", "#E69F00", "#56B4E9", "#009E73",
                     "#330066", "#CC79A7", "#0072B2", "#D55E00"))
  ggplot(output, aes(x = Order.q, y = Estimate.SC, colour = Assemblage))+
    geom_line(size = 1.2) +
    scale_colour_manual(values = cbPalette) +
    geom_ribbon(aes(ymin = SC.LCL, ymax = SC.UCL, fill = Assemblage), alpha = 0.2, linetype=0) +
    scale_fill_manual(values = cbPalette) +
    labs(x = "Order q", y = "Sample completeness") +
    theme(text = element_text(size=18)) +
    theme(legend.position = "bottom", legend.box = "vertical",
          legend.key.width = unit(1.2,"cm"),
          legend.title = element_blank())
}

#
####
# Empirical Diversity
#
# \code{Diversity_emp} Empirical Diversity with order q
#
# @param x a vector of abundances-based/incidences-based species data.\cr
# @param q a integer vector for the order of Hill number\cr
# @param datatype data type of input data: individual-based abundance data (\code{datatype = "abundance"}),
# sampling-unit-based incidence frequencies data (\code{datatype = "incidence_freq"}) or species by sampling-units incidence matrix (\code{datatype = "incidence_raw"}).\cr
# @return a vector of empirical Diversity with order q: \cr\cr

Diversity_emp = function (x, q, datatype) {
  if (datatype == "abundance") {
    p <- x[x>0]/sum(x)
    n = sum(x)

    Sub <- function(q, p){
      p = p[p>0]

      if(q==0) sum(p>0)
      else if(q==1) exp(-sum(p*log(p)))
      else sum(p^q)^(1/(1-q))
    }
    sapply(q, Sub, p=p)
  } else if (datatype == "incidence_freq") {
    T = x[1]; x = x[-1]

    Sub <- function(q, x){
      pi <- x[x>0]/T
      p = pi/sum(pi)

      if(q==0) sum(p>0)
      else if(q==1) exp(-sum(p*log(p)))
      else sum(p^q)^(1/(1-q))
    }
    sapply(q, Sub, x=x)
  }
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

even.class = function(q, qD, S, E.class) {
  tmp = c()
  if (E.class == 1)
    tmp = ifelse(q!=1, (1-qD^(1-q))/(1-S^(1-q)), log(qD)/log(S))
  if (E.class == 2)
    tmp = ifelse(q!=1, (1-qD^(q-1))/(1-S^(q-1)), log(qD)/log(S))
  if (E.class == 3)
    tmp = (qD-1)/(S-1)
  if (E.class == 4)
    tmp = (1-1/qD)/(1-1/S)
  if (E.class == 5)
    tmp = log(qD)/log(S)

  return(tmp)
}

#
####
# Evenness profile
#
# \code{Evenness.profile} Estimation or Empirical of Evenness with order q
#
# @param x a data.frame, a vector, or a list for data.
# @param q a integer vector for the order of Hill number.
# @param datatype a binary choose with 'abundance' or 'incidence_freq'
# @param method a binary calculation method with 'Estimated' or 'Empirical'
# @param E.class a integer vector between 1 to 6
# @param C a standardized coverage for calculating evenness index
# @return a list of estimated(empirical) evenness with order q, each list is combined with a matrix

Evenness.profile <- function(x, q, datatype=c("abundance","incidence_freq"), method, E.class, C=NULL) {
  if (method == "Estimated") {
    estqD = estimateD(x, q, datatype, base="coverage", level=C, nboot=0)
    estS = estimateD(x, 0, datatype, base="coverage", level=C, nboot=0)

    out = lapply(E.class, function(i) {
      tmp = sapply(names(x), function(k) even.class(q, estqD[estqD$Assemblage==k, "qD"], estS[estS$Assemblage==k, "qD"], i))
      if(class(tmp)[1] %in% c("numeric","integer")) {tmp = t(as.matrix(tmp, nrow=1))}
      rownames(tmp) = q
      tmp
      })
  } else if (method == "Empirical") {

    if (datatype == "abundance") {
      empqD = sapply(x, function(k) iNEXT:::Diversity_profile_MLE(k, q))
      empS = sapply(x, function(k) iNEXT:::Diversity_profile_MLE(k, 0))
    } else if (datatype == "incidence_freq") {
      empqD = sapply(x, function(k) iNEXT:::Diversity_profile_MLE.inc(k, q))
      empS = sapply(x, function(k) iNEXT:::Diversity_profile_MLE.inc(k, 0))
    }

    out = lapply(E.class, function(i) {
      if ((is.vector(empqD)==TRUE) & (length(empqD)==1)) {
        name = names(empqD); empqD = matrix(empqD); colnames(empqD) = name}

      tmp = sapply(names(x), function(k) even.class(q, empqD[,k], empS[k], i))
      if(class(tmp)[1] %in% c("numeric","integer")) {tmp = t(as.matrix(tmp, nrow=1))}
      rownames(tmp) = q
      tmp
    })
  }

  names(out) = paste("E", E.class, sep="")
  return(out)
}

#
####
#' Evenness main function
#'
#' \code{Evenness} Estimation (Empirical) of Evenness with order q
#'
#' R scipts "Evenness" for Chao and Ricotta (2019) Ecology paper.
#' This R code is for computing Figures 2, 3 and 4 of Chao and Ricotta (2019) paper.
#' installed and loaded before running the scripts.

#' @param x a matrix/data.frame/list/vector of abundances-based/incidences-based species data.\cr
#' @param q a integer vector of the order of Hill number\cr
#' @param datatype data type of input data: individual-based abundance data (\code{datatype = "abundance"}),
#' sampling-unit-based incidence frequencies data (\code{datatype = "incidence_freq"}) or species by sampling-units incidence matrix (\code{datatype = "incidence_raw"}).\cr
#' @param method a binary calculation method with 'Estimated' or 'Empirical'\cr
#' @param nboot an integer specifying the number of bootstrap replications, default is 30.\cr
#' @param conf a positive number < 1 specifying the level of confidence interval, default is 0.95.\cr
#' @param E.class a integer vector between 1 to 6
#' @param C a standardized coverage for calculating evenness index
#' @return A list of estimated(empirical) evenness with order q.\cr
#'         Different lists represents different classes of Evenness.\cr
#'         Each list is combined with order.q and sites.\cr
#'         If "method" is estimated, then fist list will be named "C" which means the
#'         maximum standardized coverage between all double reference sample size.\cr\cr
#' \code{$summary} individual summary of 4 steps of data. \cr\cr
#'
#' @examples
#' ## Type (1) example for abundance based data (data.frame)
#' ## Ex.1
#' data(Spider)
#' out1 <- Evenness(x = Spider, datatype = "abundance")
#' out1
#'
#' ## Type (2) example for incidence based data (list of data.frame)
#' ## Ex.2
#' data(woody_incid)
#' out2 <- Evenness(x = woody_incid[,c(1,4)], datatype = "incidence_freq")
#' out2
#'
#' @references
#' Chao,A.and Ricotta,C.(2019).Quantifying evenness and linking it to diversity, beta diversity, and similarity.
#' @export

Evenness <- function (x, q = seq(0, 2, 0.2), datatype = "abundance", method = "Estimated",
                      nboot = 30, conf = 0.95, E.class = c(1:5), C=NULL)
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

  kind <- c("Estimated", "Empirical")
  if (length(method) > 1)
    stop("only one calculation method")
  if (is.na(pmatch(method, kind)))
    stop("invalid method")
  if (pmatch(method, kind) == -1)
    stop("ambiguous method")

  class <- c(1:5)
  if (sum(E.class %in% class) != length(E.class))
    stop("invalid E.class")

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
    qD <- Evenness.profile(x, q, "abundance", method, E.class, C)
    qD <- map(qD, as.vector)
    
    if (nboot > 1) {
      Prob.hat <- lapply(1:length(x), function(i) iNEXT:::EstiBootComm.Ind(x[[i]]))
      Abun.Mat <- lapply(1:length(x), function(i) rmultinom(nboot, sum(x[[i]]), Prob.hat[[i]]))

      error = apply( matrix(sapply(1:nboot, function(b) {
                    dat = lapply(1:length(Abun.Mat),function(j) Abun.Mat[[j]][,b])
                    names(dat) = paste("Site", 1:length(dat), sep="")
                    dat.qD = Evenness.profile(dat, q, "abundance", method, E.class, C)
                    unlist(dat.qD)
                    }), nrow=length(q)*length(E.class)*length(Abun.Mat))
        , 1, sd, na.rm = TRUE)

      error = matrix(error, ncol=length(E.class))
      se = split(error, col(error))

    } else {
      se = lapply(1:length(E.class), function(x) NA)
    }

    out <- lapply(1:length(E.class), function(k) {
      tmp = data.frame(Order.q = rep(q, length(x)), Evenness = as.vector(qD[[k]]), s.e. = as.vector(se[[k]]),
                      Even.LCL = as.vector(qD[[k]] - qnorm(1-(1-conf)/2)*se[[k]]), Even.UCL = as.vector(qD[[k]] + qnorm(1-(1-conf)/2)*se[[k]]),
                      Assemblage = rep(names(x), each=length(q)), method = rep( method, length(q)*length(x))
      )
      tmp$Even.LCL[tmp$Even.LCL < 0] <- 0
      tmp
    })
    if (is.null(C) == TRUE) C = unique(estimateD(x, q = 0, datatype = datatype, base = "coverage", nboot=0)$goalSC)
    if (method=="Estimated") {out <- append(C, out)}

  } else if (datatype == "incidence") {
    qD <- Evenness.profile(x, q, "incidence_freq", method, E.class, C)
    qD <- map(qD, as.vector)
    
    if (nboot > 1) {
      nT <- lapply(1:length(x), function(i) x[[i]][1])
      Prob.hat <- lapply(1:length(x), function(i) iNEXT:::EstiBootComm.Sam(x[[i]]))
      Incid.Mat <- lapply(1:length(x), function(i) t(sapply(Prob.hat[[i]], function(p) rbinom(nboot, nT[[i]], p))))
      Incid.Mat <- lapply(1:length(x), function(i) matrix(c(rbind(nT[[i]], Incid.Mat[[i]])), ncol = nboot))

      for (i in 1:length(Incid.Mat)) {
        tmp = which(colSums(Incid.Mat[[i]][-1,]) == nT[[i]])
        if (length(tmp) > 0)
          Incid.Mat <- lapply(Incid.Mat, function(k) k[, -tmp])
      }

      if (ncol(Incid.Mat[[1]]) == 0) {
        error = as.list(rep(0, length(E.class)))
        warning("Insufficient data to compute bootstrap s.e.")
      } else {
        error = apply(  matrix(sapply(1:nboot, function(b) {
          dat = lapply(1:length(Incid.Mat),function(j) Incid.Mat[[j]][,b])
          names(dat) = paste("Site", 1:length(dat), sep="")
          dat.qD = Evenness.profile(dat, q, "incidence_freq", method, E.class, C)
          unlist(dat.qD)  }), nrow=length(q)*length(E.class)*length(Incid.Mat))
          , 1, sd, na.rm = TRUE)

        error = matrix(error, ncol=length(E.class))
        se = split(error, col(error))
      }

    } else {
      se = lapply(1:length(E.class), function(x) NA)
    }

    out <- lapply(1:length(E.class), function(k) {
      tmp = data.frame(Order.q = rep(q, length(x)), Evenness = as.vector(qD[[k]]), s.e. = as.vector(se[[k]]),
                       Even.LCL = as.vector(qD[[k]] - qnorm(1-(1-conf)/2)*se[[k]]), Even.UCL = as.vector(qD[[k]] + qnorm(1-(1-conf)/2)*se[[k]]),
                       Assemblage = rep(names(x), each=length(q)), method = rep(method, length(q)*length(x))
      )
      tmp$Even.LCL[tmp$Even.LCL < 0] <- 0
      tmp
    })
    if (is.null(C) == TRUE) C = unique(estimateD(x, q = 0, datatype = datatype, base = "coverage", nboot=0)$goalSC)
    if (method=="Estimated") {out <- append(C, out)}
  }

  if (method=="Estimated") {
    names(out) = c("Coverage", paste("E", E.class, sep = ""))
  } else if (method=="Empirical") {names(out) = paste("E", E.class, sep = "")}

  return(out)
}


#
####
#' ggplot for Evenness
#
#' \code{ggEven} The figure for estimation of Evenness with order q\cr
#'
#' @param output a table generated from Evenness function\cr
#' @return a figure of estimated sample completeness with order q\cr
#'
#' @examples
#' ## Type (1) example for abundance based data (data.frame)
#' ## Ex.1
#' data(Spider)
#' out1 <- Evenness(x = Spider, datatype = "abundance")
#' ggEven(out1)
#'
#' ## Type (2) example for incidence based data (list of data.frame)
#' ## Ex.2
#' data(woody_incid)
#' out2 <- Evenness(x = woody_incid[,c(1,4)], datatype = "incidence_freq")
#' ggEven(out2)
#'
#' @references
#' Chao,A.and Ricotta,C.(2019).Quantifying evenness and linking it to diversity, beta diversity, and similarity.
#' @export

ggEven <- function(output) {
  if (names(output[1]) == "Coverage")  output = output[-1]
  cbPalette <- rev(c("#999999", "#E69F00", "#56B4E9", "#009E73",
                     "#330066", "#CC79A7", "#0072B2", "#D55E00"))
  classdata = cbind(do.call(rbind, output),
                    class = rep(names(output), each=nrow(output[[1]])))

  fig = ggplot(classdata, aes(x=Order.q, y=Evenness, colour=Assemblage, lty = method)) +
    geom_line(size=1.2) +
    scale_colour_manual(values = cbPalette) +
    geom_ribbon(data = classdata %>% filter(method=="Estimated"),
                aes(ymin=Even.LCL, ymax=Even.UCL, fill=Assemblage),
                alpha=0.2, linetype=0) +
    geom_ribbon(data = classdata %>% filter(method=="Empirical"),
                aes(ymin=Even.LCL, ymax=Even.UCL, fill=Assemblage),
                alpha=0.2, linetype=0) +
    scale_fill_manual(values = cbPalette) +
    labs(x="Order q", y="Evenness") +
    # theme_bw(base_size = 18) +
    theme(text=element_text(size=18)) +
    theme(legend.position = "bottom", legend.box = "vertical",
          legend.key.width = unit(1.2,"cm"),
          # plot.margin = unit(c(1.5,0.3,1.2,0.3), "lines"),
          legend.title = element_blank(),
          legend.margin = margin(0,0,0,0),
          legend.box.margin = margin(-10,-10,-5,-10),
          text = element_text(size=12),
          plot.margin = unit(c(5.5,5.5,5.5,5.5), "pt")
    )
  
  if (length(output) != 1) fig = fig +
    facet_wrap(~class) +
    theme(strip.text.x = element_text(size=12, colour = "purple", face="bold"))
  
  return(fig)
}


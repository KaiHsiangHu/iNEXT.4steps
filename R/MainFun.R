#' Main function for complete 4-step analysis
#'
#' \code{iNEXT4steps} computes all statistics in the complete 4-step analysis and visualizes the output. It computes sample completeness, observed and asymptotic diversity, size-based and coverage-based standardized diversity, and evenness. 
#' 
#' 
#' @param data (a) For \code{datatype = "abundance"}, data can be input as a vector of species abundances (for a single assemblage), matrix/data.frame (species by assemblages), or a list of species abundance vectors. \cr
#' (b) For \code{datatype = "incidence_raw"}, data can be input as a list of matrix/data.frame (species by sampling units); data can also be input as a matrix/data.frame by merging all sampling units across assemblages based on species identity; in this case, the number of sampling units (nT, see below) must be input.
#' @param q a numerical vector specifying the orders of q that will be used to compute sample completeness and evenness as well as plot the relevant profiles. Default is \code{seq(0, 2, by = 0.2)}.
#' @param datatype data type of input data: individual-based abundance data (\code{datatype = "abundance"}) or species by sampling-units incidence matrix (\code{datatype = "incidence_raw"}) with all entries being 0 (non-detection) or 1 (detection).
#' @param nboot a positive integer specifying the number of bootstrap replications when assessing sampling uncertainty and constructing confidence intervals. Enter 0 to skip the bootstrap procedures. Default is 30.
#' @param conf a positive number < 1 specifying the level of confidence interval. Default is 0.95.
#' @param nT (required only when \code{datatype = "incidence_raw"} and input data in a single matrix/data.frame) a vector of positive integers specifying the number of sampling units in each assemblage. If assemblage names are not specified (i.e., \code{names(nT) = NULL}), then assemblages are automatically named as "Assemblage1", "Assemblage2",..., etc. 
#' @param details a logical variable to indicate whether the detailed numerical values for each step are displayed. Default is \code{FALSE}.
#' 
#' @import ggplot2
#' @import reshape2
#' @import dplyr
#' @import ggpubr
#' @import purrr
#' @import iNEXT.3D
#' @importFrom stats qnorm
#' @importFrom stats rbinom
#' @importFrom stats rmultinom
#' @importFrom stats sd
#' @importFrom grDevices hcl
#' 
#' @return a list of three of objects: \cr\cr
#' \code{$summary} Numerical table for each individual step.
#' \item{Assemblage}{the assemblage names.}
#' \item{qTD}{'Species richness' represents the taxonomic diversity of order q=0; 'Shannon diversity' represents the taxonomic diversity of order q=1, 'Simpson diversity' represents the taxonomic diversity of order q=2.}
#' \item{TD_obs}{the observed taxonomic diversity value of order q.}
#' \item{TD_asy}{the estimated asymptotic diversity value of order q.}
#' \item{s.e.}{the bootstrap standard error of the estimated asymptotic diversity of order q.}
#' \item{qTD.LCL, qTD.UCL}{the bootstrap lower and upper confidence limits for the estimated asymptotic diversity of order q at the specified level in the setting (with a default value of 0.95).}
#' \item{Pielou J'}{a widely used evenness measure based on Shannon entropy.} \cr\cr
#' \code{$figure} six figures including five individual figures (for STEPS 1, 2a, 2b, 3 and 4 respectively) and a complete set of five plots. \cr\cr
#' \code{$details} (only when \code{details = TRUE}). The numerical output for plotting all figures. \cr\cr
#' 
#' 
#' @examples
#' \donttest{
#' ## Complete 4-step analysis for abundance data
#' data(Data_spider)
#' Four_Steps_out1 <- iNEXT4steps(data = Data_spider, datatype = "abundance")
#' Four_Steps_out1
#' 
#' 
#' ## Complete 4-step analysis for incidence data
#' data(Data_woody_plant)
#' Four_Steps_out2 <- iNEXT4steps(data = Data_woody_plant, datatype = "incidence_raw")
#' Four_Steps_out2
#' }
#' 
#' 
#' @export

iNEXT4steps <- function(data, q = seq(0, 2, 0.2), datatype = "abundance", 
                        nboot = 30, conf = 0.95, nT = NULL, details = FALSE) 
{
  if ((length(data) == 1) && (inherits(data, c("numeric", "integer"))))
    stop("Error: Your data does not have enough information.")
  
  logic = c("TRUE", "FALSE")
  if (is.na(pmatch(details, logic)))
    stop("invalid details setting")
  if (pmatch(details, logic) == -1)
    stop("ambiguous details setting")
  
  ## SC ##
  SC.table <- Completeness(data, q = q, datatype = datatype, nboot = nboot, conf = conf, nT = nT)
  
  ## iNEXT ##
  iNEXT.table <- iNEXT3D(data, diversity = "TD", q = c(0, 1, 2), datatype = datatype, nboot = nboot, conf = conf, nT = nT)
  
  ## Asymptotic ##
  qD.table <- ObsAsy3D(data, diversity = "TD", q = q, datatype = datatype, nboot = nboot, conf = conf, nT = nT)
  
  ## Evenness ##
  Even.table <- Evenness(data, q = q, datatype = datatype, method = "Estimated", nboot = nboot, conf = conf, nT = nT, E.class = 3)
  Cmax = unique(Even.table$E3$SC)
  
  
  plot.names = c("(a) STEP 1.\n Sample completeness profiles",
                 "(b) STEP 2a.\n Size-based rarefaction/extrapolation",
                 "(c) STEP 2b.\n Asymptotic and empirical diversity profiles",
                 "(d) STEP 3.\n Coverage-based rarefaction/extrapolation",
                 "(e) STEP 4.\n Evenness profiles")
  
  table.names = c("STEP 1. Sample completeness profiles",
                  "STEP 2b. Observed diversity values and asymptotic estimates",
                  "STEP 3. Non-asymptotic coverage-based rarefaction and extrapolation analysis",
                  "STEP 4. Evenness among species abundances of orders q = 1 and 2 at Cmax based on the normalized slope of a diversity profile")
  
  
  if (length(unique((Even.table[[1]]$Assemblage))) > 1) {
    iNEXT.table[[2]]$size_based$Assemblage = factor(iNEXT.table[[2]]$size_based$Assemblage)
    
    level = levels(iNEXT.table[[2]]$size_based$Assemblage)
    
    SC.table$Assemblage = factor(SC.table$Assemblage, level)
    
    qD.table$Assemblage = factor(qD.table$Assemblage, level)
    
    Even.table[[1]]$Assemblage = factor(Even.table[[1]]$Assemblage, level)
  }

  ## 5 figures ##
  if (length(unique(SC.table$Assemblage)) <= 8) {
    
    SC.plot <- ggCompleteness(SC.table) +
      labs(title = plot.names[1]) +
      theme(text = element_text(size = 12),
            plot.margin = unit(c(5.5, 5.5, 5.5, 5.5), "pt"),
            plot.title = element_text(size = 11, colour = 'blue', face = "bold", hjust = 0))

    size.RE.plot <- ggiNEXT3D(iNEXT.table, type = 1, facet.var = "Order.q", color.var = "Order.q") +
      labs(title = plot.names[2]) +
      theme(text = element_text(size = 12),
            plot.margin = unit(c(5.5, 5.5, 5.5, 5.5), "pt"),
            legend.margin = margin(0, 0, 0, 0),
            plot.title = element_text(size = 11, colour = 'blue', face = "bold", hjust = 0))

    cover.RE.plot <- ggiNEXT3D(iNEXT.table, type = 3, facet.var = "Order.q", color.var = "Order.q") +
      labs(title = plot.names[4]) +
      theme(text = element_text(size = 12),
            plot.margin = unit(c(5.5, 5.5, 5.5, 5.5), "pt"),
            legend.margin = margin(0, 0, 0, 0),
            plot.title = element_text(size = 11, colour = 'blue', face = "bold", hjust = 0))

    AO.plot <- ggObsAsy3D(qD.table) +
      labs(title = plot.names[3]) +
      theme(text = element_text(size = 12),
            plot.margin = unit(c(5.5, 5.5, 5.5, 5.5), "pt"),
            plot.title = element_text(size = 11, colour = 'blue', face = "bold", hjust = 0))

    even.plot <- ggEvenness(Even.table) +
      labs(title = plot.names[5]) +
      theme(text = element_text(size = 12),
            plot.margin = unit(c(5.5, 5.5, 5.5, 5.5), "pt"),
            plot.title = element_text(size = 11, colour = 'blue', face = "bold", hjust = 0))

    legend.p = get_legend(SC.plot + theme(legend.direction = "vertical"))
    
    steps.plot = ggarrange(SC.plot       + guides(color = "none", fill = "none"),
                           size.RE.plot  + guides(color = "none", fill = "none", shape = "none"),
                           AO.plot       + guides(color = "none", fill = "none"),
                           cover.RE.plot + guides(color = "none", fill = "none", shape = "none"),
                           even.plot     + guides(color = "none", fill = "none"),
                           legend.p, nrow = 3, ncol = 2
    )
  } else { warning("The number of communities exceeds eight. We don't show the figures.") }
  
  estD = estimate3D(data, diversity = 'TD', q = c(0, 1, 2), datatype, base = "coverage", level = NULL, nboot = 0, nT = nT)
  
  
  ##  Outpue_summary ##
  if (sum(q %in% c(0,1,2)) == 3) 
    summary = list(summary.deal(SC.table, 1),
                   iNEXT.table[[3]] %>%  
                     lapply(FUN = function(x) if(is.numeric(x)) round(x,2) else x) %>% data.frame,
                   summary.deal(estD, 3),
                   summary.deal(Even.table, 4, estD)
    ) else {
      
      SC.tableq012 <- Completeness(data, q = c(0,1,2), datatype = datatype, nboot = 0, nT = nT)
      
      Even.tableq012 <- Evenness(data, q = c(0,1,2), datatype = datatype, method = "Estimated", nboot = 0, nT = nT, E.class = 3)
      
      summary = list(summary.deal(SC.tableq012, 1),
                     
                     iNEXT.table[[3]] %>%  
                       lapply(FUN = function(x) if(is.numeric(x)) round(x,2) else x) %>% data.frame,
                     
                     summary.deal(estD, 3),
                     
                     summary.deal(Even.tableq012, 4, estD)
                     )
    }
  
  
  names(summary) = table.names
  
  
  ##  Output ##
  if (details == FALSE) {

    if (length(unique(SC.table$Assemblage)) <= 8) {
      
      ans <- list(summary = summary,
                  
                  figure = list(SC.plot, size.RE.plot, AO.plot, cover.RE.plot, even.plot, steps.plot))
      
    } else { ans <- list(summary = summary) }

  } else if (details == TRUE) {
    
    tab = list("Sample completeness" = SC.table, 
               "iNEXT" = iNEXT.table[[2]],
               "Observed and asymptotic diversity" = qD.table, 
               "Evenness" = Even.table)

    if (length(unique(SC.table$Assemblage)) <= 8) {
      
      ans <- list(summary = summary,
                  
                  figure = list(SC.plot, size.RE.plot, AO.plot, cover.RE.plot, even.plot, steps.plot),
                  
                  details = tab)
      
    } else { ans <- list(summary = summary, details = tab)}
    
  }
  
  return(ans)
}



# Generate the summary table from the deatils of each function.
#
# \code{summary.deal} Generate the summary table of four lists from the deatils of each function.
#
# @param table a table
# @param step a numerical value of which step summary do you need
# @param Pielou a logical variable whether you want to generate step 4 summary
# @return a matrix or a list according different steps

summary.deal <- function(table, step, Pielou = NULL) {
  
  if (step == 1) {
    
    tmp = (table %>% filter(Order.q %in% c(0, 1, 2)))[, c("Order.q", "Estimate.SC", "Assemblage")]
    
    out = dcast(tmp, Assemblage ~ Order.q, value.var = "Estimate.SC") %>% 
      lapply(FUN = function(x) if(is.numeric(x)) round(x,2)
                     else x) %>% data.frame()
    colnames(out)[-1] = paste("q = ", c(0, 1, 2), sep="")
  }
  
  if (step == 3){
    
    tmp = table[,c(1, 2, 6)]
    C = round(unique(table$SC), 3)
    
    out = dcast(tmp, Assemblage ~ Order.q, value.var = colnames(tmp)[3]) %>% 
      lapply(FUN = function(x) if(is.numeric(x)) round(x,2)
             else x) %>% data.frame()
    
    colnames(out) = c(paste("Cmax = ", C, sep = ""),
                      paste("q = ", c(0, 1, 2), sep = ""))
  }
  
  if (step == 4){
    
    tmp = (table[[1]] %>%
             filter(Order.q %in% c(0, 1, 2)))[, c("Order.q", "Evenness", "Assemblage")]
    out = dcast(tmp, Assemblage ~ Order.q, value.var = "Evenness")
    
    C = round(unique(table$E3$SC), 3)
    D = (Pielou %>% filter(Order.q == 1))[,c("Assemblage", "qTD")]
    S = (Pielou %>% filter(Order.q == 0))[,c("Assemblage", "qTD")]
    out[,2] = sapply(out$Assemblage, function(x) log(D[D$Assemblage == x, "qTD"]) / log(S[S$Assemblage == x, "qTD"]))
    
    colnames(out) = c(paste("Cmax = ", C, sep = ""),
                      "Pielou J'", paste("q = ", c(1, 2), sep=""))
    out[,-1] <- round(out[,-1],2)
  }
  
  return(out)
}


# \code{sample_completeness} Estimation of Sample Completeness with order q
#
# @param x a vector of data
# @param q a integer vector of the order of Hill number
# @param datatype a binary choose with 'abundance' or 'incidence_freq'
# @return a vector of estimated sample completeness with order q

sample_completeness = function(x, q, datatype = c("abundance", "incidence_freq")){
  
  if(datatype == "abundance"){
    x = x[x > 0]
    n = sum(x)
    f1 = sum(x == 1); f2 = sum(x == 2)
    A = ifelse(f2 > 0, 2 * f2 / ((n - 1) * f1 + 2 * f2), ifelse(f1 > 0, 2 / ((n - 1) * (f1 - 1) + 2),1))
    
    sc.abund = function(q){
      
      if (q == 0){
        
        S_obs = length(x)
        f0_hat = ifelse(f2 == 0, ((n - 1) / n) * (f1 * (f1 - 1) / 2), ((n - 1) / n) * ((f1^2) / (2 * f2)))
        c_hat = S_obs / (S_obs + f0_hat)
        return(c_hat)
        
      } else if (q == 1){
        
        c_hat = 1 - (f1 / n) * (1 - A)
        return(c_hat)
        
      } else if (q == 2){
        
        x = x[x >= 2]
        c_hat = 1 - (f1 / n) * ((A * (1 - A)) / sum(x * (x - 1) / (n * (n - 1))))
        return(c_hat)
        
      } else {
        
        r <- 0:(n-1)
        sort.data = sort(unique(x))
        tab = table(x)
        
        term = sapply(sort.data, function(z){
          k = 0:(n - z)
          sum(choose(k - q, k) * exp(lchoose(n - k - 1, z - 1) - lchoose(n, z)))
        })
        
        lambda_hat = sum(tab * term) + ifelse(f1 == 0 | A == 1, 0, f1 / n * (1 - A)^(1 - n) * (A^(q - 1) - sum(choose(q - 1, r) * (A - 1)^r)))
        c_hat = 1 - ((f1 / n) * (A^(q - 1)) * (1 - A) / lambda_hat)
        return(c_hat)
        
      }
    }
    
    est = sapply(q, sc.abund)
    
  } else if (datatype == "incidence_freq") {
    
    t = x[1]
    x = x[-1]; x = x[x>0]
    u = sum(x)
    Q1 = sum(x == 1); Q2 = sum(x == 2)
    B = ifelse(Q2 > 0, 2 * Q2 / ((t - 1) * Q1 + 2 * Q2), ifelse(Q1 > 0, 2 / ((t - 1) * (Q1 - 1) + 2), 1))
    
    sc.incid = function(q){
      
      if (q == 0){
        
        S_obs = length(x)
        f0_hat = ifelse(Q2 == 0, ((t - 1) / t) * (Q1 * (Q1 - 1) / 2), ((t - 1) / t) * ((Q1^2) / (2 * Q2)))
        c_hat = S_obs / (S_obs + f0_hat)
        return(c_hat)
        
      } else if (q == 1){
        
        c_hat = 1 - (Q1 / u) * (1 - B)
        return(c_hat)
        
      } else if (q == 2){
        
        x = x[x >= 2]
        c_hat = 1 - (t - 1) * Q1 * ((B * (1 - B)) / sum(x * (x - 1)))
        return(c_hat)
        
      } else {
        
        r <- 0:(t - 1)
        sort.data = sort(unique(x))
        tab = table(x)
        
        term = sapply(sort.data,function(z){
          k = 0:(t - z)
          sum(choose(k - q, k) * exp(lchoose(t - k - 1, z - 1) - lchoose(t, z)))
        })
        
        phi_hat = sum(tab * term) + ifelse(Q1 == 0 | B == 1, 0, Q1 / t * (1 - B)^(1 - t) * (B^(q - 1) - sum(choose(q - 1,r) * (B - 1)^r)))
        c_hat = 1 - ((Q1 / t) * (B^(q - 1)) * (1 - B) / phi_hat)
        return(c_hat)
        
      }
    }
    
    est = sapply(q, sc.incid)
  }
  
  return(est)
}


#' Main function for STEP 1: Assessment of sample completeness 
#'
#' \code{Completeness} computes sample completeness estimates of orders q = 0 to 2 in increments of 0.2 (by default).
#'
#' @param data (a) For \code{datatype = "abundance"}, data can be input as a vector of species abundances (for a single assemblage), matrix/data.frame (species by assemblages), or a list of species abundance vectors. \cr
#' (b) For \code{datatype = "incidence_raw"}, data can be input as a list of matrix/data.frame (species by sampling units); data can also be input as a matrix/data.frame by merging all sampling units across assemblages based on species identity; in this case, the number of sampling units (nT, see below for this argument) must be input. 
#' @param q a numerical vector specifying the orders of sample completeness. Default is \code{seq(0, 2, by = 0.2)}.
#' @param datatype data type of input data: individual-based abundance data (\code{datatype = "abundance"}) or species by sampling-units incidence matrix (\code{datatype = "incidence_raw"}) with all entries being 0 (non-detection) or 1 (detection).
#' @param nboot a positive integer specifying the number of bootstrap replications when assessing sampling uncertainty and constructing confidence intervals. Enter 0 to skip the bootstrap procedures. Default is 30.
#' @param conf a positive number < 1 specifying the level of confidence interval. Default is 0.95.
#' @param nT (required only when \code{datatype = "incidence_raw"} and input data in a single matrix/data.frame) a vector of positive integers specifying the number of sampling units in each assemblage. If assemblage names are not specified (i.e., \code{names(nT) = NULL}), then assemblages are automatically named as "Assemblage1", "Assemblage2",..., etc. 
#' @return a matrix of estimated sample completeness of order q: 
#'         \item{Order.q}{the order of sample completeness.}
#'         \item{Estimate.SC}{the estimated sample completeness of order q.}
#'         \item{s.e.}{standard error of sample completeness estimate.}
#'         \item{SC.LCL, SC.UCL}{the bootstrap lower and upper confidence limits for the sample completeness of order q at the specified level (with a default value of \code{0.95}).}
#'         \item{Assemblage}{the assemblage name.}
#' 
#'
#' @examples
#' \donttest{
#' ## Sample completeness for abundance data
#' data(Data_spider)
#' SC_out1 <- Completeness(data = Data_spider, datatype = "abundance")
#' SC_out1
#' }
#' 
#' ## Sample completeness for incidence raw data
#' data(Data_woody_plant)
#' SC_out2 <- Completeness(data = Data_woody_plant, datatype = "incidence_raw")
#' SC_out2
#' 
#' 
#' @export

Completeness <- function (data, q = seq(0, 2, 0.2), datatype = "abundance", nboot = 30,
                          conf = 0.95, nT = NULL)
{
  TYPE <- c("abundance", "incidence", "incidence_freq", "incidence_raw")
  if (is.na(pmatch(datatype, TYPE)))
    stop("invalid datatype")
  if (pmatch(datatype, TYPE) == -1)
    stop("ambiguous datatype")
  datatype <- match.arg(datatype, TYPE)
  
  if(datatype == "incidence") stop("Completeness can only accept datatype = 'incidence_freq' or 'datatype = incidence_raw'.") 
  
  if (datatype == "incidence_raw") {
    
    if (length(data) != 1) {
      data = iNEXT.3D:::as.incfreq(data, nT = nT)
      datatype = "incidence_freq"
    }
    
    if (length(data) == 1) {
      name = names(data)
      data = list(iNEXT.3D:::as.incfreq(data, nT = nT))
      names(data) = name
      datatype = "incidence_freq"
    }
  }
  
  if (inherits(data, c("numeric", "integer"))) {
    data <- list("Assemblage_1" = data)
  }
  
  if (inherits(data, c("data.frame", "matrix"))) {
    
    datalist <- lapply(1:ncol(data), function(i) data[, i])
    
    if (is.null(colnames(data)))
      names(datalist) <- paste0("Assemblage_", 1:ncol(data))
    else names(datalist) <- colnames(data)
    
    data <- datalist
  }
  
  if (inherits(data, "list")) {
    if (is.null(names(data))) names(data) = paste0("Assemblage_", 1:length(data))
  }
  
  if (datatype == "abundance") {
    
    out <- lapply(1:length(data), function(i) {
      dq <- sample_completeness(data[[i]], q, "abundance")
      
      if (nboot > 1) {
        Prob.hat <- iNEXT.3D:::EstiBootComm.Ind(data[[i]])
        Abun.Mat <- rmultinom(nboot, sum(data[[i]]), Prob.hat)
        
        se <- apply( matrix(apply(Abun.Mat,  2, function(xb) sample_completeness(xb, q, "abundance")), nrow=length(q)),
                     1, sd, na.rm = TRUE)
      }
      else {
        se = NA
      }
      out <- data.frame(Order.q = q, 
                        Estimate.SC = dq, 
                        s.e. = se,
                        SC.LCL = dq - qnorm(1 - (1 - conf) / 2) * se, 
                        SC.UCL = dq + qnorm(1 - (1 - conf) / 2) * se,
                        Assemblage = names(data)[i])
      out$SC.LCL[out$SC.LCL < 0] <- 0
      out$SC.UCL[out$SC.UCL > 1] <- 1
      out
    })
    out <- do.call(rbind, out)
    
  } else if (datatype == "incidence_freq") {
    
    out <- lapply(1:length(data), function(i) {
      dq <- sample_completeness(data[[i]], q, "incidence_freq")
      
      if (nboot > 1) {
        nT <- data[[i]][1]
        Prob.hat <- iNEXT.3D:::EstiBootComm.Sam(data[[i]])
        Incid.Mat <- t(sapply(Prob.hat, function(p) rbinom(nboot, nT, p)))
        Incid.Mat <- matrix(c(rbind(nT, Incid.Mat)), ncol = nboot)
        
        se <- apply( matrix(apply(Incid.Mat, 2, function(yb) sample_completeness(yb, q, "incidence_freq")), nrow = length(q)),
                     1, sd, na.rm = TRUE)
      }
      else {
        se = NA
      }
      out <- data.frame(Order.q = q, 
                        Estimate.SC = dq, 
                        s.e. = se,
                        SC.LCL = dq - qnorm(1 - (1 - conf) / 2) * se, 
                        SC.UCL = dq + qnorm(1 - (1 - conf) / 2) * se,
                        Assemblage = names(data)[i])
      out$SC.LCL[out$SC.LCL < 0] <- 0
      out$SC.UCL[out$SC.UCL > 1] <- 1
      out
    })
    
    out <- do.call(rbind, out)
  }
  return(out)
}


#' ggplot for depicting sample completeness profiles
#'
#' \code{ggCompleteness} is a \code{ggplot2} extension for \code{Completeness} object to plot sample completeness with order q between 0 and 2.
#'
#' @param output output obtained from the function \code{Completeness}.
#' @return a figure depicting the estimated sample completeness with respect to the order q.
#' 
#' 
#' @examples
#' \donttest{
#' ## Sample completeness profile for abundance data
#' data(Data_spider)
#' SC_out1 <- Completeness(data = Data_spider, datatype = "abundance")
#' ggCompleteness(SC_out1)
#' }
#' 
#' ## Sample completeness profile for incidence raw data
#' data(Data_woody_plant)
#' SC_out2 <- Completeness(data = Data_woody_plant, datatype = "incidence_raw")
#' ggCompleteness(SC_out2)
#' 
#' 
#' @export

ggCompleteness <- function(output) {
  
  # Check if the number of unique 'Assemblage' is 8 or less
  if (length(unique(output$Assemblage)) <= 8){
    cbPalette <- rev(c("#999999", "#E69F00", "#56B4E9", "#009E73", 
                       "#330066", "#CC79A7", "#0072B2", "#D55E00"))
  }else{
    # If there are more than 8 assemblages, start with the same predefined color palette
    # Then extend the palette by generating additional colors using the 'ggplotColors' function
    cbPalette <- rev(c("#999999", "#E69F00", "#56B4E9", "#009E73", 
                       "#330066", "#CC79A7", "#0072B2", "#D55E00"))
    cbPalette <- c(cbPalette, ggplotColors(length(unique(output$Assemblage))-8))
  }
  
  
  ggplot(output, aes(x = Order.q, y = Estimate.SC, colour = Assemblage))+
    geom_line(size = 1.2) +
    scale_colour_manual(values = cbPalette) +
    geom_ribbon(aes(ymin = SC.LCL, ymax = SC.UCL, fill = Assemblage), alpha = 0.2, linetype=0) +
    scale_fill_manual(values = cbPalette) +
    labs(x = "Order q", y = "Sample completeness") +
    theme_bw() +
    theme(legend.position = "bottom", legend.box = "vertical",
          legend.key.width = unit(1.2, "cm"),
          legend.title = element_blank(),
          legend.margin = margin(0, 0, 0, 0),
          legend.box.margin = margin(-10, -10, -5, -10),
          text = element_text(size = 16),
          plot.margin = unit(c(5.5, 5.5, 5.5, 5.5), "pt")) +
    guides(linetype = guide_legend(keywidth = 2.5))
}


# Calculate six classes for Evenness
#
# @param q a integer vector for the order of Hill number
# @param qTD a vector for diversity
# @param S a integer value for species
# @param E a integer value between 1 to 6
# @return a vector for evenness value

even.class = function(q, qTD, S, E.class, pi) {
  tmp = c()
  
  if (E.class == 1)
    tmp = ifelse(q != 1, (1 - qTD^(1 - q)) / (1 - S^(1 - q)), log(qTD) / log(S))
  
  if (E.class == 2)
    tmp = ifelse(q != 1, (1 - qTD^(q - 1)) / (1 - S^(q - 1)), log(qTD) / log(S))
  
  if (E.class == 3)
    tmp = (qTD - 1) / (S - 1)
  
  if (E.class == 4)
    tmp = (1 - 1 / qTD) / (1 - 1 / S)
  
  if (E.class == 5)
    tmp = log(qTD) / log(S)
  
  if (E.class == 6) 
    tmp = sapply(q, function(qi) {
      
      if(qi == 0) {
        
        pi <- pi[pi > 0]
        nu <- abs(pi - (1/S))
        nu <- nu[nu > 0]
        sub <- (sum(log(abs(nu))) / sum(nu > 0) - (log(1 - 1/S) + (1 - S) * log(S)) / S)
        1 - exp(sub)
        
      }else{
        
        pi <- pi[pi > 0]
        1 - (sum(abs(pi - 1/S)^qi) / ((1 - 1/S)^qi + (S - 1) * S^(-qi)))^(1/qi)
      }
    })
  
  return(tmp)
}


# \code{Evenness.profile} Estimation or Empirical of Evenness with order q
#
# @param x a data.frame, a vector, or a list for data.
# @param q a integer vector for the order of Hill number.
# @param datatype a binary choose with 'abundance' or 'incidence_freq'
# @param method a binary calculation method with 'Estimated' or 'Observed'
# @param E.class a integer vector between 1 to 6
# @param C a standardized coverage for calculating evenness index
# @return a list of estimated(Observed) evenness with order q, each list is combined with a matrix

Evenness.profile <- function(x, q, datatype = c("abundance","incidence_freq"), method, E.class, C = NULL) {
  
  if (method == "Estimated") {
    
    estqD = estimate3D(x, diversity = 'TD', q, datatype, base = "coverage", level = C, nboot = 0)
    estS = estimate3D(x, diversity = 'TD', 0, datatype, base = "coverage", level = C, nboot = 0)
    
    out = lapply(E.class, function(i) {
      tmp = sapply(1:length(x), function(k) even.class(q, estqD[estqD$Assemblage == names(x)[k], "qTD"], 
                                                       estS[estS$Assemblage == names(x)[k], "qTD"], i, x[[k]]/sum(x[[k]])))
      
      if (inherits(tmp, c("numeric","integer"))) {tmp = t(as.matrix(tmp, nrow = 1))}
      rownames(tmp) = q
      
      tmp
    })
    
  } else if (method == "Observed") {
    
    empqD = ObsAsy3D(x, diversity = 'TD', q = q, datatype = datatype, nboot = 0, method = 'Observed')
    empS = ObsAsy3D(x, diversity = 'TD', q = 0, datatype = datatype, nboot = 0, method = 'Observed')
    
    out = lapply(E.class, function(i) {
      tmp = sapply(1:length(x), function(k) even.class(q, empqD[empqD$Assemblage == names(x)[k], "qTD"], 
                                                       empS[empS$Assemblage == names(x)[k], "qTD"], i, x[[k]]/sum(x[[k]])))
      
      if (inherits(tmp, c("numeric","integer"))) {tmp = t(as.matrix(tmp, nrow = 1))}
      rownames(tmp) = q
      
      tmp
    })
    
  }
  
  names(out) = paste("E", E.class, sep="")
  
  return(out)
}



#' Main function for STEP 4: Assessment of evenness
#'
#' \code{Evenness} computes standardized and observed evenness of order q = 0 to q = 2 in increments of 0.2 (by default) and depicts evenness profiles based on five classes of evenness measures developed 
#' in Chao and Ricotta (2019).
#' 
#' @param data (a) For \code{datatype = "abundance"}, data can be input as a vector of species abundances (for a single assemblage), matrix/data.frame (species by assemblages), or a list of species abundance vectors. \cr
#' (b) For \code{datatype = "incidence_raw"}, data can be input as a list of matrix/data.frame (species by sampling units); data can also be input as a matrix/data.frame by merging all sampling units across assemblages based on species identity; in this case, the number of sampling units (nT, see below) must be input. 
#' @param q a numerical vector specifying the orders of evenness. Default is \code{seq(0, 2, by = 0.2)}.
#' @param datatype data type of input data: individual-based abundance data (\code{datatype = "abundance"}) or species by sampling-units incidence matrix (\code{datatype = "incidence_raw"}) with all entries being 0 (non-detection) or 1 (detection).
#' @param method a binary selection of method with \code{"Estimated"} (evenness is computed under a standardized coverage value) or \code{"Observed"} (evenness is computed for the observed data).\cr
#' @param nboot a positive integer specifying the number of bootstrap replications when assessing sampling uncertainty and constructing confidence intervals. Enter 0 to skip the bootstrap procedures. Default is \code{30}.
#' @param conf a positive number < \code{1} specifying the level of confidence interval. Default is \code{0.95}.
#' @param nT (required only when \code{datatype = "incidence_raw"} and input data is matrix/data.frame) a vector of nonnegative integers specifying the number of sampling units in each assemblage. If assemblage names are not specified, then assemblages are automatically named as "Assemblage1", "Assemblage2",..., etc. 
#' @param E.class an integer vector between 1 to 5 specifying which class(es) of evenness measures are selected; default is 1:5 (select all five classes). 
#' @param SC (required only when \code{method = "Estimated"}) a standardized coverage value for calculating estimated evenness. If \code{SC = NULL}, then this function computes the diversity estimates for the minimum sample coverage among all samples extrapolated to double reference sizes (Cmax).
#' @return A list of several tables containing estimated (or observed) evenness with order q.\cr
#'         Each tables represents a class of evenness.
#'         \item{Order.q}{the order of evenness}
#'         \item{Evenness}{the computed evenness value of order q.}
#'         \item{s.e.}{standard error of evenness value.}
#'         \item{Even.LCL, Even.UCL}{the bootstrap lower and upper confidence limits for the evenness of order q at the specified level (with a default value of \code{0.95}).}
#'         \item{Assemblage}{the assemblage name.}
#'         \item{Method}{\code{"Estimated"} or \code{"Observed"}.}
#'         \item{SC}{the standardized coverage value under which evenness values are computed (only for \code{method = "Estimated"})}
#'         
#' 
#' @examples
#' ## Evenness for abundance data
#' # The observed evenness values for abundance data
#' data(Data_spider)
#' Even_out1_obs <- Evenness(data = Data_spider, datatype = "abundance", 
#'                           method = "Observed", E.class = 1:5)
#' Even_out1_obs
#' 
#' \donttest{
#' # Estimated evenness for abundance data with default sample coverage value
#' data(Data_spider)
#' Even_out1_est <- Evenness(data = Data_spider, datatype = "abundance", 
#'                           method = "Estimated", SC = NULL, E.class = 1:5)
#' Even_out1_est
#' }
#' 
#' ## Evenness for incidence raw data
#' # The observed evenness values for incidence raw data
#' data(Data_woody_plant)
#' Even_out2_obs <- Evenness(data = Data_woody_plant, datatype = "incidence_raw", 
#'                           method = "Observed", E.class = 1:5)
#' Even_out2_obs
#' 
#' \donttest{
#' # Estimated evenness for incidence data with user's specified coverage value of 0.98
#' data(Data_woody_plant)
#' Even_out2_est <- Evenness(data = Data_woody_plant, datatype = "incidence_raw", 
#'                           method = "Estimated", SC = 0.98, E.class = 1:5)
#' Even_out2_est
#' }
#' 
#' @references
#' Chao, A. and Ricotta, C. (2019). Quantifying evenness and linking it to diversity, beta diversity, and similarity. Ecology, 100(12), e02852.
#' @export

Evenness <- function (data, q = seq(0, 2, 0.2), datatype = "abundance", method = "Estimated",
                      nboot = 30, conf = 0.95, nT = NULL, E.class = 1:5, SC = NULL)
{
  TYPE <- c("abundance", "incidence", "incidence_freq", "incidence_raw")
  if (is.na(pmatch(datatype, TYPE)))
    stop("invalid datatype")
  if (pmatch(datatype, TYPE) == -1)
    stop("ambiguous datatype")
  datatype <- match.arg(datatype, TYPE)
  
  if(datatype == "incidence") stop("Completeness can only accept datatype = 'incidence_freq' or 'datatype = incidence_raw'.") 
  
  kind <- c("Estimated", "Observed")
  if (length(method) > 1)
    stop("only one calculation method")
  if (is.na(pmatch(method, kind)))
    stop("invalid method")
  if (pmatch(method, kind) == -1)
    stop("ambiguous method")
  
  class <- c(1:6)
  if (sum(E.class %in% class) != length(E.class))
    stop("invalid E.class")
  
  if (datatype == "incidence_raw") {
    
    if (length(data) != 1) {
      data = iNEXT.3D:::as.incfreq(data, nT = nT)
      datatype = "incidence_freq"
    }
    
    if (length(data) == 1) {
      name = names(data)
      data = list(iNEXT.3D:::as.incfreq(data, nT = nT))
      names(data) = name
      datatype = "incidence_freq"
    }
  }
  
  if (inherits(data, c("numeric", "integer"))) {
    data <- list("Assemblage_1" = data)
  }
  
  if (inherits(data, c("data.frame", "matrix"))) {
    
    datalist <- lapply(1:ncol(data), function(i) data[, i])
    
    if (is.null(colnames(data)))
      names(datalist) <- paste0("Assemblage_", 1:ncol(data))
    else names(datalist) <- colnames(data)
    
    data <- datalist
  }
  
  if (inherits(data, "list")) {
    if (is.null(names(data))) names(data) = paste0("Assemblage_", 1:length(data))
  }
  
  if (is.null(SC) == TRUE) {
    if (datatype == "abundance") SC = sapply(data, function(x) iNEXT.3D:::Coverage(x, "abundance", 2*sum(x))) %>% min
    if (datatype == "incidence_freq") SC = sapply(data, function(x) iNEXT.3D:::Coverage(x, "incidence_freq", 2*x[1])) %>% min
  }
  
  
  if (datatype == "abundance") {
    
    qD <- Evenness.profile(data, q, "abundance", method, E.class, SC)
    qD <- map(qD, as.vector)
    
    if (nboot > 1) {
      Prob.hat <- lapply(1:length(data), function(i) iNEXT.3D:::EstiBootComm.Ind(data[[i]]))
      
      Abun.Mat <- lapply(1:length(data), function(i) rmultinom(nboot, sum(data[[i]]), Prob.hat[[i]]))
      
      error = apply( matrix(sapply(1:nboot, function(b) {
        
        dat = lapply(1:length(Abun.Mat),function(j) Abun.Mat[[j]][,b])
        names(dat) = paste("Site", 1:length(dat), sep="")
        
        dat.qD = Evenness.profile(dat, q, "abundance", method, E.class, SC)
        
        unlist(dat.qD)
        }), nrow = length(q) * length(E.class) * length(Abun.Mat))
      , 1, sd, na.rm = TRUE)
      
      error = matrix(error, ncol=length(E.class))
      se = split(error, col(error))
      
    } else {
      se = lapply(1:length(E.class), function(x) NA)
    }
    
    out <- lapply(1:length(E.class), function(k) {
      
      tmp = data.frame(Order.q = rep(q, length(data)), 
                       Evenness = as.vector(qD[[k]]), 
                       s.e. = as.vector(se[[k]]),
                       Even.LCL = as.vector(qD[[k]] - qnorm(1-(1-conf)/2)*se[[k]]), 
                       Even.UCL = as.vector(qD[[k]] + qnorm(1-(1-conf)/2)*se[[k]]),
                       Assemblage = rep(names(data), each=length(q)), 
                       Method = rep( method, length(q)*length(data))
                       )
      tmp$Even.LCL[tmp$Even.LCL < 0] <- 0
      tmp
    })
    
  } else if (datatype == "incidence_freq") {
    
    qD <- Evenness.profile(data, q, "incidence_freq", method, E.class, SC)
    qD <- map(qD, as.vector)
    
    if (nboot > 1) {
      nT <- lapply(1:length(data), function(i) data[[i]][1])
      Prob.hat <- lapply(1:length(data), function(i) iNEXT.3D:::EstiBootComm.Sam(data[[i]]))
      
      Incid.Mat <- lapply(1:length(data), function(i) t(sapply(Prob.hat[[i]], function(p) rbinom(nboot, nT[[i]], p))))
      Incid.Mat <- lapply(1:length(data), function(i) matrix(c(rbind(nT[[i]], Incid.Mat[[i]])), ncol = nboot))
      
      error = apply(  matrix(sapply(1:nboot, function(b) {
        
        dat = lapply(1:length(Incid.Mat),function(j) Incid.Mat[[j]][,b])
        names(dat) = paste("Site", 1:length(dat), sep="")
        
        dat.qD = Evenness.profile(dat, q, "incidence_freq", method, E.class, SC)
        
        unlist(dat.qD)  
        }), nrow = length(q) * length(E.class) * length(Incid.Mat))
      , 1, sd, na.rm = TRUE)
      
      error = matrix(error, ncol=length(E.class))
      se = split(error, col(error))
      
    } else {
      se = lapply(1:length(E.class), function(x) NA)
    }
    
    out <- lapply(1:length(E.class), function(k) {
      
      tmp = data.frame(Order.q = rep(q, length(data)), 
                       Evenness = as.vector(qD[[k]]), 
                       s.e. = as.vector(se[[k]]),
                       Even.LCL = as.vector(qD[[k]] - qnorm(1-(1-conf)/2)*se[[k]]), 
                       Even.UCL = as.vector(qD[[k]] + qnorm(1-(1-conf)/2)*se[[k]]),
                       Assemblage = rep(names(data), each=length(q)), 
                       Method = rep(method, length(q)*length(data))
                       )
      tmp$Even.LCL[tmp$Even.LCL < 0] <- 0
      tmp
    })
    
  }
  
  if (method == "Estimated") {
    out <- lapply(out, function(x) x %>% mutate(SC = SC))
  }
  
  names(out) = paste("E", E.class, sep = "")
  
  return(out)
}



#' ggplot for depicting evenness profiles
#'
#' \code{ggEvenness} is a \code{ggplot2} extension for \code{Evenness} object to plot evenness with order q. 
#'
#' @param output output obtained from the function \code{Evenness}.
#' @return a figure depicting the estimated (or observed) evenness with respect to order q.
#' 
#' 
#' @examples
#' ## Evenness profiles for abundance data
#' # The observed evenness profile for abundance data
#' data(Data_spider)
#' Even_out1_obs <- Evenness(data = Data_spider, datatype = "abundance", 
#'                     method = "Observed", E.class = 1:5)
#' ggEvenness(Even_out1_obs)
#' 
#' \donttest{
#' # The estimated evenness profile for abundance data with default sample coverage value
#' data(Data_spider)
#' Even_out1_est <- Evenness(data = Data_spider, datatype = "abundance", 
#'                           method = "Estimated", SC = NULL, E.class = 1:5)
#' ggEvenness(Even_out1_est)
#' }
#' 
#' ## Evenness profiles for incidence raw data
#' # The observed evenness profile for incidence data
#' data(Data_woody_plant)
#' Even_out2_obs <- Evenness(data = Data_woody_plant, datatype = "incidence_raw", 
#'                     method = "Observed", E.class = 1:5)
#' ggEvenness(Even_out2_obs)
#' 
#' \donttest{
#' # The estimated evenness profile for incidence data with user's specified coverage value of 0.98
#' data(Data_woody_plant)
#' Even_out2_est <- Evenness(data = Data_woody_plant, datatype = "incidence_raw", 
#'                           method = "Estimated", SC = 0.98, E.class = 1:5)
#' ggEvenness(Even_out2_est)
#' }
#' 
#' @export

ggEvenness <- function(output) {
  
  classdata = cbind(do.call(rbind, output),
                    class = rep(names(output), each = nrow(output[[1]])))
  
  # Check if the number of unique 'Assemblage' is 8 or less
  if (length(unique(classdata$Assemblage)) <= 8){
    cbPalette <- rev(c("#999999", "#E69F00", "#56B4E9", "#009E73", 
                       "#330066", "#CC79A7", "#0072B2", "#D55E00"))
  }else{
    # If there are more than 8 assemblages, start with the same predefined color palette
    # Then extend the palette by generating additional colors using the 'ggplotColors' function
    cbPalette <- rev(c("#999999", "#E69F00", "#56B4E9", "#009E73", 
                       "#330066", "#CC79A7", "#0072B2", "#D55E00"))
    cbPalette <- c(cbPalette, ggplotColors(length(unique(classdata$Assemblage))-8))
  }
  
  fig = ggplot(classdata, aes(x = Order.q, y = Evenness, colour = Assemblage)) +
    geom_line(size = 1.2) +
    geom_ribbon(aes(ymin=Even.LCL, ymax=Even.UCL, fill = Assemblage),
                alpha=0.2, linetype=0) +
    scale_colour_manual(values = cbPalette) +
    scale_fill_manual(values = cbPalette) +
    labs(x = "Order q", y = "Evenness") +
    theme_bw() +
    theme(legend.position = "bottom", legend.box = "vertical",
          legend.key.width = unit(1.2, "cm"),
          legend.title = element_blank(),
          legend.margin = margin(0, 0, 0, 0),
          legend.box.margin = margin(-10, -10, -5, -10),
          text = element_text(size = 16),
          plot.margin = unit(c(5.5, 5.5, 5.5, 5.5), "pt")) +
    guides(linetype = guide_legend(keywidth = 2.5))
  
  if (length(output) != 1) fig = fig +
    facet_wrap(. ~ class) +
    theme(strip.text.x = element_text(size = 12, colour = "purple", face = "bold"))
  
  return(fig)
}


# Generate Color Palette for ggplot2
#
# This function creates a color palette suitable for ggplot2 visualizations by evenly spacing colors in the HCL color space. The function ensures that the colors are well-distributed and visually distinct, making it ideal for categorical data where each category needs to be represented by a different color.
#
# @param g An integer indicating the number of distinct colors to generate. This value should be a positive integer, with higher values resulting in a broader range of colors.
# @return A vector of color codes in hexadecimal format, suitable for use in ggplot2 charts and plots. The length of the vector will match the input parameter `g`.
# @examples
# # Generate a palette of 5 distinct colors
# ggplotColors(5)
#
# # Use the generated colors in a ggplot2 chart
# library(ggplot2)
# df <- data.frame(x = 1:5, y = rnorm(5), group = factor(1:5))
# ggplot(df, aes(x, y, color = group)) +
#   geom_point() +
#   scale_color_manual(values = ggplotColors(5))
#
ggplotColors <- function(g){
  d <- 360/g # Calculate the distance between colors in HCL color space
  h <- cumsum(c(15, rep(d,g - 1))) # Create cumulative sums to define hue values
  hcl(h = h, c = 100, l = 65) # Convert HCL values to hexadecimal color codes
}



## ========== no visible global function definition for R CMD check ========== ##
utils::globalVariables(c("Order.q", "Estimate.SC", "Assemblage", "SC.LCL", "SC.UCL", 
                         "Even.LCL", "Even.UCL"
))







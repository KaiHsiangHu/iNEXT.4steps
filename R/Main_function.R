#' iNEXT 4 steps
#'
#' \code{iNEXT.4steps}:\cr
#' A complete (random sampling) biological analysis combined with four parts:\cr
#' Step1: Sample Completeness.\cr
#' Step2: Interpolation and Extrapolation.\cr
#' Step3: Asymptotic diversity.\cr
#' Step4: Evenness.\cr
#' @param data a matrix/data.frame/list/vector of abundances-based/incidences-based species data.\cr
#' Type (1) abundance data:\cr
#' When there are N assemblages, the
#' observed species abundances should be arranged as a species (in rows) by assemblage (in columns) matrix. The first row
#' (including N entries) lists the assemblage labels or site names for the N assemblages.\cr
#' Type (2) incidence data:\cr
#' The data input format for incidence data must be raw detection/non-detection data. That is, data for each community/assemblage
#' consist of a species-by-sampling-unit matrix. Users must first merge multiple-community data by species identity to obtain a pooled
#' list of species; then the rows of the input data refer to this pooled list. \cr
#' @param datatype data type of input data: individual-based abundance data (\code{datatype = "abundance"}),
#' sampling-unit-based incidence frequencies data (\code{datatype = "incidence_freq"}) or species by sampling-units incidence matrix (\code{datatype = "incidence_raw"}).\cr
#' @param tree a structure data for phylogenetic diversity
#' @param reftime a positive value or sequence specifying the reference times for diversity computation. If NULL, then reftime is set to be the tree depth of the phylogenetic tree, which is spanned by all the observed species in the pooled assemblage. Default is NULL.
#' @param distM a distance matrix for functional diversity
#' @param threshold a proposed tau value for functional diversity
#' @param q a integer vector for the order of Hill number\cr
#' (setting except \code{step2}).\cr
#' @param size a vector of nonnegative integers specifying the sample sizes for which diversity estimates will be calculated. If \code{NULL}, the diversity estimates will
#' be calculated for those sample sizes determined by the specified/default \code{endpoint} and \code{knot}. \cr
#' (setting only for \code{step2}).\cr
#' @param endpoint an interger specifying the endpoint for rarefaction and extrapolation range. If \code{NULL}, \code{endpoint} = double of the maximum
#' reference sample size. It will be ignored if \code{size} is given. \cr
#' (setting only for \code{step2}).\cr
#' @param knots an integer specifying the number of knot between 1 and the \code{endpoint}, default is 40.\cr
#' (setting only for \code{step2}).\cr
#' @param conf a positive number < 1 specifying the level of confidence interval, default is 0.95.\cr
#' @param nboot an integer specifying the number of bootstrap replications, default is 30.\cr
#' @param details a logical variable to determine whether do you want to print out the detailed value of 4 plots, default is \code{FALSE}.\cr
#' @import devtools
#' @import ggplot2
#' @import reshape2
#' @import dplyr
#' @import ggpubr
#' @import purrr
#' @import iNEXT
#' @import iNEXTPD2
#' @import FunD
#' @importFrom stats qnorm
#' @importFrom stats rbinom
#' @importFrom stats rmultinom
#' @importFrom stats sd
#' @return a list of two of objects: \cr\cr
#' \code{$summary} individual summary of 4 steps of data. \cr\cr
#' \code{$figure} 5 figures of analysis process. \cr\cr
#' \code{$details} the information for generating \code{figure}. \cr
#' If you nees it, you should key in \code{details = TRUE}. \cr\cr
#' @examples
#' \dontrun{
#' ## Type (1) example for abundance based data (data.frame)
#' ## Ex.1
#' data(Spider)
#' out1 <- iNEXT.4steps(data = Spider, datatype = "abundance")
#' out1
#' ## Type (2) example for incidence based data (list of data.frame)
#' ## Ex.2
#' data(woody_incid)
#' out2 <- iNEXT.4steps(data = woody_incid[,c(1,4)], datatype = "incidence_freq")
#' out2
#' }
#' @references
#' Chao,A., Y.Kubota, D.Zeleny, C.-H.Chiu.
#' Quantifying sample completeness and comparing diversities among assemblages. Ecological Research.
#' @export

iNEXT.4steps <- function(data, datatype = "abundance", tree = NULL, reftime = NULL, distM = NULL, threshold = NULL,
                         q = seq(0, 2, 0.25), size = NULL, endpoint = NULL,
                         knots = 30, conf = 0.95, nboot = 30, details = FALSE) {
  qD = "TD"
  if ((length(data)==1) && (class(data) %in% c("numeric", "integer")))
    stop("Error: Your data does not have enough information.")
  logic = c("TRUE", "FALSE")
  if (is.na(pmatch(details, logic)))
    stop("invalid details setting")
  if (pmatch(details, logic) == -1)
    stop("ambiguous details setting")
  if (sum(c(0,1,2) %in% q) != 3) {q = c(q, c(0,1,2)[(c(0,1,2) %in% q) == FALSE]); q = sort(q)}
  if ((qD == "PhD") && is.null(tree))
    stop("You should input tree data for Phylogenetic diversity.")
  if ((qD == "FunD") && is.null(distM))
    stop("You should input distance data for Functional diversity.")

  plot.names = c("(a) Sample completeness profiles",
                 "(b) Size-based rarefaction/extrapolation",
                 "(c) Asymptotic and empirical diversity profiles",
                 "(d) Coverage-based rarefaction/extrapolation",
                 "(e) Evenness profiles")
  table.names = c("STEP1. Sample completeness profiles",
                  "STEP2. Asymptotic analysis",
                  "STEP3. Non-asymptotic coverage-based rarefaction and extrapolation analysis",
                  "STEP4. Evenness among species abundances")
  ## SC ##
  SC.table <- SC(data, q=q, datatype, nboot, 0.95)

  ## iNEXT ##
  if (qD == "TD") {
    inextTD.table <- iNEXT(data, q=c(0, 1, 2), datatype, size, endpoint, knots, se=TRUE, conf, nboot)
  } else if (qD == "PhD") {
    inextPD.table <- iNEXTPD(data, datatype=datatype, tree=tree, q=c(0, 1, 2), reftime=NULL,
                        endpoint=endpoint, knots=knots, size=size, nboot=nboot, conf=conf)
  } else if (qD == "FunD") {
    # inextFD = iNEXTFD(data, distM, "abundance", q=c(0,1,2), nboot=0)
    # m = inextFD$inext[inextFD$inext$Order.q==0,]$m
    # m = lapply(1:ncol(data), function(i) m[(length(m)/3*(i-1)+1):(i*length(m)/3)])
    # RE.table = FunD:::AUCtable_iNextFD(data, distM, datatype="abundance", m=m, nboot=nboot)
  }

  ## qD ##
  if (qD == "TD") {
    asy.TD <- TdAsy(data, q=q, datatype=datatype, nboot=nboot, conf=conf)
    obs.TD <- TdObs(data, q=q, datatype=datatype, nboot=0, conf=conf)

    qTD.table = rbind(asy.TD, obs.TD)
    qTD.table$s.e. = (qTD.table$qD.UCL-qTD.table$qD)/qnorm(1-(1-conf)/2)
  } else if (qD == "PhD") {
    asy.PD <- PhdAsy(data, datatype=datatype, tree=tree, q=q, reftime=reftime, type="PD", conf=conf, nboot=nboot)
    obs.PD <- PhdObs(data, datatype=datatype, tree=tree, q=q, reftime=reftime, type="PD", conf=conf, nboot=0)

    qPD.table = rbind(asy.PD, obs.PD)
    qPD.table$s.e. = (qPD.table$qPD.UCL-qPD.table$qPD)/qnorm(1-(1-conf)/2)
  } else if (qD == "FunD") {

    # estAUC = FunD:::AUCtable_est(data, distM, q=q, datatype=datatype, nboot=nboot)
    # empAUC = FunD:::AUCtable_mle(data, distM, q=q, datatype=datatype, nboot=0)
    # asy.table = cbind(rbind(data.frame(estAUC), data.frame(empAUC)),
    #                 method = rep(c("Estimated", "Empirical"), each=nrow(estAUC)))
    # colnames(asy.table) = c("Site", "order", "qD", "qD.LCL", "qD.UCL", "method")
    # asy.table$method = factor(asy.table$method, levels=c("Estimated", "Empirical"))
    # asy.table$s.e. = (asy.table$qD.UCL-asy.table$qD)/qnorm(1-(1-conf)/2)
  }

  even.table <- Evenness(data, q=q, datatype, "Estimated", nboot, conf, E.class=3)
  Cmax = even.table[1]

  if (length(unique((even.table[[2]]$Community)))>1) {
    if (qD=="TD") {
      inextTD.table$iNextEst$size_based$Assemblage = factor(inextTD.table$iNextEst$size_based$Assemblage)
      level = levels(inextTD.table$iNextEst$size_based$Assemblage)

      SC.table$Community = factor(SC.table$Community, level)
      qTD.table$Site = factor(qTD.table$Site, level)
      even.table[[2]]$Community = factor(even.table[[2]]$Community, level)

    } else if (qD=="PhD") {

      # RE.table$inext$site = as.factor(RE.table$inext$site)
      # level = levels(RE.table$inext$site)
      #
      # SC.table$Community = factor(SC.table$Community, level)
      # asy.table$Site = factor(asy.table$Site, level)
      # even.table[[2]]$Community = factor(even.table[[2]]$Community, level)

    } else if (qD=="FunD") {
      # RE.table$site = as.factor(RE.table$site)
      # level = levels(RE.table$site)
      #
      # SC.table$Community = factor(SC.table$Community, level)
      # asy.table$Site = factor(asy.table$Site, level)
      # even.table[[2]]$Community = factor(even.table[[2]]$Community, level)

    }
  }

  ## 5 figures ##
  if (length(unique(SC.table$Community)) <= 8) {

    SC.plot <- ggSC(SC.table) +
      labs(title=plot.names[1]) +
      theme(text=element_text(size=12),
            plot.margin = unit(c(5.5,5.5,5.5,5.5), "pt"),
            plot.title = element_text(size=11, colour='blue', face="bold", hjust=0))

    if (qD=="TD") {
      size.RE.plot <- ggiNEXT(inextTD.table, type=1, facet.var="Order.q", color.var="Order.q") +
        labs(title=plot.names[2]) +
        theme(text=element_text(size=12),
              plot.margin = unit(c(5.5,5.5,5.5,5.5), "pt"),
              legend.margin = margin(0, 0, 0, 0),
              plot.title = element_text(size=11, colour='blue', face="bold", hjust=0))
      size.plot <- ggplot_build(size.RE.plot + guides(color=FALSE, fill=FALSE, shape=FALSE))
      size.plot$data[[1]]$size <- 3
      size.plot <- ggplot_gtable(size.plot)
    } else if (qD=="PhD") {
      cbPalette <- rev(c("#999999", "#E69F00", "#56B4E9", "#009E73",
                         "#330066", "#CC79A7", "#0072B2", "#D55E00"))

      size.RE.plot = ggplot(RE.table$inext, aes(x = m, y = qPD, color = site)) +
        geom_line(size = 1.2, aes(x = m, y = qPD, color = site,
                                  linetype = lty)) +
        scale_colour_manual(values = cbPalette) +
        geom_ribbon(aes(ymin = qPD.LCL, ymax = qPD.UCL, fill = site),
                    alpha = 0.3, colour = NA) +
        scale_fill_manual(values = cbPalette) +
        geom_point(size = 3, data = subset(RE.table$inext, method ==
                                             "Observed")) +
        xlab("Number of individuals") + ylab("Phylogenetic diversity") +
        labs(title=plot.names[2]) +
        scale_linetype_manual(values = c("solid", "dashed"),
                              name = "lty", breaks = c("Rarefaction", "Extrapolation"),
                              labels = c("Interpolation", "Extrapolation")) +
        theme(legend.position = "bottom", legend.box = "vertical",
              legend.key.width = unit(1.2,"cm"),
              # plot.margin = unit(c(1.5,0.3,1.2,0.3), "lines"),
              legend.title = element_blank(),
              legend.margin = margin(0,0,0,0),
              legend.box.margin = margin(-10,-10,-5,-10),
              text = element_text(size=8),
              plot.margin = unit(c(5.5,5.5,5.5,5.5), "pt"),
              plot.title = element_text(size=11, colour='blue', face="bold", hjust=0)) +
        guides(linetype = guide_legend(keywidth = 2.5)) +
        facet_wrap(.~order)

      size.plot <- ggplot_build(size.RE.plot + guides(color=FALSE, fill=FALSE, shape=FALSE))
      size.plot <- ggplot_gtable(size.plot)
    } else if (qD=="FunD") {

    }

    if (qD=="TD") {
      cover.RE.plot <- ggiNEXT(inextTD.table, type=3, facet.var="Order.q", color.var="Order.q") +
        labs(title=plot.names[4]) +
        theme(text=element_text(size=12),
              plot.margin = unit(c(5.5,5.5,5.5,5.5), "pt"),
              legend.margin = margin(0, 0, 0, 0),
              plot.title = element_text(size=11, colour='blue', face="bold", hjust=0))
      cover.plot <- ggplot_build(cover.RE.plot + guides(color=FALSE, fill=FALSE, shape=FALSE))
      cover.plot$data[[1]]$size <- 3
      cover.plot <- ggplot_gtable(cover.plot)

    } else if (qD=="PhD") {
      cbPalette <- rev(c("#999999", "#E69F00", "#56B4E9", "#009E73",
                         "#330066", "#CC79A7", "#0072B2", "#D55E00"))

      cover.RE.plot = ggplot(RE.table$inext, aes(x = SC, y = qPD, color = site)) +
        geom_line(size = 1.2, aes(x = SC, y = qPD, color = site,
                                  linetype = lty)) +
        scale_colour_manual(values = cbPalette) +
        geom_ribbon(aes(ymin = qPD.LCL, ymax = qPD.UCL, fill = site),
                    alpha = 0.3, colour = NA) +
        scale_fill_manual(values = cbPalette) +
        geom_point(size = 3, data = subset(RE.table$inext, method ==
                                             "Observed")) +
        xlab("Sample Coverage") + ylab("Functional diversity") +
        labs(title=plot.names[4]) +
        scale_linetype_manual(values = c("solid", "dashed"),
                              name = "lty", breaks = c("Rarefaction", "Extrapolation"),
                              labels = c("Interpolation", "Extrapolation")) +
        theme(legend.position = "bottom", legend.box = "vertical",
              legend.key.width = unit(1.2,"cm"),
              # plot.margin = unit(c(1.5,0.3,1.2,0.3), "lines"),
              legend.title = element_blank(),
              legend.margin = margin(0,0,0,0),
              legend.box.margin = margin(-10,-10,-5,-10),
              text = element_text(size=8),
              plot.margin = unit(c(5.5,5.5,5.5,5.5), "pt"),
              plot.title = element_text(size=11, colour='blue', face="bold", hjust=0)) +
        guides(linetype = guide_legend(keywidth = 2.5)) +
        facet_wrap(.~order)

      cover.plot <- ggplot_build(cover.RE.plot + guides(color=FALSE, fill=FALSE, shape=FALSE))
      cover.plot <- ggplot_gtable(cover.plot)
    } else if (qD == "FunD") {

    }

    if (qD=="TD") {
      qTD.table$method = factor(qTD.table$method, level=c("Estimated", "Empirical"))
      asy.plot <- ggtqplotD(qTD.table) +
        labs(title=plot.names[3]) +
        theme(text=element_text(size=12),
              plot.margin = unit(c(5.5,5.5,5.5,5.5), "pt"),
              plot.title = element_text(size=11, colour='blue', face="bold", hjust=0))
    } else if (qD=="PhD") {
      asy.plot <- ggAsymDiv(asy.table) +
        labs(x = "Order q", y = "Phylogenetic Asymptotic") +
        labs(title=plot.names[3]) +
        theme(legend.position = "bottom", legend.box = "vertical",
              legend.key.width = unit(1.2,"cm"),
              # plot.margin = unit(c(1.5,0.3,1.2,0.3), "lines"),
              legend.title = element_blank(),
              legend.margin = margin(0,0,0,0),
              legend.box.margin = margin(-10,-10,-5,-10),
              text = element_text(size=8),
              plot.margin = unit(c(5.5,5.5,5.5,5.5), "pt"),
              plot.title = element_text(size=11, colour='blue', face="bold", hjust=0)
        )
    } else if (qD=="FunD") {
      asy.plot <- ggAsymDiv(asy.table) +
        labs(x = "Order q", y = "Functional Asymptotic") +
        labs(title=plot.names[3]) +
        theme(legend.position = "bottom", legend.box = "vertical",
              legend.key.width = unit(1.2,"cm"),
              # plot.margin = unit(c(1.5,0.3,1.2,0.3), "lines"),
              legend.title = element_blank(),
              legend.margin = margin(0,0,0,0),
              legend.box.margin = margin(-10,-10,-5,-10),
              text = element_text(size=8),
              plot.margin = unit(c(5.5,5.5,5.5,5.5), "pt"),
              plot.title = element_text(size=11, colour='blue', face="bold", hjust=0)
        )
    }

    even.plot <- ggEven(even.table) +
      labs(title=plot.names[5]) +
      theme(text=element_text(size=12),
            plot.margin = unit(c(5.5,5.5,5.5,5.5), "pt"),
            plot.title = element_text(size=11, colour='blue', face="bold", hjust=0)) +
      guides(linetype = FALSE)

    if (qD=="PhD") {
      SC.plot = SC.plot + theme(text=element_text(size=8))
      even.plot = even.plot + theme(text=element_text(size=8))}

    legend.p = get_legend(SC.plot + theme(legend.direction = "vertical"))
    steps.plot = ggarrange(SC.plot + guides(color=FALSE, fill=FALSE),
                           size.plot,
                           asy.plot + guides(color=FALSE, fill=FALSE),
                           cover.plot,
                           even.plot + guides(color=FALSE, fill=FALSE),
                           legend.p, nrow=3, ncol=2
    )
  } else { warning("The number of communities exceed eight. We don't show the figures.") }

  if (qD=="TD") {
    estD = estimateD(data, q=c(0,1,2), datatype, base="coverage", level=NULL, nboot=0)
  } else if (qD=="PhD") {
    estD = EstimatePD(data, tree, datatype, q=c(0,1,2), level=NULL, nboot=0)[[1]]
    colnames(estD)[6:8] = c("qD", "qD.LCL", "qD.UCL")
  } else if (qD=="FunD") {
    estD = FunD:::EstimateFD(data, distM, datatype, q=c(0,1,2), nboot=0)
    colnames(estD)[6:8] = c("qD", "qD.LCL", "qD.UCL")
  }
  ##  Outpue_summary ##
  summary = list(summary.deal(SC.table, 1),
                 summary.deal(qTD.table, 2),
                 summary.deal(estD, 3),
                 summary.deal(even.table, 4, estD)
  )
  names(summary) = table.names

  ##  Output ##

  if (details==FALSE) {

    if (length(unique(SC.table$Community)) <= 8) {
      ans <- list(summary = summary,
                  figure = list(SC.plot, size.RE.plot, asy.plot,
                                cover.RE.plot, even.plot, steps.plot))
    } else { ans <- list(summary = summary) }

  } else if (details==TRUE) {
    tab = list("Sample Completeness" = SC.table, "iNEXT" = inextTD.table$iNextEst,
               "Asymptotic Diversity" = qTD.table, "Evenness" = even.table)

    if (length(unique(SC.table$Community)) <= 8) {
      ans <- list(summary = summary,
                  figure = list(SC.plot, size.RE.plot, asy.plot,
                                cover.RE.plot, even.plot, steps.plot),
                  details = tab)
    } else { ans <- list(summary = summary, details = tab)}

  }
  return(ans)
}


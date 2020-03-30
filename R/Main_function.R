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
#' @param size a vector of nonnegative integers specifying the sample sizes for which diversity estimates will be calculated. If \code{NULL}, the diversity estimates will
#' be calculated for those sample sizes determined by the specified/default \code{endpoint} and \code{knot}. \cr
#' (setting only for \code{step2}).\cr
#' @param endpoint an interger specifying the endpoint for rarefaction and extrapolation range. If \code{NULL}, \code{endpoint} = double of the maximum
#' reference sample size. It will be ignored if \code{size} is given. \cr
#' (setting only for \code{step2}).\cr
#' @param knots an integer specifying the number of knot between 1 and the \code{endpoint}, default is 40.\cr
#' (setting only for \code{step2}).\cr
#' @param se a logical variable to calculate the bootstrap standard error and confidence interval of a level specified by conf, default is \code{TRUE}.\cr
#' (setting only for \code{step2}).\cr
#' @param conf a positive number < 1 specifying the level of confidence interval, default is 0.95.\cr
#' @param nboot an integer specifying the number of bootstrap replications, default is 30.\cr
#' @param details a logical variable to determine whether do you want to print out the detailed value of 4 plots, default is \code{FALSE}.\cr
#' @import devtools
#' @import iNEXT
#' @import ggplot2
#' @import reshape2
#' @import dplyr
#' @import ggpubr
#' @import purrr
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
#' Chao,A., Y.Kubota, D.Zeleny, C.-H.Chiu, C.-F.Li, B.Kusumoto, M.Yasuhara, S.Thorn, C.-L.Wei, M.J.Costello, and R.K.olwell(2020).
#' Quantifying sample completeness and comparing diversities among assemblages. Ecological Research.
#' @export

iNEXT.4steps <- function(data, datatype="abundance", size=NULL, endpoint=NULL,
                         knots=30, se=TRUE, conf=0.95, nboot=30, details=FALSE) {
  logic = c("TRUE", "FALSE")
  if (is.na(pmatch(details, logic)))
    stop("invalid details setting")
  if (pmatch(details, logic) == -1)
    stop("ambiguous details setting")

  plot.names = c("(a) Sample completeness profiles",
                 "(b) Size-based rarefaction/extrapolation",
                 "(c) Asymptotic and empirical diversity profiles",
                 "(d) Coverage-based rarefaction/extrapolation",
                 "(e) Evenness profiles")
  table.names = c("STEP1. Sample completeness profiles",
                  "STEP2. Asymptotic analysis",
                  "STEP3. Non-asymptotic coverage-based rarefaction and extrapolation analysis",
                  "STEP4. Evenness among species abundances")
  ## 4 Details ##
  SC.table <- SC(data, q=seq(0, 2, 0.2), datatype, nboot, conf)
  RE.table <- iNEXT(data, q=c(0, 1, 2), datatype, size, endpoint, knots, se, conf, nboot)
  asy.table <- rbind(iNEXT:::AsymDiv(data, q=seq(0, 2, 0.2), datatype, nboot, conf, method="Estimated"),
                     iNEXT:::AsymDiv(data, q=seq(0, 2, 0.2), datatype, 0, conf, method="Empirical"))
  asy.table$s.e. = (asy.table$qD.UCL-asy.table$qD)/qnorm(1-(1-conf)/2)
  even.table <- Evenness(data, q=seq(0, 2, 0.2), datatype, "Estimated", nboot, conf, E.type=3)
  Cmax = even.table[1]

  if (length(RE.table$DataInfo$site)>1) {
    level = levels(RE.table$DataInfo$site)
    SC.table$Community = factor(SC.table$Community, level)
    asy.table$Site = factor(asy.table$Site, level)
    even.table[[2]]$Community = factor(even.table[[2]]$Community, level)
  }

  ## 5 figures ##
  SC.plot <- ggSC(SC.table) +
    labs(title=plot.names[1]) +
    theme(text=element_text(size=12),
          plot.margin = unit(c(5.5,5.5,5.5,5.5), "pt"),
          plot.title = element_text(size=11, colour='blue', face="bold", hjust=0))

  size.RE.plot <- ggiNEXT(RE.table, type=1, facet.var="order", color.var="order") +
    labs(title=plot.names[2]) +
    theme(text=element_text(size=12),
          plot.margin = unit(c(5.5,5.5,5.5,5.5), "pt"),
          plot.title = element_text(size=11, colour='blue', face="bold", hjust=0))
  size.plot <- ggplot_build(size.RE.plot)
  size.plot$data[[1]]$size <- 3
  size.plot <- ggplot_gtable(size.plot)

  cover.RE.plot <- ggiNEXT(RE.table, type=3, facet.var="order", color.var="order") +
    labs(title=plot.names[4]) +
    theme(text=element_text(size=12),
          plot.margin = unit(c(5.5,5.5,5.5,5.5), "pt"),
          plot.title = element_text(size=11, colour='blue', face="bold", hjust=0))
  cover.plot <- ggplot_build(cover.RE.plot)
  cover.plot$data[[1]]$size <- 3
  cover.plot <- ggplot_gtable(cover.plot)

  asy.plot <- ggAsymDiv(asy.table) +
    labs(title=plot.names[3]) +
    theme(text=element_text(size=12),
          plot.margin = unit(c(5.5,5.5,5.5,5.5), "pt"),
          plot.title = element_text(size=11, colour='blue', face="bold", hjust=0))

  even.plot <- ggEven(even.table)[[1]] +
    labs(title=plot.names[5]) +
    theme(text=element_text(size=12),
          plot.margin = unit(c(5.5,5.5,5.5,5.5), "pt"),
          plot.title = element_text(size=11, colour='blue', face="bold", hjust=0))

  estD = estimateD(data, q=c(0,1,2), datatype, base="coverage", level=NULL, nboot=0)
  ##  Outpue_summary ##
  summary = list(summary.deal(SC.table, 1),
                 summary.deal(asy.table, 2),
                 summary.deal(estD, 3),
                 summary.deal(even.table, 4, estD)
  )
  names(summary) = table.names

  ##  Output_figures ##
  # steps.plot = grid.arrange(SC.plot, size.RE.plot, asy.plot,
  #                          cover.RE.plot, even.plot, nrow=2)
  steps.plot = ggarrange(SC.plot, size.plot, asy.plot,
                         cover.plot, even.plot
  )
  if (details==FALSE) {
    ans <- list(summary = summary,
                figure = list(SC.plot, size.RE.plot, asy.plot,
                              cover.RE.plot, even.plot, steps.plot))
  } else if (details==TRUE) {
    tab = list("Sample Completeness" = SC.table, "iNEXT" = RE.table,
               "Asymptotic Diversity" = asy.table, "Evenness" = even.table)
    ans <- list(summary = summary,
                figure = list(SC.plot, size.RE.plot, asy.plot,
                              cover.RE.plot, even.plot, steps.plot),
                details = tab)
  }
  return(ans)
}

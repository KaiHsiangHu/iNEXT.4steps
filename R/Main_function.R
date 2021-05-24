#' iNEXT 4 steps
#'
#' \code{iNEXT4steps}:\cr
#' A complete (random sampling) biological analysis combined with four parts:\cr
#' Step1: Sample Completeness.\cr
#' Step2: Interpolation and Extrapolation.\cr
#' Step3: Asymptotic diversity.\cr
#' Step4: Evenness.\cr
#' 
#' @param data a matrix/data.frame/list/vector of abundances-based/incidences-based species data.\cr
#' Type (1) abundance data:\cr
#' When there are N assemblages, the
#' observed species abundances should be arranged as a species (in rows) by assemblage (in columns) matrix. The first row
#' (including N entries) lists the assemblage labels or site names for the N assemblages.\cr
#' Type (2) incidence data:\cr
#' The data input format for incidence data must be raw detection/non-detection data. That is, data for each community/assemblage
#' consist of a species-by-sampling-unit matrix. Users must first merge multiple-community data by species identity to obtain a pooled
#' list of species; then the rows of the input data refer to this pooled list. \cr
#' @param diversity selection of diversity type: 'TD' = Taxonomic diversity, 'PD' = Phylogenetic diversity, and 'FD' = Functional diversity.
#' @param datatype data type of input data: individual-based abundance data (\code{datatype = "abundance"}), sampling-unit-based incidence frequencies data (\code{datatype = "incidence_freq"}) or species by sampling-units incidence matrix (\code{datatype = "incidence_raw"}).
#' @param nboot an integer specifying the number of replications. If confidence interval is not desired, use 0 to skip this time-consuming step.
#' @param nT needed only when \code{datatype = "incidence_raw"}, a sequence of named nonnegative integers specifying the number of sampling units in each assemblage. \cr
#' If \code{names(nT) = NULL}, then assemblage are automatically named as "assemblage1", "assemblage2",..., etc. \cr
#' It is necessary when \code{diversity = 'PD'} and \code{datatype = "incidence_raw"}.
#' @param PDtree a phylo object describing the phylogenetic tree in Newick format for all observed species in the pooled assemblage. It is necessary when \code{diversity = 'PD'}.
#' @param PDreftime Select a reference time point for \code{diversity = 'PD'}. Default is NULL.  
#' @param PDtype Select phylogenetic diversity type: \code{PDtype = "PD"} for Chao et al. (2010) phylogenetic diversity and \code{PDtype = "meanPD"} for mean phylogenetic diversity (phylogenetic Hill number). It will be use when \code{diversity = 'PD'}. Default is \code{"meanPD"}.
#' @param FDdistM a pair wise distance matrix for all pairs of observed species in the pooled assemblage. It will be use when \code{diversity = 'FD'}.
#' @param FDtype a binary selection for FD. \code{FDtype = "tau_values"} computes diversity under certain threshold values. \code{FDtype = "AUC"} computes an overall FD which integrates all threshold values between zero and one. Default is \code{"AUC"}.
#' @param FDtau a sequence between 0 and 1 specifying tau. If \code{NULL}, \code{threshold = } dmean. Default is \code{NULL}. It will be use when \code{diversity = 'FD'} and \code{FDtype = "tau_values"}.
#' @param p_row number of row for 4steps figure, default = 2.
#' @param p_col number of column for 4steps figure, default = 3.
#' @param details a logical variable to determine whether do you want to print out the detailed value of 4 plots, default is \code{FALSE}.\cr
#' 
#' @import devtools
#' @import ggplot2
#' @import reshape2
#' @import dplyr
#' @import ggpubr
#' @import purrr
#' @import ape
#' @import tidytree
#' @import phyclust
#' @import iNEXT.3D
#' @importFrom stats qnorm
#' @importFrom stats rbinom
#' @importFrom stats rmultinom
#' @importFrom stats sd
#' @return a list of two of objects: \cr\cr
#' \code{$summary} individual summary of 4 steps of data. \cr\cr
#' \code{$figure} 5 figures of analysis process. \cr\cr
#' \code{$details} the information for generating \code{figure}. \cr
#' If you need it, you should key in \code{details = TRUE}. \cr\cr
#' 
#' @examples
#' \dontrun{
#' ## Type (1) example for abundance based data (data.frame)
#' ## Ex.1
#' data(Spider)
#' out1 <- iNEXT4steps(data = Spider, diversity = "TD", datatype = "abundance")
#' out1
#' 
#' ## Ex.2
#' data(brazil)
#' data(brazil_tree)
#' out2 <- iNEXT4steps(data = brazil, diversity = "PD", datatype = "abundance", nboot = 0, PDtree = tree)
#' out2
#' 
#' ## Ex.3
#' data(brazil)
#' data(brazil_distM)
#' out3 <- iNEXT4steps(data = brazil, diversity = "FD", datatype = "abundance", nboot = 0, FDdistM = distM, FDtype = 'tau_values')
#' out3
#' 
#' ## Type (2) example for incidence based data (list of data.frame)
#' ## Ex.1
#' data(woody_incid)
#' out <- iNEXT4steps(data = woody_incid[,c(1,4)], diversity = "TD", datatype = "incidence_freq")
#' out
#' 
#' }
#' 
#' @references
#' Chao,A., Y.Kubota, D.Zeleny, C.-H.Chiu.
#' Quantifying sample completeness and comparing diversities among assemblages. Ecological Research.
#' @export

iNEXT4steps <- function(data, diversity = c("TD", "PD", "FD"), datatype = "abundance", nT = NULL,
                         PDtree = NULL, PDreftime = NULL, PDtype = 'meanPD', FDdistM = NULL, FDtype = 'AUC', FDtau = NULL,
                         nboot = 30, p_row = 2, p_col = 3, details = FALSE) {
  q = seq(0, 2, 0.25)
  if ((length(data) == 1) && (class(data) %in% c("numeric", "integer")))
    stop("Error: Your data does not have enough information.")
  
  logic = c("TRUE", "FALSE")
  if (is.na(pmatch(details, logic)))
    stop("invalid details setting")
  if (pmatch(details, logic) == -1)
    stop("ambiguous details setting")
  if ((diversity == "PD") && is.null(PDtree))
    stop("You should input tree data for Phylogenetic diversity.")
  if ((diversity == "FD") && is.null(FDdistM))
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
  SC.table <- SC(data, q = q, datatype, nboot, 0.95)

  ## iNEXT ##
  iNEXT.table <- iNEXT3D(data, diversity = diversity, q = c(0, 1, 2), datatype = datatype, nboot = nboot, nT = nT, 
                         PDtree = PDtree, PDreftime = PDreftime, PDtype = PDtype, FDdistM = FDdistM, FDtype = FDtype, FDtau = FDtau)
  
  ## Asymptotic ##
  qD.table <- rbind(asy3D(data, diversity = diversity, q = q, datatype = datatype, nboot = nboot, nT = nT, PDtree = PDtree, 
                          PDreftime = PDreftime, PDtype = PDtype, FDdistM = FDdistM, FDtype = FDtype, FDtau = FDtau),
                    obs3D(data, diversity = diversity, q = q, datatype = datatype, nboot = nboot, nT = nT, PDtree = PDtree, 
                          PDreftime = PDreftime, PDtype = PDtype, FDdistM = FDdistM, FDtype = FDtype, FDtau = FDtau))
  qD.table$s.e. = (qD.table[,4] - qD.table[,2])/qnorm(1 - (1 - 0.95)/2)
  
  ## Evenness ##
  Even.table <- Evenness(data, q = q, datatype, "Estimated", nboot, E.class = 3)
  Cmax = Even.table[1]

  if (length(unique((Even.table[[2]]$Community)))>1) {
    iNEXT.table[[2]]$size_based$Assemblage = factor(iNEXT.table[[2]]$size_based$Assemblage)
    level = levels(iNEXT.table[[2]]$size_based$Assemblage)
    
    SC.table$Assemblage = factor(SC.table$Assemblage, level)
    qD.table$Assemblage = factor(qD.table$Assemblage, level)
    Even.table[[2]]$Assemblage = factor(Even.table[[2]]$Assemblage, level)
  }

  ## 5 figures ##
  if (length(unique(SC.table$Community)) <= 8) {
    SC.plot <- ggSC(SC.table) +
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

    asy.plot <- ggasy3D(qD.table) +
      labs(title = plot.names[3]) +
      theme(text = element_text(size = 12),
            plot.margin = unit(c(5.5, 5.5, 5.5, 5.5), "pt"),
            plot.title = element_text(size = 11, colour = 'blue', face = "bold", hjust = 0))

    even.plot <- ggEven(Even.table) +
      labs(title = plot.names[5]) +
      theme(text = element_text(size = 12),
            plot.margin = unit(c(5.5, 5.5, 5.5, 5.5), "pt"),
            plot.title = element_text(size = 11, colour = 'blue', face = "bold", hjust = 0))

    legend.p = get_legend(SC.plot + theme(legend.direction = "vertical"))
    steps.plot = ggarrange(SC.plot + guides(color = FALSE, fill = FALSE),
                           size.RE.plot,
                           asy.plot + guides(color = FALSE, fill = FALSE),
                           cover.RE.plot,
                           even.plot + guides(color = FALSE, fill = FALSE),
                           legend.p, nrow = p_row, ncol = p_col
    )
  } else { warning("The number of communities exceed eight. We don't show the figures.") }
  
  estD = estimate3D(data, diversity = 'TD', q = c(0, 1, 2), datatype, base = "coverage", level = NULL, nboot = 0, nT = nT)
  est3D = estimate3D(data, diversity = diversity, q = c(0, 1, 2), datatype, base = "coverage", level = NULL, nboot = 0, nT = nT, PDtree = PDtree, PDtype = PDtype, FDdistM = FDdistM, FDtype = FDtype)
  
  if (diversity == 'FD' & FDtype == 'AUC') summary_step2 = iNEXT.table[[2]] else summary_step2 = iNEXT.table[[3]]
  ##  Outpue_summary ##
  summary = list(summary.deal(SC.table, 1),
                 summary_step2,
                 summary.deal(est3D, 3),
                 summary.deal(Even.table, 4, estD)
  )
  names(summary) = table.names
  
  ##  Output ##
  if (details == FALSE) {

    if (length(unique(SC.table$Community)) <= 8) {
      ans <- list(summary = summary,
                  figure = list(SC.plot, size.RE.plot, asy.plot,
                                cover.RE.plot, even.plot, steps.plot))
    } else { ans <- list(summary = summary) }

  } else if (details == TRUE) {
    if (diversity == 'FD' & FDtype == 'AUC')  RE = iNEXT.table[[1]] else RE = iNEXT.table[[2]]
    tab = list("Sample Completeness" = SC.table, "iNEXT" = RE,
               "Asymptotic Diversity" = qD.table, "Evenness" = Even.table)

    if (length(unique(SC.table$Community)) <= 8) {
      ans <- list(summary = summary,
                  figure = list(SC.plot, size.RE.plot, asy.plot,
                                cover.RE.plot, even.plot, steps.plot),
                  details = tab)
    } else { ans <- list(summary = summary, details = tab)}

  }
  return(ans)
}


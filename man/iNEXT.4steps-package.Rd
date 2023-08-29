\name{iNEXT.4steps-package}
\alias{iNEXT.4steps-package}
\alias{iNEXT.4steps}
\docType{package}
\title{
iNterpolation and EXTrapolation for four steps of biodiversity
}
\description{
\code{iNEXT.4steps} (iNterpolation and EXTrapolation for four steps of biodiversity) is an original R package which provide an easy complete biological analysis computation. In Chao et al. (2020) paper, they proposed a complete biological analysis process:

* \code{STEP1}. Sample completeness profiles. 

* \code{STEP2}. Size-based rarefaction and extrapolation analysis and the asymptotic diversity profile. 

* \code{STEP3}. Non-asymptotic coverage-based rarefaction and extrapolation analysis. 

* \code{STEP4}. Evenness among species abundances. 


These are the foundation of \code{iNEXT.4steps}. Here we will introduce functions about \code{STEP1} and \code{STEP4}, particularly. If you want to grasp the functions of \code{STEP2}, \code{STEP3} individually, then please search the related paper Hsieh et al. (2016) or the package \code{iNEXT.3D} from Anne Chao's github to get more details. \code{iNEXT.3D} contains two major parts:

1. Non-asymptotic diversity:

(1a) Sample-size-based (or size-based) R/E sampling curves: \code{iNEXT3D} computes rarefied and extrapolated \code{taxonomic}, \code{phylogenetic}, or \code{functional} diversity estimates under a specified sample size. This type of diversity curve plots the estimated diversity with respect to sample size.     

(1b) Sample-coverage-based (or coverage-based) R/E sampling curves: \code{iNEXT3D} computes rarefied and extrapolated \code{taxonomic}, \code{phylogenetic}, or \code{functional} diversity estimates under a specified sample coverage. This type of diversity curve plots the estimated diversity with respect to sample coverage. 

2. Asymptotic diversity: asymptotic \code{taxonomic}, \code{phylogenetic} or \code{functional} diversity estimate with q-profile, time-profile, and tau-profile.


In \code{iNEXT.4steps} package, we provide other four functions for users to calculate and visualize the related biological statistics: \code{Completeness} and \code{Evenness} for calculation, and \code{ggCompleteness} and \code{ggEvenness} for visualization.

This package contains five main functions: \cr\cr
1. \code{iNEXT4steps} Computing sample completeness, empirical and asymptotic diversity, size-based and coverage-based standardized diversity, and evenness. Besides, visualizing the output above. \cr\cr
2. \code{Completeness} Compute estimated sample completeness with order q. \cr\cr
3. \code{ggCompleteness} Visualizing the output from the function \code{Completeness}. \cr\cr
4. \code{Evenness} Estimation (Observed) of Evenness with order q \cr\cr
5. \code{ggEvenness} Visualizing the output from the function \code{Evenness}. \cr\cr

An Online version of \code{iNEXT.4steps} is also available for users without an R background: \cr
https://chao.shinyapps.io/iNEXT_4steps/. \cr
In the detailed Online \code{iNEXT.4steps} User's Guide, we illustrate all the running procedures in an easily
accessible way through numerical examples with proper interpretations of portions of the output.
All the data of those illustrative examples are included in this package. \cr\cr
}
\author{
Anne Chao, K.H. Hu \cr

Maintainer: K.H. Hu <river19951224@gmail.com>
}
\references{
Chao, A., Y. Kubota, D. Zeleny, C.-H. Chiu, C.-F. Li, B. Kusumoto, M. Yasuhara, S. Thorn, C.-L. Wei, M. J. Costello, and R. K. Colwell (2020). Quantifying sample completeness and comparing diversities among assemblages. Ecological Research, 35, 292-314. \cr\cr
T.C. Hsieh, K. H. Ma, and Chao, A. (2016). iNEXT: An R package for rarefaction and extrapolation of species diversity (Hill numbers). Methods in Ecology and Evolution, 7, 1451-1456. \cr\cr
}
\keyword{ package }

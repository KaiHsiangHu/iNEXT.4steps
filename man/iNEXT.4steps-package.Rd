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

In addition, iNEXT.4steps also introduce functions about \code{STEP1} and \code{STEP4}, particularly. If user want to grasp the output from \code{STEP2} and \code{STEP3} only, please refer to \code{iNEXT.3D} from CRAN to get more details. 

There are five main functions in \code{iNEXT.4steps}: 

1. \code{iNEXT4steps} computes sample completeness, observed and asymptotic diversity, size-based and coverage-based standardized diversity, and evenness. In addition, \code{iNEXT4steps} visualize the output above. \cr

2. \code{Completeness} computes estimated sample completeness with order q. \cr

3. \code{ggCompleteness} visualizes the output from the function \code{Completeness}. \cr

4. \code{Evenness} computes estimated (observed) evenness with order q \cr

5. \code{ggEvenness} visualizes the output from the function \code{Evenness}. \cr

An Online version of \code{iNEXT.4steps} is also available for users without an R background: \cr
https://chao.shinyapps.io/iNEXT_4steps/. \cr
In the detailed Online \code{iNEXT.4steps} User's Guide, we illustrate all the running procedures in an easily
accessible way through numerical examples with proper interpretations of portions of the output.
All the data of those illustrative examples are included in this package. \cr\cr
}
\author{
Anne Chao, Kai-Hsiang Hu \cr

Maintainer: Anne Chao <chao@stat.nthu.edu.tw>
}
\references{
Chao, A., Y. Kubota, D. Zeleny, C.-H. Chiu, C.-F. Li, B. Kusumoto, M. Yasuhara, S. Thorn, C.-L. Wei, M. J. Costello, and R. K. Colwell (2020). Quantifying sample completeness and comparing diversities among assemblages. Ecological Research, 35, 292-314. \cr\cr
}
\keyword{ package }

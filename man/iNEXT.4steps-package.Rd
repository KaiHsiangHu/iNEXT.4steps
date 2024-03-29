\name{iNEXT.4steps-package}
\alias{iNEXT.4steps-package}
\alias{iNEXT.4steps}
\docType{package}
\title{
Four-step biodiversity analysis based on iNEXT
}
\description{
This package expands \pkg{iNEXT} (Chao et al. 2014)  to include the estimation of sample completeness and evenness under a unified framework of Hill numbers. 
\pkg{iNEXT.4steps} links sample completeness, diversity estimation, interpolation and extrapolation (\pkg{iNEXT}), and evenness in a fully integrated approach. 
An Online version of \pkg{iNEXT.4steps} is also available for users without an R background: 
https://chao.shinyapps.io/iNEXT_4steps/. \cr
The pertinent background for the four-step methodology is provided in Chao et al. (2020). The four-step procedures are described in the following: 

\itemize{
      \item \code{STEP 1}. Assessment of sample completeness profile 
      
      Before performing biodiversity analysis, it is important to first quantify the sample completeness of a biological survey. Chao et al. (2020) generalized the conventional sample completeness to a class of measures parametrized by an order \eqn{q \geq 0}. When \eqn{q = 0}, sample completeness reduces to the conventional measure of completeness, i.e., the ratio of the observed species richness to the true richness (observed plus undetected). When \eqn{q = 1}, the measure reduces to the sample coverage (the proportion of the total number of individuals in the entire assemblage that belong to detected species), a concept original developed by Alan Turing in his cryptographic analysis during WWII.  When \eqn{q = 2}, it represents a generalized sample coverage with each species being proportionally weighted by its squared species abundance (i.e., each individual being proportionally weighted by its species abundance); this measure thus is disproportionally sensitive to highly abundant species. For a general order \eqn{q \geq 0} (not necessarily to be an integer) , the sample completeness of order q quantifies the proportion of the assemblage's individuals belonging to detected species, with each individual being proportionally weighted by the (q-1)th power of its abundance. Sample completeness profile depicts its estimate with respect to order \eqn{q \geq 0}; this profile fully characterizes the sample completeness of a biological survey. \cr

\pkg{iNEXT.4steps} features the estimated profile for all orders of \eqn{q \geq 0} based on the methodology developed in Chao et al. (2020). All estimates are theoretically between 0 and 1. If the estimated sample completeness profile is a horizontal line at the level of unity for all orders of \eqn{q \geq 0}, then the survey is complete, implying there is no undetected diversity. In most applications, the estimated profile increases with order q, revealing the existence of undetected diversity. The sample completeness estimate for \eqn{q = 0} provides an upper bound for the proportion of observed species; its complement represents a lower bound for the proportion of undetected species. This interpretation is mainly because data typically do not contain sufficient information to accurately estimate species richness and only a lower bound of species richness can be well estimated. By contrast, for \eqn{q \geq 1}, when data are not sparse, the sample completeness value for \eqn{q \geq 1} can be very accurately estimated measures. The values for \eqn{q \geq 2} typically are very close to unity, signifying that almost all highly abundant species (for abundance data) or highly frequent species (for incidence data) had been detected in the reference sample. \cr

  \item \code{STEP 2}. Examination of the size-based rarefaction and extrapolation sampling curves, and the asymptotic diversity profile for \eqn{0 \leq q \leq 2}  \cr

  \item \code{(STEP 2a)}. For each dataset, first examine the pattern of the size-based rarefaction and extrapolation sampling curve up to double the reference sample size for \eqn{q = 0,1} and \eqn{2}.  If the curve stays at a fixed level (this often occurs for the measures of \eqn{q = 1} and \eqn{2}), then our asymptotic estimate presented in Step 2b can be used to accurately infer the true diversity of the entire assemblage. Otherwise, our asymptotic diversity estimate represents only a lower bound (this often occurs for the measures of \eqn{q = 0}).  \cr

  \item \code{(STEP 2b)}. When the true diversity can be accurately inferred, the extent of undetected diversity within each dataset is obtained by comparing the estimated asymptotic diversity profile and empirical profile; the difference in diversity between any two assemblages can be evaluated and tested for significance.  \cr

  \item \code{STEP 3}. Examination of non-asymptotic coverage-based rarefaction and extrapolation analysis for orders \eqn{q = 0}, \eqn{1} and \eqn{2} 

When sampling data do not contain sufficient information to accurately infer true diversity, fair comparisons of diversity across multiple assemblages should be made by standardizing the sample coverage (i.e., comparing diversity for a standardized fraction of an assemblage's individuals). This comparison can be done based on seamless integration of coverage-based rarefaction and extrapolation sampling curves up to a maximum coverage (\eqn{C_{max}} = the minimum sample coverage among all samples extrapolated to double reference sizes).  \cr

  \item \code{STEP 4}. Assessment of evenness profiles 

Chao and Ricotta (2019) developed five classes of evenness measures parameterized by an order \eqn{q \geq 0}, the same order that is used to index sample completeness. All classes of evenness measures are functions of diversity and species richness, and all are standardized to the range of [0, 1] to adjust for the effect of differing species richness. Evenness profile depicts evenness estimate with respect to order \eqn{q \geq 0}. Because true species richness typically cannot be accurately estimated, evenness profile typically can only be accurately measured when both diversity and richness are computed at a fixed level of sample coverage up to a maximum coverage \eqn{C_{max}} defined in Step 3. \pkg{iNEXT.4steps} shows, by default, the relevant statistics and plot for only one class of evenness measure (based on the normalized slope of a diversity profile), but all the five classes are optionally featured.  \cr
    }
  
NOTE 1: Sufficient data are required to perform the 4-step analysis. If there are only a few species in users' data, it is likely that data are too sparse to use \pkg{iNEXT.4steps}.  \cr

NOTE 2: The analyses in STEP 2 and STEP 3 are mainly based on package \pkg{iNEXT} available from CRAN. Thus, \pkg{iNEXT.4steps} expands \pkg{iNEXT} to include the estimation of sample completeness and evenness.  \cr

NOTE 3: As with \pkg{iNEXT}, \pkg{iNEXT.4steps} only deals with taxonomic/species diversity. Researchers who are interested in phylogenetic diversity and functional diversity should use package \pkg{iNEXT.3D} available from CRAN and see the relevant paper (Chao et al. 2021) for methodology.  \cr \cr


There are five main functions in \pkg{iNEXT.4steps}: 

1. \code{iNEXT4steps} computes all statistics in the complete 4-step analysis and visualizes the output. It computes sample completeness, observed and asymptotic diversity, size-based and coverage-based standardized diversity, and evenness. 
   \cr

2. \code{Completeness} computes sample completeness estimates of order q = 0 to q = 2 in increments of 0.2 (by default). This function is specifically for users who only require sample completeness estimates. \cr

3. \code{ggCompleteness} visualizes the output obtained from the function \code{Completeness}. \cr

4. \code{Evenness} computes standardized (or observed) evenness of order q = 0 to q = 2 in increments of 0.2 (by default) based on five classes of evenness measures. This function is specifically for users who only require evenness estimates. \cr

5. \code{ggEvenness} visualizes the output obtained from the function \code{Evenness}. 

}
\author{
Anne Chao, Kai-Hsiang Hu \cr

Maintainer: Anne Chao <chao@stat.nthu.edu.tw>
}
\references{
Chao, A., Gotelli, N. G., Hsieh, T. C., Sander, E. L., Ma, K. H., Colwell, R. K. and Ellison, A. M. (2014). Rarefaction and extrapolation with Hill numbers: a framework for sampling and estimation in species biodiversity studies. Ecological Monographs, 84, 45-67. \cr\cr
Chao, A., Henderson, P. A., Chiu, C.-H., Moyes, F., Hu, K.-H., Dornelas, M and Magurran, A. E. (2021). Measuring temporal change in alpha diversity: a framework integrating taxonomic, phylogenetic and functional diversity and the iNEXT.3D standardization. Methods in Ecology and Evolution, 12, 1926-1940. \cr\cr
Chao, A. and Ricotta, C. (2019). Quantifying evenness and linking it to diversity, beta diversity, and similarity. Ecology, 100(12), e02852. \cr\cr
Chao, A., Y. Kubota, D. Zeleny, C.-H. Chiu, C.-F. Li, B. Kusumoto, M. Yasuhara, S. Thorn, C.-L. Wei, M. J. Costello, and R. K. Colwell (2020). Quantifying sample completeness and comparing diversities among assemblages. Ecological Research, 35, 292-314. \cr\cr
}
\keyword{ package }

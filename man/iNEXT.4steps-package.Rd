\name{iNEXT.4steps-package}
\alias{iNEXT.4steps-package}
\alias{iNEXT.4steps}
\docType{package}
\title{
iNterpolation and EXTrapolation for four steps of biodiversity
}
\description{
\code{iNEXT.4steps} (iNterpolation and EXTrapolation for four steps of biodiversity) expands iNEXT (Chao et al. 2014) to include the estimation of sample completeness and evenness under a unified framework of Hill numbers. \code{iNEXT.4steps} links sample completeness, diversity estimation, interpolation and extrapolation (iNEXT), and evenness in a fully integrated approach. The pertinent background for the four-step methodology is provided in Chao et al. (2020). The four-step procedures are described in the following:

* \code{STEP1}. Assessment of sample completeness profile 

The sample completeness of order \eqn{q \geq 0} quantifies the proportion of the assemblage's individuals belonging to detected species, with each individual being proportionally weighted by the (q-1)th power of its abundance. When \eqn{q = 0}, sample completeness reduces to the conventional measure of completeness, i.e., the ratio of the observed species richness to the true richness (observed plus undetected). When q=1, the measure reduces to the sample coverage (the proportion of the total number of individuals in the entire assemblage that belong to detected species), a concept original developed by Alan Turing in his cryptographic analysis during WWII.  Sample completeness profile depicts its estimate with respect to order \eqn{q \geq 0}; this profile fully characterizes the sample completeness of a biological survey.

iNEXT.link features the estimated profile based on the theory in Chao et al. (2020). If the estimated sample completeness profile is a horizontal line at the level of unity for all orders of \eqn{q \geq 0}, then the survey is complete, implying there is no undetected diversity. In most applications, the estimated profile increases with order q, revealing the existence of undetected diversity. The sample completeness estimate for \eqn{q = 0} provides an upper bound for the proportion of observed species; its complement represents a lower bound for the proportion of undetected species. If data are not sparse, then the sample completeness value for \eqn{q = 1} accurately measures the proportion of an assemblage's individuals belonging to detected species; and the values for \eqn{q \geq 2} typically are very close to unity, signifying that almost all highly abundant species (for abundance data) or highly frequent species (for incidence data) had been detected in the reference sample. 

* \code{STEP2}. Examination of the size-based rarefaction and extrapolation sampling curves, and the asymptotic diversity profile for \eqn{0 \leq q \leq 2}. 

\code{(STEP 2a)}. For each dataset, first examine the pattern of the size-based rarefaction and extrapolation sampling curve up to double the reference sample size for \eqn{q = 0}, \eqn{1} and \eqn{2}.  If the curve stays at a fixed level (this often occurs for the measures of \eqn{q = 1} and \eqn{2}), then our asymptotic estimate presented in Step 2b can be used to accurately infer the true diversity of the entire assemblage. Otherwise, our asymptotic diversity estimate represents only a lower bound. 

\code{(STEP 2b)}. When the true diversity can be accurately inferred, the extent of undetected diversity within each dataset is obtained by comparing the estimated asymptotic diversity profile and empirical profile; the difference in diversity between any two assemblages can be evaluated and tested for significance.

* \code{STEP3}. Examination of non-asymptotic coverage-based rarefaction and extrapolation analysis for orders \eqn{q = 0}, \eqn{1} and \eqn{2}. 

When sampling data do not contain sufficient information to accurately infer true diversity, fair comparisons of diversity across multiple assemblages should be made by standardizing the sample coverage (i.e., comparing diversity for a standardized fraction of an assemblage's individuals). This comparison can be done based on seamless integration of coverage-based rarefaction and extrapolation sampling curves up to a maximum coverage (C_max = the minimum sample coverage among all samples extrapolated to double reference sizes). 

** NOTE **: In STEP 2 and STEP 3, \code{iNEXT.4steps} only deals with taxonomic/species diversity. Researchers who are interested in phylogenetic diversity and functional diversity should use iNEXT.3D package from CRAN and see the relevant paper (Chao et al. 2021) for methodology.


* \code{STEP4}. Assessment of evenness profiles 

Chao and Ricotta (2019) developed five classes of evenness measures parameterized by order \eqn{q \geq 0}. All classes of evenness measures are functions of diversity and species richness, and all are standardized to the range of [0, 1] to adjust for the effect of differing species richness. Evenness profile depicts evenness estimate with respect to order \eqn{q \geq 0}. Because true species richness typically cannot be accurately estimated, evenness profile typically can only be accurately measured when both diversity and richness are computed at a fixed level of sample coverage up to a maximum coverage C_max defined in Step 3. \code{iNEXT.4steps} features one class of evenness measure (as well as Pielou's J'), but the other four classes are also provided.


There are five main functions in \code{iNEXT.4steps}: 

1. \code{iNEXT4steps} computes sample completeness, observed and asymptotic diversity, size-based and coverage-based standardized diversity, and evenness. In addition, \code{iNEXT4steps} visualize the output above. \cr

2. \code{Completeness} computes estimated sample completeness with order q. \cr

3. \code{ggCompleteness} visualizes the output from the function \code{Completeness}. \cr

4. \code{Evenness} computes estimated (or observed) evenness with order q \cr

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
Chao, A., Gotelli, N. G., Hsieh, T. C., Sander, E. L., Ma, K. H., Colwell, R. K. and Ellison, A. M. (2014). Rarefaction and extrapolation with Hill numbers: a framework for sampling and estimation in species biodiversity studies. Ecological Monographs, 84, 45-67. \cr\cr
Chao, A. and Ricotta, C. (2019). Quantifying evenness and linking it to diversity, beta diversity, and similarity. Ecology, 100(12), e02852. \cr\cr
Chao, A., Henderson, P. A., Chiu, C.-H., Moyes, F., Hu, K.-H., Dornelas, M and Magurran, A. E. (2021). Measuring temporal change in alpha diversity: a framework integrating taxonomic, phylogenetic and functional diversity and the iNEXT.3D standardization. Methods in Ecology and Evolution, 12, 1926-1940. \cr\cr
}
\keyword{ package }

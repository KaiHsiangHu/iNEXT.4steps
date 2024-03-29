% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MainFun.R
\name{Evenness}
\alias{Evenness}
\title{Evenness main function}
\usage{
Evenness(
  data,
  q = seq(0, 2, 0.2),
  datatype = "abundance",
  method = "Estimated",
  nboot = 50,
  conf = 0.95,
  nT = NULL,
  E.class = 1:5,
  SC = NULL
)
}
\arguments{
\item{data}{(a) For \code{datatype = "abundance"}, data can be input as a vector of species abundances (for a single assemblage), matrix/data.frame (species by assemblages), or a list of species abundance vectors. \cr
(b) For \code{datatype = "incidence_raw"}, data can be input as a list of matrix/data.frame (species by sampling units); data can also be input as a matrix/data.frame by merging all sampling units across assemblages based on species identity; in this case, the number of sampling units (nT, see below) must be input.}

\item{q}{a numerical vector specifying the diversity orders. Default is \code{(0, 0.2, 0.4,...,2)}.}

\item{datatype}{data type of input data: individual-based abundance data (\code{datatype = "abundance"}) or species by sampling-units incidence matrix (\code{datatype = "incidence_raw"}) with all entries being 0 (non-detection) or 1 (detection).}

\item{method}{a binary calculation method with \code{"Estimated"} or \code{"Observed"}.\cr}

\item{nboot}{a positive integer specifying the number of bootstrap replications when assessing sampling uncertainty and constructing confidence intervals. Enter 0 to skip the bootstrap procedures. Default is \code{50}.}

\item{conf}{a positive number < \code{1} specifying the level of confidence interval. Default is \code{0.95}.}

\item{nT}{(required only when \code{datatype = "incidence_raw"} and input data is matrix/data.frame) a vector of nonnegative integers specifying the number of sampling units in each assemblage. If assemblage names are not specified, then assemblages are automatically named as "Assemblage1", "Assemblage2",..., etc.}

\item{E.class}{an integer vector between 1 to 5}

\item{SC}{(required only when method = "Estimated") a standardized coverage for calculating estimated evenness. If \code{SC = NULL}, then this function computes the diversity estimates for the minimum sample coverage among all samples extrapolated to double reference sizes (Cmax).}
}
\value{
A list of several tables containing estimated (or observed) evenness with order q.\cr
        Each tables represents a class of evenness.
        \item{Order.q}{the diversity order of q.}
        \item{Evenness}{the evenness of order q.}
        \item{s.e.}{standard error of evenness.}
        \item{Even.LCL, Even.UCL}{the bootstrap lower and upper confidence limits for the evenness of order q at the specified level (with a default value of \code{0.95}).}
        \item{Assemblage}{the assemblage name.}
        \item{Method}{\code{"Estimated"} or \code{"Observed"}.}
        \item{SC}{the target standardized coverage value. (only when \code{method = "Estimated"})}
}
\description{
\code{Evenness} Estimation (or Observed) of Evenness with order q
}
\details{
R scipts "Evenness" for Chao and Ricotta (2019) Ecology paper.
This R code is for computing Figures 2 of Chao and Ricotta (2019) paper.
}
\examples{
## Type (1) example for abundance based data
# Example 1. observed evenness for abundance data
data(Data_spider)
output1 <- Evenness(data = Data_spider, datatype = "abundance", 
                    method = "Observed", E.class = 1:5)
output1

\donttest{
# Example 2. estimated evenness for abundance data with default SC = Cmax
data(Data_spider)
output2 <- Evenness(data = Data_spider, datatype = "abundance", 
                    method = "Estimated", SC = NULL, E.class = 1:5)
output2
}

## Type (2) example for incidence based data
# Example 3. observed evenness for incidence data
data(Data_woody_plant)
output3 <- Evenness(data = Data_woody_plant, datatype = "incidence_raw", 
                    method = "Observed", E.class = 1:5)
output3

\donttest{
# Example 4. estimated evenness for incidence data with default SC = Cmax
data(Data_woody_plant)
output4 <- Evenness(data = Data_woody_plant, datatype = "incidence_raw", 
                    method = "Estimated", SC = NULL, E.class = 1:5)
output4
}

}
\references{
Chao, A. and Ricotta, C. (2019). Quantifying evenness and linking it to diversity, beta diversity, and similarity. Ecology, 100(12), e02852.
}

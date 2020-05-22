% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Sub_function.R
\name{Evenness}
\alias{Evenness}
\title{Evenness main function}
\usage{
Evenness(
  x,
  q = seq(0, 2, 0.2),
  datatype = "abundance",
  method = "Estimated",
  nboot = 30,
  conf = 0.95,
  E.type = c(1:5)
)
}
\arguments{
\item{x}{a matrix/data.frame/list/vector of abundances-based/incidences-based species data.\cr}

\item{q}{a integer vector of the order of Hill number\cr}

\item{datatype}{data type of input data: individual-based abundance data (\code{datatype = "abundance"}),
sampling-unit-based incidence frequencies data (\code{datatype = "incidence_freq"}) or species by sampling-units incidence matrix (\code{datatype = "incidence_raw"}).\cr}

\item{method}{a binary calculation method with 'Estimated' or 'Empirical'\cr}

\item{nboot}{an integer specifying the number of bootstrap replications, default is 30.\cr}

\item{conf}{a positive number < 1 specifying the level of confidence interval, default is 0.95.\cr}

\item{E.type}{a integer vector between 1 to 6}
}
\value{
A list of estimated(empirical) evenness with order q.\cr
        Different lists represents different classes of Evenness.\cr
        Each list is combined with order.q and sites.\cr
        If "method" is estimated, then fist list will be named "Cmax" which means the
        maximum standardized coverage between all double reference sample size.\cr\cr
\code{$summary} individual summary of 4 steps of data. \cr\cr
}
\description{
\code{Evenness} Estimation (Empirical) of Evenness with order q
}
\details{
R scipts "Evenness" for Chao and Ricotta (2019) Ecology paper.
This R code is for computing Figures 2, 3 and 4 of Chao and Ricotta (2019) paper.
installed and loaded before running the scripts.
}
\examples{
\dontrun{
## Type (1) example for abundance based data (data.frame)
## Ex.1
data(Spider)
out1 <- Evenness(x = Spider, datatype = "abundance")
out1

## Type (2) example for incidence based data (list of data.frame)
## Ex.2
data(woody_incid)
out2 <- Evenness(x = woody_incid[,c(1,4)], datatype = "incidence_freq")
out2
}
}
\references{
Chao,A.and Ricotta,C.(2019).Quantifying evenness and linking it to diversity, beta diversity, and similarity.
}
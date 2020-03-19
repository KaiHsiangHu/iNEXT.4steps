\name{coral}
\alias{coral}
\docType{data}
\title{
Coral incidence-based data of two regions: Coral triangle and Madagascar.
}
\description{
The data for Scleractinia species considered in this example were taken from two databases:
Global Biodiversity Information Facility (GBIF; https://www.gbif.org/) and Ocean Biogeographic Information System (OBIS; http://www.iobis.org/).
Kusumoto et al. (2020) provided a more detailed description of how this dataset was compiled and coral taxonomy was scrutinized.
The global area was first divided into grids of 5 x 5 degrees.
We selected 30 such grids covering the Coral Triangle (Figure 5), which is empirically regarded as a diversity hot spot of corals and other marine organisms (e.g. Briggs 2003; Veron et al. 2009; Asaad et al. 2018).
We also selected the same number of grids covering the Madagascar area (Figure 5) where another diversity peak was detected by a global-scale analysis by Kusumoto et al. (2020).
Each defined area was further divided into sub-grids of 0.01 x 0.01 degrees, and species incidence data were recorded in each sub-grid.
}
\usage{data("coral")}
\format{
  A data frame with 536 species on the following 2 communities.
  \describe{
    \item{\code{Madagascar}}{a numeric vector}
    \item{\code{Coral.Triangle}}{a numeric vector}
  }
}
\source{
Global Biodiversity Information Facility (GBIF; https://www.gbif.org/).
Ocean Biogeographic Information System (OBIS; http://www.iobis.org/).
}
\examples{
data(coral)
\dontrun{
coral.out <- iNEXT.4steps(data = coral, datatype = "incidence_freq")
coral.out
}
}
\keyword{datasets}

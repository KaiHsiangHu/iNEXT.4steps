\name{woody_incid}
\alias{woody_incid}
\docType{data}
\title{
Woody incidence-based data
}
\description{
The datasets considered here are a subset of The National Vegetation Database of Taiwan (AS-TW-001), sampled between 20032007 within the first national vegetation inventory project (Chiou et al. 2009).
Over 3600 vegetation plots, each 20×20-m in area, were set up in various locations in Taiwan, and all woody plant individuals taller than 2 meters were recorded in each plot.
For illustration here, we selected only plots belonging to two vegetation types (according to Li et al. 2013): Pyrenaria-Machilus subtropical winter monsoon forest and Chamaecyparis montane mixed cloud forest, sampled in the northern part of Taiwan (in ecoregions 7 and 8 according to Su 1985).

}
\usage{data("woody_incid")}
\format{
  A data frame with 494 species on the following 4 sites.
  \describe{
    \item{\code{Upper_cloud}}{a numeric vector}
    \item{\code{Lower_cloud}}{a numeric vector}
    \item{\code{Lowland}}{a numeric vector}
    \item{\code{Monsoon}}{a numeric vector}
  }
}
\source{
The National Vegetation Database of Taiwan (AS-TW-001).
}
\examples{
data(woody_incid)
\dontrun{
woody_incid.out <- iNEXT.4steps(data = woody_incid, datatype = "incidence_freq")
woody_incid.out
}
}
\keyword{datasets}

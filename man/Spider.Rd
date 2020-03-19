\name{Spider}
\alias{Spider}
\docType{data}
\title{
Spider abundance-based data
}
\description{
These data were sampled in a mountain forest ecosystem in the Bavarian Forest National Park, Germany (Thorn et al. 2016, 2017).
A total of 12 experimental plots were established in "closed forest" stands (6 plots) and "open forest" stands with naturally occurring gaps and edges (6 plots) to assess the effects of microclimate on communities of epigeal (ground-dwelling) spiders.
Epigeal spiders were sampled over three years with four pitfall traps in each plot, yielding a total of 3171 individuals belonging to 85 species recorded in the pooled habitat.

}
\usage{data("Spider")}
\format{
  A data frame with 85 species on the following 2 sites.
  \describe{
    \item{\code{Open}}{a numeric vector}
    \item{\code{Closed}}{a numeric vector}
  }
}
\source{
A mountain forest ecosystem in the Bavarian Forest National Park, Germany (Thorn et al. 2016, 2017).
}
\examples{
data(Spider)
\dontrun{
Spider.out <- iNEXT.4steps(data = Spider, datatype = "abundance")
Spider.out
}
}
\keyword{datasets}

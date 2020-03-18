\name{Spider}
\alias{Spider}
\docType{data}
\title{
Spider abundance-based data
}
\description{
This data set is combined with two sites. Each column means a generation. Each row is a species.
}
\usage{data("Spider")}
\format{
  A data frame with 85 observations on the following 2 variables.
  \describe{
    \item{\code{Open}}{a numeric vector}
    \item{\code{Closed}}{a numeric vector}
  }
}
\examples{
data(Spider)
\dontrun{
Spider.out <- iNEXT.4steps(data = Spider, datatype = "abundance")
Spider.out
}
}
\keyword{datasets}

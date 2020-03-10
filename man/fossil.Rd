\name{fossil}
\alias{fossil}
\docType{data}
\title{
fossil abundance-based data
}
\description{
This data set is combined with four sites. Each column means a generation. Each row is a species.
You can find several ancients in each generation.
}
\usage{data("fossil")}
\format{
  A data frame with 178 observations on the following 2 variables.
  \describe{
    \item{\code{Miocene}}{a numeric vector}
    \item{\code{Pliocene}}{a numeric vector}
  }
}
\examples{
data(fossil)
\dontrun{
fossil.out <- iNEXT.4steps(data = fossil, datatype = "abundance")
fossil.out
}
}
\keyword{datasets}

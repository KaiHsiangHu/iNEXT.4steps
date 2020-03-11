\name{woody_incid}
\alias{woody_incid}
\docType{data}
\title{
Woody incidence-based data
}
\description{
This data set is combined with four sites. Each column means a site.
The first row is total plots(quardrats). Other rows are species.
}
\usage{data("woody_incid")}
\format{
  A data frame with 494 observations on the following 4 variables.
  \describe{
    \item{\code{Upper_cloud}}{a numeric vector}
    \item{\code{Lower_cloud}}{a numeric vector}
    \item{\code{Lowland}}{a numeric vector}
    \item{\code{Monsoon}}{a numeric vector}
  }
}
\examples{
data(woody_incid)
\dontrun{
woody_incid.out <- iNEXT.4steps(data = woody_incid, datatype = "incidence_freq")
woody_incid.out
}
}
\keyword{datasets}

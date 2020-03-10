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
    \item{\code{upper_cloud}}{a numeric vector}
    \item{\code{lower_cloud}}{a numeric vector}
    \item{\code{lowland}}{a numeric vector}
    \item{\code{monsoon}}{a numeric vector}
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

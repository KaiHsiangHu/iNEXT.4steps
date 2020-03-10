\name{coral}
\alias{coral}
\docType{data}
\title{
Coral incidence-based data of two regions: Coral triangle and Madagascar.
}
\description{
This data set is combined with two sites. Each column means a site.
The first row is total plots(quardrats). Other rows are species.
Coral triangle is a triangle region from Astralia to eastern-south Asia in the Pacific Ocean.
Each plot is constituted with size: 0.01 degree x 0.01 degree.
Madagascar is at the corner of Africa.
The range of Coral triangle is larger then Madagascar.
}
\usage{data("coral")}
\format{
  A data frame with 717 observations on the following 2 variables.
  \describe{
    \item{\code{Madagascar}}{a numeric vector}
    \item{\code{Coral.Triangle}}{a numeric vector}
  }
}
\examples{
data(coral)
\dontrun{
coral.out <- iNEXT.4steps(data = coral, datatype = "incidence_freq")
coral.out
}
}
\keyword{datasets}

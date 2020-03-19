\name{fossil}
\alias{fossil}
\docType{data}
\title{
fossil abundance-based data
}
\description{
Shin et al.(2019) investigated species richness patterns in Neogene fossil marine ostracods from Java, Indonesia.
In their data, a total of 171 species were identified from 9 sediment samples from three different geological ages:
"Middle Miocene" (1 sample), "Late Miocene" (3 samples) and "Pliocene" (5 samples).
They concluded that there was a significance increase in species richness from the Late Miocene to the Pliocene.
For illustration purposes, we pool the data of Middle Miocene and Late Miocene and compare the pooled data (referred to as “Miocene”, 4 samples pooled) with those of the Pliocene (5 samples).

}
\usage{data("fossil")}
\format{
  A data frame with 178 species on the following 2 generations.
  You can find several ancients in each generation.
  \describe{
    \item{\code{Miocene}}{a numeric vector}
    \item{\code{Pliocene}}{a numeric vector}
  }
}
\source{
Shin et al.(2019), Java, Indonesia.
}
\examples{
data(fossil)
\dontrun{
fossil.out <- iNEXT.4steps(data = fossil, datatype = "abundance")
fossil.out
}
}
\keyword{datasets}

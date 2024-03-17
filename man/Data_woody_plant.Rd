\name{Data_woody_plant}
\alias{Data_woody_plant}
\docType{data}
\title{
Woody plants incidence-based frequency data
}
\description{
The datasets considered here are a subset of The National Vegetation Database of Taiwan (AS-TW-001), sampled between 2003 and 2007 within the first national vegetation inventory project (Chiou et al. 2009).
Over 3600 vegetation plots, each 20x20-m in area, were set up in various locations in Taiwan, and all woody plant individuals taller than 2 meters were recorded in each plot.
For illustration here, we selected only plots belonging to two vegetation types (according to Li et al. 2013): Pyrenaria-Machilus subtropical winter monsoon forest and Chamaecyparis montane mixed cloud forest, sampled in the northern part of Taiwan (in ecoregions 7 and 8 according to Su 1985).

}
\usage{data(Data_woody_plant)}
\format{
  \code{Data_woody_plant} a list of two species-by-sampling-unit data frames. Each element in the data frame is 1 for a detection, and 0 for a non-detection.\cr
  A list of 2: \cr
  $ Monsoon    : num [1:329, 1:191] 0 0 0 0 0 0 0 0 0 0 ... \cr
  $ Upper_cloud: num [1:239, 1:153] 0 0 0 0 0 0 0 0 0 0 ... \cr
}
\source{
The National Vegetation Database of Taiwan (AS-TW-001).
}
\examples{
data(Data_woody_plant)
}
\references{
Chiou, C.-R., Hsieh, C.-F., Wang, J.-C., Chen, M.-Y., Liu, H.-Y., Yeh, C.-L., ... Song, M. G.-Z. (2009). The first national vegetation inventory in Taiwan. Taiwan Journal of Forest Science, 24, 295-302. \cr\cr
Li, C.-F., Chytry, M., Zeleny, D., Chen, M. -Y., Chen, T.-Y., Chiou, C.-R., ... Hsieh, C.-F. (2013). Classification of Taiwan forest vegetation. Applied Vegetation Science, 16, 698-719. https://doi.org/10.1111/avsc.12025 \cr\cr
Su, H. -J. (1985). Studies on the climate and vegetation types of the natural forests in Taiwan (III) A Scheme of Geographical Climatic Regions. Quarterly Journal of Chinese Forestry, 18, 33-44. \cr\cr
}
\keyword{datasets}

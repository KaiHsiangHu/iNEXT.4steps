\name{Data_woody_plant}
\alias{Data_woody_plant}
\docType{data}
\title{
  Incidence raw data
}
\description{
  This dataset was taken from The National Vegetation Database of Taiwan, sampled between 2003 and 2007 (Chiou et al. 2009). Only data in plots (each 20x20-m in area) belonging to two 
vegetation types (monsoon forest and upper_cloud forest) were used (Li et al. 2013). All woody plant individuals taller than 2 meters were recorded in each plot. 
In the monsoon forest, 329 species and 6814 incidences were recorded in 191 plots. In the upper cloud forest, 239 species and 3371 incidences were recorded in 153 plots (each plot is regarded as a sampling unit). 

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
  Li, C.-F., Chytry, M., Zeleny, D., Chen, M. -Y., Chen, T.-Y., Chiou, C.-R., ... Hsieh, C.-F. (2013). Classification of Taiwan forest vegetation. Applied Vegetation Science, 16, 698-719. \cr
  https://doi.org/10.1111/avsc.12025 \cr\cr
  
}
\keyword{datasets}

\name{hyperDataSubset}
\alias{hyperDataSubset}
\docType{data}
\title{
Hypericum Perennial Dataset
}
\description{
Demographic data of Hypericum cumulicola in Archbold Biological Station (FL, USA). Life cycle, experimental design and data are described in Quintana-Ascencio & Menges (2003). Data contains a subset of individuals from population "bald 1" and annual period "1997-1998". Full dataset can be obtained upon request to the authors (pquintan@mail.ucf.edu). For experimental design see references below.
}
\usage{data(hyperDataSubset)}
\format{
  The format is:
 chr "HyperDataSubset"
}
\details{
data-frame with headings:
	- id: unique plant id (this file contains only a subset of all individuals)
	- bald: population (this subset contains only one population)
	- year: transition from t to t+1 (this subset contains only data for 1997-1998)
	- surv: survival (0 = dead, 1= alive, NAs if not yet recruited)
	- size: length of longest stem in individual (cm) in time t
	- stage: recruits vs established individuals in time t (0 = individual was recruited in time t, 1 = already established individual prior to time t, NA = individual not yet recruited in time t)
	- fec0Flowering: probability of reproduction (0= no flowering, 1 = individual was flowering in time t)
	- fec1Fruits: number of fruits per plant (NA if not Fec0Flowering = 0)
	- sizeNext: length of longest stem in individual (cm) in time t+1
	- stageNext: recruits vs established individuals in time t+1 (0 = individual was recruited in time t+1, 1 = already established individual prior to time t+1, NA = individual not yet recruited or dead in t+1)
}

\references{
Quintana-Ascencio & Menges. 2003. A fire-explicit population viability analyses of Hypericum cumulicola in Florida Rosemary scrub. Conservation Biology 17, p433-449.

}
\author{
Pedro Quintana Ascencio & Eric Menges
}

\examples{
data(hyperDataSubset)
print(head(hyperDataSubset))

%WOULD YOU LIKE ME TO ADD THE ONLINE APPENDIX OF THE MEE PAPER HERE TOO?

}
\keyword{datasets}


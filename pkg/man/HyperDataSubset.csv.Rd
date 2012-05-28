\name{hyperDataSubset}
\alias{hyperDataSubset}
\docType{data}
\title{
Hypericum Perennial Dataset
}
\description{
Demographic data of Hypericum cumulicola in "Florida rosemary scrub at" Archbold Biological Station (FL, USA). Life cycle, experimental design and data are described in Quintana-Ascencio & Menges (2003). Data contains a subset of individuals from population "bald 1" and annual period "1997-1998". Full dataset can be obtained upon request to the authors (pedro.quintana-ascencio@ucf.edu).
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
	
	- ontogeny: recruits vs established individuals in time t (0 = individual was recruited in time t, 1 = already established individual prior to time t, NA = individual not yet recruited in time t)
	
	- fec0Flowering: probability of reproduction (0= no flowering, 1 = individual was flowering in time t)
	
	- fec1Fruits: number of fruits per plant (NA if not Fec0Flowering = 0)
	
	- sizeNext: length of longest stem in individual (cm) in time t+1
	
	- ontogenyNext: recruits vs established individuals in time t+1 (0 = individual was recruited in time t+1, 1 = already established individual prior to time t+1, NA = individual not yet recruited or dead in t+1)
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


#Due to the sampling design described above, here we consider only individual with certain recruit origin:
d <- hyperDataSubset
d <- subset(d,is.na(d$size)==FALSE | d$ontogenyNext==1)

#Side-experiments revealed that the following vital rates are size-independent and equal to:
  #Number of seeds produced per fruit
    fec2Seeds<-13.78
  #Probability of seedling establishment
    fec3Establish<-0.001336
  #Probability of seedling survival half a year after germinating, corresponding to the next annual census
    fec4HalfYearSeedlingSurvival<-0.14
  #Probability of a seed going into the seed bank
    goSB<-0.08234528
  #Probability of a seed staying in the seed bank
    staySB<-0.672
#Note that the aforementioned vital rates are function of time since last fire. See Quintana-Ascencio et al (2003) for more information.

#A simple re-organization of the data, getting rid of non-critical information
d<-d[,c("stage","stageNext","surv","size","sizeNext","fec0Flowering","fec1Fruits")]

#The following states the continuous (max height of individual plant) part of the IPM. Note that the IPM to be constructed here contains a discrete stage: seedbank.
d$stage<-d$stageNext<-"continuous"
d$stage[is.na(d$size)]<-NA
#If individual did not survive, it is labelled as dead to t+1.
d$stageNext[d$surv==0]<-"dead"
#Adds probability of seeds going into (continuous -> seedbank), staying (seedbank -> seedbank) and leaving (continuous -> seedbank) the discrete stage.
d$number<-1
sb<-data.frame(stage=c("seedbank","seedbank","continuous"),stageNext=c("seedbank","continuous","seedbank"),surv=1,size=NA,sizeNext=NA,fec0Flowering=NA,fec1Fruits=NA,number=c(staySB,(1-staySB)*fec3Establish*fec4HalfYearSeedlingSurvival,1))
#d<-rbind(d,sb)
#d$stage<-as.factor(d$stage)
#d$stageNext<-as.factor(d$stageNext)




}
\keyword{datasets}


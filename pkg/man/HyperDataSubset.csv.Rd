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


#Due to the sampling design described above, here we consider only individual with certain recruit origin:
d <- hyperDataSubset
d <- subset(d,is.na(d$size)==FALSE | d$stageNext==1)

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
d<-rbind(d,sb)
d$stage<-as.factor(d$stage)
d$stageNext<-as.factor(d$stageNext)

#Create survival object with a square term 
so<-makeSurvObj(d,explanatoryVariables="size+size2")
picSurv(d,so)

#Create growth object
go<-makeGrowthObj(d,explanatoryVariables="size", responseType = "sizeNext")
picGrow(d,go)
abline(a=0,b=1,lty=2)

#Create fecundity object

 #First, info on mean size of offspring is necessary
meanOffspringSize=mean(d[is.na(d$size)==TRUE&d$stageNext=="continuous","sizeNext"])

fo<-makeFecObj(d,Family=c("binomial","poisson"),
               explanatoryVariables=c("size","size"),
               meanOffspringSize=mean(d[is.na(d$size)==TRUE&is.na(d$sizeNext)==FALSE,"sizeNext"]),
               sdOffspringSize=sd(d[is.na(d$size)==TRUE&is.na(d$sizeNext)==FALSE,"sizeNext"]),
               fecConstants=data.frame(fec2Seeds=fec2Seeds,fec3Establish=fec3Establish,fec4HalfYearSeedlingSurvival=fec4HalfYearSeedlingSurvival),   
               offspringSplitter=data.frame(seedbank=goSB,continuous=(1-goSB)),
               offspringTypeRates=data.frame(seedbank=c(1,1,1,0,0),continuous=rep(1,5),row.names=c("fec0Flowering","fec1Fruits","fec2Seeds","fec3Establish","fec4HalfYearSeedlingSurvival")),
               Transform=c("none","none"))

#Define discrete transition matrix
dto<-makeDiscreteTrans(d)
dummy<-as.matrix(fo@meanOffspringSize)
dimnames(dummy)<-list(1,"seedbank")
dto@meanToCont<-as.matrix(dummy,dimnames=c(1,"seedbank"))
dummy<-as.matrix(fo@sdOffspringSize)
dimnames(dummy)<-list(1,"seedbank")
dto@sdToCont<-as.matrix(dummy,dimnames=c(1,"seedbank"))
dto@discreteSurv[1,1]<-staySB+(1-staySB)*fec3Establish*fec4HalfYearSeedlingSurvival

#Create the T matrix
Tmatrix<-createIPMTmatrix(growObj=go,survObj=so,discreteTrans=dto,minSize=0,maxSize=80,nBigMatrix=81,correction="constant")

#Create the F matrix
Fmatrix<-createIPMFmatrix(fecObj=fo,minSize=0,maxSize=80,nBigMatrix=81,correction="constant")

#Build a T matrix reflecting only the continuous part of the model and check that binning, etc is adequate
Tmatrix.continuous.only <- createIPMTmatrix(growObj=go,survObj=so,minSize=0,maxSize=80,nBigMatrix=81,correction="constant")
diagnosticsTmatrix(Tmatrix.continuous.only,growObj=go,survObj=so,dff=d, correction="constant")

#Forming the IPM as a result of adding the T and F matrices
IPM <- Tmatrix + Fmatrix

#Population growth rate for the whole life cycle of Hypericum is
eigen(IPM)$value[1]
#Population growth rate excluding the seed bank stage is
eigen(IPM[2:82,2:82])$value[1]



}
\keyword{datasets}


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
	
	- size: length of longest stem in individual (cm) in time t
	
	- ontogeny: recruits vs established individuals in time t (1 = individual was recruited in time t, 0 = already established individual prior to time t, NA = individual not yet recruited in time t)
	
	- fec0: probability of reproduction (0= no flowering, 1 = individual was flowering in time t, NA = individual not alive in year t)
	
	- fec1: number of fruits per plant (NA if fec0 = 0)
	
	- surv: survival (0 = dead, 1= alive, NAs if not yet recruited or past dead)
	
	- sizeNext: length of longest stem in individual (cm) in time t+1
	
	- ontogenyNext: recruits vs established individuals in time t+1 (1 = individual was recruited in time t+1, 0 = already established individual prior to time t+1, NA = individual not yet recruited or dead in t+1)
}

\references{
Quintana-Ascencio, Menges & Weekley. 2003. A fire-explicit population viability analyses of Hypericum cumulicola in Florida Rosemary scrub. Conservation Biology 17, p433-449.

}
\author{
Pedro Quintana Ascencio & Eric Menges
}

\examples{

#Access data from the long-term censuses on Hypericum cumulicula carried out by Eric Menges, Pedro Quintana-Ascencio and coworkers at Archbold Biological Station. Here only a subset of individuals from population 'bald 1' and for the annual transition '1997-1998' are shown.
data("hyperDataSubset")
d<-hyperDataSubset

#Variables are:
  #id: unique identifier for each individual
  #bald: population. Here only bald 1
  #year: annual transition of the long-term data. Here only 1997-1998
  #surv: survival (1) or not (0) of individuals between 1997 and 1998
  #size: maximum height of the stems of each individual
  #ontogeny: because the demography of Hypericum is very dynamic (turnover is very high) the experimental design described in Quintana-Ascencio et al. (2003) consists on establishing new permanent plots every year at each population, in addition to censusing old plots. Here we differentiate between individuals that appear for the first time in time t because they were recruits (1) and those that, not being new recruits, where measured for the first time in t because they were in a new permanent plot.
  #fec0: probability of flowering (1) or not (0)
  #fec1: number of fruits per individual
  #sizeNext: same as "size" above, for t+1
  #stageNext: same as "stage" above, for t+1

#Due to the sampling design described above, here we consider only individual with a certain recruit origin:
d <- subset(d,is.na(d$size)==FALSE | d$ontogenyNext==1)

#Side-experiments revealed that the following vital rates are size-independent and equal to:
  #Number of seeds produced per fruit
    fec2<-13.78
  #Probability of seedling establishment
    fec3<-0.001336
  #Probability of seedling survival half a year after germinating, corresponding to the next annual census
    fec4<-0.14
  #Probability of a seed going into the seed bank
    goSB<-0.08234528
  #Probability of a seed staying in the seed bank
    staySB<-0.672
#Note that the aforementioned vital rates are function of time since last fire, but because here we are only dealing with one population and one year transition, we treat them as constants. See Quintana-Ascencio et al (2003) for more information.

#A simple re-organization of the data, getting rid of non-critical information
d<-d[,c("surv","size","sizeNext","fec0","fec1")]

#The following states the continuous (max height of individual plant) part of the IPM. Note that the IPM to be constructed here contains a discrete stage: seedbank.
d$stageNext<-d$stage<-"continuous"
d$stage[is.na(d$size)]<-NA
#If individual did not survive, it is labelled as dead to t+1.
d$stageNext[d$surv==0]<-"dead"
#Adds probability of seeds going into (continuous -> seedbank), staying (seedbank -> seedbank) and leaving (continuous -> seedbank) the discrete stage.
d$number<-1
#sb<-data.frame(stage=c("seedbank","seedbank","continuous"),stageNext=c("seedbank","continuous","seedbank"),surv=1,size=NA,sizeNext=NA,fec0=NA,fec1=NA,number=c(staySB,(1-staySB)*fec3*fec4,1)) # WHY THE CONTINUOUS TO SEEDBANK TRANSITION HERE? PLANTS CAN'T SURVIVE INTO THE SEEDBANK, RIGHT? TRANSITION ONLY POSSIBLE THROUGH FECUNDITY, RIGHT? EELKE HAS TAKEN IT OUT BELOW.
#sb<-data.frame(stage=c("seedbank","seedbank"),stageNext=c("seedbank","continuous"),surv=1,size=NA,sizeNext=NA,fec0=NA,fec1=NA,number=c(staySB,(1-staySB)*fec3*fec4))
#d<-rbind(d,sb)
d$stage<-as.factor(d$stage)
d$stageNext<-as.factor(d$stageNext)

#Carry out comparisons to establish the best survival model
testSurv <- survModelComp(d, expVars = c(surv~1, surv~size, surv~size + size2), testType = "AIC",makePlot = TRUE,legendPos = "bottomleft")

#Carry out comparisons to establish the best growth model
testGrow <- growthModelComp(d,expVars = c(sizeNext~1, sizeNext~size, sizeNext~size + size2), regressionType = "constantVar", testType = "AIC", makePlot = TRUE, legendPos = "bottomright")

#Create survival object using regression model indicated by testSurv
so <- makeSurvObj(d, Formula = surv~size + size2)
picSurv(d,so)

#Create growth object using regression model indicated by testGrown
go<-makeGrowthObj(d, Formula = sizeNext~size) 
picGrow(d,go)
abline(a=0,b=1,lty=2)
#JESS/EELKE: THE 1:1 LINE GOES BEYOND THE PLOT IN MY SCREEN AND IT DOES NOT DEPART FROM (0,0)

#Create fecundity object using regression models
fo <- makeFecObj(d, Formula=c(fec0~size, fec1~size),
                 Family=c("binomial","poisson"),
                 Transform=c("none", "none"),
                 meanOffspringSize=mean(d[is.na(d$size)==TRUE & is.na(d$sizeNext)==FALSE,"sizeNext"]),
                 sdOffspringSize=sd(d[is.na(d$size)==TRUE & is.na(d$sizeNext)==FALSE,"sizeNext"]),
                 fecConstants=data.frame(fec2=fec2,fec3=fec3,fec4=fec4),
                 offspringSplitter=data.frame(seedbank=goSB,continuous=(1-goSB)),
                 vitalRatesPerOffspringType=data.frame(seedbank=c(1,1,1,0,0),
                                               continuous=c(1,1,1,1,1),
                                               row.names=c("fec0","fec1","fec2","fec3","fec4")))


#Define discrete transition matrix - THIS IS WEIRD - UNDER RE-WORKING BY EELKE
dto<-makeDiscreteTrans(d, 
		discreteTrans = matrix(c(staySB,(1-staySB)*fec3*fec4,(1-staySB)*(1-fec3*fec4),0,sum(d$number[d$stage=="continuous"&d$stageNext=="continuous"],na.rm=TRUE),sum(d$number[d$stage=="continuous"&d$stageNext=="dead"],na.rm=TRUE)),ncol=2,nrow=3,dimnames=list(c("seedbank","continuous","dead"),c("seedbank","continuous"))),
		meanToCont = matrix(mean(d$sizeNext[is.na(d$stage)&d$stageNext=="continuous"]),ncol=1,nrow=1,dimnames=list(c("mean"),c("seedbank"))),
		sdToCont = matrix(sd(d$sizeNext[is.na(d$stage)&d$stageNext=="continuous"]),ncol=1,nrow=1,dimnames=list(c(""),c("seedbank"))))

#choose number of bins for discretization in the IPM
nBigMatrix <- 100

#Create the P matrix describing growth-survival transitions
#  The argument correction="discretizeExtremes" places parts of the growth distribution that fall
#  below minSize or above maxSize into the first and last bin 
#      
Pmatrix<-createIPMPmatrix(growObj=go,survObj=so,discreteTrans=dto,
                          minSize=0,maxSize=80,nBigMatrix=nBigMatrix,
                          correction="discretizeExtremes")

#Create the F matrix descributing fecundity transitions
#  The argument correction="discretizeExtremes" places parts of the continuous offspring distribution that fall
#  below minSize or above maxSize into the first and last bin
#
Fmatrix<-createIPMFmatrix(fecObj=fo,
                          minSize=0,maxSize=80,nBigMatrix=nBigMatrix,
                          correction="discretizeExtremes")

#Build a P matrix reflecting only the continuous part of the model and check that binning, etc is adequate
PmatrixContinuousOnly <- createIPMPmatrix(growObj=go,survObj=so,minSize=0,maxSize=70,nBigMatrix=nBigMatrix,correction="discretizeExtremes")
diagnosticsPmatrix(PmatrixContinuousOnly,growObj=go,survObj=so,dff=d, correction="discretizeExtremes")

#Form the IPM as a result of adding the P and F matrices
IPM <- Pmatrix + Fmatrix

#Population growth rate for the whole life cycle of Hypericum is
eigen(IPM)$value[1]
#Population growth rate excluding the seed bank stage is
eigen(IPM[2:(nBigMatrix+1),2:(nBigMatrix+1)])$value[1]



}
\keyword{datasets}


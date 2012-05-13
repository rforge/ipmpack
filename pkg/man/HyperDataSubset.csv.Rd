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

#Carry out comparisons to establish the best survival model
testSurv <- survModelComp(dataf=d, expVars = c("1", "size", "size + size2"),
                          testType = "AIC", makePlot = TRUE)

#Carry out comparisons to establish the best growth model
testGrow <- growthModelComp(dataf=d, expVars = c("1", "size", "size + size2"),
                            regressionType = "constantVar", 
                            respType = "sizeNext", testType = "AIC", makePlot = TRUE)

#Create survival object with a square term (from the above)
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

#The following creates a 3x3 panel with some key features of the IPM for Hypericum presented in figure 2 of the manuscript to which this online appendix is associated

par(mfrow=c(3,3),mar=c(4.1,4.1,1,1))

#Panel 1: Survival
ncuts=15
os <- order(d$size)
os.size <- (d$size)[os]
psz <- tapply(os.size, as.numeric(cut(os.size, ncuts)), mean, na.rm = TRUE)
picSurv(d,so,ncuts=ncuts)
legend(-4,1.1,"a)",bty="n",cex=1.5)

#Panel 2: Growth
picGrow(d,go)
dummy<-go@fit$coef[1]+go@fit$coef[2]*os.size
lines(os.size,dummy,col=2)
abline(a=0,b=1,lty=2)
legend(-6,63,"b)",bty="n",cex=1.5)

#Panel 3: Probability of flowering
os.fec0Flowering <- (d$fec0Flowering)[os]
ps <- tapply(os.fec0Flowering, as.numeric(cut(os.size, ncuts)), mean, na.rm = TRUE)
plot(as.numeric(psz), as.numeric(ps), pch = 19, xlab = "Size at t", ylab = "Flowering at t", main = "Flowering")
dummy<-fo@fitFec[[1]]$coef[1]+fo@fitFec[[1]]$coef[2]*os.size
lines(os.size,exp(dummy)/(1+exp(dummy)),col=2)
legend(-4,1.1,"c)",bty="n",cex=1.5)

#Panel 4: Number of fruits as a function of individual size
plot(d$size,d$fec1Fruits,pch=19,xlab="Size at t",ylab="# Fruits", main="Fruit production")
dummy<-fo@fitFec[[2]]$coef[1]+fo@fitFec[[2]]$coef[2]*os.size
lines(os.size,exp(dummy),col=2)
legend(-8,700,"d)",bty="n",cex=1.5)

#Panel 5: IPM kernel
image(Tmatrix@meshpoints, Tmatrix@meshpoints, log(t(IPM[2:82,2:82])), main = "IPM kernel",xlab = "Size in t", ylab = "Size in t+1")
legend(-6,83,"e)",bty="n",cex=1.5)

#Panel 6: IPM elements elasticity
elast <- elas(IPM)
image(Tmatrix@meshpoints, Tmatrix@meshpoints, t(elast[2:82,2:82]), main = "Elasticity",xlab = "Size at t", ylab = "Size at t+1")
legend(-6,83,"f)",bty="n",cex=1.5)

#Panel 7: Age-specific trajectories
ageing=survivorship(Tmatrix.continuous.only,1,20)
plot(ageing$surv.curv,type="l",ylab="Survival (lx) & force of mortality (qx)",xlab="Age",main="Age-specific trajectories",col="dark gray")
points(ageing$mortality,type="l",col="red",lty=2)
legend(-1.5,0.57,"g)",bty="n",cex=1.5)

#Panel 8: mean and variation of passage time & mean and variation of life expectancy?
mLE<-meanLifeExpect(Tmatrix)
vLE<-varLifeExpect(Tmatrix)
plot(Tmatrix@meshpoints,mLE[2:ncol(Tmatrix)],ylim=c(0,max(mLE+vLE)),type="l",
     ylab="Years",xlab="Size",main="Life expectancy",lty=1,col="red")
polygon(c(Tmatrix@meshpoints,rev(Tmatrix@meshpoints)),c(vLE[2:ncol(Tmatrix)]+mLE[2:ncol(Tmatrix)],rev((-vLE[2:ncol(Tmatrix)]+mLE[2:ncol(Tmatrix)]))),col=rgb(1,0,0,0.5),border=10,lwd=0.5,density=115)
legend(-8,11,"h)",bty="n",cex=1.5)

#Panel 9: Mean passage time to a chosen size
chosenSize=30
plot(Tmatrix@meshpoints,passageTime(chosenSize,Tmatrix)[2:ncol(Tmatrix)],type="l",lty=1,col="blue",
      ylab="Years",xlab="Size",main="Passage time to 30cm")
abline(v=chosenSize,lty=4,col="dark gray")
legend(-8,4.3,"i)",bty="n",cex=1.5)



}
\keyword{datasets}


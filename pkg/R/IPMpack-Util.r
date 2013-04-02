


# Function to extract IPM output from a list
# of P (survival + growth) and F (fecundity) matrices
# (usually from Bayes fit) 
#
# Parameters - PmatrixList
#            - targetSize - the size you want passage time estimated for.
#            - FmatrixList
#
# Returns - a list 

getIPMoutput <- function(PmatrixList,targetSize=c(),FmatrixList=NULL){
	
	if (length(targetSize)==0)  { 
		print("no target size for passage time provided; taking meshpoint median")
		targetSize <- median(PmatrixList[[1]]@meshpoints)
	}
	nsamps <- length(PmatrixList)
	h1 <- PmatrixList[[1]]@meshpoints[2]-PmatrixList[[1]]@meshpoints[1]
	stableStage <- LE <- pTime <- matrix(NA,nsamps,length(PmatrixList[[1]]@.Data[1,]))
	lambda <- rep(NA,nsamps)
	for (k in 1:nsamps) {
		Pmatrix <- PmatrixList[[k]]
		LE[k,]<-meanLifeExpect(Pmatrix) 
		pTime[k,]<-passageTime(targetSize,Pmatrix) 
		
		if (class(FmatrixList)!="NULL") {
			IPM <- Pmatrix + FmatrixList[[k]]
			lambda[k] <- Re(eigen(IPM)$value[1])
			stableStage[k,] <- eigen(IPM)$vector[,1]
			#normalize stable size distribution
			stableStage[k,] <- stableStage[k,]/(h1*sum(stableStage[k,]))
		}
	}
	
	return(list(LE=LE,pTime=pTime,lambda=lambda,stableStage=stableStage))
	
}


# Function to extract IPM output from posteriors 
# (usually from Bayes fit)  - quicker to do all at once
# rather than build list of IPM P matrices, then list of IPM F matrices
#
# Parameters - survObjlist - list of survival objects
#            - growObjList - list of growth objects
#            - targetSize - the size you want passage time estimated for.
#            - nBigMatrix - the number of bins
#            - minSize - the minimum size
#            - maxSize - the maximum size
#            - cov - do you want to fit a discrete covariate
#            - fecObjList - list of fecundity objects
#            - envMat - matrix of env transitions (only if cov=TRUE)
#            - nsizeToAge - numeric describing how many size to age defined (0 - 100s)
# 
#
# Returns - a list 

getIPMoutputDirect <- function(survObjList,growObjList,targetSize=c(),
		nBigMatrix,minSize,maxSize,discreteTrans = 1,
		cov=FALSE,fecObjList=NULL, envMat=NULL,
		nsizeToAge=0, sizeStart = 10,
		integrateType = "midpoint", correction = "none", storePar=TRUE,
		chosenCov = data.frame(covariate = 1),
		onlyLowerTriGrowth=FALSE){
	
	# adjust the sample lengths to they are all the same
	if (length(targetSize)==0)  targetSize <- 0.2*(minSize+maxSize)
	nsamp <- max(length(growObjList),length(survObjList),length(fecObjList))
	if (length(survObjList)<nsamp)  
		survObjList <- sample(survObjList,size=nsamp,replace=TRUE)
	if (length(growObjList)<nsamp)  
		growObjList <- sample(growObjList,size=nsamp,replace=TRUE)
	if (class(fecObjList)!="NULL") {
		if (length(fecObjList)<nsamp)  
			fecObjList <- sample(fecObjList,size=nsamp,replace=TRUE)
	}
	
	# store chosen parameters
	if (storePar){
	surv.par <- matrix(NA,nsamp,length(survObjList[[1]]@fit$coefficients))
	grow.par <- matrix(NA,nsamp,length(growObjList[[1]]@fit$coefficients)+1)
	for (k in 1:nsamp) {
		surv.par[k,] <- survObjList[[k]]@fit$coefficients
		grow.par[k,] <- c(growObjList[[k]]@fit$coefficients,growObjList[[k]]@sd)
	}} else { surv.par <- grow.par <- c()}
	
	#set up storage
	if (class(discreteTrans)=="discreteTrans") nDisc <- (ncol(discreteTrans@discreteTrans)-1) else nDisc <- 0
	
	if (class(envMat)!="NULL") nEnv <- envMat@nEnvClass else nEnv <- 1
	LE <- pTime <- matrix(NA,nsamp,(nBigMatrix+nDisc)*nEnv)
	if (class(fecObjList)=="NULL") {
		lambda <- stableStage <- c()
	} else {
		stableStage <- matrix(NA,nsamp,(nBigMatrix+nDisc)*nEnv)
		lambda <- rep(NA,nsamp)
	}
	if (nsizeToAge==0) { resAge <- resSize <- c() } else { resAge <- resSize <- matrix(NA,nsamp,nsizeToAge)} 
	if (length(sizeStart)==0) { if (minSize<0) sizeStart <- 0.5*minSize else sizeStart <- 2*minSize }
	
	#go!
	for (k in 1:nsamp) {
		
		if (!cov) {
			Pmatrix <- makeIPMPmatrix(nBigMatrix = nBigMatrix, minSize = minSize, 
					maxSize = maxSize,  growObj = growObjList[[k]],
					survObj = survObjList[[k]],discreteTrans=discreteTrans,
					integrateType=integrateType, correction=correction) 
			
		} else {
			Pmatrix <- makeCompoundPmatrix(nEnvClass = nEnv,
					nBigMatrix = nBigMatrix, minSize = minSize,
					maxSize = maxSize, envMatrix=envMat,growObj = growObjList[[k]],
					survObj = survObjList[[k]],discreteTrans=discreteTrans,
					integrateType=integrateType, correction=correction)    
			
		}
		
		if (onlyLowerTriGrowth & !cov) {
			Pmatrix@.Data <- Pmatrix@.Data*lower.tri(Pmatrix@.Data, diag = TRUE)
			nvals <- colSums(Pmatrix@.Data,na.rm=TRUE)
			Pmatrix@.Data <- t((t(Pmatrix@.Data)/nvals) *
							surv(size = Pmatrix@meshpoints, 
									cov = chosenCov,
									survObj = survObjList[[k]]))                                    
		}
		
	
		LE[k,] <- meanLifeExpect(Pmatrix) 
		pTime[k,] <- passageTime(targetSize,Pmatrix) 
		if (k==1) h1 <- diff(Pmatrix@meshpoints)[1]
		
		if (class(fecObjList)!="NULL") {
			if (!cov) { 
				Fmatrix <- makeIPMFmatrix(nBigMatrix = nBigMatrix, minSize = minSize, 
						maxSize = maxSize,  
						fecObj=fecObjList[[k]],
						integrateType=integrateType, correction=correction)
			} else {
				Fmatrix <- makeCompoundFmatrix(nEnvClass = nEnv,
						nBigMatrix = nBigMatrix, minSize = minSize, 
						maxSize = maxSize, envMatrix=envMat,
						fecObj=fecObjList[[k]],integrateType=integrateType, correction=correction)
			}
			
			
			
			IPM <- Pmatrix + Fmatrix
			lambda[k] <- Re(eigen(IPM)$value[1])
			stableStage[k,] <- eigen(IPM)$vector[,1]
			#normalize stable size distribution
			stableStage[k,] <- stableStage[k,]/(h1*sum(stableStage[k,]))
			
			#print("here2")
		}
		
		# get size to age results
		if (nsizeToAge>0) { 
			res2 <- sizeToAge(Pmatrix=Pmatrix,startingSize=minSize*1.3,
					targetSize=seq(sizeStart,maxSize*0.9,length=nsizeToAge))
			resAge[k,] <- res2$timeInYears
			resSize[k,] <- res2$targetSize
		}
		
	}
	
	return(list(LE=LE,pTime=pTime,lambda=lambda,stableStage=stableStage,
					meshpoints=Pmatrix@meshpoints,resAge=resAge,resSize=resSize,
					surv.par=surv.par,grow.par=grow.par))
	
}



## Function to get passage time FROM a particular size TO a range of sizes
## (i.e. size to age) when provided with a Pmatrix, a starting size, and a list
## of target sizes
#
# Parameters - Pmatrix
#            - startingSize
#            - targetSizes
#
# Returns - list containing vector of targets, vector of corresponding times, and the startingSize
#
sizeToAge <- function(Pmatrix,startingSize,targetSize) {
	
	#locate where the first size is in the meshpoints of Pmatrix
	diffv <- abs(startingSize-Pmatrix@meshpoints)
	start.index <- median(which(diffv==min(diffv),arr.ind=TRUE))
	timeInYears <- rep(NA,length(targetSize))
	
	#loop over to see where its going
	for (k in 1:length(targetSize)) {
		pTime <- passageTime(targetSize[k],Pmatrix)
		timeInYears[k] <- pTime[start.index]
	}
	
	return(list(timeInYears=timeInYears,targetSize=targetSize,startingSize=startingSize))
	
}






## FUNCTION FOR MAKING PICTURE OF SURVIVAL DATA AND FITTED SURVIVAL OBJECT ############################
## and PICTURE GROWTH DATA AND FITTED GROWTH OBJECT  
# ** note that growth may need redefining if new growth  objects are defined


# parameters - dataf - the dataframe
#            - survObj - a survival object
#            - ncuts - the number of cuts used for binning survival
# returns - 

picSurv <- function(dataf, survObj, ncuts = 20, makeTitle = "Survival", ...) { 
	
	#organize data and plot mean of ncut successive sizes, so trend more obvious
	os<-order(dataf$size); os.surv<-(dataf$surv)[os]; os.size<-(dataf$size)[os]; 
	psz<-tapply(os.size,as.numeric(cut(os.size,ncuts)),mean,na.rm=TRUE); #print(psz)
	ps<-tapply(os.surv,as.numeric(cut(os.size,ncuts)), mean, na.rm = TRUE);#print(ps)
	
	if (length(grep("covariate",names(survObj@fit$model)))==0) {  
		#plot data
		plot(as.numeric(psz),as.numeric(ps),pch=19,
				xlab="Size at t", ylab = "Survival to t+1", main = makeTitle, ...)
		#Plot fitted models
		points(dataf$size[order(dataf$size)],surv(dataf$size[order(dataf$size)],data.frame(covariate=1),survObj),type="l",col=2)
	} else {
		plot(as.numeric(psz),as.numeric(ps),
				type="n",pch=19,xlab="Size at t", ylab="Survival to t+1",main="Survival",...)
		dataf$covariate <- as.factor(dataf$covariate)
		levels(dataf$covariate) <- 1:length(unique(dataf$covariate))
		os.cov<-(dataf$covariate)[os]
		sizes <- dataf$size[!is.na(dataf$size)]; sizes <- sizes[order(sizes)]
		ud <- unique(dataf$covariate); ud <- ud[!is.na(ud)]
		for (k in 1:length(ud)) { 
			tp <- os.cov==ud[k]
			psz<-tapply(os.size[tp], as.numeric(cut(os.size[tp], ncuts)), mean, na.rm = TRUE); #print(psz)
			ps<-tapply(os.surv[tp],as.numeric(cut(os.size[tp],ncuts)),mean,na.rm=TRUE);#print(ps)
			points(as.numeric(psz),as.numeric(ps),pch=19,col=k)
			newd <- data.frame(size=sizes,size2=sizes^2,size3=sizes^3,
					covariate=rep(as.factor(ud[k]),length(sizes)))
			if(length(grep("expsize",survObj@fit$formula))==1)
				newd$expsize=exp(sizes)
			if(length(grep("logsize",survObj@fit$formula))==1)
				newd$logsize=log(sizes)
			if(length(grep("logsize2",survObj@fit$formula))==1)
				newd$logsize=(log(sizes))^2
			
			pred.surv <- predict(survObj@fit,newd,type="response")
			points(newd$size,pred.surv,type="l",col=k)
		}
	} 
	
	
	
	
	
}


## Function defining growth using Hossfeld function
#
# Parameters - size - DBH
#            - par - three parameters
#
# Returns - growth increment (not log scale)
Hossfeld <- function(size,par) {
	deltDBH <- (par[2]*par[3]*size^(par[3]-1))/((par[2]+((size^par[3])/par[1]))^2)
	return(deltDBH)
}


## Function to fit Hossfeld function using optim 
#
# Parameters - par - three parameters
#            - dataf - a data-frame
#
# Returns the SS

wrapHossfeld <- function(par, dataf) { 
	pred <- Hossfeld(dataf$size, par[1:3]) 
	ss <- sum((pred - dataf$incr)^2, na.rm = T)
	return(ss) 
} 


# Function to make pic of growth object fit with data
#
# parameters - dataf - the dataframe
#            - growObj - a growth object
#
# returns - 
#

picGrow <- function(dataf, growObj, mainTitle = "Growth",...) {
	predVar <- attr(growObj@fit$terms,"predvars")[[2]]  #jess quick fix. this function does not work with declineVar either at the moment
	if (class(growObj)=="growthObjTruncIncr") { 
		predVar <- "incr"	
	} else {
		predVar <- attr(growObj@fit$terms,"predvars")[[2]]
	}
		
	if(predVar == "sizeNext") {
		plot(dataf$size, dataf$sizeNext,pch=19, xlab="Size at t", ylab="Size at t+1", main = mainTitle,...)
		abline(a = 0, b = 1)
	}else{
		dataf$incr <- dataf$sizeNext - dataf$size
		plot(dataf$size, dataf$incr, pch = 19, xlab = "Size at t", ylab="Size increment", main = mainTitle,...)
		abline(a = 0, b = 0)
	}
	
	if (length(grep("covariate", names(growObj@fit$model))) > 0) {  
		#convert to 1:n for indexing later and to relate to discrete
		dataf$covariate <- as.factor(dataf$covariate)
		levels(dataf$covariate) <- 1:length(unique(dataf$covariate))
		ud <- unique(dataf$covariate); ud <- ud[!is.na(ud)]
		for (k in 1:length(ud)) { 
			tp <- dataf$covariate == ud[k]
			points(dataf$size[tp], dataf$sizeNext[tp], pch = 19, col = k)            
		}
		ud <- as.factor(ud)
	} else {
		ud <- 0
	} 
	
	sizes <- dataf$size[!is.na(dataf$size)]; sizes <- sizes[order(sizes)]
	
	for (k in 1:length(ud)) { 
		newd <- data.frame(size = sizes, size2 = sizes ^ 2,size3 = sizes ^ 3,
				covariate = as.factor(rep(ud[k],length(sizes))))
					
		if(length(grep("expsize", names(growObj@fit$coefficients))) == 1)
			newd$expsize = exp(sizes)
		if(length(grep("logsize", names(growObj@fit$coefficients))) == 1)
			newd$logsize = log(sizes)
		if(length(grep("logsize2", names(growObj@fit$coefficients))) == 1)
			newd$logsize=(log(sizes))^2
		
		
		if (length(grep("decline",tolower(as.character(class(growObj)))))>0 | 
				length(grep("trunc",tolower(as.character(class(growObj)))))>0) { 
				pred.size <- .predictMuX(growObj, newd, covPred = k)
			} else  {
				pred.size <- predict(growObj@fit,newd,type = "response")	
			}
		if (length(grep("incr", tolower(as.character(class(growObj))))) == 0) {
			points(sizes, pred.size, type = "l", col = k + 1)	
		} else { 
			if (length(grep("logincr", tolower(as.character(class(growObj))))) > 0) {
				points(sizes, sizes + exp(pred.size), type = "l", col = k + 1)		
			} else { 
				lines(sizes, pred.size, col = k + 1)	
			}
		}
	}
}

## FUNCTION FOR TURNING DATA INTO MATRIX DEFINING ENVIRONMENTAL TRANSITIONS ############################
## data is vector of env level at t, and one timestep later, at t+1

makeEnvObj <- function(dataf){
	#turn into index starting at 1
	minval <-  min(c(dataf$covariate,dataf$covariateNext),na.rm=TRUE)
	startEnv <- dataf$covariate-minval+1
	nextEnv <- dataf$covariateNext-minval+1
	
	
	nEnvClass <- max(c(startEnv,nextEnv), na.rm=TRUE)
	desired.mat <- matrix(0,nEnvClass,nEnvClass) 
	mats<-table(startEnv,nextEnv)
	rx <- as.numeric(rownames(mats));#print(rx)
	cx <- as.numeric(colnames(mats))
	desired.mat[cbind(rep(rx,length(cx)),rep(cx,each=length(rx)))]=c(as.matrix(mats))
	
	rc <- new("envMatrix",
			nEnvClass = nEnvClass)
	
	rc@.Data <- t(t(desired.mat)/colSums(desired.mat))
	
	return(rc) 
	
}



## FUNCTIONS FOR SIMULATING DATA ############################
## Return a data-frame with the headings used in the 'makeObject' functions


## Generate a simple data-frame for only continuous covariates
#  with a total of 1000 measurements with columns called
# size, sizeNext, surv, covariate, covariateNext, fec,
#
#
generateData <- function(nSamp=1000, type="simple"){
	if (type=="simple"){
	covariate <- sample(0:1, size=nSamp, replace=TRUE, prob = c(0.2, 0.8))
	covariateNext <- sample(0:1, size=nSamp, replace=TRUE, prob = c(0.8, 0.2))
	size <- rnorm(nSamp,5,2)
	#size <- exp(rnorm(1000, -1, 1.1))
	sizeNext <- 1+0.8*size-0.9*covariate+rnorm(nSamp,0,1)
	seedlings <- sample(1:nSamp,size=100,replace=TRUE)
	size[seedlings] <- NA; sizeNext[seedlings] <- rnorm(100,2,0.5)
	fec <- surv <- rep(NA, length(size))
	surv[!is.na(size)] <- rbinom(sum(!is.na(size)),1,logit(-1+0.2*size[!is.na(size)]))
	fec[!is.na(size)] <- rnorm(sum(!is.na(size)),exp(-7+0.9*size[!is.na(size)]),1)
	fec[size<quantile(size,0.20,na.rm=TRUE) | fec<0] <- 0
	fec <- fec*10
	
	stage <- stageNext <- rep("continuous",nSamp)
	stage[is.na(size)] <- NA
	stageNext[is.na(sizeNext)] <- "dead"
	
	dataf <- data.frame(size=size,sizeNext=sizeNext,surv=surv,
			covariate=covariate,covariateNext=covariateNext,
			fec=fec, stage=stage,stageNext=stageNext)
	
	dataf$sizeNext[dataf$surv==0] <- NA
	
	}
	
	if (type=="discrete") dataf <-.generateDataDiscrete(nSamp=nSamp)	
	if (type=="stochastic") dataf <-.generateDataStoch(nSamp=nSamp)	
	
	if (type!="simple" & type!="stochastic" & type!="discrete"){
		stop("unsupported type: supported types include simple, stochastic, or discrete")
	}
	
	return(dataf)
}

## Generate a simple data-frame for continuous and discrete covariates
#  with a total of 1000 measurements in columns called
# size, sizeNext, surv, fec, stage, stageNext number
# Stage contains names including "continuous", and then a range
# of names for discrete stages, e.g., in this example,
#  "dormant" "seedAge1"   "seedOld" 
#
# 
.generateDataDiscrete <- function(nSamp=1000){
	size <- rnorm(nSamp,5,2)
	sizeNext <- 1+0.8*size+rnorm(nSamp,0,1)
	surv <- rbinom(nSamp,1,logit(-1+0.2*size))
	sizeNext[surv==0] <- NA
	fec <- rnorm(length(size),exp(-7+0.9*size),1)
	fec[size<quantile(size,0.20) | fec<0] <- 0
	stage <- rep("continuous",nSamp)
	stageNext <- rep("continuous",nSamp)
	sizeNext[surv==0] <- NA
	stageNext[surv==0] <- c("dead")
	number <- rep(1,nSamp)
	become.dormant <- which(rank(size)%in%sample(rank(size),50,prob=surv*fec))
	sizeNext[become.dormant] <- NA; stageNext[become.dormant] <- c("dormant")
	were.dormant <- which(rank(sizeNext)%in%sample(rank(sizeNext),50,prob=surv*fec))
	size[were.dormant] <- NA; stage[were.dormant] <- c("dormant")
	dataf <- rbind(data.frame(size=size,sizeNext=sizeNext,surv=surv,
					fec=fec,stage=stage,stageNext=stageNext,number=number),
			data.frame(size=NA,sizeNext=NA,surv=rep(c(1,0),2),fec=0,
					stage=rep(c("seedAge1","seedOld"),each=2),
					stageNext=rep(c("seedOld","dead"),2),
					number=c(202,220,115,121)),
			data.frame(size=NA,sizeNext=rnorm(113,3,2),surv=1,fec=0,
					stage=c(rep("seedAge1",33),rep("seedOld",30),rep(NA,50)),
					stageNext=c("continuous"),number=1))
	
	
	return(dataf)
}

## Generate a simple data-frame for continuous and discrete covariates
#  with a total of 1000 measurements in columns called
# size, sizeNext, surv, fec, stage, stageNext number
#
# 
.generateDataStoch <- function(nSamp=1000){
	covariate1 <- rnorm(nSamp)
	covariate2 <- rnorm(nSamp)
	covariate3 <- rnorm(nSamp)
	size <- rnorm(nSamp,5,2)
	sizeNext <- 1+0.9*size+3*covariate1+0.01*covariate2+0.2*covariate3+rnorm(nSamp,0,0.1)
	
	fec <- surv <- rep(NA, length(size))
	surv[!is.na(size)] <- rbinom(sum(!is.na(size)),1,logit(-1+0.2*size[!is.na(size)]))
	fec[!is.na(size)] <- rnorm(sum(!is.na(size)),exp(-7+0.9*size[!is.na(size)]),1)
	fec[size<quantile(size,0.20,na.rm=TRUE) | fec<0] <- 0
	fec <- 10*fec
	
	seedlings <- sample(1:nSamp,size=100,replace=TRUE)
	size[seedlings] <- NA; 
	sizeNext[seedlings] <- rnorm(100,-2,0.1)
	surv[seedlings] <- 1
	#set to flower when covariate1 is around 1.5
	pfec <- 1*(runif(length(size))<logit(size+covariate1)); #print(pfec)
	fec[pfec==0] <- 0
	#fill in stage
	stage <- stageNext <- rep("continuous",nSamp)
	stage[is.na(size)] <- NA
	stageNext[is.na(sizeNext)] <- "dead"
	
	dataf <- data.frame(size=size,sizeNext=sizeNext,surv=surv,
			covariate1=covariate1,covariate2=covariate2,covariate3=covariate3,
			fec=fec, stage=stage,stageNext=stageNext, number=rep(1,length(size)))
	
	dataf$sizeNext[dataf$surv==0] <- NA
	
	return(dataf)
}


## FUNCTION FOR MAKING A LIST OF IPMS ############################################
# to do for stoch env with a single discrete covariate. ##########################

makeListIPMs <- function(dataf, nBigMatrix=10, minSize=-2,maxSize=10, 
		integrateType="midpoint", correction="none",
		explSurv=surv~size+size2+covariate,
		explGrow=sizeNext~size+size2+covariate, 
		regType="constantVar",explFec=fec~size,Family="gaussian", 
		Transform="none",fecConstants=data.frame(NA)) {
	
	#convert to 1:n for indexing later
	dataf$covariate <- as.factor(dataf$covariate)
	levels(dataf$covariate) <- 1:length(unique(dataf$covariate))
	
	print(explSurv)
	sv1 <- makeSurvObj(dataf=dataf,
			Formula=explSurv)
	gr1 <- makeGrowthObj(dataf=dataf,
			Formula=explGrow,
			regType=regType)
	
	fv1 <- makeFecObj(dataf=dataf,Formula=explFec, Family=Family, Transform=Transform, 
			fecConstants=fecConstants) 
	
	covs <- unique(dataf$covariate)
	covs <- covs[!is.na(covs)]
	
	#print(covs)
	
	IPM.list <- list()
	for (k in 1:length(covs)) { 
		
		tpF <- makeIPMFmatrix(nBigMatrix = nBigMatrix, minSize = minSize,
				maxSize = maxSize, chosenCov = data.frame(covariate=as.factor(k)),
				fecObj = fv1,integrateType=integrateType, correction=correction)
		tpS <- makeIPMPmatrix(nBigMatrix = nBigMatrix, minSize = minSize,
				maxSize = maxSize, chosenCov = data.frame(covariate=as.factor(k)),
				growObj = gr1, survObj = sv1,
				integrateType=integrateType, correction=correction)
		IPM.list[[k]] <- tpF+tpS
	}
	return(IPM.list)
	
}



### simulation Carlina ######################################################

## Function to simulate something a bit like Carlina
##
## Parameters - nSamp - starting pop sizes
#             - nYrs - no years in simulation
#             - nSampleYrs - no years to extract for the data
#             - ... bunch of parameters
#             - meanYear - means for year effects
#             - matVarYear - variance covariances for year effects
#             - densDep - model density dependence in seedling establishment or not
#
# Returns - list including dataf - data-frame of data
#                                - various of the simulation parameters for convenience
#
#    PARAMETERS FROM TABLE 2 IN REES AND ELLNER 2009
simulateCarlina <- function(nSamp=200,nYrs=1000,nSampleYrs=15,
		m0=-1.37,ms=0.59,
		b0=-12.05,bs=3.64,
		A=-1,B=2,
		ag=1.14,bg=0.74,sig=0.29,
		mean.kids=3.16,sd.kids=0.5,
		meanYear=c(0,0,0),
		matVarYear=matrix(c(1.03,0,0,0,0.037,0.041,0,0.041,0.075),3,3), 
		densDep=TRUE,maxPerYr=1000,maxStoreSeedlingsPerYr=200,sizes=c()) {
		
	#initiate and set up year index
	if (length(sizes)==0) sizes <- rnorm(nSamp,3,0.5)
	startYr <- (nYrs-nSampleYrs)
	
	#recruits
	nrec <- c(20,42,12,17,8,19,58,45,44,2,56,25,75,92,94,6,4,34,104)
	
	#total plants
	totpl <- c(21,57,47,25,25,33,88,94,97,26,85,80,122,175,160,10,6,189)
		
	#set up dataframe
	dataf <-matrix(NA,(maxPerYr+maxStoreSeedlingsPerYr)*nYrs,11)
	maxPop <- nrow(dataf)
	n.per.yr  <-  rep(NA,nYrs)
	trueGrow <- rep(NA,nYrs)	
	
	count <- 0

	#stoch sims
	tmp <- rmvnorm(nYrs,mean=meanYear,sigma=matVarYear)
	
	for (t in 1:nYrs) {
		
		#yr effects
		nSeedlings <- sample(nrec,size=1,replace=FALSE)
		#print(tmp)
		n.per.yr[t] <- length(sizes)

		m.year <- tmp[t,1]
		cg.year <- tmp[t,2]
		b.year <- tmp[t,3]
		
		if (n.per.yr[t]>0) { 
			#survival
			sx <- 1*(logit(m0+ms*sizes+m.year)>runif(n.per.yr[t]))
			#flowering
			fx <- 1*(logit(b0+bs*sizes)>runif(n.per.yr[t]))
			#fertility
			seedsx <- exp(A+B*sizes)*fx*sx
			#seedsx[seedsx>0] <- rpois(sum(seedsx>0),seedsx[seedsx>0])
		} else {
			sx <- fx <- seedsx <- c()
			print(c("extinct in year ", t))
		}
		
		if (densDep) pEst <- min(nSeedlings/max(sum(seedsx, na.rm=TRUE),1),1) else pEst <- 1
		
		babies <- rnorm(ceiling(pEst*max(sum(seedsx,na.rm=TRUE),1)),
				mean.kids+b.year,sd.kids) 
		#if (count>0) print(c("yr",t-startYr,pEst*max(sum(seedsx,na.rm=TRUE),1),length(babies)))
		#will end up with nrec babies at least in density dependent case
		if (length(babies)<nSeedlings & densDep) nSeedlings <- length(babies)
		
		#growth
		sizeNext <- rnorm(length(sizes),ag+bg*sizes+cg.year,sig)
		
		#remove dead or flowered
		fx[sx==0] <- NA
		sizeNext[sx==0 | fx==1] <- NA
		
		#storage
		if (t > startYr) {
			#print(t)
			#don't 'do it at all if too big (so as not to do only adults, etc)
			if ((count+length(sizes))>maxPop) { print("large pop size, breaking");break()}
			if ((count+length(sizes)+length(babies))>maxPop) { print("large pop size, breaking");break()}
			
			chs <- (count+1):(count+length(sizes))
			dataf[chs,1] <- sizes
			dataf[chs,2] <- sizeNext
			dataf[chs,3] <- 1*sx
			dataf[chs,4] <- 1*fx
			dataf[chs,5] <- seedsx
			dataf[chs,6] <- t
			dataf[chs,7] <- nSeedlings
			dataf[chs,8] <- m.year
			dataf[chs,9] <- cg.year
			dataf[chs,10] <- b.year
			
			count <- count + length(sizes)
			
			print(lenth(sizes))	
			
			nbabes.store <- min(length(babies),maxStoreSeedlingsPerYr)
			chs <- (count+1):(count+nbabes.store)
			
			dataf[chs,2] <- babies[1:nbabes.store]
			dataf[chs,6] <- t
			dataf[chs,7] <- nSeedlings
			dataf[chs,8] <- m.year
			dataf[chs,9] <- cg.year
			dataf[chs,10] <- b.year
			dataf[chs,11] <- "sexual"
			
			count <- count + nbabes.store
					
			print(nbabes.store)	
			
			} 
		
						
			#new pop - i.e. those that grew, survive, and didn't flower; as well as babies
			sizes <- c(sizeNext[sx==1 & fx==0 & !is.na(fx)],babies)
			if (length(sizes)==0) print("extinct")

			#get true growth rate 
			trueGrow[t] <- log(length(sizes)/n.per.yr[t])
			
			#cull size at t t down to 1000, if not density dependent
			if (!densDep) sizes <- sample(sizes,size=maxPerYr, replace=TRUE)

	}
	
	
	list.par <- list(m0=m0,ms=ms,
			b0=b0,bs=bs,
			A=A,B=B,
			ag=ag,bg=bg,sig=sig,
			mean.kids=mean.kids,sd.kids=sd.kids,
			meanYear=meanYear,
			matVarYear=matVarYear,
			nrec=nrec)
	
	dataf <- dataf[1:count,]
	
	dataf <- data.frame(dataf,stringsAsFactors = FALSE)	
	#print(table(dataf[,6]))

	#plot(cumsum(trueGrow[startYr:nYrs])/(1:length(trueGrow[startYr:nYrs])),type="l")
	trueGrow <- mean(trueGrow[startYr:nYrs],na.rm=TRUE)  
	#abline(h=trueGrow)
	
	colnames(dataf) <- c("size","sizeNext","surv","flower","fec",
			"year","nSeedlings","m.year","cg.year","b.year","offspringNext")
	
	dataf$size <- as.numeric(dataf$size)
	dataf$sizeNext <- as.numeric(dataf$sizeNext)
	dataf$surv <- as.numeric(dataf$surv)
	dataf$flower <- as.numeric(dataf$flower)
	dataf$fec <- as.numeric(dataf$fec)
	dataf$nSeedlings <- as.numeric(dataf$nSeedlings)
	dataf$year <- as.numeric(dataf$year)
	dataf$cg.year <- as.numeric(dataf$cg.year)
	dataf$m.year <- as.numeric(dataf$m.year)
	dataf$b.year <- as.numeric(dataf$b.year)	
	
	#print(table(dataf$year))
	
	dataf$fec[dataf$fec==0] <- NA
	
	return(list(dataf=dataf,meanYear=meanYear,matVarYear=matVarYear,list.par=list.par,trueGrow=trueGrow))
	
}


#Find years where can estimate all three stochastic vital rates(survival, growth and baby size)
.identifyPossibleYearsCarlina <- function(dataf){
	
	yr1 <- table(dataf$year[!is.na(dataf$size) & 
							!is.na(dataf$sizeNext) & is.na(dataf$offspringNext)])
	yr2 <- table(dataf$year[!is.na(dataf$size) & 
							!is.na(dataf$surv) & is.na(dataf$offspringNext)])
	yr3 <- table(dataf$year[!is.na(dataf$sizeNext) & 
							!is.na(dataf$offspringNext)])
	
	good.yrs <- intersect(as.numeric(as.character(names(yr1)[yr1>1])),
			as.numeric(as.character(names(yr2))[yr2>1]))
	good.yrs <- intersect(good.yrs,as.numeric(as.character(names(yr3)[yr3>1])))
	
	return(is.element(dataf$year,good.yrs))
}



## Convert Increment - where exact dates of census vary but some multiplier of yearly increments
## are desired; this function takes a data-frame (with columns size, sizeNext,
## and, importantly exactDate, exactDatel)
## and returns a data-frame with sizeNext modified proportional to the time
## elapsed in the desired yaerly increments, adding an additional column denoted 'increment'. 
#
# Parameters - dataf - a dataframe with headings size, sizeNext, exactDate,exactDatel
#            - nYrs - the number of years between censuses desired (e.g. for CTFS data, 5 years)
#
# Returns - a dataframe with the same headings
convertIncrement <- function(dataf, nYrs=1) {
	incr <- dataf$sizeNext - dataf$size
	if (class(dataf$exactDatel) == "Date")  {
		timeElapsed <- (difftime(dataf$exactDatel,dataf$exactDate)[[1]])/(365*nYrs)
	} else {
		timeElapsed <- (dataf$exactDatel-dataf$exactDate)/(365*nYrs)
	}
	incrNew <- incr / timeElapsed
	dataf$sizeNext <- dataf$size + incrNew
	dataf$incr <- dataf$incrNew
	return(dataf)
	
}



# Function to plot the results of a stochastic simulation
# structure run
#
# Parameters - tVals - time points
#            - st - output of trackPopStructManyCovSeedBank or trackPopStructManyCov
#            - covTest - the key covariate for germination / flowering
#            - nRunIn - how many to leave off pics
# Returns - 
#
plotResultsStochStruct <- function(tVals,st,covTest,nRunIn=15,log="y",...) { 
	
	par(mfrow=c(2,2),bty="l")
	plot(tVals[nRunIn:length(tVals)],
			colSums(st$rc[,nRunIn:length(tVals)]+1),
			xlab="Time", 
			ylab="Population size",type="l",log=log,...)
	abline(v=1:max(tVals),lty=3)
	covTestplot <- exp(mean(colSums(st$rc[,nRunIn:length(tVals)])) +
					((covTest-mean(covTest))/sd(covTest))*
					sd(colSums(st$rc[,nRunIn:length(tVals)])))
	points(tVals,covTestplot+1,type="l",lty=3,col=2)
	
	if (log=="y") st$rc <- log(st$rc)
	
	image(tVals[nRunIn:length(tVals)],
			st$IPM.here@meshpoints,
			t(st$rc[,nRunIn:length(tVals)]+1),
			ylab="Size", xlab="Time")
}




# makeCovDf creates a dataframe of size variables for prediction
# TODO: make able to use 'covariates'
.makeCovDf <- function(size, explanatoryVariables) {
	sizeSorted <- sort(size)
	splitVars <- strsplit(explanatoryVariables, split = "\\+")
	expVar <- unlist(splitVars[grep("size", splitVars)])
	covDf <- data.frame(size = sizeSorted)
	for(i in 1:length(expVar)) {
		if(is.null(expVar[i])) next
		expVar[i] <- sub('[[:space:]]', '', expVar[i])
		if(expVar[i] == "size2") {
			covDf$size2 <- sizeSorted ^ 2
		}
		if(expVar[i] == "size3"){
			covDf$size3 <- sizeSorted ^ 3
		}
		if(expVar[i] == "expsize") {
			covDf$expsize <- exp(sizeSorted)
		}
		if(expVar[i] == "logsize") {
			covDf$logsize <- log(sizeSorted)
		}
		if(expVar[i] == "logsize2") {
			covDf$logsize2 <- log(sizeSorted ^ 2)
		}
	}
	return(covDf)
}


# Function to compare model fits for growth and survival objects built with different linear combinations of covariates. 
# Growth can have multiple response forms. Uses .makeCovDf.  Separate plot functions. 
#
#
# Returns - list with a summary table of covariates and scores, and a sub-list of growth or survival objects.
#
growthModelComp <- function(dataf, 
		expVars = c(sizeNext~1, sizeNext~size, sizeNext~size + size2), 
		regressionType = "constantVar",
		testType = "AIC",
		makePlot = FALSE,
		mainTitle = "",
		legendPos = "topright", ...) {
	varN <- length(expVars)
	typeN <- length(regressionType)
	treatN <- varN * typeN
	summaryTable <- data.frame()
	grObj <- vector("list", length = treatN)
	i <- 1
	for(v in 1:varN) {
		for(t in 1:typeN) {
			grObj[[i]] <- makeGrowthObj(dataf = dataf, regType = regressionType[t], Formula = expVars[[v]]) 
			if (length(grep("decline",tolower(as.character(class(grObj[[i]])))))>0) { 
				summaryTable <- rbind(summaryTable, cbind(as.character(unlist(expVars))[v], regressionType[t],  match.fun(testType)(grObj[[i]]@fit$fit)))
			} else { 
			summaryTable <- rbind(summaryTable, cbind(as.character(unlist(expVars))[v], regressionType[t],  match.fun(testType)(grObj[[i]]@fit)))
			}	
			i <- i + 1
		}
	}
	summaryTable <- as.data.frame(summaryTable)
	names(summaryTable) <- c("Exp. Vars", "Reg. Type", testType)
	outputList <- list(summaryTable = summaryTable, growthObjects = grObj)
	
	# PLOT SECTION #
	if(makePlot == TRUE) {
		plotGrowthModelComp(grObj = grObj, summaryTable = summaryTable, dataf = dataf, expVars = expVars,  testType = "AIC",  plotLegend = TRUE, mainTitle = mainTitle, legendPos, ...)
	}
	return(outputList)
}


survModelComp <- function(dataf, 
		expVars = c(surv~1, surv~size, surv~size + size2), 
		testType = "AIC",
		makePlot = FALSE,
		mainTitle = "", ncuts = 20,
		legendPos = "bottomleft", ...) {
	varN <- length(expVars)
	treatN <- varN
	summaryTable <- data.frame()
	svObj <- vector("list", length = treatN)
	i <- 1
	for(v in 1:varN) {
		svObj[[i]] <- makeSurvObj(dataf = dataf, Formula = expVars[[v]]) 
		summaryTable <- rbind(summaryTable, cbind(as.character(unlist(expVars))[v], match.fun(testType)(svObj[[i]]@fit)))
		i <- i + 1
	}
	summaryTable <- as.data.frame(summaryTable)
	names(summaryTable) <- c("Exp. Vars", testType)
	outputList <- list(summaryTable = summaryTable, survObjects = svObj)
	
	# PLOT SECTION #
	if(makePlot == TRUE) {
		## this is the surv picture    
		plotSurvModelComp(svObj = svObj, summaryTable = summaryTable, dataf = dataf, expVars = expVars, testType = "AIC", plotLegend = TRUE, mainTitle = mainTitle, ncuts = ncuts, legendPos = legendPos, ...)
	}
	return(outputList)
}	

# Plot functions for model comparison.  Plots the series of fitted models for growth and survival objects.  
# Can plot a legend with the model covariates and model test criterion scores (defaults to AIC).

plotGrowthModelComp <- function(grObj, summaryTable, dataf, expVars,  testType = "AIC", plotLegend = TRUE, mainTitle = "", legendPos = "topright", ...) {
	treatN <- length(grObj)
	sizeSorted <- unique(sort(dataf$size))
	if(length(grep("sizeNext", unlist(as.character(expVars))))) {
		y.lab <- "Size at t + 1"
		dataSizeNext <- dataf$sizeNext
	}
	if(length(grep("incr", unlist(as.character(expVars))))) {
		y.lab <- "Growth"
		dataSizeNext <- dataf$sizeNext - dataf$size
	}
	if(length(grep("logincr", unlist(as.character(expVars))))){
		y.lab <- "log(growth)"
		dataSizeNext <- log(dataf$sizeNext - dataf$size)
	}
	plot(dataf$size, dataSizeNext, pch = 19, xlab = "Size at t", ylab = y.lab, main = mainTitle, cex = 0.8,...)
	for(p in 1:treatN) {
		
		#PROBLEM HERE 
		expVarsHere <- paste(attr(terms((expVars[[p]])),"term.labels"),collapse="+")
	
		newd <- .makeCovDf(sizeSorted, expVarsHere)
		if (length(grep("decline",tolower(as.character(class(grObj[[p]])))))>0) {
			pred.size <- .predictMuX(grObj[[p]], newd)
		} else { pred.size <- predict(grObj[[p]]@fit, newd, type = "response")}
		lines(sizeSorted, pred.size, type = "l", col = (p + 1))
	}
	if(plotLegend) {
		legend(legendPos, legend = sprintf("%s: %s = %.1f", expVars, testType, as.numeric(as.character(summaryTable[,3]))), col = c(2:(p + 1)), lty = 1, xjust = 1, bg = "white")
	}
}

plotSurvModelComp <- function(svObj, summaryTable, dataf,  expVars, testType = "AIC", plotLegend = TRUE, mainTitle = "", ncuts = 20, legendPos = "bottomleft", ...) {
	treatN <- length(svObj)
	#ncuts <- 20  # survival bins
	os <- order(dataf$size)  # order size
	osSurv <- (dataf$surv)[os] # order survival data according to size
	osSize<-(dataf$size)[os] # ordered size data
	binnedSize <- tapply(osSize, as.numeric(cut(osSize, breaks=ncuts)), mean, na.rm = TRUE); # bin Size data
	binnedSurv <- tapply(osSurv, as.numeric(cut(osSize, breaks=ncuts)), mean, na.rm = TRUE) #bin Survival probabilities
	plot(binnedSize, binnedSurv, pch = 19, xlab = "Size at t", ylab = "Survival to t + 1", main = mainTitle, cex = 0.8,...)
	for(p in 1:treatN) {
		expVarsHere <- paste(attr(terms(expVars[[p]]),"term.labels"),collapse="+")
		newd <- .makeCovDf(osSize, expVarsHere[p])
		lines(dataf$size[order(dataf$size)], surv(dataf$size[os], data.frame(covariate=1), svObj[[p]]), col = (p + 1))           
	}
	if(plotLegend) {
		legend(legendPos, legend = sprintf("%s: %s = %.1f", expVars, testType, as.numeric(as.character(summaryTable[,2]))), col = c(2:(p + 1)), lty = 1, xjust = 1, bg = "white")
	}
}




## Function to add pdf to model comparison pics
addPdfGrowthPic <- function(respType = "sizeNext",
		sizesPlotAt=c(20,50,60),
		sizeRange=c(20,400),
		incrRange=c(-10,50), 
		scalar=100,
		growthObjList,
		cols=1:5,
		cov=data.frame(covariate=1),
		minShow=1e-2,
		jitt=2,  #how far apart should pdfs be if you are plotting several
		...){
	
	nval <- length(growthObjList)
	sizes <- seq(sizeRange[1],sizeRange[2],length=500)
	incr <- seq(incrRange[1],incrRange[2],length=500)
	logincr <- log(incr); logincr[!is.finite(logincr)] <- NA
	
	for (j in 1:nval) {
		#print(j)
		for (k in 1:length(sizesPlotAt)) {           
			if (respType=="sizeNext") {
				pred <- growth(sizesPlotAt[k],sizes,cov,growthObjList[[j]])*scalar
				pred[pred<minShow] <- NA
				points(sizesPlotAt[k]+pred+jitt*(j-1),
						sizes,type="l", col=cols[j],...)                
			}
			if (respType=="incr") {
				pred <- growth(sizesPlotAt[k],sizesPlotAt[k]+incr,cov,growthObjList[[j]])*scalar
				pred[pred<minShow] <- NA
				points(sizesPlotAt[k]+pred+jitt*(j-1),
						incr,type="l", col=cols[j],...)
			}
			if (respType=="logincr") {
				pred <- growth(sizesPlotAt[k],sizesPlotAt[k]+logincr,cov,growthObjList[[j]])*scalar
				pred[pred<minShow] <- NA
				#print(range(pred,na.rm=TRUE))
				points(sizesPlotAt[k]+pred+jitt*(j-1),
						logincr,type="l", col=cols[j],...)
			}
		}}
}

## Function to take fit of these and output a list of growth objects
getListRegObjects <- function(Obj,nsamp=1000) {
	
	require(mvtnorm)
	require(MASS)
	
	#generate new set parameters from mvn
	npar <- length(Obj@fit$coefficients)
	newpar <- rmvnorm(nsamp, mean = Obj@fit$coefficients, sigma = vcov(Obj@fit))
	
	objList <- list()
	
	for (j in 1:nsamp) {
		objList[[j]] <- Obj
		objList[[j]]@fit$coefficients <- newpar[j,]
	}
	
	return(objList)
}



## Function to take fit of these and output a list of growth objects
getListRegObjectsFec <- function(Obj,nsamp=1000) {
	
	require(mvtnorm)
	require(MASS)

	objList <- list()
	
	#generate new set parameters from mvn
	for (j in 1:nsamp) {
		for (k in 1:length(Obj@fitFec)) { 
		npar <- length(Obj@fitFec[[k]]$coefficients)
		newpar <- rmvnorm(nsamp, mean = Obj@fitFec[[k]]$coefficients, 
					sigma = vcov(Obj@fitFec[[k]]))
		objList[[j]] <- Obj
		objList[[j]]@fitFec[[k]]$coefficients <- newpar[j,]
		}
	}	
		
	return(objList)
}



## Function to coerce Growth object to parameters and variance desired
coerceGrowthObj <- function(growthObj, coeff, sd){

	if (length(growthObj@fit$coefficients) !=length(coeff)) {
		print("warning: number of desired coefficients to not match number of growth object coefficients")
	} 	
	growthObj@fit$coefficients[] <- as.numeric(coeff)
	growthObj@sd[] <- as.numeric(sd)
	return(growthObj)
}


## Function to coerce Survival object to parameters desired
coerceSurvObj <- function(survObj, coeff){
	
	if (length(survObj@fit$coefficients) !=length(coeff)) print("warning: number of desired coefficients to not match number of growth object coefficients")
	survObj@fit$coefficients[] <- as.numeric(coeff)
	
	return(survObj)
}





# Function to extract IPM output from a list
# of T (survival + growth) and F (fecundity) matrices
# (usually from Bayes fit) 
#
# Parameters - Tmatrixlist
#            - target.size - the size you want passage time estimated for.
#            - Fmatrixlist
#
# Returns - a list 

getIPMoutput <- function(Tmatrixlist,target.size=c(),Fmatrixlist=NULL){
	
	if (length(target.size)==0)  { 
		print("no target size for passage time provided; taking meshpoint median")
		target.size <- median(TmatrixList[[1]]@meshpoints)
	}
	nsamps <- length(Tmatrixlist)
	h1 <- Tmatrixlist[[1]]@meshpoints[2]-Tmatrixlist[[1]]@meshpoints[1]
	stable.size <- LE <- ptime <- matrix(NA,nsamps,length(Tmatrixlist[[1]]@.Data[1,]))
	lambda <- rep(NA,nsamps)
	for (k in 1:nsamps) {
		Tmatrix <- Tmatrixlist[[k]]
		LE[k,]<-MeanLifeExpect(Tmatrix) 
		ptime[k,]<-PassageTime(target.size,Tmatrix) 
		
		if (class(Fmatrixlist)!="NULL") {
			IPM <- Tmatrix + Fmatrixlist[[k]]
			lambda[k] <- Re(eigen(IPM)$value[1])
			stable.size[k,] <- eigen(IPM)$vector[,1]
			#normalize stable size distribution
			stable.size[k,] <- stable.size[k,]/(h1*sum(stable.size[k,]))
		}
	}
	
	return(list(LE=LE,ptime=ptime,lambda=lambda,stable.size=stable.size))
	
}


# Function to extract IPM output from posteriors 
# (usually from Bayes fit)  - quicker to do all at once
# rather than build list of IPM T matrices, then list of IPM F matrices
#
# Parameters - survObjlist - list of survival objects
#            - growObjList - list of growth objects
#            - target.size - the size you want passage time estimated for.
#            - nBigMatrix - the number of bins
#            - minSize - the minimum size
#            - maxSize - the maximum size
#            - cov - do you want to fit a discrete covariate
#            - fecObjList - list of fecundity objects
#            - env.mat - matrix of env transitions (only if cov=TRUE)
#            - n.size.to.age - numeric describing how many size to age defined (0 - 100s)
# 
#
# Returns - a list 

getIPMoutputDirect <- function(survObjList,growObjList,target.size=c(),
		nBigMatrix,minSize,maxSize,DiscreteTrans = 1,
		cov=FALSE,fecObjList=NULL, env.mat=NULL,
		n.size.to.age=0, size.start=10,
		integrate.type="midpoint", correction="none"){
	
	# adjust the sample lengths to they are all the same
	if (length(target.size)==0)  target.size <- 0.2*(minSize+maxSize)
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
	surv.par <- matrix(NA,nsamp,length(survObjList[[1]]@fit$coefficients))
	grow.par <- matrix(NA,nsamp,length(growObjList[[1]]@fit$coefficients)+1)
	for (k in 1:nsamp) {
		surv.par[k,] <- survObjList[[k]]@fit$coefficients
		grow.par[k,] <- c(growObjList[[k]]@fit$coefficients,sd(growObjList[[k]]@fit$residuals))
	}
	
	#set up storage
	if (is.data.frame(DiscreteTrans)) ndisc <- ncol(DiscreteTrans) else ndisc <- 0
	if (class(env.mat)!="NULL") n.env <- env.mat@n.env.class else n.env <- 1
	LE <- ptime <- matrix(NA,nsamp,(nBigMatrix+ndisc)*n.env)
	if (class(fecObjList)=="NULL") {
		lambda <- stable.size <- c()
	} else {
		stable.size <- matrix(NA,nsamp,(nBigMatrix+ndisc)*n.env)
		lambda <- rep(NA,nsamp)
	}
	if (n.size.to.age==0) { resAge <- resSize <- c() } else { resAge <- resSize <- matrix(NA,nsamp,n.size.to.age)} 
	if (length(size.start)==0) { if (minSize<0) size.start <- 0.5*minSize else size.start <- 2*minSize }
	
	#go!
	for (k in 1:nsamp) {
		
		if (!cov) {
			Tmatrix <- create.IPM.Tmatrix(nBigMatrix = nBigMatrix, minSize = minSize, 
					maxSize = maxSize, growObj = growObjList[[k]],
					survObj = survObjList[[k]],DiscreteTrans=DiscreteTrans,
					integrate.type=integrate.type, correction=correction) 
			
		} else {
			Tmatrix <- create.compound.Tmatrix(n.env.class = n.env,
					nBigMatrix = nBigMatrix, minSize = minSize, 
					maxSize = maxSize, envMatrix=env.mat,growObj = growObjList[[k]],
					survObj = survObjList[[k]],DiscreteTrans=DiscreteTrans,
					integrate.type=integrate.type, correction=correction)    
	
		}
		
		LE[k,] <- MeanLifeExpect(Tmatrix) 
		ptime[k,] <- PassageTime(target.size,Tmatrix) 
		if (k==1) h1 <- diff(Tmatrix@meshpoints)[1]
		
		if (class(fecObjList)!="NULL") {
			if (!cov) { 
				Fmatrix <- create.IPM.Fmatrix(nBigMatrix = nBigMatrix, minSize = minSize, 
						maxSize = maxSize, 
						fecObj=fecObjList[[k]],
						integrate.type=integrate.type, correction=correction)
			} else {
				Fmatrix <- create.compound.Fmatrix(n.env.class = n.env,
						nBigMatrix = nBigMatrix, minSize = minSize, 
						maxSize = maxSize, envMatrix=env.mat,
						fecObj=fecObjList[[k]],integrate.type=integrate.type, correction=correction)
			}

	
			
			IPM <- Tmatrix + Fmatrix
			lambda[k] <- Re(eigen(IPM)$value[1])
			stable.size[k,] <- eigen(IPM)$vector[,1]
			#normalize stable size distribution
			stable.size[k,] <- stable.size[k,]/(h1*sum(stable.size[k,]))
			
			#print("here2")
		}
		
		# get size to age results
		if (n.size.to.age>0) { 
			res2 <- SizeToAge(Tmatrix=Tmatrix,starting.size=minSize*1.3,
					target.size=seq(size.start,maxSize*0.9,length=n.size.to.age))
			resAge[k,] <- res2$timeInYears
			resSize[k,] <- res2$target.size
		}
		
	}
	
	return(list(LE=LE,ptime=ptime,lambda=lambda,stable.size=stable.size,
					meshpoints=Tmatrix@meshpoints,resAge=resAge,resSize=resSize,
					surv.par=surv.par,grow.par=grow.par))
	
}



## Function to get passage time FROM a particular size TO a range of sizes
## (i.e. size to age) when provided with a Tmatrix, a starting size, and a list
## of target sizes
#
# Parameters - Tmatrix
#            - starting.size
#            - target.sizes
#
# Returns - list containing vector of targets, vector of corresponding times, and the starting.size
#
SizeToAge <- function(Tmatrix,starting.size,target.size) {
	
	#locate where the first size is in the meshpoints of Tmatrix
	diffv <- abs(starting.size-Tmatrix@meshpoints)
	start.index <- median(which(diffv==min(diffv),arr.ind=TRUE))
	timeInYears <- rep(NA,length(target.size))
	
	#loop over to see where its going
	for (k in 1:length(target.size)) {
		ptime <- PassageTime(target.size[k],Tmatrix)
		timeInYears[k] <- ptime[start.index]
	}
	
	return(list(timeInYears=timeInYears,target.size=target.size,starting.size=starting.size))
	
}






## FUNCTION FOR MAKING PICTURE OF SURVIVAL DATA AND FITTED SURVIVAL OBJECT ############################
## and PICTURE GROWTH DATA AND FITTED GROWTH OBJECT  
# ** note that growth may need redefining if new growth  objects are defined


# parameters - dataf - the dataframe
#            - survObj - a survival object
#            - ncuts - the number of cuts used for binning survival
# returns - 

picSurv <- function(dataf,survObj,ncuts=20,...) { 
	
	#organize data and plot mean of ncut successive sizes, so trend more obvious
	os<-order(dataf$size); os.surv<-(dataf$surv)[os]; os.size<-(dataf$size)[os]; 
	psz<-tapply(os.size,as.numeric(cut(os.size,ncuts)),mean,na.rm=TRUE); #print(psz)
	ps<-tapply(os.surv,as.numeric(cut(os.size,ncuts)),mean,na.rm=TRUE);#print(ps)
	
	if (length(grep("covariate",names(survObj@fit$model)))==0) {  
		#plot data
		plot(as.numeric(psz),as.numeric(ps),pch=19,
				xlab="Size at t", ylab="Survival to t+1",main="Survival",...)
		#Plot fitted models
		points(dataf$size[order(dataf$size)],surv(dataf$size[order(dataf$size)],1,survObj),type="l",col=2)
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
			psz<-tapply(os.size[tp],as.numeric(cut(os.size[tp],ncuts)),mean,na.rm=TRUE); #print(psz)
			ps<-tapply(os.surv[tp],as.numeric(cut(os.size[tp],ncuts)),mean,na.rm=TRUE);#print(ps)
			points(as.numeric(psz),as.numeric(ps),pch=19,col=k)
			newd <- data.frame(size=sizes,size2=sizes^2,size3=sizes^3,
					covariate=rep(as.factor(ud[k]),length(sizes)))
			if(length(grep("logsize",survObj@fit$formula)))
				newd$logsize=log(sizes)
			if(length(grep("logsize2",survObj@fit$formula)))
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

picGrow <- function(dataf,growObj) {
	
	
	plot(dataf$size,dataf$sizeNext,pch=19,xlab="Size at t", ylab="Size at t+1",main="Growth")
	
	if (length(grep("covariate",names(growObj@fit$model)))>0) {  
		#convert to 1:n for indexing later and to relate to discrete
		dataf$covariate <- as.factor(dataf$covariate)
		levels(dataf$covariate) <- 1:length(unique(dataf$covariate))
		ud <- unique(dataf$covariate); ud <- ud[!is.na(ud)]
		for (k in 1:length(ud)) { 
			tp <- dataf$covariate==ud[k]
			points(dataf$size[tp], dataf$sizeNext[tp],pch=19,col=k)            
		}
		ud <- as.factor(ud)
	} else {
		ud <- 0
	} 
	
	sizes <- dataf$size[!is.na(dataf$size)]; sizes <- sizes[order(sizes)]
	
	for (k in 1:length(ud)) { 
		newd <- data.frame(size=sizes,size2=sizes^2,
				covariate=as.factor(rep(ud[k],length(sizes))))
		#print(ud[k])
		
		if(length(grep("logsize",growObj@fit$formula)))
			newd$logsize=log(sizes)
		if(length(grep("logsize2",growObj@fit$formula)))
			newd$logsize=(log(sizes))^2
		
		pred.size <- predict(growObj@fit,newd,type="response")
		if (class(growObj)=="growthObj" | class(growObj)=="growthObj.declinevar") {
			points(sizes,pred.size,type="l",col=k)
		} else { if (class(growObj)=="growthObj.incr") {
				points(sizes,sizes+pred.size,type="l",col=k)
			}}
	}
	
}

## FUNCTION FOR TURNING DATA INTO MATRIX DEFINING ENVIRONMENTAL TRANSITIONS ############################
## data is vector of env level at t, and one timestep later, at t+1

makeEnvObj <- function(dataf){
	#turn into index starting at 1
	minval <-  min(c(dataf$covariate,dataf$covariateNext))
	startEnv <- dataf$covariate-minval+1
	nextEnv <- dataf$covariateNext-minval+1
	
	
	n.env.class <- max(c(startEnv,nextEnv))
	desired.mat <- matrix(0,n.env.class,n.env.class) 
	mats<-table(startEnv,nextEnv)
	rx <- as.numeric(rownames(mats));#print(rx)
	cx <- as.numeric(colnames(mats))
	desired.mat[cbind(rep(rx,length(cx)),rep(cx,each=length(rx)))]=c(as.matrix(mats))
	
	rc <- new("env.matrix",
			n.env.class = n.env.class)
	
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
generateData <- function(){
	covariate <- sample(0:1, size=1000, replace=TRUE, prob = c(0.2, 0.8))
	covariateNext <- sample(0:1, size=1000, replace=TRUE, prob = c(0.8, 0.2))
	size <- rnorm(1000,5,2)
	#size <- exp(rnorm(1000, -1, 1.1))
	sizeNext <- 1+0.8*size-0.9*covariate+rnorm(1000,0,1)
	seedlings <- sample(1:1000,size=100,replace=TRUE)
	size[seedlings] <- NA; sizeNext[seedlings] <- rnorm(100,-2,0.1)
	fec <- surv <- rep(NA, length(size))
	surv[!is.na(size)] <- rbinom(sum(!is.na(size)),1,logit(-1+0.2*size[!is.na(size)]))
	fec[!is.na(size)] <- rnorm(sum(!is.na(size)),exp(-7+0.9*size[!is.na(size)]),1)
	fec[size<quantile(size,0.20,na.rm=TRUE) | fec<0] <- 0
	fec <- fec*10
	
	stage <- stageNext <- rep("continuous",1000)
	stage[is.na(size)] <- NA
	stageNext[is.na(sizeNext)] <- "dead"
	
	dataf <- data.frame(size=size,sizeNext=sizeNext,surv=surv,
			covariate=covariate,covariateNext=covariateNext,
			fec=fec, stage=stage,stageNext=stageNext)
	return(dataf)
}

## Generate a simple data-frame for continuous and discrete covariates
#  with a total of 1000 measurements in columns called
# size, sizeNext, surv, fec, stage, stageNext number
# Stage contains names including "continuous", and then a range
# of names for discrete stages, e.g., in this example,
#  "dormant" "seed.age.1"   "seed.old" 
#
# 
generateDataDiscrete <- function(){
	size <- rnorm(1000,5,2)
	sizeNext <- 1+0.8*size+rnorm(1000,0,1)
	surv <- rbinom(1000,1,logit(-1+0.2*size))
	sizeNext[surv==0] <- NA
	fec <- rnorm(length(size),exp(-7+0.9*size),1)
	fec[size<quantile(size,0.20) | fec<0] <- 0
	stage <- rep("continuous",1000)
	stageNext <- rep("continuous",1000)
	sizeNext[surv==0] <- NA
	stageNext[surv==0] <- c("dead")
	number <- rep(1,1000)
	become.dormant <- which(rank(size)%in%sample(rank(size),50,prob=surv*fec))
	sizeNext[become.dormant] <- NA; stageNext[become.dormant] <- c("dormant")
	were.dormant <- which(rank(sizeNext)%in%sample(rank(sizeNext),50,prob=surv*fec))
	size[were.dormant] <- NA; stage[were.dormant] <- c("dormant")
	dataf <- rbind(data.frame(size=size,sizeNext=sizeNext,surv=surv,
					fec=fec,stage=stage,stageNext=stageNext,number=number),
			data.frame(size=NA,sizeNext=NA,surv=rep(c(1,0),2),fec=0,
					stage=rep(c("seed.age.1","seed.old"),each=2),stageNext=rep(c("seed.old","dead"),2),
					number=c(202,220,115,121)),
			data.frame(size=NA,sizeNext=rnorm(113,3,2),surv=1,fec=0,
					stage=c(rep("seed.age.1",33),rep("seed.old",30),rep(NA,50)),
					stageNext=c("continuous"),number=1))
	
	return(dataf)
}

## Generate a simple data-frame for continuous and discrete covariates
#  with a total of 1000 measurements in columns called
# size, sizeNext, surv, fec, stage, stageNext number
#
# 
generateDataStoch <- function(){
	covariate1 <- rnorm(1000)
	covariate2 <- rnorm(1000)
	covariate3 <- rnorm(1000)
	size <- rnorm(1000,5,2)
	sizeNext <- 1+0.9*size+3*covariate1+0.01*covariate2+0.2*covariate3+rnorm(1000,0,0.1)
	
	fec <- surv <- rep(NA, length(size))
	surv[!is.na(size)] <- rbinom(sum(!is.na(size)),1,logit(-1+0.2*size[!is.na(size)]))
	fec[!is.na(size)] <- rnorm(sum(!is.na(size)),exp(-7+0.9*size[!is.na(size)]),1)
	fec[size<quantile(size,0.20,na.rm=TRUE) | fec<0] <- 0
	fec <- 10*fec
	
	seedlings <- sample(1:1000,size=100,replace=TRUE)
	size[seedlings] <- NA; 
	sizeNext[seedlings] <- rnorm(100,-2,0.1)
	surv[seedlings] <- 1
	#set to flower when covariate1 is around 1.5
	pfec <- 1*(runif(length(size))<logit(size+covariate1)); #print(pfec)
	fec[pfec==0] <- 0
	#fill in stage
	stage <- stageNext <- rep("continuous",1000)
	stage[is.na(size)] <- NA
	stageNext[is.na(sizeNext)] <- "dead"
	
	dataf <- data.frame(size=size,sizeNext=sizeNext,surv=surv,
			covariate1=covariate1,covariate2=covariate2,covariate3=covariate3,
			fec=fec, stage=stage,stageNext=stageNext, number=rep(1,length(size)))
	return(dataf)
}


## FUNCTION FOR MAKING A LIST OF IPMS ############################################
# to do for stoch env with a single discrete covariate. ##########################

makeListIPMs <- function(dataf, nBigMatrix=10, minSize=-2,maxSize=10, 
		integrate.type="midpoint", correction="none",
		explSurv="size+size2+covariate",explGrow="size+size2+covariate", 
		regType="constantVar",responseType="sizeNext",explFec="size",fec.constants=1) {
	
	#convert to 1:n for indexing later
	dataf$covariate <- as.factor(dataf$covariate)
	levels(dataf$covariate) <- 1:length(unique(dataf$covariate))
	
	sv1 <- makeSurvObj(dataf,
			explanatoryVariables=explSurv)
	gr1 <- makeGrowthObj(dataf,
			explanatoryVariables=explGrow,
			responseType=responseType,
			regType=regType)
	fv1 <- makeFecObj(dataf,explanatoryVariables=explFec, fec.constants=fec.constants) 
	
	covs <- unique(dataf$covariate)
	covs <- covs[!is.na(covs)]
	
	#print(covs)
	
	IPM.list <- list()
	for (k in 1:length(covs)) { 
		
		tpF <- create.IPM.Fmatrix(nBigMatrix = nBigMatrix, minSize = minSize,
				maxSize = maxSize, chosen.cov = k,
				fecObj = fv1,integrate.type=integrate.type, correction=correction)
		tpS <- create.IPM.Tmatrix(nBigMatrix = nBigMatrix, minSize = minSize,
				maxSize = maxSize, chosen.cov = k,growObj = gr1, survObj = sv1,
				integrate.type=integrate.type, correction=correction)
		IPM.list[[k]] <- tpF+tpS
	}
	return(IPM.list)
	
}




## Convert Increment - where exact dates of census vary but some multiplier of yearly increments
## are desired; this function takes a data-frame (with columns size, sizeNext,
## and, importantly exactDate, exactDatel)
## and returns a data-frame with sizeNext modified proportional to the time
## elapsed in the desired yaerly increments, adding an additional column denoted 'increment'. 
#
# Parameters - dataf - a dataframe with headings size, sizeNext, exactDate,exactDatel
#            - nyrs - the number of years between censuses desired (e.g. for CTFS data, 5 years)
#
# Returns - a dataframe with the same headings
convertIncrement <- function(dataf, nyrs=1) {
	incr <- dataf$sizeNext - dataf$size
	if (class(dataf$exactDatel) == "Date")  {
		timeElapsed <- (difftime(dataf$exactDatel,dataf$exactDate)[[1]])/(365*nyrs)
	} else {
		timeElapsed <- (dataf$exactDatel-dataf$exactDate)/(365*nyrs)
	}
	incrNew <- incr / timeElapsed
	dataf$sizeNext <- dataf$size + incrNew
	dataf$incr <- dataf$incrNew
	return(dataf)
	
}

## Function to run all analyses with simplest model
## fits for a dataf object
#
# Parameters - dataf - dataframe with right headings, i.e. size, sizeNext, surv
#            - chosen.size - size for which passage time desired
#            - minSize - lower limit for IPM - defaults to fraction of smallest observed 
#            - maxSize - upper limit for IPM - default to produce of largest observed
#            - nBigMatrix - numerical resolution of IPM - defaults to 500
#            - do.log - is data on a log scale? (for plotting) - default is TRUE
#            - do.plot - figures desired? - default is TRUE
#
# Returns - list including LE - Life Expectancy
#                          ptime - passage time to chosen size
#                          Tmatrix - Tmatrix
#                          growObj -
#                          survObj - survival object
#
run.Simple.Model <- function(dataf,
		chosen.size,
		minSize=c(),
		maxSize=c(),
		nBigMatrix=500,
		do.plot=TRUE,
		is.log=TRUE,
		integrate.type="midpoint", correction="none") {
	
	# set up IPM limits if not defined
	if (length(minSize)==0) {
		if (min(dataf$size,na.rm=TRUE)<0) minSize <- 2*min(dataf$size,na.rm=TRUE) else minSize <- 0.5*min(dataf$size,na.rm=TRUE)}
	if (length(maxSize)==0) maxSize <- 1.5*max(dataf$size,na.rm=TRUE)
	
	# Make necessary objects - here using simplest polynomials
	sv1 <- makeSurvObj(dataf)
	gr1 <- makeGrowthObj(dataf)
	
	#print(gr1)
	
	# Establish function for plotting 
	if (is.log) {
		conv <- function(x) return(exp(x)) 
		xrange <- range(exp(dataf$size),na.rm=TRUE)
		xrange[1] <- max(xrange[1],exp(minSize))
		xrange[2] <- min(xrange[2],exp(maxSize))
		axes <- "x"
	} else {
		conv <- function(x) return(x)
		xrange <- range(dataf$size,na.rm=TRUE)
		xrange[1] <- max(xrange[1],minSize)
		xrange[2] <- min(xrange[2],maxSize)
		axes <- ""
	}
	
	# Plot data and fitted models corresponding to these objects, if do.plot is TRUE
	if (do.plot) { 
		par(mfrow=c(3,2),bty="l")    
		#growth
		p1 <- picGrow(dataf,gr1)
		#survival
		p2 <- picSurv(dataf,sv1,ncuts=50)
	}
	
	# Make IPM Tmatrix with these objects, and chosen size range, and resolution (nBigMatrix)
	tmp <- create.IPM.Tmatrix(nBigMatrix = nBigMatrix, minSize = minSize, maxSize = maxSize,
			growObj = gr1, survObj = sv1,integrate.type=integrate.type, correction=correction)
	
	# Get the mean life expect from every size value in IPM
	LE <- MeanLifeExpect(tmp)
	varLE <- VarLifeExpect(tmp); #print(varLE)
	if (do.plot) { 
		plot(conv(tmp@meshpoints), LE,type = "l",xlab = "Size", ylab = "Mean life expectancy",log=axes,
				xlim=xrange, main="Life expectancy", ylim=range(c(pmax((LE-(sqrt(varLE))),0), (LE+(sqrt(varLE)))))) 
		points(conv(tmp@meshpoints), pmax((LE-(sqrt(varLE))),0),type="l",lty=3)
		points(conv(tmp@meshpoints), (LE+(sqrt(varLE))),type="l",lty=3)
	}
	
	# Get the passage time to the targeted size class
	ptime <- PassageTime(chosen.size,tmp); #print(ptime)
	if (do.plot) { 
		plot(conv(tmp@meshpoints),ptime,type = "l",xlab = "Size at start",log=axes,
				ylab = "Time to reach chosen size",ylim=range(ptime[tmp@meshpoints<chosen.size],na.rm=TRUE),
				xlim=conv(range(tmp@meshpoints[tmp@meshpoints<chosen.size],na.rm=TRUE)), main="Time to reach particular size")
		abline(v = conv(chosen.size),col = 2) #show the target size in red
	}
	
	
	
	
	return(list(ptime=ptime,LE=LE,Tmatrix=tmp,growObj=gr1, survObj=sv1))
}




# Function to plot the results of a stochastic simulation
# structure run
#
# Parameters - tvals - time points
#            - st - output of TrackPopStructManyCovSeedBank or TrackPopStructManyCov
#            - covtest - the key covariate for germination / flowering
#            - n.runin - how many to leave off pics
# Returns - 
#
plotResultsStochStruct <- function(tvals,st,covtest,n.runin=15,log="y",...) { 
		
	par(mfrow=c(2,2),bty="l")
	plot(tvals[n.runin:length(tvals)],
			colSums(st$rc[,n.runin:length(tvals)]+1),
				xlab="Time", 
			ylab="Population size",type="l",log=log,...)
	abline(v=1:max(tvals),lty=3)
	covtestplot <- exp(mean(colSums(st$rc[,n.runin:length(tvals)])) +
					((covtest-mean(covtest))/sd(covtest))*
					sd(colSums(st$rc[,n.runin:length(tvals)])))
	points(tvals,covtestplot+1,type="l",lty=3,col=2)

	if (log=="y") st$rc <- log(st$rc)
	
	image(tvals[n.runin:length(tvals)],
			st$IPM.here@meshpoints,
			t(st$rc[,n.runin:length(tvals)]+1),
			ylab="Size", xlab="Time")
}


# Function to coerce a general matrix pop model
# into my IPM form
#
# Parameters - a matrix
#
# Returns an object of class "IPM.matrix"
#
coerceMatrixIPM <- function(amat) {
	
	newmat <- as(amat,"IPM.matrix")
	newmat@n.env.class <- 1
	newmat@nBigMatrix <- length(newmat[1,])
	newmat@meshpoints <- 1:newmat@nBigMatrix
	newmat@env.index <- 1
	
	return(newmat)
	
}


# Function to build a discrete Tmatrix, with the same slots as an IPM.matrix
# provided with bins and the usual type of data-frame (columns size, sizeNext, surv)
#
# Parameters - dataf - a dataframe
#            - bins - the lower and upper edge of desired bins
#            - n.env - the environment level (currently just defaults)
#
# Returns - an object of class IPM.matrix with dim length(bins)*length(bins) containing
#         - survival transitions
#
create.MPM.Tmatrix <- function(dataf, bins, n.env=1) {
	
	loc.now <- findInterval(dataf$size[!is.na(dataf$size) & !is.na(dataf$sizeNext)],bins)+1
	loc.next <- findInterval(dataf$sizeNext[!is.na(dataf$size) & !is.na(dataf$sizeNext)],bins)+1    
	
	nbins <- max(c(loc.next,loc.now))
	
	MPM <- matrix(0,nbins,nbins)
	MPM[,] <- table(loc.next,loc.now)
	MPM <- t(t(MPM)/as.numeric(table(loc.now)))
	
	rc <- new("IPM.matrix",
			n.env.class = 1, 
			nBigMatrix = nbins,
			nrow = 1*nbins,
			ncol =1*nbins,
			meshpoints = 1:nbins,
			env.index = rep(1:n.env, each=nbins))
	
	rc[,] <- MPM  
	
	return(rc)
}

# Function to build a usual discrete Fmatrix, with the same slots as an IPM.matrix
# provided with bins and the usual type of data-frame
#
# Parameters - dataf - a dataframe
#            - bins - the lower and upper edge of desired bins
#            - p.est - probability of seed establishment
#            - n.env - the environment level (currently just defaults)
#
# Returns - an object of class IPM.matrix with dim length(bins)*length(bins) containing
#         - fecundity transitions
#
# ! assumes no relationship between adult size class and their baby's size class
#

create.MPM.Fmatrix <- function(dataf, bins,offspring.classes=1, offspring.prop=1, n.env=1) {
	
	loc.now <- findInterval(dataf$size[dataf$fec>0 & !is.na(dataf$size) & !is.na(dataf$fec)],bins)+1
	n.now <- sapply(split(dataf$fec[dataf$fec>0 & !is.na(dataf$size) & !is.na(dataf$fec)],loc.now),median)
	
	nbins <- max(loc.now); 
	#print(nbins)
	
	offspring.prop <- offspring.prop/sum(offspring.prop)
	
	MPM <- matrix(0,nbins,nbins)
	for (j in 1:length(offspring.classes)) MPM[offspring.classes[j],as.numeric(names(n.now))] <-  offspring.prop[j]*n.now
	
	rc <- new("IPM.matrix",
			n.env.class = 1, 
			nBigMatrix = nbins,
			nrow = 1*nbins,
			ncol =1*nbins,
			meshpoints = 1:nbins,
			env.index = rep(1:n.env, each=nbins))
	
	rc[,] <- MPM  
	
	return(rc)
}

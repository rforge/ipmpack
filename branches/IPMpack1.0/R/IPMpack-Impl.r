
## FUNCTIONS FOR TURNING DATA INTO GROWTH AND SURVIVAL OBJECTS ############################
## uses linear and logistic regressions, with various polynomials
# dataf must have columns "size", "sizeNext", "surv",
# and facultatively, "covariate" and "covariatel" (if single discrete covariate, like light)
# or covariate1, covariate2, covariate3, covariate4,... (if continuous and discrete)
# and "fec", "age";  age is for picking out seedling sizes for full IPM



## 1. Growth  models  #############################################


#More general function for growth than older versions listed below
#
# Parameters - dataf - the data-frame (which must contain the headings found in the formulas)
#              formula= - a model formula that requires
#                    "sizeNext" or "incr" as a reponse variable
#                    "size" as a possible covariate possibly with
#                        other combinations including size2 (size^2), size3(size^3), logsize(log(size)
#                        and potentially a discrete covariate (called covariate)
#              responseType - the response variable desired, crucial for building
#                             the right kind of object. Possible levels are "sizeNext", "incr", "logincr"
#              regType - options are "constantVar" (i.e. use lm), declineVar (use gls with model for variance)
# 
# Returns - a growth object                  
#
#
makeGrowthObj <- function(dataf,
		explanatoryVariables="size",
		responseType="sizeNext",
		regType="constantVar",
		Family="gaussian") {
	
	if (responseType=="incr" & length(dataf$incr) == 0) {
		print("building incr as sizeNext - size")
		dataf$incr <- dataf$sizeNext - dataf$size
	}
	if (responseType=="logincr" & length(dataf$logincr) == 0) {
		print("building logincr as log(sizeNext - size) - pre-build if this is not appropriate")
		dataf$logincr <- log(dataf$sizeNext - dataf$size)
	}
	
	Formula <- as.formula(paste(responseType, '~', explanatoryVariables, sep = ''))
	#create appropriate size based covariates
	dataf$size2 <- dataf$size ^ 2
	dataf$size3 <- dataf$size ^ 3
	if (length(grep("logsize", Formula)) > 0) dataf$logsize <- log(dataf$size)
	
	#setup for discrete covariates if data suggests may be implemented by the
	#presence of "covariate" and "covariateNext"
		if ("covariate" %in% strsplit(as.character(explanatoryVariables), "[+-\\*]")[[1]] & length(dataf$covariate) > 0) { 
		dataf$covariate <- as.factor(dataf$covariate)
		levels(dataf$covariate) <- 1:length(unique(dataf$covariate))
	}
	if ("covariateNext" %in% strsplit(as.character(explanatoryVariables), "[+-\\*]")[[1]] & length(dataf$covariateNext) > 0) { 
		dataf$covariateNext <- as.factor(dataf$covariateNext)
		levels(dataf$covariateNext) <- 1:length(unique(dataf$covariateNext))
	}
	#eval fit and make the objects
	if (Family == "poisson") {
		fit <- glm(Formula, data=dataf, family = "poisson")
		gr1 <- new("growthObjPois")
		gr1@fit <- fit
		if (regType != "constantVar") print("Warning: your regType is ignored because a poisson model is fitted")
	} else {
	if (regType == "constantVar")  {
		fit <- lm(Formula, data=dataf)
	} else { 
		if (regType == "declineVar"){
			require(nlme)
			Formula <- as.formula(Formula)	
			fit.here <- gls(Formula, na.action = na.omit, weights = varExp(form =  ~fitted(.)), data = dataf)
			fit <- list(coefficients = fit.here$coefficients,
						sigmax2 = summary(fit.here)$sigma^2,
						var.exp.coef = as.numeric(fit.here$modelStruct$varStruct[1]), 
						fit = fit.here)
		}
	}
    #make the objects
	#with sizeNext as response
	if (responseType == "sizeNext") { 
		
		if (class(fit) == "lm") { 
			gr1 <- new("growthObj")
			gr1@fit <- fit
		} else {
			if (class(fit.here) == "gls") { 
				gr1 <- new("growthObjDeclineVar")
				gr1@fit <- fit
			} else {
				print("unknown formula;
								please use lm or gls for declining variance models")
			}
		}
	} else {
		if (responseType == "incr") { 
			
			if (class(fit) == "lm") { 
				gr1 <- new("growthObjIncr")
				gr1@fit <- fit
			} else {
				if (class(fit.here) == "gls") { 
					gr1 <- new("growthObjIncrDeclineVar")
					gr1@fit <- fit
				} else {
					print("unknown formula;
									please use lm or gls for declining variance models")
				}
			}
		} else {
			if (responseType == "logincr") {
				
				if (class(fit) == "lm") { 
					gr1 <- new("growthObjLogIncr")
					gr1@fit <- fit
				} else {
					if (class(fit.here) == "gls") { 
						gr1 <- new("growthObjLogIncrDeclineVar")
						gr1@fit <- fit
					} else {
						print("unknown formula;
										please use lm or gls for declining variance models")
					}
				}
			}
		}
		}
	}
	return(gr1)
}




#More general function for growth than older versions listed below but now for multiple covariates
#
# Parameters - dataf - the data-frame (which must contain the headings found in the formulas)
#              formula= - a model formula that requires
#                    "sizeNext" or "incr" as a reponse variable
#                    "size" as a possible covariate possibly with
#                        other combinations including size2 (size^2), size3(size^3), logsize(log(size)
#                        and potentially a discrete covariate (called covariate)
#              responseType - the response variable desired, crucial for building
#                             the right kind of object. Possible levels are "sizeNext", "incr", "logincr"
#              regType - options are "constantVar" (i.e. use lm), declineVar (use gls with model for variance)
# 
# Returns - a growth object                  
#
#
makeGrowthObjManyCov <- function(dataf,
		explanatoryVariables = "size+size2+covariate1",
		responseType = "sizeNext",
		regType = "constantVar"){
	
	if (responseType == "incr" & length(dataf$incr) == 0) {
		print("building incr as sizeNext-size")
		dataf$incr <- dataf$sizeNext-dataf$size
	}
	
	if (responseType == "logincr" & length(dataf$logincr) == 0) {
		print("building logincr as log(sizeNext-size) - pre-build if this is not appropriate")
		dataf$logincr <- log(dataf$sizeNext-dataf$size)
	}
	
	Formula <- as.formula(paste(responseType, '~', explanatoryVariables, sep = ''))
	
	#create appropriate size based covariates
	dataf$size2 <- dataf$size ^ 2
	dataf$size3 <- dataf$size ^ 3
	if (length(grep("logsize", Formula)) > 0) dataf$logsize <- log(dataf$size)
	
	#eval fit
	if (regType == "constantVar")  {
		fit <- lm(Formula, data = dataf)
	} else { 
		if (regType == "declineVar"){
			require(nlme)
			Formula <- as.formula(Formula)	
			fit.here <- gls(Formula,	na.action = na.omit, weights = varExp(form =  ~fitted(.)), data = dataf)
			fit <- list(coefficients=fit.here$coefficients,
						sigmax2=summary(fit.here)$sigma^2,
						var.exp.coef=as.numeric(fit.here$modelStruct$varStruct[1]),
						fit=fit.here)
			
		}
	}
	
	#make the objects
	#with sizeNext as response
	if (responseType == "sizeNext") { 
		
		if (class(fit) == "lm") { 
			gr1 <- new("growthObjMultiCov")
			gr1@fit <- fit
		} else {
			if (class(fit.here) == "gls") { 
				gr1 <- new("growthObjMultiCovDeclineVar")
				gr1@fit <- fit
			} else {
				print("unknown formula;
								please use lm or gls for declining variance models")
			}
		}
	} else {
		if (responseType == "incr") { 
			
			if (class(fit) == "lm") { 
				gr1 <- new("growthObjMultiCovIncr")
				gr1@fit <- fit
			} else {
				if (class(fit.here) == "gls") { 
					gr1 <- new("growthObjMultiCovIncrDeclineVar")
					gr1@fit <- fit
				} else {
					print("unknown formula;
									please use lm or gls for declining variance models")
				}
			}
			
		} else {
			if (responseType == "logincr") {
				
				if (class(fit) == "lm") { 
					gr1 <- new("growthObjMultiCovLogIncr")
					gr1@fit <- fit
				} else {
					if (class(fit.here) == "gls") { 
						gr1 <- new("growthObjMultiCovLogIncrDeclineVar")
						gr1@fit <- fit
					} else {
						print("unknown formula;
										please use lm or gls for declining variance models")
					}
				}
			}
		}
	}
	return(gr1)
}




## Function to create a new Hossfeld growth object
#
# Parameters - dataf - a dataframe
#
# Returns - a Hossfeld growth object

#no covariate, and one polynom, linear regression
makegrowthObjHossfeld <- function(dataf) {  
	if (length(dataf$incr)==0) dataf$incr <- dataf$sizeNext-dataf$size
	dataf$incr[dataf$incr<0] <- 0
	tmp <- optim(c(1, 1, 1), wrapHossfeld, dataf = dataf, method = "Nelder-Mead")
	print(tmp$convergence)
	gr1 <- new("growthObjHossfeld")
	gr1@paras <- tmp$par
	resids <- Hossfeld(dataf$size, tmp$par) - dataf$incr 
	gr1@sd <- sd(resids, na.rm = T)
	return(gr1)
}


# Function to create a truncated increment growth model 
# currently only defined for increment, leftVal is min
#
# Paramteres - dataf - a dataframe
#
#
# Retrurns - growth object with truncated increment
#
#no covariate, and one polynom, linear regression on increment
#makeGrowthObjIncrTrunc <- function(dataf,
#		explanatoryVariables = "size",
#		responseType = "incr",
#		leftVal = 0) {
#	require(censReg)
#	
#	dataf$size2 <- dataf$size ^ 2
#	
#	Formula <- paste(responseType, '~', explanatoryVariables, sep = '')
#
#	if (length(grep("logsize", Formula)) > 0) dataf$logsize <- log(dataf$size)
#			
#	if (responseType == "incr") { 
#		if (length(dataf$incr) == 0) {
#			print("building incr as sizeNext-size")
#			dataf$incr <- dataf$sizeNext - dataf$size
#			dataf$incr[dataf$incr<leftVal] <- NA
#		}}
#				
#	fit <- censReg(as.formula(Formula),data=dataf, left=leftVal)
#	gr1 <- new("growthObjTruncIncr")
#	gr1@fit <- fit$estimate
#	gr1@varcov <- vcov(fit)
#	return(gr1)
#}




# 2. Survival models  #######################################################################################################




#General function for survival 
#
# Parameters - dataf - the data-frame (which must contain the headings found in the formulas)
#              formula - a model formula that requires
#                    "surv" as a reponse variable
#                    "size" as a possible covariate possibly with
#                        other combinations including size2 (size^2), size3(size^3), logsize(log(size)
#                        and potentially a discrete covariate (called covariate)
# 
# Returns - a survival object                   
#
#
makeSurvObj <- function(dataf,
		explanatoryVariables="size+size2"){
	
	
	#build appropriate size based covariates
	dataf$size2 <- dataf$size^2
	dataf$size3 <- dataf$size^3
	if (length(grep("logsize",explanatoryVariables))>0) dataf$logsize <- log(dataf$size)
	
	#build formula
	formula<-paste('surv','~',explanatoryVariables,sep='')
	
	#print(formula)
	
	#setup for discrete covariates if data suggests may be implemented by the
	#presence of "covariate" and "covariateNext"
	if ("covariate"%in%strsplit(as.character(explanatoryVariables),"[+-\\*]")[[1]]&length(dataf$covariate)>0) { 
		dataf$covariate <- as.factor(dataf$covariate)
		levels(dataf$covariate) <- 1:length(unique(dataf$covariate))
		
	}
	if ("covariateNext"%in%strsplit(as.character(explanatoryVariables),"[+-\\*]")[[1]]&length(dataf$covariateNext)>0) { 
		dataf$covariateNext <- as.factor(dataf$covariateNext)
		levels(dataf$covariateNext) <- 1:length(unique(dataf$covariateNext))
	}
	
	
	#  print(formula)
	
	fit <- glm(formula,family=binomial,data=dataf)
	sv1 <- new("survObj")
	sv1@fit <- fit
	return(sv1)
}


#More general function for survival than older versions (listed below) 
#
# Parameters - dataf - the data-frame (which must contain the headings found in the formulas)
#              formula - a model formula that requires
#                    "surv" as a reponse variable
#                    "size" as a possible covariate possibly with
#                        other combinations including size2 (size^2), size3(size^3), logsize(log(size)
#                        and potentially a discrete covariate (called covariate)
# 
# Returns - a survival object                   
#
#
makeSurvObjManyCov <- function(dataf,
		explanatoryVariables="size+size2+covariate1"){
	
	#build appropriate size based covariates
	dataf$size2 <- dataf$size^2
	dataf$size3 <- dataf$size^3
	if (length(grep("logsize",explanatoryVariables))>0) dataf$logsize <- log(dataf$size)
	
	#build formula
	formula<-paste('surv','~',explanatoryVariables,sep='')
	
	fit <- glm(formula,family=binomial,data=dataf)
	sv1 <- new("survObjMultiCov")
	sv1@fit <- fit
	return(sv1)
}


# 3. fecundity models  #######################################################################################################

makeFecObj <- function(dataf,
		fecConstants=as.numeric(NA),
		explanatoryVariables="size",
		Family="gaussian",
		Transform="none",
		fecNames=NA,
		meanOffspringSize=NA,
		varOffspringSize=NA,
		offspringSplitter=data.frame(continuous=1),
		fecByDiscrete=data.frame(NA)){
	
	#order stage names from discrete to continuous
	stages <- names(tapply(c(levels(dataf$stage),levels(dataf$stageNext)),c(levels(dataf$stage),levels(dataf$stageNext)),length))
	stages <- stages[stages!="dead"] 
	stages <- c(stages[stages!="continuous"],"continuous") 
	if ((sum(names(offspringSplitter)%in%stages)/length(offspringSplitter))<1) {
		stop("Error - the variable names in your offspringSplitter data.frame are not all part of the levels of stage or stageNext in your data file. Please fix this by adjusting your offspringSplitter entry to include the correct variable names, e.g. offspringSplitter=data.frame(continuous=.7,seedAge1=.3)")
	}
	dummy<-rep(0,length(stages));names(dummy)<-stages;dummy<-as.data.frame(t(as.matrix(dummy)))
	for (i in names(offspringSplitter)) dummy[i]<-offspringSplitter[i]
	offspringSplitter <- dummy
	
	##warnings
	if (length(dataf$stage)==0) {
		print("Warning - no column named stage - assuming all continuous")
		dataf$stageNext <- dataf$stage <- rep("continuous", length(dataf[,1]))
		dataf$stage[is.na(dataf$size)] <- NA
		dataf$stageNext[is.na(dataf$sizeNext)] <- "dead"
	}
	
	if (ncol(offspringSplitter)>1 & (ncol(offspringSplitter)-1)!=ncol(fecByDiscrete)) {
		print("Warning - offspring splitter indicates more than just continuous stages. No fecundity by the discrete stages supplied in fecByDiscrete; assumed that is 0")
		#fecByDiscrete <- matrix(0,col(offspringSplitter)-1,col(offspringSplitter)-1)
		fecByDiscrete <- offspringSplitter[,1:(ncol(offspringSplitter)-1)]
		fecByDiscrete[] <- 0
	}
	
	if (sum(offspringSplitter)!=1) {
		print("Warning - offspring splitter does not sum to 1. It is now rescaled to sum to 1.")
		offspringSplitter <- offspringSplitter / sum(offspringSplitter) 
	}
	
	if ("covariate"%in%strsplit(as.character(explanatoryVariables),"[+-\\*]")[[1]]&length(dataf$covariate)>0) { 
		dataf$covariate <- as.factor(dataf$covariate)
		levels(dataf$covariate) <- 1:length(unique(dataf$covariate))
		
	}
	if ("covariateNext"%in%strsplit(as.character(explanatoryVariables),"[+-\\*]")[[1]]&length(dataf$covariateNext)>0) { 
		dataf$covariateNext <- as.factor(dataf$covariateNext)
		levels(dataf$covariateNext) <- 1:length(unique(dataf$covariateNext))
	}
	#print(table(dataf$covariate))
	
	f1 <- new("fecObj")
	dataf$size2 <- dataf$size^2
	if (length(grep("logsize",as.character(explanatoryVariables)))>0) dataf$logsize <- log(dataf$size)
	
	if (is.na(fecNames)) fecNames <- names(dataf)[grep("fec",names(dataf))]
	if (length(fecNames)>length(explanatoryVariables)) {
		misE <- (length(explanatoryVariables)+1):length(fecNames)
		print(c("number in explanatoryVariables not the same as the number of fecundity columns in the data file, using default of `size' for missing ones which are:",fecNames[misE],". (which might be exactly what you want)"))
		explanatoryVariables <- c(explanatoryVariables,rep("size",length(fecNames)-length(explanatoryVariables)))
	}
	if (length(fecNames)>length(Family)) {
		misE <- (length(Family)+1):length(fecNames)
		print(c("number of families not the same as the number of fecundity columns in the data file, using default of `gaussian' for missing ones which are:",fecNames[misE],". (which might be exactly what you want)"))
		Family <- c(Family,rep("gaussian",length(fecNames)-length(Family)))
	}
	if (length(fecNames)>length(Transform)) {
		misE <- (length(Transform)+1):length(fecNames)
		print(c("number of transforms not the same as the number of fecundity columns in the data file, using default of `none' for missing ones which are:",fecNames[misE],". (which might be exactly what you want)"))
		Transform <- c(Transform,rep("none",length(fecNames)-length(Transform)))
	}
	
	for (i in 1:length(fecNames)) {
		if (Transform[i]=="log") dataf[,fecNames[i]] <- log(dataf[,fecNames[i]])
		if (Transform[i]=="sqrt") dataf[,fecNames[i]] <- sqrt(dataf[,fecNames[i]])
		if (Transform[i]=="-1") dataf[,fecNames[i]] <- dataf[,fecNames[i]]-1
		dataf[!is.finite(dataf[,fecNames[i]]),fecNames[i]] <- NA
		f1@fitFec[[i]] <- glm(paste(fecNames[i],'~',explanatoryVariables[i],sep=''),family=Family[i],data=dataf)
	}
	
	if (is.na(meanOffspringSize)) {
		offspringdata<-subset(dataf,is.na(dataf$stage)&dataf$stageNext=="continuous")
		meanOffspringSize <- mean(offspringdata$sizeNext)
		varOffspringSize <- var(offspringdata$sizeNext) }
	f1@fecConstants <- fecConstants
	f1@meanOffspringSize <- meanOffspringSize
	f1@varOffspringSize <- varOffspringSize
	f1@offspringSplitter <- offspringSplitter 
	f1@fecByDiscrete <- fecByDiscrete
	f1@Transform <- Transform
	return(f1)
}



## NO different from the above yet, except in what it produces

makeFecObjManyCov <- function(dataf,
		fecConstants=as.numeric(NA),
		explanatoryVariables="size",
		Family="gaussian",
		Transform="none",
		meanOffspringSize=NA,
		varOffspringSize=NA,
		offspringSplitter=data.frame(continuous=1),
		fecByDiscrete=data.frame(NA)){
	
	##warnings
	if (length(dataf$stage)==0) {
		print("Warning - no column named stage - assuming all continuous")
		dataf$stageNext <- dataf$stage <- rep("continuous", length(dataf[,1]))
		dataf$stage[is.na(dataf$size)] <- NA
		dataf$stageNext[is.na(dataf$sizeNext)] <- "dead"
	}
	
	if(ncol(offspringSplitter)>1 & (ncol(offspringSplitter)-1)!=ncol(fecByDiscrete)) {
		print("Warning - offspring splitter indicates more than just continuous stages. No fecundity by the discrete stages supplied in fecByDiscrete; assumed that is 0")
		fecByDiscrete <- matrix(0,col(offspringSplitter)-1,col(offspringSplitter)-1)
	}
	
	if(sum(offspringSplitter)!=1) {
		print("Warning - offspring splitter does not sum to 1. It is now rescaled to sum to 1.")
		
	}

	# the covariate covariatenext transform removed here, because presumably if you are doing ManyCov you 
	# have something else in mind!
	
	f1 <- new("fecObjMultiCov")
	dataf$size2 <- dataf$size^2
	if (length(grep("logsize",explanatoryVariables))>0) dataf$logsize <- log(dataf$size)
	
	fecNames <- names(dataf)[grep("fec",names(dataf))]
	if (length(fecNames)>length(explanatoryVariables)) {
		misE <- length(explanatoryVariables):length(fecNames)
		print(c("number in explanatoryVariables not the same as the number of fecundity columns in the data file, using default of `size' for missing ones which are:",fecNames[misE]))
		explanatoryVariables <- c(explanatoryVariables,rep("size",length(fecNames)-length(explanatoryVariables)))
	}
	if (length(fecNames)>length(Family)) {
		misE <- length(Family):length(fecNames)
		print(c("number of families not the same as the number of fecundity columns in the data file, using default of `gaussian' for missing ones which are:",fecNames[misE]))
		Family <- c(Family,rep("gaussian",length(fecNames)-length(Family)))
	}
	if (length(fecNames)>length(Transform)) {
		misE <- length(Transform):length(fecNames)
		print(c("number of transforms not the same as the number of fecundity columns in the data file, using default of `none' for missing ones which are:",fecNames[misE]))
		Transform <- c(Transform,rep("none",length(fecNames)-length(Transform)))
	}
	
	for (i in 1:length(fecNames)) {
		
		if (Transform[i]=="log") dataf[,fecNames[i]] <- log(dataf[,fecNames[i]])
		if (Transform[i]=="sqrt") dataf[,fecNames[i]] <- sqrt(dataf[,fecNames[i]])
		if (Transform[i]=="-1") dataf[,fecNames[i]] <- dataf[,fecNames[i]]-1
		dataf[!is.finite(dataf[,fecNames[i]]),fecNames[i]] <- NA
		
		#print(range(dataf[,fecNames[i]]))
		#print(range(dataf[,"size"], na.rm=TRUE))
		#print(range(dataf[,"covariate"]))
		
		f1@fitFec[[i]] <- glm(paste(fecNames[i],'~',explanatoryVariables[i],sep=''),family=Family[i],data=dataf)
	}
	
	if (is.na(meanOffspringSize)) {
		offspringdata<-subset(dataf,is.na(dataf$stage)&dataf$stageNext=="continuous")
		meanOffspringSize <- mean(offspringdata$sizeNext)
		varOffspringSize <- var(offspringdata$sizeNext) }
	f1@fecConstants <- fecConstants
	f1@meanOffspringSize <- meanOffspringSize
	f1@varOffspringSize <- varOffspringSize
	f1@offspringSplitter <- offspringSplitter 
	f1@fecByDiscrete <- fecByDiscrete
	f1@Transform <- Transform
	return(f1)
}





# 4. Discrete Transition models  #######################################################################################################

## Function to take a data-frame and make a discrete transition object
## for combining with a continuous T matrix
#
# Parameters - dataf - dataframe with headings of at least
#                      size, sizeNext, surv, fec, stage, stageNext, number
#
# Returns - an object of class discreteTrans
#
makeDiscreteTrans <- function(dataf) {
	
	#order stage names from discrete to continuous
	stages <- names(tapply(c(levels(dataf$stage),levels(dataf$stageNext)),c(levels(dataf$stage),levels(dataf$stageNext)),length))
	stages <- stages[stages!="dead"] 
	stages <- c(stages[stages!="continuous"],"continuous") 
	#define the number of classes
	nclasses <- length(stages)
	#define matrices to hold the transition between all classes
	discrete.trans <- matrix(0,nrow=nclasses,ncol=nclasses, dimnames=list(stages,stages))
	#define matrix to hold sd and mean of re-entry into continous + matrix of  survival for all discrete stages
	sd.to.cont <- mean.to.cont <- discrete.surv <- matrix(NA,nrow=1,ncol=nclasses-1,dimnames=list(1,stages[1:length(stages)-1]))
	# define matrix to hold transitions from the continuous to the discrete
	distrib.to.discrete <- matrix(NA,ncol=1,nrow=nclasses-1,dimnames=list(stages[1:length(stages)-1],"continuous"))
	
	#loop over discrete stages and fill 
	for (j in stages[1:(length(stages)-1)]) {
		for (i in stages) discrete.trans[i,j] <- sum(dataf[dataf$stage==j & dataf$stageNext==i,]$number,na.rm=TRUE)
		discrete.surv[,j] <- sum(discrete.trans[,j],na.rm=T) / sum(dataf[dataf$stage == j,]$number, na.rm = TRUE)
		discrete.trans[,j] <- discrete.trans[,j] / sum(discrete.trans[,j], na.rm = TRUE)
		mean.to.cont[,j] <- mean(dataf[dataf$stage == j & dataf$stageNext == i,]$sizeNext, na.rm = TRUE)
		sd.to.cont[,j] <- sd(dataf[dataf$stage == j & dataf$stageNext == i, ]$sizeNext,na.rm = TRUE)
	}
	
	for (i in stages[1:(length(stages) - 1)])
		distrib.to.discrete[i,] <- sum(dataf[dataf$stage == "continuous" & dataf$stageNext == i,]$number, na.rm = TRUE)
	distrib.to.discrete <- distrib.to.discrete/sum(distrib.to.discrete,na.rm=TRUE)
	
	
	subdata <- subset(dataf, dataf$stage == "continuous" & dataf$surv == 1)
	subdata$cont.to.discrete <- 1
	subdata$cont.to.discrete[subdata$stageNext == "continuous"] <- 0
	subdata$size2 <- subdata$size ^ 2
	surv.to.discrete <- glm(cont.to.discrete ~ size + size2, family = binomial, data = subdata)
			
	rownames(discrete.trans) <- stages	
	colnames(discrete.trans) <- stages	
	
	#define new object
	disTrans <- new("discreteTrans")
	disTrans@nclasses <- nclasses
	disTrans@discrete.trans <- discrete.trans
	disTrans@discrete.surv <- discrete.surv
	disTrans@mean.to.cont <- mean.to.cont
	disTrans@sd.to.cont <- sd.to.cont
	disTrans@distrib.to.discrete <- distrib.to.discrete
	disTrans@surv.to.discrete <- surv.to.discrete
	
	return(disTrans)
	
}


## 5. Models including Data Augmentation models #############################################################


# Function to augment mortality data for large trees
# Parameters dataf - existing data frame
#            size.thresh - the size above which tree death is being augmented
#            prop.dead - the proportion of these expected to be dead

.deathDataAugment <- function (dataf, size.thresh, prop.dead) { 
	
	n.now <- sum(dataf$size > size.thresh)
	n.new.dead <- ceiling(prop.dead * n.now / (1 - prop.dead))
	new.size <- rnorm(n.new.dead,size.thresh, sd(dataf$size[dataf$size > size.thresh]))
	
	datanew <- data.frame(size = new.size, sizeNext = rep(NA, n.new.dead), surv = rep(0, n.new.dead), 
			covariate = rep(0, n.new.dead), covariateNext = rep(0, n.new.dead),
			fec = rep(NA, n.new.dead)) 
	
	dataf.new <- rbind(dataf,datanew)
	
	return(dataf.new)
	
}






## 6. Models building Bayes posteriors and corresponding growth / survival and Matrices  ##############

# ! these all need better coding up of priors ###

## Using MCMCglmm to get a list of growth posteriors on size
# with either no covariate, or just 1 discrete covariate
#
# Parameters - dataf - dataframe
#              formula= - a model formula that requires
#                    "sizeNext" or "incr" as a reponse variable
#                    "size" as a possible covariate possibly with
#                        other combinations including size2 (size^2), size3(size^3), logsize(log(size)
#                        and potentially a discrete covariate (called covariate)
#              responseType - the response variable desired, crucial for building
#                             the right kind of object. Possible levels are "sizeNext", "incr", "logincr"
#            - meanB the mean of the prior of the coefficients for survival (should be the same length as desired coeff)
#            - varB the variance of the prior of the coeff for survival (note could add for growth also)
#            - nitt - the total number of iterations
#
# Returns - list including list of growth objects, + list of survival objects
makePostGrowthObjs <- function(dataf,
		explanatoryVariables = "size+size2+covariate",
		responseType = "sizeNext",
		meanB=rep(0,3), varB = rep(1e10), burnin = 3000, nitt = 50000) {
	
	require(MCMCglmm)
	
	if (responseType == "incr" & length(dataf$incr) == 0) {
		print("building incr as sizeNext-size")
		dataf$incr <- dataf$sizeNext - dataf$size
	}
	
	if (responseType == "logincr" & length(dataf$logincr) == 0) {
		print("building logincr as log(sizeNext-size) - pre-build if this is not appropriate")
		dataf$logincr <- log(dataf$sizeNext - dataf$size)
	}
	
	dataf$size2 <- dataf$size ^ 2
	dataf$size3 <- dataf$size ^ 3
	if (length(grep("logsize", explanatoryVariables)) > 0) dataf$logsize <- log(dataf$size)
	
	#setup for discrete covariates if data suggests may be implemented by the
	#presence of "covariate" and "covariateNext"
	if ("covariate" %in% strsplit(as.character(explanatoryVariables), "[+-\\*]")[[1]] & length(dataf$covariate) > 0) { 
		dataf$covariate <- as.factor(dataf$covariate)
		levels(dataf$covariate) <- 1:length(unique(dataf$covariate))
		
	}
	if ("covariateNext"%in%strsplit(as.character(explanatoryVariables), "[+-\\*]")[[1]]&length(dataf$covariateNext) > 0) { 
		dataf$covariateNext <- as.factor(dataf$covariateNext)
		levels(dataf$covariateNext) <- 1:length(unique(dataf$covariateNext))
	}
	
		
	#get rid of NAs
	dataf <- dataf[!is.na(dataf$size) & !is.na(dataf$sizeNext),]
	
	#fit growth model
	Formula <- as.formula(paste(responseType, '~', explanatoryVariables, sep = ''))
	fit <- MCMCglmm(Formula, data = dataf, verbose = FALSE, burnin = burnin, nitt = nitt)
	dummyFit <- lm(Formula, data = dataf)
	
	#create list of growth models reflecting posterior
	gr <- list()
	for (k in 1:length(fit$Sol[,1])) {
		dummyFit$coefficients <- fit$Sol[k,]
		dummyFit <- alteredFit(dummyFit = dummyFit, newCoef = fit$Sol[k,],  desiredSd = sqrt(fit$VCV[k, 1]))
		if (responseType=="sizeNext") gr[[k]] <-  new("growthObj")
		if (responseType=="incr") gr[[k]] <-  new("growthObjIncr")
		if (responseType=="logincr") gr[[k]] <-  new("growthObjLogIncr")
		gr[[k]]@fit <- dummyFit       
	}
	
	return(gr)
	
}

# replace the growth object fit with a new, desired variance for predict
alteredFit <- function(dummyFit = dummyFit, 
		newCoef = dummyFit$coefficients, 
		desiredSd = 1) {
	dummyFit$coefficients[] <- newCoef
	dummyFit$residuals <- rnorm(length(dummyFit$residuals), mean = 0, sd = desiredSd)	
	# need to use qr here to assign dummyFit so that there is no warning when decomposed for n
	# Error in rnorm(residDf, 0, sd = desiredSd) : object 'residDf' not found
	return(dummyFit)	
}


## Using MCMCglmm to get a list of  survival posterios
# with either no covariate, or just 1 discrete covariate
#
# Parameters - dataf - dataframe
#              formula - a model formula that requires
#                    "surv" as a reponse variable
#                    "size" as a possible covariate possibly with
#                        other combinations including size2 (size^2), size3(size^3), logsize(log(size)
#                        and potentially a discrete covariate (called covariate)
# #            - meanB the mean of the prior of the coefficients for survival (should be the same length as desired coeff)
#            - varB the variance of the prior of the coeff for survival (note could add for growth also)
#            - nitt - the total number of iterations
#
# Returns - list including list of growth objects, + list of survival objects
makePostSurvivalObjs <- function(dataf,
		explanatoryVariables="size+size2",
		meanB = rep(0, 3), varB=rep(1e10),burnin=3000, nitt = 50000) {
	
	require(MCMCglmm)
	#build appropriate size based covariates
	dataf$size2 <- dataf$size ^ 2
	dataf$size3 <- dataf$size ^ 3
	if (length(grep("logsize",explanatoryVariables)) > 0) dataf$logsize <- log(dataf$size)
	
	#setup for discrete covariates if data suggests may be implemented by the
	#presence of "covariate" and "covariateNext"
	if ("covariate"%in%strsplit(as.character(explanatoryVariables),"[+-\\*]")[[1]]&length(dataf$covariate)>0) { 
		dataf$covariate <- as.factor(dataf$covariate)
		levels(dataf$covariate) <- 1:length(unique(dataf$covariate))
		
	}
	if ("covariateNext"%in%strsplit(as.character(explanatoryVariables),"[+-\\*]")[[1]]&length(dataf$covariateNext)>0) { 
		dataf$covariateNext <- as.factor(dataf$covariateNext)
		levels(dataf$covariateNext) <- 1:length(unique(dataf$covariateNext))
	}
			
	#get rid of NAs
	dataf <- dataf[!is.na(dataf$size)& !is.na(dataf$surv),]
	
	
	#build formula
	Formula<-paste('surv','~',explanatoryVariables,sep='')
	Formula <- as.formula(Formula)
	
	#fit survival model 
	#avoid fitting a prior unless you need to (takes longer with)
    fit<-MCMCglmm(Formula, data=dataf[!is.na(dataf$surv),], 
			verbose=FALSE, prior = list(R = list(V = 1, fix = 1)),
			family="categorical", burnin=burnin,nitt = nitt)        
	
	dummy.fit <- glm(Formula, data=dataf,family=binomial)
	
	#create list survival models reflecting posterior
	sv <- list()
	for (k in 1:length(fit$Sol[,1])) {
		dummy.fit$coefficients <- fit$Sol[k,]
		sv[[k]] <-  new("survObjOverDisp")
		sv[[k]]@fit <- dummy.fit       
	}
	
	return(sv)
	
}




## Using MCMCglmm to get a list of  fecundity posteriors - this assumes that
# only applies to the first three functions
#
# Parameters - dataf - dataframe
#            - covfec - is there a discrete covariate in prob of reproduction
#            - covseeds - is there a discrete covariate in fecundity production
#            - meanB the mean of the prior of the coefficients for XX
#            - varB the variance of the prior of the coeff for XX
#            - nitt - the total number of iterations
#
# Returns - list including list of growth objects, + list of survival objects
makePostFecObjs <- function(dataf,
		fecConstants=as.numeric(NA),
		explanatoryVariables="size",
		Family="gaussian",
		Transform="none",
		meanOffspringSize=NA,
		varOffspringSize=NA,
		offspringSplitter=data.frame(continuous=1),
		fecByDiscrete=data.frame(NA),burnin=3000,nitt=50000) {
	
	require(MCMCglmm)
	
	##warnings
	if (length(dataf$stage)==0) {
		print("Warning - no column named stage - assuming all continuous")
		dataf$stageNext <- dataf$stage <- rep("continuous", length(dataf[,1]))
		dataf$stage[is.na(dataf$size)] <- NA
		dataf$stageNext[is.na(dataf$sizeNext)] <- "dead"
	}
	
	if(ncol(offspringSplitter)>1 & (ncol(offspringSplitter)-1)!=ncol(fecByDiscrete)) {
		print("Warning - offspring splitter indicates more than just continuous stages. No fecundity by the discrete stages supplied in fecByDiscrete; assumed that is 0")
		fecByDiscrete <- matrix(0,col(offspringSplitter)-1,col(offspringSplitter)-1)
	}
	
    #setup for discrete covariates if data suggests may be implemented by the
	#presence of "covariate" and "covariateNext"
	if ("covariate"%in%strsplit(as.character(explanatoryVariables),"[+-\\*]")[[1]]&length(dataf$covariate)>0) { 
		dataf$covariate <- as.factor(dataf$covariate)
		levels(dataf$covariate) <- 1:length(unique(dataf$covariate))
		
	}
	if ("covariateNext"%in%strsplit(as.character(explanatoryVariables),"[+-\\*]")[[1]]&length(dataf$covariateNext)>0) { 
		dataf$covariateNext <- as.factor(dataf$covariateNext)
		levels(dataf$covariateNext) <- 1:length(unique(dataf$covariateNext))
	}
			
	dataf$size2 <- dataf$size^2
	if (length(grep("logsize",explanatoryVariables))>0) dataf$logsize <- log(dataf$size)
	fecNames <- names(dataf)[grep("fec",names(dataf))]
	if (length(fecNames)>length(explanatoryVariables)) {
		print("number in explanatoryVariables not the same as the number of fecundity columns in the data file, taking first for all")
		explanatoryVariables <- rep(explanatoryVariables[1],length(fecNames))
	}
	if (length(fecNames)>length(Family)) {
		print("number of families not the same as the number of fecundity columns in the data file, taking first for all")
		Family <- rep(Family[1],length(fecNames))
	}
	if (length(fecNames)>length(Transform)) {
		print("number of transforms not the same as the number of fecundity columns in the data file, taking first for all")
		Transform <- rep(Transform[1],length(fecNames))
	}
	
	if (is.na(meanOffspringSize)) {
		offspringdata<-subset(dataf,is.na(dataf$stage)&dataf$stageNext=="continuous")
		meanOffspringSize <- mean(offspringdata$sizeNext)
		varOffspringSize <- var(offspringdata$sizeNext)
	}
	
	#get rid of NAs
	dataf <- dataf[!is.na(dataf$size) & !is.na(dataf$sizeNext),]
	
	fit <- dummy.fit <- list()
	
	for (i in 1:length(fecNames)) {
		
		if (Transform[i]=="log") dataf[,fecNames[i]] <- log(dataf[,fecNames[i]])
		if (Transform[i]=="sqrt") dataf[,fecNames[i]] <- sqrt(dataf[,fecNames[i]])
		if (Transform[i]=="-1") dataf[,fecNames[i]] <- dataf[,fecNames[i]]-1
		dataf[!is.finite(dataf[,fecNames[i]]),fecNames[i]] <- NA
		
		Formula <- paste(fecNames[i],'~',explanatoryVariables[i],sep='')
		Formula <- as.formula(Formula)
		
		fit[[i]] <- MCMCglmm(Formula,
				data=dataf[!is.na(dataf[,fecNames[i]]),], 
				verbose=FALSE, burnin=burnin,nitt=nitt,family=Family[i])
		
		dummy.fit[[i]] <- glm(Formula,
				data=dataf[!is.na(dataf[,fecNames[i]]),],family=Family[i])
		
	}
	
	#print(length(fit[[1]]$Sol[,1]))	
	
	#create list of growth models reflecing posterior
	fv <- list()
	for (k in 1:length(fit[[1]]$Sol[,1])) {
		fv[[k]] <-  new("fecObj")
		
		for (i in 1:length(fecNames)) { 
			dummy.fit[[i]]$coefficients <- fit[[i]]$Sol[k,]
			##TODO check over-ride
			dummy.fit[[i]]$residuals <- rnorm(length(dummy.fit[[i]]$residuals),0,sqrt(fit[[i]]$VCV[k,1]))
			fv[[k]]@fitFec[[i]] <- dummy.fit[[i]]
					
		}
		
		fv[[k]]@fecConstants <- fecConstants
		fv[[k]]@meanOffspringSize <- meanOffspringSize
		fv[[k]]@varOffspringSize <- varOffspringSize
		fv[[k]]@offspringSplitter <- offspringSplitter
		fv[[k]]@fecByDiscrete <- fecByDiscrete
		fv[[k]]@Transform <- Transform 
	}
	print(k)
	
	return(fv)
	
}


# Function to take a list of growth and survival objects and make a list of Tmatrices
#
# Parameters - growObjList - a list of growth objects
#            - survObjList - a list of survival objects
#            - nBigMatrix - the number of bins
#            - minSize - the minimum size
#            - maxSize - the maximum size
#            - cov - is a discrete covariate considered
#            - envMat - enviromental matrix for transition between
# 
# Returns    - a list of Tmatrices
makeListTmatrix <- function(growObjList,survObjList,
		nBigMatrix,minSize,maxSize, cov=FALSE, envMat=NULL,
		integrateType="midpoint",correction="none") {
	
	if (length(growObjList)>length(survObjList)) { 
		survObjList <- sample(survObjList,size=length(growObjList),replace=TRUE)
	} else { 
		if (length(growObjList)<length(survObjList))  
			growObjList <- sample(growObjList,size=length(survObjList),replace=TRUE)
	}
	
	nsamp <- length(growObjList)
	
	TmatrixList <- list()
	for ( k in 1:length(growObjList)) { 
		if (!cov) {
			TmatrixList[[k]] <- createIPMTmatrix(nBigMatrix = nBigMatrix, minSize = minSize, 
					maxSize = maxSize, growObj = growObjList[[k]],
					survObj = survObjList[[k]],integrateType=integrateType, correction=correction) 
		} else {
			TmatrixList[[k]] <- createCompoundTmatrix(nEnvClass = length(envMat[1,]),
					nBigMatrix = nBigMatrix, minSize = minSize, 
					maxSize = maxSize, envMatrix=envMat,
					growObj = growObjList[[k]],
					survObj = survObjList[[k]],integrateType=integrateType, correction=correction)    
		}
	}
	
	return(TmatrixList)
}

# Function to take a list of growth and survival objects and make a list of Fmatrices

makeListFmatrix <- function(fecObjList,nBigMatrix,minSize,maxSize, cov=FALSE, 
		envMat=NULL,integrateType="midpoint",correction="none") {
	
	nsamp <- max(length(fecObjList))
	if (length(fecObjList)<nsamp)  
		fecObjList <- sample(fecObjList,size=nsamp,replace=TRUE)
	
	FmatrixList <- list()
	for ( k in 1:nsamp) {
		if (!cov) { 
			FmatrixList[[k]] <- createIPMFmatrix(nBigMatrix = nBigMatrix, minSize = minSize, 
					maxSize = maxSize, 
					fecObj=fecObjList[[k]],integrateType=integrateType, correction=correction)
		} else {
			FmatrixList[[k]] <- createCompoundFmatrix(nEnvClass = length(envMat[1,]),
					nBigMatrix = nBigMatrix, minSize = minSize, 
					maxSize = maxSize, envMatrix=envMat,
					fecObj=fecObjList[[k]],integrateType=integrateType, correction=correction)
		}
		
		FmatrixList[[k]] <-  FmatrixList[[k]]
	}
	return(FmatrixList)
}



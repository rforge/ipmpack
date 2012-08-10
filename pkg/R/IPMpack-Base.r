

##can you do this to get rid of the warnings? check 
#require(nlme) #apparently doesn't help!

logit <- function(x) { u<-exp(pmin(x,50)); return(u/(1+u))}


### FUNCTIONS FOR PREDICTING GROWTH AND SURVIVAL ###########################
### AS A FUNCTION OF  SIZE, COV, AND GROWTH AND SURV OBJECTS ################

# Register a growth  generic 2
#
# Parameters - (testing update)
#   size - size measurement this time-step
#   sizeNext - size measurement the next time-step (midpoints)
#   cov - the covariate (light, etc) this time-step
#   grow.obj - the growth object 
#
#   Note "growth" takes a SINGLE value for covariate although can take a VECTOR for size,sizeNext
#
#
# Returns -
#   the growth transition from size to sizeNext
setGeneric("growth",
		function(size, sizeNext, cov, growthObj) standardGeneric("growth"))


# Register a survival  generic
#
#Parameters -
#   size - size measurement this time-step
#   cov - the covariate (light env, etc)
#   surv.obj - the survival object
#
#   Note "surv" takes a SINGLE value for covariate although can take a VECTOR for size
#
#Returns -
#  The survival probability at size
setGeneric("surv",
		function(size, cov, survObj) standardGeneric("surv"))


# Register a cumulative growth  generic
#
# Parameters -
#   size - size measurement this time-step
#   sizeNext - size measurement the next time-step (midpoints)
#   cov - the covariate (light, etc) this time-step
#   grow.obj - the growth object 
#
#   Note "growth" takes a SINGLE value for covariate although can take a VECTOR for size,sizeNext
#
#
# Returns -
#   the growth transition from size to sizeNext
setGeneric("growthCum",
		function(size, sizeNext, cov, growthObj) standardGeneric("growthCum"))


### CLASSES UNDERLYING THE GROWTH / SURVIVAL FUNCTIONS ######

## GROWTH OBJECTS ##
# Create a generic growth object containing a lm
setClass("growthObj",
		representation(fit = "lm",sd="numeric"))

setClass("growthObjPois",
		representation(fit = "glm"))

#setClass("growthObjNegBin",
#		representation(fit = "list"))


# Create a generic growth object with normal errors on increment
setClass("growthObjIncr",
		representation(fit = "lm",sd="numeric"))

# Create a generic growth object with truncated normal errors on increment
setClass("growthObjTruncIncr",
		representation(fit = "list", varcov="matrix"))

# Create a generic growth object with log normal errors on increment
setClass("growthObjLogIncr",
		representation(fit = "lm",sd="numeric"))

# Create a generic growth object with declining errors 
setClass("growthObjDeclineVar",
		representation(fit = "list"))

# Create a generic growth object with declining errors for increment
setClass("growthObjIncrDeclineVar",
		representation(fit = "list"))

# Create a generic growth object with declining errors for logincrement
setClass("growthObjLogIncrDeclineVar",
		representation(fit = "list"))

# Create a generic growth object containing the Hossfeld parameters 
setClass("growthObjHossfeld",
		representation(paras="numeric",
				sd="numeric", 
				logLik="numeric", 
				hessian="matrix"))

# Create a generic growth object with a poisson model 
setClass("growthObjPois",
		representation(fit = "glm"))

## SURVIVAL OBJECTS ##
# Create a generic survival object
setClass("survObj",
		representation(fit = "glm"))

# Create a generic survival object for use where over-dispersion
# modeled, using Diggles approximate correction for the transform
setClass("survObjOverDisp",
		representation(fit = "glm"))



## FECUNDITY OBJECTS ##
# Create a generic fecundity object
setClass("fecObj",
		representation(fitFec = "list",
				fecNames = "character",
				fecConstants = "data.frame",
				offspringSplitter = "data.frame",
				fecByDiscrete = "data.frame",
				offspringTypeRates = "data.frame",
				Transform = "character",
				offspringRel = "lm",
				sdOffspringSize = "numeric")
)

setClass("fecObjMultiCov",
		representation(fitFec = "list",
				fecNames = "character",
				fecConstants = "data.frame",
				offspringSplitter = "data.frame",
				fecByDiscrete = "data.frame",
				offspringTypeRates = "data.frame",
				Transform = "character",
				offspringRel = "lm",
				sdOffspringSize = "numeric")
)




## DISCRETE TRANSITION MATRIX OBJECTS ##
# Create a generic discrete transition matrix object
setClass("discreteTrans",
		representation(nclasses = "numeric",
				discreteTrans = "matrix",
				discreteSurv = "matrix",
				meanToCont = "matrix",
				sdToCont = "matrix",
				distribToDiscrete = "matrix",
				survToDiscrete = "glm"))





######## DEFINE METHODS ##########################################################################################

#Method to obtain probability of survival using
# logistic regression on size with a single covariate
#
#Parameters -
#   size = current size (vector)
#   cov = current discrete covariate (.e.g. light..., single value)
#   survObj = a survival object, containig e.g. a glm fitted to data
#
#Returns -
#  survival probability for given sizes and covariate level
setMethod("surv", 
		c("numeric","data.frame","survObj"),
		function(size,cov,survObj){

			newd <- data.frame(cbind(cov,size=size),
						stringsAsFactors = FALSE)
		
			newd$size2 <- size^2
			newd$size3 <- size^3
			
			if (length(grep("logsize",
							survObj@fit$formula))>0) newd$logsize=log(size)
			if (length(grep("logsize2",
							survObj@fit$formula))>0) newd$logsize2=(log(size))^2
			
#				print(head(newd))
			
			u <- predict(survObj@fit,newd,type="response")
			return(u);
		})

#Method to obtain probability of survival using
# logistic regression on size with a single covariate
#  where the logistic regression was modeled with over-dispersion
#  (e.g., using MCMCglmm) - !over-dispersion assumed to be set to 1
#
#Parameters -
#   size = current size (vector)
#   cov = current discrete covariate (.e.g. light..., single value)
#   survObj = a survival object, containig e.g. a glm fitted to data
#
#Returns -
#  survival probability for given sizes and covariate level
setMethod("surv", 
		c("numeric","data.frame","survObjOverDisp"),
		function(size,cov,survObj){
			newd <- data.frame(cbind(cov,size=size),
					stringsAsFactors = FALSE)
			newd$size2 <- size^2
			newd$size3 <- size^3
			
			if (length(grep("logsize",
							survObj@fit$formula))>0) newd$logsize=log(size)
			if (length(grep("logsize2",
							survObj@fit$formula))>0) newd$logsize2=(log(size))^2
			
			u <- predict(survObj@fit,newd,type="link")
			c2 <- ((16 * sqrt(3))/(15 * pi))^2  #from MCMCglmm course notes, search for c2
			u <- logit(u/sqrt(1+c2)) 
			return(u);
		})



# Method to obtain growth transitions -
# here linear regression to size (with powers) with a single covariate
#
#Parameters -
#   size = size now (vector)
#   sizeNext = size going to  (vector)
#   cov = the covariate (.e.g. light..., single value)
#   growthObj = a growth object
#
#Returns -
#  growth transition probability from size to sizeNext at that covariate level 
setMethod("growth", 
		c("numeric","numeric","data.frame","growthObj"),
		function(size,sizeNext,cov,growthObj){
			newd <- data.frame(cbind(cov,size=size),
					stringsAsFactors = FALSE)
			newd$size2 <- size^2
			newd$size3 <- size^3
			if (length(grep("logsize",
							growthObj@fit$formula))>0) { newd$logsize <- log(size)}
			
			mux <- predict(growthObj@fit,newd,type="response")
			sigmax <-growthObj@sd
			u <- dnorm(sizeNext,mux,sigmax,log=FALSE)  
			return(u);
		})

#  growth transition (poisson model) probability from size to sizeNext at that covariate level 
#	NOTE DO NOT USE THIS WITH AN IPM!!
setMethod("growth", 
		c("numeric","numeric","data.frame","growthObjPois"),
		function(size,sizeNext,cov,growthObj){
			newd <- data.frame(cbind(cov,size=size),
					stringsAsFactors = FALSE)
			newd$size2 <- size^2
			newd$size3 <- size^3
			
			if (length(grep("logsize",
							growthObj@fit$formula))>0) { newd$logsize <- log(size)}
			
			mux <- predict(growthObj@fit,newd,type="response")
			u <- dpois(sizeNext,mux,log=FALSE)  
			return(u);
		})



# growth for predicting next incr 
setMethod("growth", 
		c("numeric","numeric","data.frame","growthObjIncr"),
		function(size,sizeNext,cov,growthObj){
			newd <- data.frame(cbind(cov,size=size),
					stringsAsFactors = FALSE)
			newd$size2 <- size^2
			newd$size3 <- size^3
			
			if (length(grep("logsize",
							growthObj@fit$formula))>0) newd$logsize=log(size)
			
			mux <- predict(growthObj@fit,newd,type="response")
			
			#print(mux)
			
			sigmax <-growthObj@sd
			u <- dnorm(sizeNext,size+mux,sigmax,log=FALSE)  
			return(u); 
		})



# growth for predicting next truncated incr 
setMethod("growth", 
		c("numeric","numeric","data.frame","growthObjTruncIncr"),
		function(size,sizeNext,cov,growthObj){
			require(truncnorm)
			
			newd <- data.frame(cbind(cov,size=size),
					stringsAsFactors = FALSE)
			newd$size2 <- size^2
			newd$size3 <- size^3
			
			if (length(grep("logsize",
							names(growthObj@fit$coefficients)))>0) newd$logsize=log(size)
						
			mux <- .predictMuX(grObj=growthObj,newData=newd,covPred=cov)
			sigmax <- growthObj@fit$sigma
			u <- dtruncnorm(sizeNext,a=size,b=Inf,mean=size+mux,sd=sigmax)  
			return(u); 
		})



# growth for predicting next logincr 
setMethod("growth", 
		c("numeric","numeric","data.frame","growthObjLogIncr"),
		function(size,sizeNext,cov,growthObj){
		
			newd <- data.frame(cbind(cov,size=size),
					stringsAsFactors = FALSE)
			newd$size2 <- size^2
			newd$size3 <- size^3
			
			if (length(grep("logsize",
							growthObj@fit$formula))>0) newd$logsize=log(size)
			
			mux <- predict(growthObj@fit,newd,type="response")
			sigmax <-growthObj@sd
			u <- dlnorm(sizeNext-size,mux,sigmax,log=FALSE)  
			return(u);
		})


# growth for predicting next logincr with decline var
setMethod("growth", 
		c("numeric","numeric","data.frame","growthObjLogIncrDeclineVar"),
		function(size,sizeNext,cov,growthObj){
			newd <- data.frame(cbind(cov,size=size),
					stringsAsFactors = FALSE)
			newd$size2 <- size^2
			newd$size3 <- size^3
			
			if (length(grep("logsize",
							names(growthObj@fit$coefficients)))>0) newd$logsize=log(size)
			
			mux <- .predictMuX(grObj=growthObj,newData=newd,covPred=cov)
			sigmax2 <- growthObj@fit$sigmax2
			var.exp.coef<-growthObj@fit$var.exp.coef
			sigmax2<-sigmax2*exp(2*(var.exp.coef*mux));
		
			u <- dlnorm(sizeNext-size,mux,sqrt(sigmax2),log=FALSE)  
			return(u);
		})

# growth for predicting next logincr with decline var
setMethod("growthCum", 
		c("numeric","numeric","data.frame","growthObjLogIncrDeclineVar"),
		function(size,sizeNext,cov,growthObj){
	
			newd <- data.frame(cbind(cov,size=size),
					stringsAsFactors = FALSE)
			newd$size2 <- size^2
			newd$size3 <- size^3
			
			if (length(grep("logsize",
							names(growthObj@fit$coefficients)))>0) newd$logsize=log(size)
			
			mux <- .predictMuX(grObj=growthObj,newData=newd,covPred=cov)
			sigmax2 <- growthObj@fit$sigmax2
			var.exp.coef<-growthObj@fit$var.exp.coef
			sigmax2<-sigmax2*exp(2*(var.exp.coef*mux));
			
			u <- plnorm(sizeNext-size,mux,sqrt(sigmax2),log=FALSE)  
			return(u);
		})








# Slightly alternative approach, using cumulative probs ##

# growth for predicting next size with a polynomial or log
# using pnorm (i.e. getting cumulative at boundary points and doing difference)
setMethod("growthCum", 
		c("numeric","numeric","data.frame","growthObj"),
		function(size,sizeNext,cov,growthObj){
	
			newd <- data.frame(cbind(cov,size=size),
					stringsAsFactors = FALSE)
			newd$size2 <- size^2
			newd$size3 <- size^3
			if (length(grep("logsize",
							growthObj@fit$formula))>0) newd$logsize=log(size)
			mux <- predict(growthObj@fit,newd,type="response")
			sigmax <-growthObj@sd
			u <- pnorm(sizeNext,mux,sigmax,log=FALSE)  
			return(u);
		})

# growth for predicting next incr with a polynomial or log
# using pnorm (i.e. getting cumulative at boundary points and doing difference)
setMethod("growthCum", 
		c("numeric","numeric","data.frame","growthObjIncr"),
		function(size,sizeNext,cov,growthObj){
	
			newd <- data.frame(cbind(cov,size=size),
					stringsAsFactors = FALSE)
			newd$size2 <- size^2
			newd$size3 <- size^3
			if (length(grep("logsize",
							growthObj@fit$formula))>0) newd$logsize=log(size)
			
			mux <- predict(growthObj@fit,newd,type="response")
			sigmax <-growthObj@sd
			u <- pnorm(sizeNext,size+mux,sigmax,log=FALSE)  
			return(u); 
		})



# growth for predicting next truncated incr with cumulative 
setMethod("growthCum", 
		c("numeric","numeric","data.frame","growthObjTruncIncr"),
		function(size,sizeNext,cov,growthObj){
			require(truncnorm)
		
			newd <- data.frame(cbind(cov,size=size),
					stringsAsFactors = FALSE)
			newd$size2 <- size^2
			newd$size3 <- size^3
			
			if (length(grep("logsize",
							names(growthObj@fit$coefficients)))>0) newd$logsize=log(size)
			
			mux <- .predictMuX(grObj=growthObj,newData=newd,covPred=cov)
			sigmax <- sqrt(growthObj@fit$sigmax)
			
			u <- ptruncnorm(sizeNext,a=size,b=Inf,mean=size+mux,sd=sigmax)  
			return(u); 
		})


# growth for predicting next logincr with a polynomial or log
# using pnorm (i.e. getting cumulative at boundary points and doing difference)
setMethod("growthCum", 
		c("numeric","numeric","data.frame","growthObjLogIncr"),
		function(size,sizeNext,cov,growthObj){
		
			newd <- data.frame(cbind(cov,size=size),
					stringsAsFactors = FALSE)
			newd$size2 <- size^2
			newd$size3 <- size^3
			if (length(grep("logsize",
							growthObj@fit$formula))>0) newd$logsize=log(size)
			
			mux <- predict(growthObj@fit,newd,type="response")
			sigmax <-growthObj@sd
			u <- plnorm(sizeNext-size,mux,sigmax,log=FALSE)  
			return(u);
		})


#Simple growth methods, using  declining variance in growth
# using pnorm (i.e. getting cumulative at boundary points and doing difference)
setMethod("growthCum", 
		c("numeric","numeric","data.frame","growthObjDeclineVar"),
		function(size,sizeNext,cov,growthObj){
		
			newd <- data.frame(cbind(cov,size=size),
					stringsAsFactors = FALSE)
			newd$size2 <- size^2
			newd$size3 <- size^3
			
			if (length(grep("logsize",
							names(growthObj@fit$coefficients))) > 0) newd$logsize=log(size)
			
			mux <- .predictMuX(grObj=growthObj,newData=newd,covPred=cov)
			sigmax2 <- growthObj@fit$sigmax2
			var.exp.coef<-growthObj@fit$var.exp.coef
			sigmax2<-sigmax2*exp(2*(var.exp.coef*mux));
			u <- pnorm(sizeNext,mux,sqrt(sigmax2),log=FALSE)  
			return(u);
		})




#Simple growth methods, using  declining variance in growth
setMethod("growth", 
		c("numeric","numeric","data.frame","growthObjDeclineVar"),
		function(size,sizeNext,cov,growthObj){
			newd <- data.frame(cbind(cov,size=size),
					stringsAsFactors = FALSE)
			newd$size2 <- size^2
			newd$size3 <- size^3
			
			if (length(grep("logsize",
							names(growthObj@fit$coefficients)))>0) newd$logsize=log(size)
			
			mux <- .predictMuX(grObj=growthObj,newData=newd)
			sigmax2 <- growthObj@fit$sigmax2
			var.exp.coef<-growthObj@fit$var.exp.coef
			sigmax2<-sigmax2*exp(2*(var.exp.coef*mux));
			u <- dnorm(sizeNext,mux,sqrt(sigmax2),log=FALSE)  
			return(u);
		})


# function to predict growth object - NOTE THAT THIS will not work with factors / contrasts
.predictMuX <- function(grObj, newData, covPred = 1) {
	dataBase <- newData
	coefNames <- attr(grObj@fit$coefficients, "names")
	coefValues <- as.matrix(grObj@fit$coefficients)
	covNames <- coefNames[grep("covariate", coefNames)]
	covPos <- as.numeric(unlist(strsplit(covNames, "covariate")))
	covPos <- as.numeric(covPos[!is.na(covPos)])
	covDf <- as.data.frame(matrix(0, nrow = dim(newData)[1], ncol = length(covPos)))
	names(covDf) <- covNames
	# Turn "on" intended covariate and build newDataFrame
	if(covPred != 1) {
		covDf[, (covPred - 1)] <- 1
	}
	newData <- cbind(dataBase, covDf)
	newDataSubset <- as.matrix(cbind(1, newData[, (names(newData) %in% coefNames)]))
	predValues <- as.matrix(newDataSubset) %*% matrix(coefValues, length(coefValues), 1)
	return(as.numeric(predValues))
}


#same but with declining variance in growth on incrment
setMethod("growth", 
	c("numeric","numeric","data.frame","growthObjIncrDeclineVar"),
		function(size, sizeNext, cov, growthObj){
			newd <- data.frame(cbind(cov,size=size),
					stringsAsFactors = FALSE)
			newd$size2 <- size^2
			newd$size3 <- size^3
			
			if (length(grep("logsize",
							names(growthObj@fit$coefficients))) > 0) newd$logsize = log(size)
							
			mux <- .predictMuX(grObj=growthObj,newData=newd,covPred=cov)
			sigmax2 <- growthObj@fit$sigmax2
			var.exp.coef<-growthObj@fit$var.exp.coef
			sigmax2<-sigmax2*exp(2*(var.exp.coef*mux));
			u <- dnorm(sizeNext, size + mux, sqrt(sigmax2), log = FALSE)  
			return(u);
		})


## Define a new growth method for Hossfeld growth (classic midpoint rule approach)
setMethod("growth", c("numeric", "numeric", "data.frame", "growthObjHossfeld"), 
		function(size, sizeNext, cov, growthObj) { 
			mux <- size+Hossfeld(size, growthObj@paras) 
			sigmax <- growthObj@sd 
			u <- dnorm(sizeNext, mux, sigmax, log = FALSE) 
			return(u)
		}) 


## Define a new growth method for Hossfeld growth (classic midpoint rule approach)
setMethod("growthCum", c("numeric", "numeric", "data.frame", "growthObjHossfeld"), 
		function(size, sizeNext, cov, growthObj) { 
			mux <- size+Hossfeld(size, growthObj@paras) 
			sigmax <- growthObj@sd 
			u <- pnorm(sizeNext, mux, sigmax, log = FALSE) 
			return(u)
		}) 




# Method combining growth and survival for doing outer (not a generic, so that don't have to
# define for all the different classes / growth and surv take care of that)
growSurv <- function(size,sizeNext,cov,growthObj,survObj){
	growth(size, sizeNext, cov, growthObj) * surv(size,cov,survObj)
}






### CLASSES AND FUNCTIONS FOR MATRICES (ENV, TMATRIX [compound or not], FMATRIX) #################

#Class for the matrix that holds the env matrix 
setClass("envMatrix",
		representation(nEnvClass = "numeric"), #number of covariate levels
		contains="matrix")


#Class for the matrix that holds the IPM
setClass("IPMmatrix",
		representation(nDiscrete = "numeric", #number of discrete classes
				nEnvClass = "numeric", #number of covariate levels
				nBigMatrix = "numeric", #the resolution of the IPM
				meshpoints = "numeric",
				env.index = "numeric",
				names.discrete = "character"),
		contains="matrix")


#Function creates a single T.IPM (survival and growth
#transitions only) for a chosen size range, env category
#(single one at a time) and growth and survival objects
#
#Parameters -
#   nEnvClass - the number of env classes, cannot be !=1, defaults to 1 
#   nBigMatrix - the number of size bins in the model
#   minSize - lower end of desired size range
#   maxSize - upper end of desired size range
#   chosenCov - current level of covariate - can be vector for continuous
#                 temporal stochasticity case
#   growObj - a growth object
#   survObj - a survival object
#   discreteTrans - discrete transition object, defaults to 1
#   integrateType - what type of integration?
#                    current options are "midpoint" (using pdf)
#                    or "cumul" (using cdf)
#   correction - do you not want to correct for integration errors ('none')
#                or correct by multiplying each column by a constant ('constant')

#
#Returns -
#  an IPM object (with or without discrete classes)
createIPMTmatrix <- function (nEnvClass = 1, nBigMatrix = 50, minSize = -1, maxSize = 50, 
		chosenCov = data.frame(covariate = 1), growObj, survObj, 
		discreteTrans = 1, integrateType = "midpoint", correction = "none") {
	if (class(growObj) == "growthObjPois") 
		print("warning: IPMs not appropriate with discrete growth processes")
	b <- minSize + c(0:nBigMatrix) * (maxSize - minSize)/nBigMatrix
	y <- 0.5 * (b[1:nBigMatrix] + b[2:(nBigMatrix + 1)])
	h <- y[2] - y[1]
	if (integrateType == "midpoint") {
		get.matrix <- t(outer(y, y, growSurv, cov = chosenCov, 
						growthObj = growObj, survObj = survObj)) * h
	}
	if (integrateType == "cumul") {
		get.matrix.cum <- t(outer(y, b, growthCum, cov = chosenCov, 
						growthObj = growObj))
		get.matrix <- get.matrix.cum[2:(nBigMatrix + 1), ] - 
				get.matrix.cum[1:nBigMatrix, ]
		get.matrix[nBigMatrix, nBigMatrix] <- get.matrix[nBigMatrix, 
				nBigMatrix] + (1 - sum(get.matrix[, nBigMatrix]))
		get.matrix <- t(t(get.matrix) * surv(size = y, cov = chosenCov, 
						survObj = survObj))
	}
	if (correction == "constant") {
		nvals <- colSums(get.matrix,na.rm=TRUE)
		loc0 <- which(nvals == 0 , arr.ind = TRUE)
		if (length(loc0) > 0) {
			print("warnings - columns that sum to 0 or that have NAs - assuming survival is along the diagonal; plot your Tmatrix to check it")
			get.matrix[,loc0] <- 0
			get.matrix[cbind(loc0, loc0)] <- surv(size = y[loc0], 
					cov = chosenCov, survObj = survObj)
		}
		nvals <- colSums(get.matrix,na.rm=TRUE)
		get.matrix <- t((t(get.matrix)/nvals) * surv(size = y, 
						cov = chosenCov, survObj = survObj))
	}
	rc <- new("IPMmatrix", nDiscrete = 0, nEnvClass = 1, nBigMatrix = nBigMatrix, 
			nrow = 1 * nBigMatrix, ncol = 1 * nBigMatrix, meshpoints = y, 
			env.index = rep(1:nEnvClass, each = nBigMatrix), names.discrete = "")
	rc[, ] <- get.matrix
	if (class(discreteTrans) == "discreteTrans") {
		nDisc <- ncol(discreteTrans@discreteSurv)
		survToDiscrete <- predict(discreteTrans@survToDiscrete, 
				data.frame(size = y, size2 = (y * y)), type = "response")
		cont.to.cont <- get.matrix * matrix(1 - survToDiscrete, 
				nrow = nBigMatrix, ncol = nBigMatrix, byrow = T)
		disc.to.disc <- discreteTrans@discreteTrans[1:nDisc, 
						1:nDisc] * matrix(c(discreteTrans@discreteSurv), 
						nrow = nDisc, ncol = nDisc, byrow = T)
		disc.to.cont <- matrix(NA, ncol = nDisc, nrow = nBigMatrix)
		cont.to.disc <- matrix(NA, nrow = nDisc, ncol = nBigMatrix)
		for (j in 1:nDisc) {
			tmp <- dnorm(y, discreteTrans@meanToCont[j], discreteTrans@sdToCont[j]) * 
					h
			if (correction == "constant") 
				tmp <- tmp/sum(tmp)
			disc.to.cont[, j] <- discreteTrans@discreteSurv[,j] * 
					discreteTrans@discreteTrans["continuous", j] * tmp
			cont.to.disc[j, ] <- discreteTrans@distribToDiscrete[j, ] * surv(y, chosenCov, survObj) * survToDiscrete
		}
		get.disc.matrix <- rbind(cbind(disc.to.disc, cont.to.disc), 
				cbind(disc.to.cont, cont.to.cont))
		rc <- new("IPMmatrix", nDiscrete = nDisc, nEnvClass = 1, 
				nBigMatrix = nBigMatrix, nrow = 1 * nBigMatrix + 
						nDisc, ncol = 1 * nBigMatrix + nDisc, meshpoints = y, 
				env.index = rep(1:nEnvClass, each = nBigMatrix), 
				names.discrete = rownames(discreteTrans@discreteTrans)[1:nDisc])
		rc[, ] <- get.disc.matrix
	}
	return(rc)
}


#Function creates a combo of T.IPMs for a chosen
#size range, with range env categories and transitions 
#between them
#and growth and survival objects
#
#Parameters -
#   nEnvClass - the number of env classes, defaults to 2, should match dim envMatrix
#   nBigMatrix - the number of size bins in the model
#   minSize - lower end of desired size range
#   maxSize - upper end of desired size range
#   envMatrix - a matrix describing transiions between env
#   growObj - a growth object
#   survObj - a survival object
#   discreteTrans - discrete transition object, defaults to 1
#   integrateType - what type of integration?
#                    current options are "midpoint" (using pdf)
#                    or "cumul" (using cdf)
#   correction - do you not want to correct for integration errors ('none')
#                or correct by multiplying each column by a constant ('constant')
#
#Returns -
#  an IPM object


createCompoundTmatrix <- function(nEnvClass = 2,
		nBigMatrix = 50,
		minSize = -1,
		maxSize = 50,
		envMatrix,
		growObj,
		survObj,
		discreteTrans=1,
		integrateType="midpoint",
		correction="none") {
	
	
	#warnings...
	if (nEnvClass!=nrow(envMatrix)) {
		print(paste("Dim of envMatrix not equal to nEnvClass. Adjusted to",nrow(envMatrix)))
		nEnvClass <- nrow(envMatrix)
	}
	
	# boundary points b and mesh points y
	b<-minSize+c(0:nBigMatrix)*(maxSize-minSize)/nBigMatrix;
	y<-0.5*(b[1:nBigMatrix]+b[2:(nBigMatrix+1)]);
	
	# step size for mid point rule, see equations 4 and 5
	h<-y[2]-y[1]
	
	#establish how how many discrete classes there are
	if (class(discreteTrans)=="discreteTrans") nDisc <- ncol(discreteTrans@discreteSurv) else nDisc <- 0
	
	#indexes for slotting in IPMs
	indexes <- rep(1:nEnvClass,each=(nBigMatrix+nDisc))
	
	#megamatrix
	megamatrix <- matrix(0,(nBigMatrix+nDisc)*nEnvClass,(nBigMatrix+nDisc)*nEnvClass) 
	
	#print(indexes)
	
	#loop over habitats / environments
	for (k in 1:nEnvClass) { 
		#IPM for individuals starting in env k
		
		#print(k)
		
		if (integrateType=="midpoint") { 
			get.matrix <- (maxSize-minSize)*
					t(outer(y,y,growSurv,cov=data.frame(covariate=as.factor(k)),
									growthObj=growObj,survObj=survObj))/nBigMatrix  
		}
		if (integrateType=="cumul") {
			get.matrix.cum <- 
					t(outer(y,b,growthCum,cov=data.frame(covariate=as.factor(k)),
									growthObj=growObj))
			get.matrix <- get.matrix.cum[2:(nBigMatrix+1),]-get.matrix.cum[1:nBigMatrix,]
			get.matrix <- t(t(get.matrix)*surv(size=y,cov=data.frame(covariate=as.factor(k)),survObj=survObj))
			
		}
		
		#fix any integration issues reducing survival by dividing by col sums and multiply by survival
	if (correction == "constant") {
		nvals <- colSums(get.matrix,na.rm=TRUE)
		loc0 <- which(nvals == 0 , arr.ind = TRUE)
		if (length(loc0) > 0) {
			print("warnings - columns that sum to 0 or that have NAs - assuming survival is along the diagonal; plot your Tmatrix to check it")
			get.matrix[,loc0] <- 0
			get.matrix[cbind(loc0, loc0)] <- surv(size = y[loc0], 
					cov = chosenCov, survObj = survObj)
		}
		nvals <- colSums(get.matrix,na.rm=TRUE)
		get.matrix <- t((t(get.matrix)/nvals) * surv(size = y, 
						cov = chosenCov, survObj = survObj))
	}
	
		#names of discrete classes default
		nmes <- ""
		
		
		# In case of discrete classes, take the IPM constructed above and add discrete classes defined in discreteTrans
		if (class(discreteTrans)=="discreteTrans") {
			nmes <- rownames(discreteTrans@discreteTrans)
			survToDiscrete <- predict(discreteTrans@survToDiscrete,data.frame(size=y,size2=(y*y)),type="response")
			cont.to.cont <- get.matrix*matrix(1-survToDiscrete,nrow=nBigMatrix,ncol=nBigMatrix,byrow=T)
			disc.to.disc <- discreteTrans@discreteTrans[1:nDisc,1:nDisc]*matrix(c(discreteTrans@discreteSurv),
					nrow=nDisc,ncol=nDisc,byrow=T)
			disc.to.cont <- matrix(NA,ncol=nDisc,nrow=nBigMatrix)
			cont.to.disc <- matrix(NA,nrow=nDisc,ncol=nBigMatrix)
			
			#print(discreteTrans@meanToCont)
			#print(discreteTrans@sdToCont)
			
			for (j in 1:nDisc) {
				tmp<-dnorm(y,discreteTrans@meanToCont[j],discreteTrans@sdToCont[j])*h
				if (correction=="constant") tmp<-tmp/sum(tmp)
				
				#print(tmp)
				
				disc.to.cont[,j] <- discreteTrans@discreteSurv[,j]*discreteTrans@discreteTrans["continuous",j]*tmp
				cont.to.disc[j,] <- discreteTrans@distribToDiscrete[j,]*surv(y,cov=data.frame(covariate=as.factor(k)),survObj)*survToDiscrete
			}
			
			get.matrix <- rbind(cbind(disc.to.disc,cont.to.disc),cbind(disc.to.cont,cont.to.cont))
		}
		
		# transit them
		subset <- c(1:nEnvClass)[envMatrix[,k]>0]
		for (j in subset) { 
			megamatrix[indexes==j,indexes==k] <- get.matrix*envMatrix[j,k]; 
		}
		
	}
	
	rc <- new("IPMmatrix",
			nEnvClass = nEnvClass, 
			nBigMatrix = nBigMatrix,
			nrow = nEnvClass*(nBigMatrix + nDisc),
			ncol = nEnvClass*(nBigMatrix + nDisc),
			meshpoints = y,
			env.index = rep(1:nEnvClass, each = nBigMatrix,
					names.discrete = nmes))
	
	rc[,] <- megamatrix
	
	return(rc) 
}


## Extra functions for use with outer in building the IPMFmatrix

## Get raw numbers of offspring produced by every size class by multiplying up the constants,
## and doing all the "predict: values needed; and taking out only the babies that go to the continuous classes

.fecRaw <- function(x,cov=data.frame(covariate=1),fecObj) { 
	
	newd <- data.frame(cbind(cov,size=x),
			stringsAsFactors = FALSE)
	
	newd$size2 <- x^2
	newd$size3 <- x^3
	
	if (length(fecObj@offspringRel)>1) {
		if (length(grep("logsize", fecObj@offspringRel$formula))>0) { newd$logsize <- log(x)}
	}
	
	fecObj@fecConstants[is.na(fecObj@fecConstants)] <- 1
	
	#fecundity rates
	fecValues <- matrix(c(rep(1,length(fecObj@fitFec)),unlist(fecObj@fecConstants)),
			ncol=length(x),nrow=length(fecObj@fitFec)+
					length(fecObj@fecConstants))
	#rownames(fecValues) <- c(fecObj@fecNames,names(fecObj@fecConstants))
	for (i in 1:length(fecObj@fitFec)) fecValues[i,] <- predict(fecObj@fitFec[[i]],newd,type="response")
	if (length(grep("log",fecObj@Transform))>0) for (i in grep("log",fecObj@Transform)) fecValues[i,]<-exp(fecValues[i,])
	if (length(grep("sqrt",fecObj@Transform))>0) for (i in grep("sqrt",fecObj@Transform)) fecValues[i,]<-(fecValues[i,])^2
	if (length(grep("-1",fecObj@Transform))>0) for (i in grep("-1",fecObj@Transform)) fecValues[i,]<-fecValues[i,]+1
	prodFecValues <- apply(fecValues[which(fecObj@offspringTypeRates[,"continuous"]==1),],2,prod)*unlist(fecObj@offspringSplitter["continuous"])
	return(list(prodFecValues,fecValues))
}


## A function that outer can use showing numbers from x to y via production and distribution offspring
.fecPreCensus <- function(x,y,cov=data.frame(covariate=1),fecObj) {
	newd <- data.frame(cbind(cov,size=x),
			stringsAsFactors = FALSE)	
	newd$size2 <- x^2
	newd$size3 <- x^3
		
	if (length(grep("logsize",
					fecObj@offspringRel$formula))>0) { newd$logsize <- log(x)}
	u <- .fecRaw(x=x,cov=cov,fecObj=fecObj)[[1]]*
			dnorm(y,predict(fecObj@offspringRel,newdata=newd, type="response"),
					fecObj@sdOffspringSize)
	
	#print(cbind(y,predict(fecObj@offspringRel)))
	
	
	return(u)
}


# TO DO .fecPostCensus properly (i.e. include size changes between the census and reproduction event) the following steps have to be made:
# 1. use growSurv to determine what the distribution of size is at the reproduction event given initial size x
# 2. over this distribution of sizes x2 at the reproduction event, what are the expected number of offspring: per x2: .fecRaw(x=x2,...)
# 3. multiply the expected number of offspring per x2 with the probability that an offspring is of size y, using dnorm(y,predict(..., newdata=newd2, ...) where newd2 is calculated for all levels of x2
# all in all, Eelke wonders if the outer-solution is still useful in the .fecPostCensus case, since x2 needs to have a certain distribution, which is not passed down to .fecPostCensus... 

## REMOVE GROWTH FROM THIS - note that this means 
#### growth obj generally not needed down below.....
## A function that outer can use showing numbers from x to y via production, growth, survival and distribution offspring
.fecPostCensus <- function(x,y,cov=data.frame(covariate=1),fecObj, growObj,survObj) {
	newd <- data.frame(cbind(cov,size=x),
			stringsAsFactors = FALSE)
	
	newd$size2 <- x^2
	newd$size3 <- x^3
	if (length(grep("logsize",fecObj@offspringRel$formula))>0 |
			length(grep("logsize",growObj@fit$formula))>0) { newd$logsize <- log(x)}            
	u <- .fecRaw(x=x,cov=cov,fecObj=fecObj)[[1]]*
			dnorm(y,predict(fecObj@offspringRel,newdata=newd, type="response"),fecObj@sdOffspringSize)*
   			surv(size=x, cov=cov, survObj=survObj)
	return(u)
}


## A function that outer can use giving pnorm for offspring reprod
.offspringCum <- function(x,y,cov=data.frame(covariate=1),fecObj) {
	newd <- data.frame(cbind(cov,size=x),
			stringsAsFactors = FALSE)
	
	newd$size2 <- x^2
	newd$size3 <- x^3
	if (length(grep("logsize",fecObj@offspringRel$formula))>0) { newd$logsize <- log(x)}            
	u <- pnorm(y,predict(fecObj@offspringRel,newdata=newd, type="response"),fecObj@sdOffspringSize)
	return(u)
}


#Function creates a single F.IPM (fecundity transitions only)
#for a chosen size range, env category (single one at a time)
#and survival and fecundity objects (assume survival 
#precedes growth; could do otherwise...)
#
#Parameters -
#   nEnvClass - the number of env classes, cannot be !=1, defaults to 1 
#   nBigMatrix - the number of size bins in the model
#   minSize - lower end of desired size range
#   maxSize - upper end of desired size range
#   chosenCov - current level of covariate
#   fecObj - a fecundity object
#   integrateType - options include "midpoint" "cumul" 
#   etc...
#Returns -
#  an IPM object

createIPMFmatrix <- function(fecObj,
		nEnvClass = 1,
		nBigMatrix = 50,
		minSize = -1,
		maxSize = 50,
		chosenCov = data.frame(covariate=1),
		integrateType="midpoint",
		correction="none",
		preCensus=TRUE,
		survObj=NULL,
		growObj=NULL) {
	
	# boundary points b and mesh points y
	b<-minSize+c(0:nBigMatrix)*(maxSize-minSize)/nBigMatrix;
	y<-0.5*(b[1:nBigMatrix]+b[2:(nBigMatrix+1)]);
	
	# step size for mid point rule, see equations 4 and 5
	h<-y[2]-y[1]
	
	fecObj@fecConstants[is.na(fecObj@fecConstants)] <- 1
	
	tmp <- matrix(0,ncol=length(y),nrow=length(y))
	
    # 1. pre-census
	if (preCensus) { 
		#print("here")
		if (integrateType=="midpoint"&fecObj@offspringSplitter$continuous>0) { 
			tmp <- t(outer(X=y,Y=y,.fecPreCensus,cov=chosenCov,fecObj=fecObj))*h 
		}
		if (integrateType=="cumul"&fecObj@offspringSplitter$continuous>0) {
			#offspring extremes (pnorm) 
			tmp.cum <- t(outer(X=y,Y=b,.offspringCum,cov=chosenCov,
							fecObj=fecObj))
			tmp <- tmp.cum[2:(nBigMatrix+1),]-tmp.cum[1:nBigMatrix,]
			#put in seed production
			tmp <- t(t(tmp)*.fecRaw(x=y,cov=chosenCov,fecObj=fecObj)[[1]])      
		}
		
		if (correction=="constant"&fecObj@offspringSplitter$continuous>0) {
			# in this case, column sums should equal raw fecundity
			correction <-.fecRaw(x=y,cov=chosenCov,fecObj=fecObj)[[1]]/colSums(tmp)
			tmp <- t(t(tmp)*correction)
		}
		
	# 2. post-census
	} else {
		#print ("Warning: in the current version of IPMpack, createIPMFmatrix still ignores the growObj you provided for your post-breeding F matrix. This will be included in a later version. Survival until breeding is already included in this version.")
		if (integrateType=="midpoint"&fecObj@offspringSplitter$continuous>0) { 
			tmp <- t(outer(X=y,Y=y,.fecPostCensus,
							cov=chosenCov,fecObj=fecObj, growObj=growObj,
							survObj=survObj))*h 
		}
		if (integrateType=="cumul"&fecObj@offspringSplitter$continuous>0) {
			#make the extreme bins offspring matrix
			tmp.cum <- t(outer(X=y,Y=b,.offspringCum,cov=chosenCov,
							fecObj=fecObj))
			tmpBabies <- tmp.cum[2:(nBigMatrix+1),]-tmp.cum[1:nBigMatrix,]
			
			#make the extreme bins growth matrix
			tmp.cum <- t(outer(X=y,Y=b,growthCum,cov=chosenCov,
							growObj=growObj))
			tmpGrowth <- tmp.cum[2:(nBigMatrix+1),]-tmp.cum[1:nBigMatrix,]
			tmpGrowth[nBigMatrix,nBigMatrix] <- tmpGrowth[nBigMatrix,nBigMatrix]+
					(1-sum(tmpGrowth[,nBigMatrix]))
			
			#put in survival and seed production
			tmp <- t(t(tmpBabies*tmpGrowth)*surv(size=y,cov=chosenCov,survObj=survObj)*.fecRaw(x=y,cov=chosenCov,fecObj=fecObj)[[1]])
			
		}
		
		if (correction=="constant"&fecObj@offspringSplitter$continuous>0) {
			# in this case, column sums should equal raw fecundity * survival
			correction <-.fecRaw(x=y,cov=chosenCov,fecObj=fecObj)[[1]]*surv(size=y,cov=chosenCov,survObj=survObj)/colSums(tmp)
			tmp <- t(t(tmp)*correction)
		}
		
	}
	
	get.matrix <- to.cont <- tmp
	
	#discrete classes
	nDisc <- length(fecObj@offspringSplitter)-1
	namesDiscrete <- "NA"
	if (nDisc>0) {
		namesDiscrete <- colnames(fecObj@offspringSplitter[1:nDisc])
		to.discrete <- matrix(0,nrow=nDisc,ncol=nBigMatrix)
		for (i in 1:nDisc) to.discrete[i,] <- apply(.fecRaw(x=y,cov=chosenCov,fecObj=fecObj)[[2]][which(fecObj@offspringTypeRates[,namesDiscrete[i]]==1),],2,prod)*unlist(fecObj@offspringSplitter[namesDiscrete[i]])
		from.discrete <- matrix(0,ncol=nDisc,nrow=nDisc+nBigMatrix)
		if (names(fecObj@fecByDiscrete)[1]!="NA.") {
			if (sum(names(fecObj@fecByDiscrete)!=namesDiscrete)>0) stop ("Error - the names of the discrete classes as you provided for the data.frame fecByDiscrete are not 100% the same discrete class names in your data.frame offspringSplitter. They should also be in alphabetical order.")
			if (sum(fecObj@fecByDiscrete)>0) {
				print ("Warning - number and sizes of offspring produced by individuals in discrete classes cannot be calculated when offspring size is a function of parent size. The Fmatrix contains zeros instead. Only solution at this point: change the F matrix yourself afterwards.")
			}
		}
		get.matrix <- cbind(from.discrete,rbind(to.discrete,to.cont))
	}
	
	#print(colSums(get.matrix))
	
	#warning about negative numbers
	if (min(get.matrix)<0) { 
		print("Warning: fertility values < 0 exist in matrix, consider transforms. Negative values set to zero") 
		get.matrix[get.matrix<0] <- 0
	}
	
	rc <- new("IPMmatrix",
			nDiscrete = nDisc,
			nEnvClass = 1, 
			nBigMatrix = nBigMatrix,
			nrow = 1*nBigMatrix+nDisc,
			ncol =1*nBigMatrix+nDisc,
			meshpoints = y,
			env.index = rep(1:nEnvClass,each=nBigMatrix),
			names.discrete=namesDiscrete)
	rc[,] <-get.matrix   
	
	return(rc)
}





#Function creates a combo of F.IPMs for a chosen
#size range, with range env categories and transitions 
#between them
#and growth and survival objects
#
#Parameters -
#   nEnvClass - the number of env classes, defaults to 2, should match dim envMatrix
#   nBigMatrix - the number of size bins in the model
#   minSize - lower end of desired size range
#   maxSize - upper end of desired size range
#   envMatrix - a matrix describing transiions between env
#   fecObj - a fecundity object
#   integrateType - NOT YET IMPLEMENTED
#
#Returns -
#  an IPM object


createCompoundFmatrix <- function(nEnvClass = 2,
		nBigMatrix = 50,
		minSize = -1,
		maxSize = 50,
		envMatrix,
		fecObj,
		integrateType="midpoint",
		correction="none",
		preCensus=TRUE,
		survObj=NULL,
		growObj=NULL) {
	
	#warnings...
	if (nEnvClass!=nrow(envMatrix)) {
		print(paste("Dim of envMatrix not equal to nEnvClass. Adjusted to",nrow(envMatrix)))
		nEnvClass <- nrow(envMatrix)
	}
	
	# boundary points b and mesh points y
	b<-minSize+c(0:nBigMatrix)*(maxSize-minSize)/nBigMatrix;
	y<-0.5*(b[1:nBigMatrix]+b[2:(nBigMatrix+1)]);
	
	# step size for mid point rule, see equations 4 and 5
	h<-y[2]-y[1]
	
	#establish how how many discrete classes there are
	if (ncol(fecObj@offspringSplitter)>1) nDisc <- ncol(fecObj@offspringSplitter)-1 else nDisc <- 0
	
	#indexes for slotting in IPMs
	indexes <- rep(1:nEnvClass,each=(nBigMatrix+nDisc))
	#print(indexes)
	
	
	#megamatrix
	megamatrix <- matrix(0,(nBigMatrix+nDisc)*nEnvClass,(nBigMatrix+nDisc)*nEnvClass) 
	
	#loop over habitats / environments
	for (k in 1:nEnvClass) {
		
		get.matrix <- createIPMFmatrix(nEnvClass =1,
				nBigMatrix = nBigMatrix,
				minSize = minSize,
				maxSize = maxSize,
				chosenCov = data.frame(covariate=as.factor(k)),
				fecObj=fecObj,
				integrateType=integrateType,
				correction=correction,preCensus=preCensus,
				survObj=survObj,
				growObj=growObj)
		
		#print(range(get.matrix))
		#print(get.matrix[1:5,1:5])
		
		# transit them
		subset <- c(1:nEnvClass)[envMatrix[,k]>0]                 
		for (j in subset) {
			megamatrix[indexes==j,indexes==k] <- get.matrix@.Data*envMatrix[j,k]; 
		}
		
	}
	
	#warning about negative numbers should appear from createIPMFmatrix
	
	
	rc <- new("IPMmatrix",
			nEnvClass = nEnvClass, 
			nBigMatrix = nBigMatrix,
			nrow = nEnvClass*(nBigMatrix+nDisc),
			ncol =nEnvClass*(nBigMatrix+nDisc),
			meshpoints = y,
			env.index = rep(1:nEnvClass,each=nBigMatrix),
			names.discrete=get.matrix@names.discrete)
	
	
	rc[,] <- megamatrix
	
	return(rc) 
}




#Function creates a single F.IPM (fecundity transitions only)
#for a chosen size range, env category (single one at a time)
#and survival and fecundity objects (assume survival 
#precedes growth; could do otherwise...)
#
#Parameters -
#   nEnvClass - the number of env classes, cannot be !=1, defaults to 1 
#   nBigMatrix - the number of size bins in the model
#   minSize - lower end of desired size range
#   maxSize - upper end of desired size range
#   chosenCov - current level of covariate
#   fecObj - a fecundity object
#   integrateType - options include "midpoint" "cumul" 
#   etc...
#Returns -
#  an IPM object


createIPMCmatrix <- function(clonalObj,
		nEnvClass = 1,
		nBigMatrix = 50,
		minSize = -1,
		maxSize = 50,
		chosenCov = data.frame(covariate=1),
		integrateType="midpoint",
		correction="none") {
	
	rc <- createIPMFmatrix(fecObj=clonalObj,
						  nEnvClass=nEnvClass,
						  nBigMatrix=nBigMatrix,
						  minSize=minSize,
						  maxSize=maxSize,
						  chosenCov=chosenCov,
						  integrateType=integrateType,
						  correction=correction)
	
	return(rc)
}





#Function creates a combo of F.IPMs for a chosen
#size range, with range env categories and transitions 
#between them
#and growth and survival objects
#
#Parameters -
#   nEnvClass - the number of env classes, defaults to 2, should match dim envMatrix
#   nBigMatrix - the number of size bins in the model
#   minSize - lower end of desired size range
#   maxSize - upper end of desired size range
#   envMatrix - a matrix describing transiions between env
#   fecObj - a fecundity object
#   integrateType - NOT YET IMPLEMENTED
#
#Returns -
#  an IPM object


createCompoundCmatrix <- function(nEnvClass = 2,
		nBigMatrix = 50,
		minSize = -1,
		maxSize = 50,
		envMatrix,
		clonalObj,
		integrateType="midpoint",
		correction="none") {
	
		rc <- createCompoundFmatrix(nEnvClass = nEnvClass,
				nBigMatrix = nBigMatrix,
				minSize = minSize,
				maxSize = maxSize,
				envMatrix = envMatrix,
				fecObj = clonalObj,
				integrateType=integrateType,
				correction=correction)
	
	return(rc) 
}





# For a single Tmatrix (!not compound and no discrete stages!), this functions makes a series
# of diagnostic plots - this is defined for growthObj,
# growthObjIncr - modification required
# if other objects used
#
# Parameters - the Tmatrix
#            - growObj - growth object used to build it
#            - survObj - survival object used to build it
#            - dff - the data from which it was built
#            - integrateType - "midpoint", or "cumul" - should
#                 correspond to what the IPM was built with
#            - do you want to implement the corrections for integration? 
# Returns - 
#

diagnosticsTmatrix <- function (Tmatrix, growObj, survObj, dff, 
		integrateType = "midpoint", 
		correction = "none", cov = data.frame(covariate = 1), 
		sizesToPlot=c()) {
	print("Range of Tmatrix is ")
	print(range(c(Tmatrix)))
	if (Tmatrix@meshpoints[1] > 0) 
		new.min <- Tmatrix@meshpoints[1]/2
	else new.min <- Tmatrix@meshpoints[1] * 1.5
	
	#colours for 1) current; 2) bigger size range; 3) bigger no bins
	cols <- c("black","tomato","darkblue")
	ltys <- c(1,1,3)
	
	#matrix with bigger size range
	Tmatrix1 <- createIPMTmatrix(nEnvClass = 1, nBigMatrix = length(Tmatrix@meshpoints), 
			minSize = new.min, maxSize = 1.5 * max(Tmatrix@meshpoints), 
			chosenCov = cov, growObj = growObj, survObj = survObj, 
			integrateType = integrateType, correction = correction)
	
	if (sum(is.na(Tmatrix1)) > 0) {
		print("Tmatrix with extended size range returns NAs; changing these to 0, and putting survival value onto diagonal for columns that sum to zero")
		Tmatrix1[is.na(Tmatrix1)] <- 0
		bad <- which(colSums(Tmatrix1) == 0, arr.ind = TRUE)
		if (length(bad) > 0) 
			Tmatrix1[cbind(bad, bad)] <- surv(size = Tmatrix1@meshpoints[bad], 
					cov = cov, survObj = survObj)
	}
	
	
	#matrix with bigger number of bins	  
	Tmatrix2 <- createIPMTmatrix(nEnvClass = 1, nBigMatrix = floor(length(Tmatrix@meshpoints) * 
							1.5), minSize =Tmatrix@meshpoints[1], maxSize = max(Tmatrix@meshpoints), 
			chosenCov = cov, growObj = growObj, survObj = survObj, 
			integrateType = integrateType, correction = correction)
	
	if (sum(is.na(Tmatrix2)) > 0) {
		print("Tmatrix with extended number of bins returns NAs; changing these to 0, and putting survival value onto diagonal for columns that sum to zero")
		Tmatrix2[is.na(Tmatrix2)] <- 0
		bad <- which(colSums(Tmatrix2) == 0, arr.ind = TRUE)
		if (length(bad) > 0) 
			Tmatrix2[cbind(bad, bad)] <- surv(size = Tmatrix2@meshpoints[bad], 
					cov = cov, survObj = survObj)
	}
	
	#start plots - put original Tmatrix in black
	#par(mfrow = c(3, 3), bty = "l")
	xlims <- range(c(Tmatrix@meshpoints, Tmatrix1@meshpoints,
					dff$size, dff$sizeNext), na.rm = TRUE)
	
	par(mfrow = c(1, 3), bty = "l",pty="s", mar=c(5,4,4,1))
	a1 <- hist(c(dff$size, dff$sizeNext), xlab = "Sizes observed", axes=FALSE,
			ylab = "Frequency", main = "", xlim = xlims, col="lightgrey", border="lightgrey",plot=TRUE)
	axis(1); axis(2)
	
	lcs <- c(0.8,0.6,0.4)	
	lcs.x <- mean(xlims) 
	
	text(lcs.x ,max(a1$counts)*lcs[1],"Current",pos=3,col=cols[1])
	arrows(Tmatrix@meshpoints[1], max(a1$counts)*lcs[1],Tmatrix@meshpoints[length(Tmatrix@meshpoints)], max(a1$counts)*lcs[1],col=cols[1], length=0.1, code=3,lty=ltys[1])
	text(lcs.x ,max(a1$counts)*lcs[2],"Extended range",pos=3,col=cols[2])
	arrows(Tmatrix1@meshpoints[1], max(a1$counts)*lcs[2],Tmatrix1@meshpoints[length(Tmatrix1@meshpoints)], max(a1$counts)*lcs[2],col=cols[2], length=0.1, code=3,lty=ltys[1])
	text(lcs.x ,max(a1$counts)*lcs[3],"Increased bins",pos=3,col=cols[3])
	arrows(Tmatrix2@meshpoints[1], max(a1$counts)*lcs[3],Tmatrix2@meshpoints[length(Tmatrix2@meshpoints)], max(a1$counts)*lcs[3],col=cols[3], length=0.1, code=3,lty=ltys[1])
	
	#legend("topright", legend = "fitted range", col = "black", lty = 2, bty = "n") ##!!! change this
	title("Size range")
	
	#survival sums
	lims <- range(c(colSums(Tmatrix), surv(Tmatrix@meshpoints, cov, survObj)))
	plot(colSums(Tmatrix), surv(Tmatrix@meshpoints, cov, survObj), 
			type = "n", xlab = "Surviving in Tmatrix", ylab = "Should be surviving",col=cols[1],lty=ltys[1], 
			xlim=lims,ylim=lims)
	abline(0, 1, lwd=2,col="grey")
	points(colSums(Tmatrix), surv(Tmatrix@meshpoints, cov, survObj), type = "l", col = cols[1],lty=ltys[1])
	points(colSums(Tmatrix1), surv(Tmatrix1@meshpoints, cov, survObj), type = "l", col = cols[2],lty=ltys[2])
	points(colSums(Tmatrix2), surv(Tmatrix2@meshpoints, cov, survObj), type = "l", col = cols[3],lty=ltys[3])
	title("Survival")
	#text(lims[2],lims[2],"(0,1)",pos=1)	
	
	# mtext("Points should lie along the 0,1 line, shown in grey", side=1,outer=TRUE,line=-1)
	
	
	
	LE <- meanLifeExpect(Tmatrix)
	LE1 <- meanLifeExpect(Tmatrix1)
	LE2 <- meanLifeExpect(Tmatrix2)
	
	plot(Tmatrix@meshpoints, LE, type = "l", xlim = range(Tmatrix1@meshpoints), 
			ylim = range(c(LE, LE1)), xlab = "Sizes", ylab = "Life expectancy",col=cols[1],lty=ltys[1])
	points(Tmatrix1@meshpoints, LE1, type = "l", col = cols[2],lty=ltys[2])
	points(Tmatrix2@meshpoints, LE2, type = "l", col = cols[3],lty=ltys[3])
	legend("topleft", legend = c("Current", "Extended range", "Increased bins"), col = cols, lty = ltys, bty = "n")
	title("Life expectancy")
	
	#mtext("Increasing size range (red) or number of bins (blue) should not alter life expectancy estimates", side=1,outer=TRUE,line=-1)
	
	
	
	if (length(sizesToPlot)==0) sizesToPlot <- as.numeric(quantile(dff$size,c(0.25, 0.5,0.75),na.rm=TRUE))
	
	loctest <- rep(NA,length(sizesToPlot))
	
	print("Please hit any key for the next plot")	
	scan(quiet="TRUE")
	
	par(mfrow = c(2, 3), bty = "l",pty="s")
	
	
	for (kk in c(1,3)) {
		if (kk==1) Tmat <- Tmatrix	
		if (kk==3) Tmat <- Tmatrix2	
		
		h <- diff(Tmat@meshpoints)[1]
		testSizes <- seq(min(Tmat@meshpoints), max(Tmat@meshpoints), 
				length = 5000)
		for (j in 1:3) {
			
			loctest[j] <- which(abs(Tmat@meshpoints-sizesToPlot[j])==min(abs(Tmat@meshpoints-sizesToPlot[j])),arr.ind=TRUE)[1]
			
			#print(loctest[j])
			
			ps <- surv(Tmat@meshpoints[loctest[j]], cov, survObj)
			newd <- data.frame(size = Tmat@meshpoints[loctest[j]], 
					size2 = Tmat@meshpoints[loctest[j]]^2, covariate = Tmat@env.index[1])
			if (length(grep("logsize", names(growObj@fit$coefficients)))) 
				newd$logsize = log(Tmat@meshpoints[loctest[j]])
			if (length(growObj@fit$model$covariate) > 0) 
				if (is.factor(growObj@fit$model$covariate)) 
					newd$covariate <- as.factor(newd$covariate)
			if (length(grep("decline", tolower(as.character(class(growObj))))) > 0 | 
					length(grep("trunc", tolower(as.character(class(growObj))))) > 0) {
				mux <- .predictMuX(growObj, newd)
			} else {
				mux <- predict(growObj@fit, newd, type = "response")
			}
			if (length(grep("incr", tolower(as.character(class(growObj))))) >  0 & 
					length(grep("logincr", tolower(as.character(class(growObj))))) == 0) {
				mux <- Tmat@meshpoints[loctest[j]] + mux
			}
			if (length(grep("decline", tolower(as.character(class(growObj))))) ==  0 &
					length(grep("trunc", tolower(as.character(class(growObj))))) == 0) {
				sigmax2 <- growObj@sd^2
			} else {
				sigmax2 <- growObj@fit$sigmax2
				var.exp.coef <- growObj@fit$var.exp.coef
				sigmax2 <- sigmax2 * exp(2 * (var.exp.coef * mux))
				if (length(grep("trunc", tolower(as.character(class(growObj))))) > 0) 
					sigmax2 <- growObj@fit$sigmax2
			}
			
			#  range.x <- range(Tmat@meshpoints[loctest[j]] + c(-3.5 * sqrt(sigmax2), +3.5 * sqrt(sigmax2)))
			range.x <- range(mux + c(-3.5 * sqrt(sigmax2), +3.5 * sqrt(sigmax2)))
			
			plot(Tmat@meshpoints, Tmat@.Data[, loctest[j]]/h/ps, 
					type = "n", xlim = range.x, xlab = "Size next", ylab = "Kernel")
			
			#if (j == 1 & kk==1) title("Numerical resolution and growth")
			if (j == 2 & kk==1) mtext("Numerical resolution and growth",3,line=4,font=2)
			if (j ==1 & kk==1) title("Current Tmatrix")
			if (j ==1 & kk==3) title("Increased bins")
			
			
			
			for (k in 1:length(Tmat@meshpoints)) {
				points(c(Tmat@meshpoints[k]) + c(-h/2, h/2), rep(Tmat@.Data[k,loctest[j]], 2)/h/ps, type = "l",col=cols[kk])
				points(rep(Tmat@meshpoints[k] + h/2, 2), c(0, Tmat@.Data[k, loctest[j]]/h/ps), type = "l", 
						lty = 1,col=cols[kk])
				points(rep(Tmat@meshpoints[k] - h/2, 2), c(0,Tmat@.Data[k, loctest[j]]/h/ps), type = "l", 
						lty = 1,col=cols[kk])
			}
			if (length(grep("logincr", tolower(as.character(class(growObj))))) == 0 &
					length(grep("trunc", tolower(as.character(class(growObj))))) == 0) {
				points(testSizes, dnorm(testSizes, mux, sqrt(sigmax2)), type = "l", col = 2)
			} else {
				if (length(grep("trunc", tolower(as.character(class(growObj))))) > 0) {
					require(truncnorm)
					points(testSizes, dtruncnorm(testSizes, a = Tmat@meshpoints[loctest[j]], 
									b = Inf, mean = mux, sd = sqrt(sigmax2)), type = "l",col = 2)
				} else {
					points(testSizes, dlnorm(testSizes - Tmat@meshpoints[loctest[j]], 
									mux, sqrt(sigmax2)), type = "l", col = 2)
				}
			}
			
			legend("topright",legend=paste("size=",round(Tmat@meshpoints[loctest[j]],2)), col = "white", 
					lty = 1, bty = "n")
			
		}
	}
	
}






##### Functions to identify sensible numbers of bins - help file on desktop  with data-frame setup
#####


convergeLambda<-function(growObj, survObj, fecObj, nBigMatrix, minSize, maxSize, 
		discreteTrans = 1, integrateType = "midpoint", correction = "none", preCensus = TRUE, tol=1e-4,binIncrease=5){
	
	lambda.new<-1000
	delta<-1000
	while(delta>tol) {
		lambda.old <-lambda.new
		nBigMatrix <- nBigMatrix + binIncrease
		Tmatrix <- createIPMTmatrix(nBigMatrix = nBigMatrix, minSize = minSize, 
				maxSize = maxSize, growObj = growObj, survObj = survObj, 
				discreteTrans = discreteTrans, integrateType = integrateType, 
				correction = correction)
		Fmatrix <- createIPMFmatrix(nBigMatrix = nBigMatrix, minSize = minSize, 
				maxSize = maxSize, fecObj = fecObj, integrateType = integrateType, 
				correction = correction, preCensus = preCensus,survObj=survObj,growObj=growObj)
		
		IPM <- Tmatrix + Fmatrix
		lambda.new <- Re(eigen(IPM)$value[1])
		
		delta<-abs(lambda.new-lambda.old)
		print(delta)
	}
	
	print(c("Final lambda from iteration:",lambda.new))
	print(c("Number of bins:",nBigMatrix))
	
	output<-list(binIncrease=binIncrease,IPM=IPM,lambda=lambda.new)
	
	return(output)
}



convergeR0<-function(growObj, survObj, fecObj, nBigMatrix, minSize, maxSize, 
		discreteTrans = 1, integrateType = "midpoint", correction = "none", preCensus = TRUE, tol=1e-4,binIncrease=5){
	
	R0.new<-1000
	delta<-1000
	while(delta>tol) {
		R0.old <-R0.new
		nBigMatrix <- nBigMatrix + binIncrease
		Tmatrix <- createIPMTmatrix(nBigMatrix = nBigMatrix, minSize = minSize, 
				maxSize = maxSize, growObj = growObj, survObj = survObj, 
				discreteTrans = discreteTrans, integrateType = integrateType, 
				correction = correction)
		Fmatrix <- createIPMFmatrix(nBigMatrix = nBigMatrix, minSize = minSize, 
				maxSize = maxSize, fecObj = fecObj, integrateType = integrateType, 
				correction = correction, preCensus = preCensus,survObj=survObj,growObj=growObj)
		
		R0.new <- R0Calc(Tmatrix,Fmatrix)
		
		delta<-abs(R0.new-R0.old)
		print(delta)
	}
	
	IPM <- Tmatrix + Fmatrix
	
	print(c("Final R0 from iteration:",R0.new))
	print(c("Number of bins:",nBigMatrix))
	
	output<-list(binIncrease=binIncrease,IPM=IPM,R0=R0.new)
	
	return(output)
}




convergeLifeExpectancyFirstBin<-function(growObj, survObj,nBigMatrix, minSize, maxSize, 
		discreteTrans = 1, integrateType = "midpoint", correction = "none", 
		tol=1e-1,binIncrease=5){
	
	LE.new<-1000
	delta<-1000
	while(delta>tol) {
		LE.old <- LE.new
		nBigMatrix <- nBigMatrix + binIncrease
		Tmatrix <- createIPMTmatrix(nBigMatrix = nBigMatrix, minSize = minSize, 
				maxSize = maxSize, growObj = growObj, survObj = survObj, 
				discreteTrans = discreteTrans, integrateType = integrateType, 
				correction = correction)
		
		LE.new <- meanLifeExpect(Tmatrix)    
		
		delta <- abs(LE.new[1]-LE.old[1])
		print(delta)
	}
	
	print(c("Final life expectancy of first bin from iteration:",LE.new[1]))
	print(c("Number of bins:",nBigMatrix))
	
	output<-list(binIncrease=binIncrease,Tmatrix=Tmatrix,LE=LE.new)
	
	return(output)
}


convergeLifeExpectancyLastBin<-function(growObj, survObj,nBigMatrix, minSize, maxSize, 
		discreteTrans = 1, integrateType = "midpoint", correction = "none", 
		tol=1e-1,binIncrease=5){
	
	LE.new<-1000
	delta<-1000
	while(delta>tol) {
		LE.old <- LE.new
		nBigMatrix <- nBigMatrix + binIncrease
		Tmatrix <- createIPMTmatrix(nBigMatrix = nBigMatrix, minSize = minSize, 
				maxSize = maxSize, growObj = growObj, survObj = survObj, 
				discreteTrans = discreteTrans, integrateType = integrateType, 
				correction = correction)
		
		LE.new <- meanLifeExpect(Tmatrix)    
		
		delta <- abs(LE.new[length(LE.new)]-LE.old[length(LE.old)])
		print(delta)
	}
	
	print(c("Final life expectancy of last bin from iteration:",LE.new[length(LE.old)]))
	print(c("Number of bins:",nBigMatrix))
	
	output<-list(binIncrease=binIncrease,Tmatrix=Tmatrix,LE=LE.new)
	
	return(output)
}




### FUNCTIONS FOR EXTRACTING INTEGRATED DEMOGRAPHIC MEASURES ########################
### i.e. life-expectancy and passage time


#Generic for mean life expectancy
#parameters - an IPM
# returns - the life expectancy for every starting size. 
meanLifeExpect <- function(IPMmatrix){
			require(MASS)
			nBigMatrix <- length(IPMmatrix@.Data[1,]) #this nBigMatrix actually contains discrete, env, etc
			#tmp <-  ginv(diag(IPMmatrix@nEnvClass*nBigMatrix)-IPMmatrix)
			tmp <-  ginv(diag(nBigMatrix)-IPMmatrix)
			lifeExpect <- colSums(tmp)
			return(lifeExpect)
		}



#Generic for variance life expectancy (p119 Caswell)
#parameters - an IPM
# returns - the variance in life expectancy for every starting size. 

varLifeExpect <- function(IPMmatrix){
			require(MASS)
			nBigMatrix <- length(IPMmatrix@.Data[1,])
			#tmp <-  ginv(diag(IPMmatrix@nEnvClass*nBigMatrix)-IPMmatrix)
			tmp <-  ginv(diag(nBigMatrix)-IPMmatrix)
			#varLifeExpect <- (2*diag(tmp)-diag(length(IPMmatrix[,1])))%*%tmp-(tmp*tmp)
			#varLifeExpect <- colSums(varLifeExpect)
			varLifeExpect <- colSums(2*(tmp%*%tmp)-tmp)-colSums(tmp)*colSums(tmp)                  
			return(varLifeExpect)
		}



		
		
#Generic for survivorship
#parameters - IPMmatrix - an IPM
#           - size1 - a size at age 1
#           - maxAge - a maxAge
# returns - a list including the survivorship up to the max age,
#                      this broken down by stage,
#                       and mortality over age 
		
		
survivorship <- function(IPMmatrix, loc, maxAge=300){
			nBigMatrix <- length(IPMmatrix@.Data[1,])
			#n <- IPMmatrix@nEnvClass*nBigMatrix
			n <- nBigMatrix
			A1 <- tmp <-  IPMmatrix
			stage.agesurv <- matrix(NA,n,maxAge)
			surv.curv <- rep (NA,maxAge)
			
			#identify the starting size you want to track - removed - specify bin directly
			#loc <- which(abs(size1-IPMmatrix@meshpoints)==min(abs(size1-IPMmatrix@meshpoints)),arr.ind=T)[1]
			popvec <- matrix(0,n,1)
			popvec[floor(loc),1] <- 1
			
			for (a in 1:maxAge) {
				surv.curv[a]<-sum(A1[,loc]); 
				stage.agesurv[c(1:n),a]<-A1[,]%*%popvec
				A1<-A1%*%tmp
			}
			
			mortality <- -log(surv.curv[2:length(surv.curv)]/surv.curv[1:(length(surv.curv)-1)])
			
			return(list(surv.curv=surv.curv,stage.agesurv=stage.agesurv, mortality = mortality))
		}
		
		
		


#Generic for first passage time 
#parameters - an IPM
#           - a size for which passage time is required            
# returns - the passage time to this size from each of the sizes in the IPM 

passageTime <- function(chosenSize,IPMmatrix){
			require(MASS)
			
			loc <- which(abs(chosenSize-IPMmatrix@meshpoints) ==
							min(abs(chosenSize - IPMmatrix@meshpoints)),arr.ind=TRUE)[1]
			matrix.dim <- length(IPMmatrix[1,])
			
			Tprime <- IPMmatrix
			Tprime[,loc] <- 0
			
			Mprime <- 1-colSums(IPMmatrix)
			Mprime[loc]<-0
			Mprime <- rbind(Mprime,rep(0,matrix.dim))
			Mprime[2,loc] <- 1
			
			Bprime <- Mprime%*% ginv(diag(matrix.dim)-Tprime)
			#print(round(Bprime[2,],2))
			#print(sum(Bprime[2,]<1e-6))
			Bprime[2,][Bprime[2,]==0] <- 1
			
			diagBprime <- diag(Bprime[2,])
			#plot(IPMmatrix@meshpoints,diag(diagBprime),type="l",log="y")
			#abline(v=chosenSize)
			Tc <- diagBprime%*%Tprime%*%ginv(diagBprime)
			eta1 <- ginv(diag(matrix.dim)-Tc)             
			
			time.to.absorb <- colSums(eta1)
			time.to.absorb[loc:length(time.to.absorb)] <- 0
			return(time.to.absorb)
		}




#Generic for first variance first passage time (not sure!!!)
#parameters - an IPM
#           - a size for which passage time is required            
# returns - the variance passage time to this size from each of the sizes in the IPM 
varPassageTime <- function(chosenSize,IPMmatrix){
			require(MASS)
			
			loc <- which(abs(chosenSize-IPMmatrix@meshpoints)==min(abs(chosenSize-IPMmatrix@meshpoints)),arr.ind=TRUE)
			matrix.dim <- length(IPMmatrix[1,])
			
			Tprime <- IPMmatrix
			Tprime[,loc] <- 0
			
			Mprime <- 1-colSums(IPMmatrix)
			Mprime[loc]<-0
			Mprime <- rbind(Mprime,rep(0,matrix.dim))
			Mprime[2,loc] <- 1
			
			Bprime <- Mprime%*% solve(diag(matrix.dim)-Tprime)
			
			Tc <- diag(Bprime[2,])%*%Tprime%*%ginv(diag(Bprime[2,]))
			eta1 <- ginv(diag(matrix.dim)-Tc)             
			
			vartimeAbsorb <- colSums(2*(eta1%*%eta1)-eta1)-colSums(eta1)*colSums(eta1)                  
			
			return(vartimeAbsorb)
		}




#Generic for life expectancy in a modeled stoch env  -NEEDS DOUBLE-CHECKING
#parameters - a compound IPM of dim nEnvClass*nBigMatrix
#           - an environmental matrix
# returns - the life expectancy for each of the sizes in the IPM (columns)
#           for each of the starting env states
#
# CURRENTLY NOT WORKING 
#
.stochLifeExpect <- function(IPMmatrix,envMatrix){
			require(MASS)
			
			matrix.dim <- length(IPMmatrix[1,])
			nstages <- IPMmatrix@nBigMatrix
			nstates <- IPMmatrix@nEnvClass
			
			pis <- Re(eigen(envMatrix)$vector[,1])
			pis <- pis/(sum(pis))
			
			#print(pis)
			
			#ckron <- kronecker(envMatrix,diag(nstages))  #doesnt work??
			m <- IPMmatrix  #%*%ckron  #eq 29 in carols paper
			
			Itilda <- diag(matrix.dim)
			
			#not used in this one
			eatildas <- array(dim=c(nstates*nstages,nstages,nstates))
			eatildas[,,] <-0
			for (i in 1:nstates){
				eatildas[,,i] <- Itilda[,((i-1)*nstages+1):(nstages*i)]
			}
			
			qatildas <- array(dim=c(nstates*nstages,nstages,nstates));
			qatildas[,,]<-0
			for (i in 1:nstates) {
				indext <-((i-1)*nstages+1):(nstages*i) 
				qatildas[indext,,i] <- IPMmatrix[indext,indext]/envMatrix[i,i]
				#print( IPMmatrix[cbind(indext,indext)]/envMatrix[i,i] )
			}                            #array of qatildas, eqn 26
			#need to remove env effect since mega-matrix pre-built 
			
			#print(qatildas)
			
			etilda <- array(dim=c(nstates*nstages,nstages));
			etilda[,]<-0
			for (i in 1:nstates) { 
				etilda[((i-1)*nstages+1):(nstages*i),] <- diag(nstages);
			}                            #etilda, eqn 27
			
			I <- diag(nstages);                 #identity matrix of dimension K x K
			
			Ns.markov <- array(dim=c(nstages,nstages,nstates)); #set up for array of Ns
			Ns.markov[,,]<-0
			for (i in 1:nstates){ 
				Ns.markov[,,i] <- I + t(etilda)%*%(solve(Itilda-m))%*%qatildas[,,i];
				#eqn 28, conditional on initial state 1
			}
			
			#print(Ns.markov)
			
			#average over all initial states %%%%%
			Nbar.markov <- rep(0,nstages);
			for (i in 1:nstates) { 
				Nbar.markov <-  Nbar.markov + pis[i] * Ns.markov[,,i]; 
			}                         #eqn 29, weight each fundamntal matrix by expctd frequency 
			
			
			
			#lifeexp, column sums
			lifeexp.markov<-matrix(0,nstates,nstages); #set up array
			for (i in 1:nstates) { 
				lifeexp.markov[i,] <- apply(Ns.markov[,,i], 2, sum);
			}
			
			return(lifeexp.markov)
		}


##Function to estimate Stochastic Passage Time

stochPassageTime <- function(chosenSize,IPMmatrix,envMatrix){
			require(MASS)
			#get the right index for the size you want
			loc <- which(abs(chosenSize-IPMmatrix@meshpoints)==min(abs(chosenSize-
													IPMmatrix@meshpoints)),arr.ind=TRUE)
			#expand out to find that size in every env
			#locs.all <- loc*c(1:IPMmatrix@nEnvClass)
			locs.all <- loc+((IPMmatrix@nBigMatrix)*(0:(envMatrix@nEnvClass-1)))
			matrix.dim <- length(IPMmatrix[1,])
			
			Tprime <- IPMmatrix
			Tprime[,locs.all] <- 0
			
			dhat <- 1-colSums(IPMmatrix)
			dhat[locs.all]<-0
			dhat <- rbind(dhat,rep(0,matrix.dim))
			dhat[2,locs.all] <- 1
			
			bhat <- dhat%*% solve(diag(matrix.dim)-Tprime)
			
			Mc <- diag(bhat[2,])%*%Tprime%*%solve(diag(bhat[2,]))
			eta1 <- solve(diag(matrix.dim)-Mc)             
			
			time.to.absorb <-colSums(eta1)
			
			return(time.to.absorb)
		}



### Short term changes / matrix iteration functions ##############################################


## TO DO - ADJUST TO ALLOW DISCRETE CLASSES ##

## Function to predict distribution in x time-steps given starting
## distribution and IPM (with a single covariate)
#
#
#  Parameters - startingSizes - vector of starting sizes
#             - IPM the IPM (Tmatrix if only intrested in grow surv; Tmatrix + Fmatrix otherwise)
#             - n.time.steps - number of time steps
#             - startingEnv - vector of starting env, same length as startingSizes, or length=1
#
# Returns - a list including starting numbers in each IPM size class (n.new.dist0) and
#                            final numbers in each IPM size class (n.new.dist)
#
#
predictFutureDistribution <- function(startingSizes,IPM, n.time.steps, startingEnv=1) {
	
	# turn starting sizes into the resolution of the IPM bins
	
	# setup slightly different for coompound or non compound dists
	if (IPM@nEnvClass>1) {
		if (length(startingEnv)==1) startingEnv <- rep(startingEnv, length(startingSizes))
		compound <- TRUE
		env.index <- IPM@env.index
		n.new.dist <- rep(0,length(IPM[1,]))
		for (ev in 1:IPM@nEnvClass) { 
			index.new.dist <- findInterval(startingSizes[startingEnv==ev],IPM@meshpoints)+1
			index.new.dist[index.new.dist>length(IPM@meshpoints)] <- length(IPM@meshpoints)
			loc.sizes <- table(index.new.dist); 
			n.new.dist[ev==IPM@env.index][as.numeric(names(loc.sizes))] <- loc.sizes
		}
		n.new.dist0 <- n.new.dist
	} else {
		compound <- FALSE
		index.new.dist <- findInterval(startingSizes,IPM@meshpoints)+1
		index.new.dist[index.new.dist>length(IPM@meshpoints)] <- length(IPM@meshpoints)
		loc.sizes <- table(index.new.dist); 
		env.index <- rep(1,length(IPM@meshpoints))
		n.new.dist <- rep(0,length(IPM@meshpoints))
		n.new.dist[as.numeric(names(loc.sizes))] <- loc.sizes
		n.new.dist0 <- n.new.dist
	}
	
	for (t in 1:n.time.steps) n.new.dist <- IPM@.Data%*%n.new.dist
	
	plot(IPM@meshpoints,n.new.dist0[env.index==1],type="l",xlab="size",
			ylab="n in each size class", ylim=range(c(n.new.dist0,n.new.dist)))
	points(IPM@meshpoints,n.new.dist[env.index==1],type="l",col=2)
	if (compound) {
		for (j in 1:max(env.index)) {
			points(IPM@meshpoints,n.new.dist0[env.index==j],type="l",col=1,lty=j)
			points(IPM@meshpoints,n.new.dist[env.index==j],type="l",col=2,lty=j)
		}
		
	}
	legend("topright",legend=c("current","future"),col=1:2,lty=1,bty="n")
	
	
	return(list(n.new.dist0=n.new.dist0,n.new.dist=n.new.dist))     
}


## TO DO - ADJUST TO ALLOW DISCRETE CLASSES ##

## Function to see how long it takes to get from a starting distribution to an final size
##
#  Parameters - startingSizes - vector of starting sizes (in any order)
#             - IPM the IPM (just a Tmatrix)
#             - endSize - the end size
#             - startingEnv - vector of starting env, same length as startingSizes, or length=1
#             - maxT - the max number of time-steps tested
#             - propReach - the proportion of the starting pop that have to be > than the endSize for it to count
#                  (plots and returned values of survivorship from preliminary runs will give a notion of how low this has to be)
#
# Returns - a list containing: ts.dist - the time-series of size distribution
#                              time.reach - the time for n.reach to be at sizes > endSize
#                              survivorship - survivorship over the course of the time elapsed for that pop
#                              
#
# THE PROBLEM WITH THIS FUNCTION IS THAT EITHER 1) YOU MAKE IT IN ONE TIME-STEP; OR 2) EVERYONE IS DEAD SO TUNING
# propReach BECOMES THE KEY - AND THE EXACT VALUE TO PROVIDE VALUES 1 < x < maxT CAN BE LUDICRIOUSLY SENSITIVE
#
timeToSize <- function(startingSizes,IPM,endSize, startingEnv=1, maxT=100, propReach=0.01) {
	
	cutoff <- which(IPM@meshpoints>endSize,arr.ind=TRUE)[1]
	n.reach <- propReach*length(startingSizes)
	
	# setup slightly different for coompound or non compound dists
	if (IPM@nEnvClass>1) {
		#if startingEnv is not a vector, assume all start in startingEnv
		if (length(startingEnv)==1) startingEnv <- rep(startingEnv, length(startingSizes))
		compound <- TRUE
		env.index <- IPM@env.index
		n.new.dist <- rep(0,length(IPM[1,]))
		for (ev in 1:IPM@nEnvClass) { 
			index.new.dist <- findInterval(startingSizes[startingEnv==ev],IPM@meshpoints)+1
			index.new.dist[index.new.dist>length(IPM@meshpoints)] <- length(IPM@meshpoints)
			loc.sizes <- table(index.new.dist); 
			n.new.dist[ev==IPM@env.index][as.numeric(names(loc.sizes))] <- loc.sizes
		}
		n.new.dist0 <- n.new.dist
	} else {
		compound <- FALSE
		index.new.dist <- findInterval(startingSizes,IPM@meshpoints)+1
		index.new.dist[index.new.dist>length(IPM@meshpoints)] <- length(IPM@meshpoints)
		loc.sizes <- table(index.new.dist); 
		env.index <- rep(1,length(IPM@meshpoints))
		n.new.dist <- rep(0,length(IPM@meshpoints))
		n.new.dist[as.numeric(names(loc.sizes))] <- loc.sizes
		n.new.dist0 <- n.new.dist
	}
	
	ts.dist <- matrix(NA,length(n.new.dist),maxT)
	
	survivorship <- rep(NA,maxT)
	for (t in 1:maxT) { 
		#print(t)
		n.new.dist <- IPM@.Data%*%n.new.dist
		#plot(n.new.dist)
		ts.dist[,t] <- n.new.dist
		
		if (!compound) {
			tot <- sum(n.new.dist[cutoff:length(IPM@meshpoints)])
			survivorship[t] <- sum(n.new.dist)/length(startingSizes)
		} else {
			tot <-sumN <- 0
			for (ev in 1:IPM@nEnvClass) {
				tot <- tot+sum(n.new.dist[env.index==ev][cutoff:length(IPM@meshpoints)])
				sumN <- sumN + sum(n.new.dist[env.index==ev])
			}
			survivorship[t] <- sumN/length(startingSizes)
		}
		
		
		if (tot>n.reach){
			time.reach <- t
			break()
		}
		
	}
	
	if (t==maxT) time.reach <- maxT
	
	par(mfrow=c(1,3),bty="l")
	plot(IPM@meshpoints,n.new.dist0[env.index==1],type="l",xlab="size", ylab="n", ylim=range(c(n.new.dist0+10,n.new.dist)))
	points(IPM@meshpoints,n.new.dist[env.index==1],type="l",col=2)
	legend("topleft",legend=c("starting distribution", "final distribution"),col=c(1,2),lty=1,bty="n")
	abline(v=IPM@meshpoints[cutoff],lty=3)
	
	plot(survivorship[1:t], xlab="Time", ylab="survivorship", type="l")
	
	if (time.reach>5) { 
		image(as.numeric(IPM@meshpoints),1:time.reach,log(ts.dist),ylab="Time steps", xlab="Size classes", main="numbers in size classes over time")
		contour(as.numeric(IPM@meshpoints),1:time.reach,log(ts.dist),add=TRUE,levels=exp(seq(0,max(log(ts.dist)),length=10)))
	}
	print(paste("Time to reach:",time.reach))
	
	return(list(ts.dist=ts.dist, time.reach=time.reach, survivorship=survivorship))     
}








# Calculate R0 
#
# Parameters - Fmatrix, Tmatrix
#
# Returns R0
R0Calc<-function(Tmatrix, Fmatrix){
	require(MASS)
	Imatrix <- matrix(0, length(Tmatrix[1,]), length(Tmatrix[1,])); 
	diag(Imatrix) <- 1
	Nmatrix <- ginv(Imatrix - Tmatrix);
	Rmatrix <- Fmatrix %*% Nmatrix
	ave.R0 <- Re(eigen(Rmatrix)$values[1])
	return(ave.R0)
}

# Calculate lambda and stable stage dist
# for constant env for really huge matrices
# using Matrix package for numerical efficiency
#
# Parameters - Amat - an IPM object
#            - tol - tolerance (i.e. precision required)
#
# Returns list containing lambda and stableStage
#
#ROB WILL MODIFY THIS CODE TO INCLUDE REPRODUCTIVE VALUES
largeMatrixCalc <- function(Tmatrix, Fmatrix, tol = 1.e-8){
	require(Matrix)
	A2 <- Matrix(Tmatrix + Fmatrix);
	nt <- Matrix(1,length(Tmatrix[1,]), 1);
	nt1 <- nt; 
	
	h1 <- diff(Tmatrix@meshpoints)[1]
	
	qmax <- 1000;
	lam <- 1; 
	while(qmax > tol) {
		nt1 <- A2 %*% nt;
		qmax <- sum(abs((nt1 - lam * nt)@x));  
		lam <- sum(nt1@x); 
		nt@x <- (nt1@x) / lam; #slight cheat  
		#cat(lam,qmax,"\n");
	} 
	nt <- matrix(nt@x, length(Tmatrix[1,]), 1); 
	stableDist <- nt / (h1 * sum(nt)); #normalize so that integral=1
	lamStable <- lam; 
	
	# Check works   
	qmax <- sum(abs(lam * nt - (Tmatrix + Fmatrix) %*% nt))
	cat("Convergence: ", qmax, " should be less than ", tol, "\n")
	
	
	return(list(lam = lam, stableDist = stableDist, h1 = h1)) 
	
}





## Sensitivity of parameters - works for an IPM built out of
## growth, survival, and fecundity objects with a fitted regression
## as one of their slots. Note that fecObj and growObj need to reflect
## the 
##
##
sensParams <- function (growObj, survObj, fecObj, nBigMatrix, minSize, maxSize, 
		discreteTrans = 1, integrateType = "midpoint", correction = "none", preCensus = TRUE, 
		delta=1e-4) {
	
	nfec <- 0
	fec.coeff.names <- c()
	for (i in 1:length(fecObj@fitFec)) {
		nfec <- nfec + length(fecObj@fitFec[[i]]$coefficients)
		fec.coeff.names <- c(fec.coeff.names, paste("reprod", 
						i, names(fecObj@fitFec[[i]]$coefficients)))
	}
	npar <- length(growObj@fit$coeff) + 1 + length(survObj@fit$coeff) + 
			length(fecObj@offspringRel$coeff) + 1 + #addedÉ
			(sum(!is.na(fecObj@fecConstants))) + nfec
	elam <- rep(0, npar)
	if ((sum(!is.na(fecObj@fecConstants))) > 0) {
		nmes <- c(paste("grow", names(growObj@fit$coeff)), "sd growth", 
				paste("surv", names(survObj@fit$coeff)), 	
				paste("offspring rel", names(fecObj@offspringRel$coeff)), "sd offspring",	#added
				paste("reprod constant",which(!is.na(fecObj@fecConstants), arr.ind = TRUE)[,2]), 
				fec.coeff.names)
	} else {
		nmes <- c(paste("grow", names(growObj@fit$coeff)), "sd growth", 
				paste("surv", names(survObj@fit$coeff)), 
				paste("offspring rel", names(fecObj@offspringRel$coeff)), "sd offspring",	#added
				fec.coeff.names)
	}
	names(elam) <- nmes[1:npar]
	slam <- elam
	
	
	Tmatrix <- createIPMTmatrix(nBigMatrix = nBigMatrix, minSize = minSize, 
			maxSize = maxSize, growObj = growObj, survObj = survObj, 
			discreteTrans = discreteTrans, integrateType = integrateType, 
			correction = correction)
	Fmatrix <- createIPMFmatrix(nBigMatrix = nBigMatrix, minSize = minSize, 
			maxSize = maxSize, fecObj = fecObj, integrateType = integrateType, 
			correction = correction,preCensus = preCensus,survObj=survObj,growObj=growObj)
	IPM <- Tmatrix + Fmatrix
	lambda1 <- Re(eigen(IPM)$value[1])
	
	#growth
	for (param.test in 1:length(growObj@fit$coeff)) {
		growObj@fit$coefficients[param.test] <- growObj@fit$coefficients[param.test] * 
				(1 + delta)
		Tmatrix <- createIPMTmatrix(nBigMatrix = nBigMatrix, 
				minSize = minSize, maxSize = maxSize, growObj = growObj, 
				survObj = survObj, discreteTrans = discreteTrans, 
				integrateType = integrateType, correction = correction)
		Fmatrix <- createIPMFmatrix(nBigMatrix = nBigMatrix, 
				minSize = minSize, maxSize = maxSize, fecObj = fecObj, 
				integrateType = integrateType, correction = correction,
				preCensus = preCensus,survObj=survObj,growObj=growObj)
		IPM <- Tmatrix + Fmatrix
		lambda2 <- Re(eigen(IPM)$value[1])
		growObj@fit$coefficients[param.test] <- growObj@fit$coefficients[param.test]/(1 + delta)
		slam[param.test] <- (lambda2 - lambda1)/(growObj@fit$coefficients[param.test] * delta)
		elam[param.test] <- (lambda2 - lambda1)/(lambda1*delta)
	}
	
	#growth sd
	param.test <- param.test + 1
	sd.store <- growObj@sd
	growObj@sd <- growObj@sd * (1 + delta)
	Tmatrix <- createIPMTmatrix(nBigMatrix = nBigMatrix, minSize = minSize, 
			maxSize = maxSize, growObj = growObj, survObj = survObj, 
			discreteTrans = discreteTrans, integrateType = integrateType, 
			correction = correction)
	Fmatrix <- createIPMFmatrix(nBigMatrix = nBigMatrix, minSize = minSize, 
			maxSize = maxSize, fecObj = fecObj, integrateType = integrateType, 
			correction = correction,preCensus = preCensus,survObj=survObj,growObj=growObj)
	IPM <- Tmatrix + Fmatrix
	lambda2 <- Re(eigen(IPM)$value[1])
	growObj@sd <- sd.store
	slam[param.test] <- (lambda2 - lambda1)/(growObj@sd*delta)
	elam[param.test] <- (lambda2 - lambda1)/(lambda1*delta)
	
	#survival
	count <- param.test
	for (param.test in 1:length(survObj@fit$coeff)) {
		survObj@fit$coefficients[param.test] <- survObj@fit$coefficients[param.test] * 
				(1 + delta)
		Tmatrix <- createIPMTmatrix(nBigMatrix = nBigMatrix, 
				minSize = minSize, maxSize = maxSize, growObj = growObj, 
				survObj = survObj, discreteTrans = discreteTrans, 
				integrateType = integrateType, correction = correction)
		Fmatrix <- createIPMFmatrix(nBigMatrix = nBigMatrix, 
				minSize = minSize, maxSize = maxSize, fecObj = fecObj, 
				integrateType = integrateType, correction = correction,
				preCensus = preCensus,survObj=survObj,growObj=growObj)
		IPM <- Tmatrix + Fmatrix
		lambda2 <- Re(eigen(IPM)$value[1])
		survObj@fit$coefficients[param.test] <- survObj@fit$coefficients[param.test]/(1 + 
					delta)
		slam[param.test + count] <- (lambda2 - lambda1)/(survObj@fit$coefficients[param.test] * delta)
		elam[param.test + count] <- (lambda2 - lambda1)/(lambda1*delta)
	}
	
	#offspring rel
	count <- count + param.test
	for (param.test in 1:length(fecObj@offspringRel$coeff)) {
		fecObj@offspringRel$coefficients[param.test] <- fecObj@offspringRel$coefficients[param.test] * 
				(1 + delta)
		Tmatrix <- createIPMTmatrix(nBigMatrix = nBigMatrix, 
				minSize = minSize, maxSize = maxSize, growObj = growObj, 
				survObj = survObj, discreteTrans = discreteTrans, 
				integrateType = integrateType, correction = correction)
		Fmatrix <- createIPMFmatrix(nBigMatrix = nBigMatrix, 
				minSize = minSize, maxSize = maxSize, fecObj = fecObj, 
				integrateType = integrateType, correction = correction,
				preCensus = preCensus,survObj=survObj,growObj=growObj)
		IPM <- Tmatrix + Fmatrix
		lambda2 <- Re(eigen(IPM)$value[1])
		fecObj@offspringRel$coefficients[param.test] <- fecObj@offspringRel$coefficients[param.test] / 
				(1 + delta)
		slam[param.test + count] <- (lambda2 - lambda1)/(fecObj@offspringRel$coefficients[param.test] * delta)
		elam[param.test + count] <- (lambda2 - lambda1)/(lambda1*delta)
	}
	
	#offspring sd
	count <- count + param.test
	fecObj@sdOffspringSize <- fecObj@sdOffspringSize *  (1 + delta)
	Tmatrix <- createIPMTmatrix(nBigMatrix = nBigMatrix, 
			minSize = minSize, maxSize = maxSize, growObj = growObj, 
			survObj = survObj, discreteTrans = discreteTrans, 
			integrateType = integrateType, correction = correction)
	Fmatrix <- createIPMFmatrix(nBigMatrix = nBigMatrix, 
			minSize = minSize, maxSize = maxSize, fecObj = fecObj, 
			integrateType = integrateType, correction = correction,
			preCensus = preCensus,survObj=survObj,growObj=growObj)
	IPM <- Tmatrix + Fmatrix
	lambda2 <- Re(eigen(IPM)$value[1])
	fecObj@sdOffspringSize <- fecObj@sdOffspringSize / (1 + delta)
	slam[count + 1] <- (lambda2 - lambda1)/(fecObj@sdOffspringSize * delta)
	elam[count + 1] <- (lambda2 - lambda1)/(lambda1*delta)
	
	
	
	
	#fecConstant	
	chs <- which(!is.na(fecObj@fecConstants), arr.ind = TRUE)[,2]
	if (length(chs) > 0) {
		count <- count + 1
		for (param.test in 1:length(chs)) {
			fecObj@fecConstants[chs[param.test]] <- fecObj@fecConstants[chs[param.test]] * 
					(1 + delta)
			Tmatrix <- createIPMTmatrix(nBigMatrix = nBigMatrix, 
					minSize = minSize, maxSize = maxSize, growObj = growObj, 
					survObj = survObj, discreteTrans = discreteTrans, 
					integrateType = integrateType, correction = correction)
			Fmatrix <- createIPMFmatrix(nBigMatrix = nBigMatrix, 
					minSize = minSize, maxSize = maxSize, fecObj = fecObj, 
					integrateType = integrateType, correction = correction,
					preCensus = preCensus,survObj=survObj,growObj=growObj)
			IPM <- Tmatrix + Fmatrix
			lambda2 <- Re(eigen(IPM)$value[1])
			fecObj@fecConstants[1, chs[param.test]] <- fecObj@fecConstants[1, 
					chs[param.test]]/(1 + delta)
			slam[param.test + count] <- (lambda2 - lambda1)/as.numeric(fecObj@fecConstants[1, 
							chs[param.test]] * delta)
			elam[param.test + count] <- (lambda2 - lambda1)/(lambda1*delta)
		}
		count <- param.test + count
	} else { count <- count+1}
	for (i in 1:length(fecObj@fitFec)) {
		for (param.test in 1:length(fecObj@fitFec[[i]]$coefficients)) {
			fecObj@fitFec[[i]]$coefficients[param.test] <- fecObj@fitFec[[i]]$coefficients[param.test] * 
					(1 + delta)
			Fmatrix1 <- createIPMFmatrix(nBigMatrix = nBigMatrix, 
					minSize = minSize, maxSize = maxSize, fecObj = fecObj, 
					integrateType = integrateType, correction = correction,
					preCensus = preCensus,survObj=survObj,growObj=growObj)
			IPM <- Tmatrix + Fmatrix1
			lambda2 <- Re(eigen(IPM)$value[1])
			fecObj@fitFec[[i]]$coefficients[param.test] <- fecObj@fitFec[[i]]$coefficients[param.test]/(1 + 
						delta)
			slam[param.test + count] <- (lambda2 - lambda1)/(fecObj@fitFec[[i]]$coefficients[param.test] * 
						delta)
			elam[param.test + count] <- (lambda2 - lambda1)/(lambda1*delta)
		}
		count <- count + param.test
	}
	return(list(slam = slam, elam = elam))
}





## identical for R0 ########################################################################

sensParamsR0 <- function (growObj, survObj, fecObj, nBigMatrix, minSize, maxSize, 
		discreteTrans = 1, integrateType = "midpoint", correction = "none", preCensus = TRUE, 
		delta=1e-4) {
	
	nfec <- 0
	fec.coeff.names <- c()
	for (i in 1:length(fecObj@fitFec)) {
		nfec <- nfec + length(fecObj@fitFec[[i]]$coefficients)
		fec.coeff.names <- c(fec.coeff.names, paste("reprod", 
						i, names(fecObj@fitFec[[i]]$coefficients)))
	}
	npar <- length(growObj@fit$coeff) + 1 + length(survObj@fit$coeff) + 
			length(fecObj@offspringRel$coeff) + 1 + #addedÉ
			(sum(!is.na(fecObj@fecConstants))) + nfec
	elam <- rep(0, npar)
	if ((sum(!is.na(fecObj@fecConstants))) > 0) {
		nmes <- c(paste("grow", names(growObj@fit$coeff)), "sd growth", 
				paste("surv", names(survObj@fit$coeff)), 	
				paste("offspring rel", names(fecObj@offspringRel$coeff)), "sd offspring",	#added
				paste("reprod constant",which(!is.na(fecObj@fecConstants), arr.ind = TRUE)[,2]), 
				fec.coeff.names)
	} else {
		nmes <- c(paste("grow", names(growObj@fit$coeff)), "sd growth", 
				paste("surv", names(survObj@fit$coeff)), 
				paste("offspring rel", names(fecObj@offspringRel$coeff)), "sd offspring",	#added
				fec.coeff.names)
	}
	names(elam) <- nmes[1:npar]
	slam <- elam
	
	
	Tmatrix <- createIPMTmatrix(nBigMatrix = nBigMatrix, minSize = minSize, 
			maxSize = maxSize, growObj = growObj, survObj = survObj, 
			discreteTrans = discreteTrans, integrateType = integrateType, 
			correction = correction)
	Fmatrix <- createIPMFmatrix(nBigMatrix = nBigMatrix, minSize = minSize, 
			maxSize = maxSize, fecObj = fecObj, integrateType = integrateType, 
			correction = correction,preCensus = preCensus,survObj=survObj,growObj=growObj)
	
	R01 <- R0Calc(Tmatrix,Fmatrix)
	
	#growth
	for (param.test in 1:length(growObj@fit$coeff)) {
		growObj@fit$coefficients[param.test] <- growObj@fit$coefficients[param.test] * 
				(1 + delta)
		Tmatrix <- createIPMTmatrix(nBigMatrix = nBigMatrix, 
				minSize = minSize, maxSize = maxSize, growObj = growObj, 
				survObj = survObj, discreteTrans = discreteTrans, 
				integrateType = integrateType, correction = correction)
		Fmatrix <- createIPMFmatrix(nBigMatrix = nBigMatrix, 
				minSize = minSize, maxSize = maxSize, fecObj = fecObj, 
				integrateType = integrateType, correction = correction,
				preCensus = preCensus,survObj=survObj,growObj=growObj)
		R02 <- R0Calc(Tmatrix,Fmatrix)
		growObj@fit$coefficients[param.test] <- growObj@fit$coefficients[param.test]/(1 + delta)
		slam[param.test] <- (R02 - R01)/(growObj@fit$coefficients[param.test] * delta)
		elam[param.test] <- (R02 - R01)/(R01*delta)
	}
	
	#growth sd
	param.test <- param.test + 1
	sd.store <- growObj@sd
	growObj@sd <- growObj@sd * (1 + delta)
	Tmatrix <- createIPMTmatrix(nBigMatrix = nBigMatrix, minSize = minSize, 
			maxSize = maxSize, growObj = growObj, survObj = survObj, 
			discreteTrans = discreteTrans, integrateType = integrateType, 
			correction = correction)
	Fmatrix <- createIPMFmatrix(nBigMatrix = nBigMatrix, minSize = minSize, 
			maxSize = maxSize, fecObj = fecObj, integrateType = integrateType, 
			correction = correction,preCensus = preCensus,survObj=survObj,growObj=growObj)
	R02 <- R0Calc(Tmatrix,Fmatrix)
	growObj@sd <- sd.store
	slam[param.test] <- (R02 - R01)/(growObj@sd *delta)
	elam[param.test] <- (R02 - R01)/(R01*delta)
	
	#survival
	count <- param.test
	for (param.test in 1:length(survObj@fit$coeff)) {
		survObj@fit$coefficients[param.test] <- survObj@fit$coefficients[param.test] * 
				(1 + delta)
		Tmatrix <- createIPMTmatrix(nBigMatrix = nBigMatrix, 
				minSize = minSize, maxSize = maxSize, growObj = growObj, 
				survObj = survObj, discreteTrans = discreteTrans, 
				integrateType = integrateType, correction = correction)
		Fmatrix <- createIPMFmatrix(nBigMatrix = nBigMatrix, 
				minSize = minSize, maxSize = maxSize, fecObj = fecObj, 
				integrateType = integrateType, correction = correction,
				preCensus = preCensus,survObj=survObj,growObj=growObj)
		R02 <- R0Calc(Tmatrix,Fmatrix)
		survObj@fit$coefficients[param.test] <- survObj@fit$coefficients[param.test]/(1 + 
					delta)
		slam[param.test + count] <- (R02 - R01)/(survObj@fit$coefficients[param.test] * 
					delta)
		elam[param.test + count] <- (R02 - R01)/(R01*delta)
	}
	
	#offspring rel
	count <- count + param.test
	for (param.test in 1:length(fecObj@offspringRel$coeff)) {
		fecObj@offspringRel$coefficients[param.test] <- fecObj@offspringRel$coefficients[param.test] * 
				(1 + delta)
		Tmatrix <- createIPMTmatrix(nBigMatrix = nBigMatrix, 
				minSize = minSize, maxSize = maxSize, growObj = growObj, 
				survObj = survObj, discreteTrans = discreteTrans, 
				integrateType = integrateType, correction = correction)
		Fmatrix <- createIPMFmatrix(nBigMatrix = nBigMatrix, 
				minSize = minSize, maxSize = maxSize, fecObj = fecObj, 
				integrateType = integrateType, correction = correction,
				preCensus = preCensus,survObj=survObj,growObj=growObj)
		R02 <- R0Calc(Tmatrix,Fmatrix)
		fecObj@offspringRel$coefficients[param.test] <- fecObj@offspringRel$coefficients[param.test] / 
				(1 + delta)
		slam[param.test + count] <- (R02 - R01)/(fecObj@offspringRel$coefficients[param.test] * delta)
		elam[param.test + count] <- (R02 - R01)/(R01*delta)
	}
	
	#offspring sd
	count <- count + param.test
	fecObj@sdOffspringSize <- fecObj@sdOffspringSize *  (1 + delta)
	Tmatrix <- createIPMTmatrix(nBigMatrix = nBigMatrix, 
			minSize = minSize, maxSize = maxSize, growObj = growObj, 
			survObj = survObj, discreteTrans = discreteTrans, 
			integrateType = integrateType, correction = correction)
	Fmatrix <- createIPMFmatrix(nBigMatrix = nBigMatrix, 
			minSize = minSize, maxSize = maxSize, fecObj = fecObj, 
			integrateType = integrateType, correction = correction,
			preCensus = preCensus,survObj=survObj,growObj=growObj)
	R02 <- R0Calc(Tmatrix,Fmatrix)
	fecObj@sdOffspringSize <- fecObj@sdOffspringSize / (1 + delta)
	slam[count + 1] <- (R02 - R01)/(fecObj@sdOffspringSize * delta)
	elam[count + 1] <- (R02 - R01)/(R01*delta)
	
	
	
	
	#fecConstant	
	chs <- which(!is.na(fecObj@fecConstants), arr.ind = TRUE)[,2]
	if (length(chs) > 0) {
		count <- count + 1
		for (param.test in 1:length(chs)) {
			fecObj@fecConstants[chs[param.test]] <- fecObj@fecConstants[chs[param.test]] * 
					(1 + delta)
			Tmatrix <- createIPMTmatrix(nBigMatrix = nBigMatrix, 
					minSize = minSize, maxSize = maxSize, growObj = growObj, 
					survObj = survObj, discreteTrans = discreteTrans, 
					integrateType = integrateType, correction = correction)
			Fmatrix <- createIPMFmatrix(nBigMatrix = nBigMatrix, 
					minSize = minSize, maxSize = maxSize, fecObj = fecObj, 
					integrateType = integrateType, correction = correction,
					preCensus = preCensus,survObj=survObj,growObj=growObj)
			R02 <- R0Calc(Tmatrix,Fmatrix)
			fecObj@fecConstants[1, chs[param.test]] <- fecObj@fecConstants[1, 
					chs[param.test]]/(1 + delta)
			slam[param.test + count] <- (R02 - R01)/as.numeric(fecObj@fecConstants[1, 
							chs[param.test]] * delta)
			elam[param.test + count] <- (R02 - R01)/(R01*delta)
		}
		count <- param.test + count
	} else { count <- count+1}
	for (i in 1:length(fecObj@fitFec)) {
		for (param.test in 1:length(fecObj@fitFec[[i]]$coefficients)) {
			fecObj@fitFec[[i]]$coefficients[param.test] <- fecObj@fitFec[[i]]$coefficients[param.test] * 
					(1 + delta)
			Fmatrix1 <- createIPMFmatrix(nBigMatrix = nBigMatrix, 
					minSize = minSize, maxSize = maxSize, fecObj = fecObj, 
					integrateType = integrateType, correction = correction,
					preCensus = preCensus,survObj=survObj,growObj=growObj)
			R02 <- R0Calc(Tmatrix,Fmatrix)
			fecObj@fitFec[[i]]$coefficients[param.test] <- fecObj@fitFec[[i]]$coefficients[param.test]/(1 + delta)
			slam[param.test + count] <- (R02 - R01)/(fecObj@fitFec[[i]]$coefficients[param.test] * delta)
			elam[param.test + count] <- (R02 - R01)/(R01*delta)
		}
		count <- count + param.test
	}
	return(list(slam = slam, elam = elam))
}





## identical for LE ########################################################################

sensParamsLifeExpect <- function (growObj, survObj, nBigMatrix, minSize, maxSize, 
		discreteTrans = 1, integrateType = "midpoint", correction = "none", preCensus = TRUE, 
		delta=1e-4,chosenBin=1) {
	
	npar <- length(growObj@fit$coeff) + 1 + length(survObj@fit$coeff)
	elam <- rep(0, npar)
	nmes <- c(paste("grow", names(growObj@fit$coeff)), "sd growth", 
			paste("surv", names(survObj@fit$coeff)))
	
	names(elam) <- nmes[1:npar]
	slam <- elam
	
	Tmatrix <- createIPMTmatrix(nBigMatrix = nBigMatrix, minSize = minSize, 
			maxSize = maxSize, growObj = growObj, survObj = survObj, 
			discreteTrans = discreteTrans, integrateType = integrateType, 
			correction = correction)
	
	LE1 <- meanLifeExpect(Tmatrix)[chosenBin]
	
	#growth
	for (param.test in 1:length(growObj@fit$coeff)) {
		growObj@fit$coefficients[param.test] <- growObj@fit$coefficients[param.test] * 
				(1 + delta)
		Tmatrix <- createIPMTmatrix(nBigMatrix = nBigMatrix, 
				minSize = minSize, maxSize = maxSize, growObj = growObj, 
				survObj = survObj, discreteTrans = discreteTrans, 
				integrateType = integrateType, correction = correction)
		LE2 <- meanLifeExpect(Tmatrix)[chosenBin]
		growObj@fit$coefficients[param.test] <- growObj@fit$coefficients[param.test]/(1 + delta)
		slam[param.test] <- (LE2 - LE1)/(growObj@fit$coefficients[param.test] * delta)
		elam[param.test] <- (LE2 - LE1)/(LE1*delta)
	}
	
	#growth sd
	param.test <- param.test + 1
	sd.store <- growObj@sd
	growObj@sd <- growObj@sd * (1 + delta)
	Tmatrix <- createIPMTmatrix(nBigMatrix = nBigMatrix, minSize = minSize, 
			maxSize = maxSize, growObj = growObj, survObj = survObj, 
			discreteTrans = discreteTrans, integrateType = integrateType, 
			correction = correction)
	LE2 <- meanLifeExpect(Tmatrix)[chosenBin]
	growObj@sd <- sd.store
	slam[param.test] <- (LE2 - LE1)/(growObj@sd*delta)
	elam[param.test] <- (LE2 - LE1)/(LE1*delta)
	
	#survival
	count <- param.test
	for (param.test in 1:length(survObj@fit$coeff)) {
		survObj@fit$coefficients[param.test] <- survObj@fit$coefficients[param.test] * 
				(1 + delta)
		Tmatrix <- createIPMTmatrix(nBigMatrix = nBigMatrix, 
				minSize = minSize, maxSize = maxSize, growObj = growObj, 
				survObj = survObj, discreteTrans = discreteTrans, 
				integrateType = integrateType, correction = correction)
		LE2 <- meanLifeExpect(Tmatrix)[chosenBin]
		survObj@fit$coefficients[param.test] <- survObj@fit$coefficients[param.test]/(1 + delta)
		slam[param.test + count] <- (LE2 - LE1)/(survObj@fit$coefficients[param.test] * delta)
		elam[param.test + count] <- (LE2 - LE1)/(LE1*delta)
	}
	return(list(slam = slam, elam = elam))
}


### sens params for discrete survival bit

sensParamsDiscrete <-  function (growObj, survObj, fecObj, nBigMatrix, minSize, maxSize, 
		discreteTrans = 1, integrateType = "midpoint", correction = "none", preCensus = TRUE, delta=1e-4) {
	
	
	#all the transitions - note that this is in the order of the columns, 
	nmes <- paste("",c(outer(colnames(discreteTrans@discreteTrans),colnames(discreteTrans@discreteTrans),paste,sep="-")),sep="")
	
	# all survival out of discrete stages
	nmes <- c(nmes,paste("survival from",colnames(discreteTrans@discreteSurv),sep=""))
	
	# all means out of discrete stages
	nmes <- c(nmes,paste("mean from ",colnames(discreteTrans@meanToCont),sep=""))
	
	# all sds out of discrete stages
	nmes <- c(nmes,paste("sd from ",colnames(discreteTrans@sdToCont),sep=""))
	
	#  not sure makes sense to do distribToDiscrete???	
	
	#coefficients linking continuous survival into discrete
	nmes <- c(nmes, paste("survival from continuous ", names(discreteTrans@survToDiscrete$coefficients),sep=""))
	
	
	elam <- rep(NA,length(nmes))
	names(elam) <- nmes
	slam <- elam
	Tmatrix <- createIPMTmatrix(nBigMatrix = nBigMatrix, minSize = minSize, 
			maxSize = maxSize, growObj = growObj, survObj = survObj, 
			discreteTrans = discreteTrans, integrateType = integrateType, 
			correction = correction)
	Fmatrix <- createIPMFmatrix(nBigMatrix = nBigMatrix, minSize = minSize, 
			maxSize = maxSize, fecObj = fecObj, integrateType = integrateType, 
			correction = correction,preCensus = preCensus,survObj=survObj,growObj=growObj)
	IPM <- Tmatrix + Fmatrix
	lambda1 <- Re(eigen(IPM)$value[1])
	
	#all the transitions
	param.test <- 0
	for (j in 1:nrow(discreteTrans@discreteTrans)) {
		for (k in 1:nrow(discreteTrans@discreteTrans)) {
			
			#print(c(rownames(discreteTrans@discreteTrans)[k],colnames(discreteTrans@discreteTrans)[j]))
			
			param.test <- param.test+1
			
			if (discreteTrans@discreteTrans[k,j]==0) next()
			
			#pick element of matrix
			adj <- rep(discreteTrans@discreteTrans[k,j]*delta,nrow(discreteTrans@discreteTrans)-1)
			discreteTrans@discreteTrans[k,j] <- discreteTrans@discreteTrans[k,j] * (1 + delta)
			#alter the other values in the columns so as continue to sum to oneÉ
			if (sum(discreteTrans@discreteTrans[k,-j]>0)>0) adj <- adj/sum(discreteTrans@discreteTrans[k,-j]>0)
			adj[discreteTrans@discreteTrans[-k,j]==0] <- 0
			discreteTrans@discreteTrans[-k,j] <- discreteTrans@discreteTrans[-k,j]-adj
			
			#print(colSums(discreteTrans@discreteTrans))
			
			
			Tmatrix <- createIPMTmatrix(nBigMatrix = nBigMatrix, 
					minSize = minSize, maxSize = maxSize, growObj = growObj, 
					survObj = survObj, discreteTrans = discreteTrans, 
					integrateType = integrateType, correction = correction)
			
			Fmatrix <- createIPMFmatrix(nBigMatrix = nBigMatrix, 
					minSize = minSize, maxSize = maxSize, fecObj = fecObj, 
					integrateType = integrateType, correction = correction,
					preCensus = preCensus,survObj=survObj,growObj=growObj)
			
			IPM <- Tmatrix + Fmatrix
			lambda2 <- Re(eigen(IPM)$value[1])
			discreteTrans@discreteTrans[k,j] <-  discreteTrans@discreteTrans[k,j]/(1 + delta)	
			discreteTrans@discreteTrans[-k,j] <- discreteTrans@discreteTrans[-k,j]+adj
			
			slam[param.test] <- (lambda2 - lambda1)/(discreteTrans@discreteTrans[k,j] * delta)
			elam[param.test] <- (lambda2 - lambda1)/(lambda1*delta)
		}}
	
	
	# all survival out of discrete stages
	count <- param.test
	for (param.test in 1:length(discreteTrans@discreteSurv)) { 
		
		discreteTrans@discreteSurv[param.test] <- discreteTrans@discreteSurv[param.test]* (1 + delta)
		
		Tmatrix <- createIPMTmatrix(nBigMatrix = nBigMatrix, 
				minSize = minSize, maxSize = maxSize, growObj = growObj, 
				survObj = survObj, discreteTrans = discreteTrans, 
				integrateType = integrateType, correction = correction)
		
		Fmatrix <- createIPMFmatrix(nBigMatrix = nBigMatrix, 
				minSize = minSize, maxSize = maxSize, fecObj = fecObj, 
				integrateType = integrateType, correction = correction,
				preCensus = preCensus,survObj=survObj,growObj=growObj)
		
		IPM <- Tmatrix + Fmatrix
		lambda2 <- Re(eigen(IPM)$value[1])
		
		discreteTrans@discreteSurv[param.test] <- discreteTrans@discreteSurv[param.test]/ (1 + delta)
		slam[param.test+count] <- (lambda2 - lambda1)/(discreteTrans@discreteSurv[param.test] * delta)
		elam[param.test+count] <- (lambda2 - lambda1)/(lambda1*delta)
		
	}
	
	
	# all mean values coming out of discrete stages
	count <- param.test+count
	for (param.test in 1:length(discreteTrans@meanToCont)) { 
		
		discreteTrans@meanToCont[param.test] <- discreteTrans@meanToCont[param.test]* (1 + delta)
		
		Tmatrix <- createIPMTmatrix(nBigMatrix = nBigMatrix, 
				minSize = minSize, maxSize = maxSize, growObj = growObj, 
				survObj = survObj, discreteTrans = discreteTrans, 
				integrateType = integrateType, correction = correction)
		
		Fmatrix <- createIPMFmatrix(nBigMatrix = nBigMatrix, 
				minSize = minSize, maxSize = maxSize, fecObj = fecObj, 
				integrateType = integrateType, correction = correction,
				preCensus = preCensus,survObj=survObj,growObj=growObj)
		
		IPM <- Tmatrix + Fmatrix
		lambda2 <- Re(eigen(IPM)$value[1])
		
		discreteTrans@meanToCont[param.test] <- discreteTrans@meanToCont[param.test]/ (1 + delta)
		slam[param.test+count] <- (lambda2 - lambda1)/(discreteTrans@meanToCont[param.test] * delta)
		elam[param.test+count] <- (lambda2 - lambda1)/(lambda1*delta)
		
	}
	
	
	# all sd values coming out of discrete stages
	count <- param.test+count
	for (param.test in 1:length(discreteTrans@sdToCont)) { 
		
		discreteTrans@sdToCont[param.test] <- discreteTrans@sdToCont[param.test]* (1 + delta)
		
		Tmatrix <- createIPMTmatrix(nBigMatrix = nBigMatrix, 
				minSize = minSize, maxSize = maxSize, growObj = growObj, 
				survObj = survObj, discreteTrans = discreteTrans, 
				integrateType = integrateType, correction = correction)
		
		Fmatrix <- createIPMFmatrix(nBigMatrix = nBigMatrix, 
				minSize = minSize, maxSize = maxSize, fecObj = fecObj, 
				integrateType = integrateType, correction = correction,
				preCensus = preCensus,survObj=survObj,growObj=growObj)
		
		IPM <- Tmatrix + Fmatrix
		lambda2 <- Re(eigen(IPM)$value[1])
		
		discreteTrans@sdToCont[param.test] <- discreteTrans@sdToCont[param.test]/ (1 + delta)
		slam[param.test+count] <- (lambda2 - lambda1)/(discreteTrans@sdToCont[param.test][param.test] * delta)
		elam[param.test+count] <- (lambda2 - lambda1)/(lambda1*delta)
		
	}
	
	# parameters linking size to survival
	count <- param.test+count
	for (param.test in 1:length(discreteTrans@survToDiscrete$coefficients)) { 
		
		discreteTrans@survToDiscrete$coefficients[param.test] <- discreteTrans@survToDiscrete$coefficients[param.test]* (1 + delta)
		
		Tmatrix <- createIPMTmatrix(nBigMatrix = nBigMatrix, 
				minSize = minSize, maxSize = maxSize, growObj = growObj, 
				survObj = survObj, discreteTrans = discreteTrans, 
				integrateType = integrateType, correction = correction)
		
		Fmatrix <- createIPMFmatrix(nBigMatrix = nBigMatrix, 
				minSize = minSize, maxSize = maxSize, fecObj = fecObj, 
				integrateType = integrateType, correction = correction,
				preCensus = preCensus,survObj=survObj,growObj=growObj)
		
		IPM <- Tmatrix + Fmatrix
		lambda2 <- Re(eigen(IPM)$value[1])
		
		discreteTrans@survToDiscrete$coefficients[param.test] <- discreteTrans@survToDiscrete$coefficients[param.test]/ (1 + delta)
		slam[param.test+count] <- (lambda2 - lambda1)/(discreteTrans@survToDiscrete$coefficient[param.test] * delta)
		elam[param.test+count] <- (lambda2 - lambda1)/(lambda1*delta)
		
	}
	
	print("Did not calculate sensitivities and elasticities for:")
	print(c(names(slam[is.na(slam)])))
	print("Values of zero")
	
	slam <- slam[!is.na(slam)]
	elam <- elam[!is.na(slam)]
	
	return(list(slam = slam, elam = elam))
	
}




sensParamsDiscreteR0 <-  function (growObj, survObj, fecObj, nBigMatrix, minSize, maxSize, 
		discreteTrans = 1, integrateType = "midpoint", correction = "none", preCensus = TRUE, delta=1e-4) {
	
	
	#all the transitions - note that this is in the order of the columns, 
	nmes <- paste("",c(outer(colnames(discreteTrans@discreteTrans),colnames(discreteTrans@discreteTrans),paste,sep="-")),sep="")
	
	# all survival out of discrete stages
	nmes <- c(nmes,paste("survival from",colnames(discreteTrans@discreteSurv),sep=""))
	
	# all means out of discrete stages
	nmes <- c(nmes,paste("mean from ",colnames(discreteTrans@meanToCont),sep=""))
	
	# all sds out of discrete stages
	nmes <- c(nmes,paste("sd from ",colnames(discreteTrans@sdToCont),sep=""))
	
	#  not sure makes sense to do distribToDiscrete???	
	
	#coefficients linking continuous survival into discrete
	nmes <- c(nmes, paste("survival from continuous ", names(discreteTrans@survToDiscrete$coefficients),sep=""))
	
	
	elam <- rep(NA,length(nmes))
	names(elam) <- nmes
	slam <- elam
	Tmatrix <- createIPMTmatrix(nBigMatrix = nBigMatrix, minSize = minSize, 
			maxSize = maxSize, growObj = growObj, survObj = survObj, 
			discreteTrans = discreteTrans, integrateType = integrateType, 
			correction = correction)
	Fmatrix <- createIPMFmatrix(nBigMatrix = nBigMatrix, minSize = minSize, 
			maxSize = maxSize, fecObj = fecObj, integrateType = integrateType, 
			correction = correction,preCensus = preCensus,survObj=survObj,growObj=growObj)
	R01 <- R0Calc(Tmatrix,Fmatrix)
	
	#all the transitions
	param.test <- 0
	for (j in 1:nrow(discreteTrans@discreteTrans)) {
		for (k in 1:nrow(discreteTrans@discreteTrans)) {
			
			#print(c(rownames(discreteTrans@discreteTrans)[k],colnames(discreteTrans@discreteTrans)[j]))
			
			param.test <- param.test+1
			
			if (discreteTrans@discreteTrans[k,j]==0) next()
			
			#pick element of matrix
			adj <- rep(discreteTrans@discreteTrans[k,j]*delta,nrow(discreteTrans@discreteTrans)-1)
			discreteTrans@discreteTrans[k,j] <- discreteTrans@discreteTrans[k,j] * (1 + delta)
			#alter the other values in the columns so as continue to sum to oneÉ
			if (sum(discreteTrans@discreteTrans[k,-j]>0)>0) adj <- adj/sum(discreteTrans@discreteTrans[k,-j]>0)
			adj[discreteTrans@discreteTrans[-k,j]==0] <- 0
			discreteTrans@discreteTrans[-k,j] <- discreteTrans@discreteTrans[-k,j]-adj
			
			#print(colSums(discreteTrans@discreteTrans))
			
			
			Tmatrix <- createIPMTmatrix(nBigMatrix = nBigMatrix, 
					minSize = minSize, maxSize = maxSize, growObj = growObj, 
					survObj = survObj, discreteTrans = discreteTrans, 
					integrateType = integrateType, correction = correction)
			
			Fmatrix <- createIPMFmatrix(nBigMatrix = nBigMatrix, 
					minSize = minSize, maxSize = maxSize, fecObj = fecObj, 
					integrateType = integrateType, correction = correction,
					preCensus = preCensus,survObj=survObj,growObj=growObj)
			
			R02 <- R0Calc(Tmatrix,Fmatrix)
			discreteTrans@discreteTrans[k,j] <-  discreteTrans@discreteTrans[k,j]/(1 + delta)	
			discreteTrans@discreteTrans[-k,j] <- discreteTrans@discreteTrans[-k,j]+adj
			
			slam[param.test] <- (R02 - R01)/(discreteTrans@discreteTrans[k,j] * delta)
			elam[param.test] <- (R02 - R01)/(R01*delta)
		}}
	
	
	# all survival out of discrete stages
	count <- param.test
	for (param.test in 1:length(discreteTrans@discreteSurv)) { 
		
		discreteTrans@discreteSurv[param.test] <- discreteTrans@discreteSurv[param.test]* (1 + delta)
		
		Tmatrix <- createIPMTmatrix(nBigMatrix = nBigMatrix, 
				minSize = minSize, maxSize = maxSize, growObj = growObj, 
				survObj = survObj, discreteTrans = discreteTrans, 
				integrateType = integrateType, correction = correction)
		
		Fmatrix <- createIPMFmatrix(nBigMatrix = nBigMatrix, 
				minSize = minSize, maxSize = maxSize, fecObj = fecObj, 
				integrateType = integrateType, correction = correction,
				preCensus = preCensus,survObj=survObj,growObj=growObj)
		
		R02 <- R0Calc(Tmatrix,Fmatrix)
		
		discreteTrans@discreteSurv[param.test] <- discreteTrans@discreteSurv[param.test]/ (1 + delta)
		slam[param.test+count] <- (R02 - R01)/(discreteTrans@discreteSurv[param.test] * delta)
		elam[param.test+count] <- (R02 - R01)/(R01*delta)
		
	}
	
	
	# all mean values coming out of discrete stages
	count <- param.test+count
	for (param.test in 1:length(discreteTrans@meanToCont)) { 
		
		discreteTrans@meanToCont[param.test] <- discreteTrans@meanToCont[param.test]* (1 + delta)
		
		Tmatrix <- createIPMTmatrix(nBigMatrix = nBigMatrix, 
				minSize = minSize, maxSize = maxSize, growObj = growObj, 
				survObj = survObj, discreteTrans = discreteTrans, 
				integrateType = integrateType, correction = correction)
		
		Fmatrix <- createIPMFmatrix(nBigMatrix = nBigMatrix, 
				minSize = minSize, maxSize = maxSize, fecObj = fecObj, 
				integrateType = integrateType, correction = correction,
				preCensus = preCensus,survObj=survObj,growObj=growObj)
		
		R02 <- R0Calc(Tmatrix,Fmatrix)
		
		discreteTrans@meanToCont[param.test] <- discreteTrans@meanToCont[param.test]/ (1 + delta)
		slam[param.test+count] <- (R02 - R01)/(discreteTrans@meanToCont[param.test] * delta)
		elam[param.test+count] <- (R02 - R01)/(R01*delta)
		
	}
	
	
	# all sd values coming out of discrete stages
	count <- param.test+count
	for (param.test in 1:length(discreteTrans@sdToCont)) { 
		
		discreteTrans@sdToCont[param.test] <- discreteTrans@sdToCont[param.test]* (1 + delta)
		
		Tmatrix <- createIPMTmatrix(nBigMatrix = nBigMatrix, 
				minSize = minSize, maxSize = maxSize, growObj = growObj, 
				survObj = survObj, discreteTrans = discreteTrans, 
				integrateType = integrateType, correction = correction)
		
		Fmatrix <- createIPMFmatrix(nBigMatrix = nBigMatrix, 
				minSize = minSize, maxSize = maxSize, fecObj = fecObj, 
				integrateType = integrateType, correction = correction,
				preCensus = preCensus,survObj=survObj,growObj=growObj)
		
		R02 <- R0Calc(Tmatrix,Fmatrix)
		
		discreteTrans@sdToCont[param.test] <- discreteTrans@sdToCont[param.test]/ (1 + delta)
		slam[param.test+count] <- (R02 - R01)/(discreteTrans@sdToCont[param.test][param.test] * delta)
		elam[param.test+count] <- (R02 - R01)/(R01*delta)
		
	}
	
	# parameters linking size to survival
	count <- param.test+count
	for (param.test in 1:length(discreteTrans@survToDiscrete$coefficients)) { 
		
		discreteTrans@survToDiscrete$coefficients[param.test] <- discreteTrans@survToDiscrete$coefficients[param.test]* (1 + delta)
		
		Tmatrix <- createIPMTmatrix(nBigMatrix = nBigMatrix, 
				minSize = minSize, maxSize = maxSize, growObj = growObj, 
				survObj = survObj, discreteTrans = discreteTrans, 
				integrateType = integrateType, correction = correction)
		
		Fmatrix <- createIPMFmatrix(nBigMatrix = nBigMatrix, 
				minSize = minSize, maxSize = maxSize, fecObj = fecObj, 
				integrateType = integrateType, correction = correction,
				preCensus = preCensus,survObj=survObj,growObj=growObj)
		
		R02 <- R0Calc(Tmatrix,Fmatrix)
		
		discreteTrans@survToDiscrete$coefficients[param.test] <- discreteTrans@survToDiscrete$coefficients[param.test]/ (1 + delta)
		slam[param.test+count] <- (R02 - R01)/(discreteTrans@survToDiscrete$coefficient[param.test] * delta)
		elam[param.test+count] <- (R02 - R01)/(R01*delta)
		
	}
	
	print("Did not calculate sensitivities and elasticities for:")
	print(c(names(slam[is.na(slam)])))
	print("Values of zero")
	
	slam <- slam[!is.na(slam)]
	elam <- elam[!is.na(slam)]
	
	return(list(slam = slam, elam = elam))
	
}





### FUNCTIONS FOR EXTRACTING STOCH GROWTH RATE ########################

# Generic approach to get stoch rate
# by sampling list IPM
#
# Parameters - listIPMmatrix - list of IPMs corresponding to different year types
#            - nRunIn - the burnin before establishing lambda_s
#            - tMax - the total time-horizon for getting lambda_s
#
# Returns lambda_s (no density dependence)

stochGrowthRateSampleList <- function(listIPMmatrix,nRunIn,tMax){
			require(MASS)
			
			nmatrices <- length(listIPMmatrix)
			
			nt<-rep(1,length(listIPMmatrix[[1]][,1]))
			Rt<-rep(NA,tMax)
			
			for (t in 1:tMax) {
				year.type <- sample(1:nmatrices,size=1,replace=FALSE)
				nt1<-listIPMmatrix[[year.type]] %*% nt	
				sum.nt1<-sum(nt1)
				Rt[t]<-log(sum.nt1)
				nt<-nt1/sum.nt1
				
			}
			
			res <- mean(Rt[nRunIn:tMax],na.rm=TRUE)
			return(res)
		}

		
		## Function to estimate stochastic growth rate
		## with density dependence in seedling establishment
#
# Parameters - listTmatrix - list of Tmatrices
#            - listFmatrix - list of Fmatrices
#            - nRunIn
#            - tMax
#            - seedList - observed recruits in each of the years
#
# Returns stoch growth rate
		
stochGrowthRateSampleListDD <- function (listTmatrix, listFmatrix,nRunIn, tMax,seedList) {
			require(MASS)
			nmatrices <- length(listTmatrix)
			nt <- rep(1, length(listTmatrix[[1]][, 1]))
			Rt <- rep(NA, tMax)
			for (t in 1:tMax) {
				year.type <- sample(1:nmatrices, size = 1, replace = FALSE)
				nseeds <- sum(listFmatrix[[year.type]]%*%nt)
				pEst <- min(seedList[min(year.type+1,nmatrices)]/nseeds,1)
				nt1 <- (listTmatrix[[year.type]]+pEst*listFmatrix[[year.type]])%*% nt
				sum.nt1 <- sum(nt1)
				Rt[t] <- log(sum.nt1)-log(sum(nt))
				nt <- nt1
			}
			res <- mean(Rt[nRunIn:tMax], na.rm = TRUE)
			return(res)
}
		
		

# Approach to get stoch rate
# with time-varying covariates
#
# Parameters - covariate - the covariates (temperature, etc)
#            - nRunIn - the burnin before establishing lambda_s
#            - tMax - the total time-horizon for getting lambda_s
#
# Returns lambda_s 


stochGrowthRateManyCov <- function(covariate,nRunIn,tMax,
		growthObj,survObj,fecObj,
		nBigMatrix,minSize,maxSize, nMicrosites,
		integrateType="midpoint",correction="none"){
	require(MASS)
	
	
	nt<-rep(1,nBigMatrix)
	Rt<-rep(NA,tMax)
	fecObj@fecConstants[is.na(fecObj@fecConstants)] <- 1 
	tmp.fecObj <- fecObj
	
	#print("fec.const start")
	#print(fecObj@fecConstants)
	
	#density dep in seedling establishment 
	if (sum(nMicrosites)>0) { dd <- TRUE; seeds <- 10000 } else { dd <- FALSE; p.est <- 1}
	
	
	for (t in 1:tMax) {
		if (dd) p.est <- min(nMicrosites[min(t,length(nMicrosites))]/seeds,1) 	
					
		#if you don't do this, rownames return errors...
		covariatet <- covariate[t,]
		row.names(covariatet) <- NULL
		
		#but if you have only one column, then it can forget its a data-frame
		if (ncol(covariate)==1) { 
			covariatet <- data.frame(covariatet)
			colnames(covariatet) <- colnames(covariate)	
		}
		
		tpS <- createIPMTmatrix(nBigMatrix = nBigMatrix, minSize = minSize,
				maxSize = maxSize, chosenCov = covariatet,
				growObj = growthObj, survObj = survObj,
				integrateType=integrateType, correction=correction)
				
		tpF <- createIPMFmatrix(nBigMatrix = nBigMatrix, minSize = minSize,
				maxSize = maxSize, chosenCov = covariatet,
				fecObj = tmp.fecObj,
				integrateType=integrateType, correction=correction)

		#total seeds for next year 
		if (dd) seeds <- sum(tpF%*%nt)
				
		IPM.here <- p.est*tpF@.Data+tpS@.Data
		nt1<-IPM.here %*% nt	
		sum.nt1<-sum(nt1)
		if (!dd) { 
			Rt[t]<-log(sum.nt1)
			nt<-nt1/sum.nt1
		} else {
			Rt[t]<-log(sum.nt1)-log(sum(nt))
			nt<-nt1	
		}
		
		
	}
	
	res <- mean(Rt[nRunIn:tMax],na.rm=TRUE)
	return(res)
}


# Approach to get track pop struct
# with time-varying covariates; assuming density
# dep in seedling establishment (i.e limited no microsites)
#
# Parameters - covariate - the covariates (temperature, etc)
#            - nRunIn - the burnin before establishing lambda_s
#            - tMax - the total time-horizon for getting lambda_s
#
# Returns matrix with time as columns, and pop struct as rows


trackPopStructManyCov<-function(covariate,nRunIn,tMax,
		growthObj,survObj,fecObj,
		nBigMatrix,minSize,maxSize,
		nMicrosites,integrateType="midpoint",correction="none"){
	require(MASS)
	
	nt <- rep(1,nBigMatrix)
	rc <- matrix(NA,nBigMatrix,tMax)
	fecObj@fecConstants[is.na(fecObj@fecConstants)] <- 1 
	tmp.fecObj <- fecObj
	#density dep? 
	if (sum(nMicrosites)>0) { dd <- TRUE; seeds <- 10000 } else { dd <- FALSE}
	
	
	
	for (t in 1:tMax) {
		if (dd) p.est <- min(nMicrosites[min(t,length(nMicrosites))]/seeds,1) 	
	
		#if you don't do this, rownames return errors...
		covariatet <- covariate[t,]
		row.names(covariatet) <- NULL
		
		#but if you have only one column, then it can forget its a data-frame
		if (ncol(covariate)==1) { 
			covariatet <- data.frame(covariatet)
			colnames(covariatet) <- colnames(covariate)	
		}
				
		tpS <- createIPMTmatrix(nBigMatrix = nBigMatrix, minSize = minSize,
				maxSize = maxSize, chosenCov = covariatet,
				growObj = growthObj, survObj = survObj,
				integrateType=integrateType, correction=correction)
		tpF <- createIPMFmatrix(nBigMatrix = nBigMatrix, minSize = minSize,
				maxSize = maxSize, chosenCov = covariatet,
				fecObj = tmp.fecObj,
				integrateType=integrateType, correction=correction)
		
		#total seeds for next year 
		if (dd) seeds <- sum(tpF%*%nt)
		
		IPM.here <- p.est*tpF+tpS
		nt1<-IPM.here %*% nt	
		rc[,t] <- nt1
		nt<-nt1
	
		
	}
	
	return(list(rc=rc,IPM.here=tpS))
}






## FUNCTIONS FOR SENSITIVITY AND ELASTICITY ############################


#parameters - an IPM (with full survival and fecundity complement)
# returns - the sensitivity of every transition for pop growth
sens<-function(A) {
	w<-Re(eigen(A)$vectors[,1]); 
	v<-Re(eigen(t(A))$vectors[,1]);
	vw<-sum(v*w);
	s<-outer(v,w)
	return(s/vw); 
}   

#parameters - an IPM (with full survival and fecundity complement)
# returns - the elasticity of every transition for pop growth
elas<-function(A) {
	s<-sens(A)
	lam<-Re(eigen(A)$values[1]);
	return((s*A)/lam);
}






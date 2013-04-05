

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

setClass("growthObjNegBin",
		representation(fit = "list"))

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
				vitalRatesPerOffspringType = "data.frame",
				Transform = "character",
				offspringRel = "lm",
				sdOffspringSize = "numeric")
)



# =============================================================================
# =============================================================================
## DISCRETE TRANSITION MATRIX OBJECTS ##
# Create a generic discrete transition matrix object
setClass("discreteTrans",
		representation(discreteTrans = "matrix",
				meanToCont = "matrix",
				sdToCont = "matrix",
				moveToDiscrete = "glm"))

# =============================================================================
# =============================================================================
## INTEGER FECUNDITY OBJECTS ##
# Create a generic fecundity object
setClass("fecObjInteger",
		representation(fitFec = "list",
				fecNames = "character",
				fecConstants = "data.frame",
				offspringSplitter = "data.frame",
				fecByDiscrete = "data.frame",
				vitalRatesPerOffspringType = "data.frame",
				Transform = "character",
				offspringRel = "glm",
				thetaOffspringSize = "numeric",
				distOffspring = "character")
)

# =============================================================================
# =============================================================================
## DISCRETE TRANSITION MATRIX INTEGER OBJECTS ##
# Create a generic discrete transition matrix object
setClass("discreteTransInteger",
		representation(discreteTrans = "matrix",
				meanToCont = "matrix",
				thetaToCont = "matrix",
				moveToDiscrete = "glm",
				distToCont = "character"))





######## DEFINE METHODS ##########################################################################################

# =============================================================================
# =============================================================================
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
			
			if (length(grep("expsize",
							survObj@fit$formula))>0) newd$expsize=exp(size)
			if (length(grep("logsize",
							survObj@fit$formula))>0) newd$logsize=log(size)
			if (length(grep("logsize2",
							survObj@fit$formula))>0) newd$logsize2=(log(size))^2
			
#				print(head(newd))
			
			u <- predict(survObj@fit,newd,type="response")
			return(u);
		})

# =============================================================================
# =============================================================================
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
			
			if (length(grep("expsize",
							survObj@fit$formula))>0) newd$expsize=exp(size)
			if (length(grep("logsize",
							survObj@fit$formula))>0) newd$logsize=log(size)
			if (length(grep("logsize2",
							survObj@fit$formula))>0) newd$logsize2=(log(size))^2
			
			u <- predict(survObj@fit,newd,type="link")
			c2 <- ((16 * sqrt(3))/(15 * pi))^2  #from MCMCglmm course notes, search for c2
			u <- logit(u/sqrt(1+c2)) 
			return(u);
		})

# =============================================================================
# =============================================================================
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
			if (length(grep("expsize",
							growthObj@fit$formula))>0) newd$expsize <- exp(size)
			if (length(grep("logsize",
							growthObj@fit$formula))>0) { newd$logsize <- log(size)}
			
			mux <- predict(growthObj@fit,newd,type="response")
			sigmax <-growthObj@sd
			u <- dnorm(sizeNext,mux,sigmax,log=FALSE)  
			return(u);
		})

# =============================================================================
# =============================================================================
#  growth transition (poisson model) probability from size to sizeNext at that covariate level 
#	NOTE DO NOT USE THIS WITH AN IPM!!
setMethod("growth", 
		c("numeric","numeric","data.frame","growthObjPois"),
		function(size,sizeNext,cov,growthObj){
			newd <- data.frame(cbind(cov,size=size),
					stringsAsFactors = FALSE)
			newd$size2 <- size^2
			newd$size3 <- size^3
			
			if (length(grep("expsize",
							growthObj@fit$formula))>0) newd$expsize <- exp(size)
			if (length(grep("logsize",
							growthObj@fit$formula))>0) { newd$logsize <- log(size)}
			
			mux <- predict(growthObj@fit,newd,type="response")
			u <- dpois(sizeNext,mux,log=FALSE)  
			return(u);
		})

# =============================================================================
# =============================================================================
#  growth transition (poisson model) probability from size to sizeNext at that covariate level 
#	NOTE DO NOT USE THIS WITH AN IPM!!
setMethod("growth", 
		c("numeric","numeric","data.frame","growthObjNegBin"),
		function(size,sizeNext,cov,growthObj){
			newd <- data.frame(cbind(cov,size=size),
					stringsAsFactors = FALSE)
			newd$size2 <- size^2
			newd$size3 <- size^3
			
			if (length(grep("expsize",
							growthObj@fit$formula))>0) newd$expsize <- exp(size)
			if (length(grep("logsize",
							growthObj@fit$formula))>0) { newd$logsize <- log(size)}
			
			mux <- predict(growthObj@fit[[1]],newd,type="response")
			u <- dnbinom(sizeNext,mu=mux,size=growthObj@fit[[2]],log=FALSE)  
			return(u);
		})

# =============================================================================
# =============================================================================
# growth for predicting next incr 
setMethod("growth", 
		c("numeric","numeric","data.frame","growthObjIncr"),
		function(size,sizeNext,cov,growthObj){
			newd <- data.frame(cbind(cov,size=size),
					stringsAsFactors = FALSE)
			newd$size2 <- size^2
			newd$size3 <- size^3
			
			if (length(grep("expsize",
							growthObj@fit$formula))>0) newd$expsize <- exp(size)
			if (length(grep("logsize",
							growthObj@fit$formula))>0) newd$logsize=log(size)
			
			mux <- predict(growthObj@fit,newd,type="response")
			
			#print(mux)
			
			sigmax <-growthObj@sd
			u <- dnorm(sizeNext,size+mux,sigmax,log=FALSE)  
			return(u); 
		})


# =============================================================================
# =============================================================================
# growth for predicting next truncated incr 
setMethod("growth", 
		c("numeric","numeric","data.frame","growthObjTruncIncr"),
		function(size,sizeNext,cov,growthObj){
			require(truncnorm)
			
			newd <- data.frame(cbind(cov,size=size),
					stringsAsFactors = FALSE)
			newd$size2 <- size^2
			newd$size3 <- size^3
			
			if (length(grep("expsize",
							growthObj@fit$formula))>0) newd$expsize <- exp(size)
			if (length(grep("logsize",
							names(growthObj@fit$coefficients)))>0) newd$logsize=log(size)
			
			mux <- .predictMuX(grObj=growthObj,newData=newd,covPred=cov)
			sigmax <- growthObj@fit$sigma
			u <- dtruncnorm(sizeNext,a=size,b=Inf,mean=size+mux,sd=sigmax)  
			return(u); 
		})


# =============================================================================
# =============================================================================
# growth for predicting next logincr 
setMethod("growth", 
		c("numeric","numeric","data.frame","growthObjLogIncr"),
		function(size,sizeNext,cov,growthObj){
			
			newd <- data.frame(cbind(cov,size=size),
					stringsAsFactors = FALSE)
			newd$size2 <- size^2
			newd$size3 <- size^3
			
			if (length(grep("expsize",
							growthObj@fit$formula))>0) newd$expsize <- exp(size)
			if (length(grep("logsize",
							growthObj@fit$formula))>0) newd$logsize=log(size)
			
			mux <- predict(growthObj@fit,newd,type="response")
			sigmax <-growthObj@sd
			u <- dlnorm(sizeNext-size,mux,sigmax,log=FALSE)  
			return(u);
		})

# =============================================================================
# =============================================================================
# growth for predicting next logincr with decline var
setMethod("growth", 
		c("numeric", "numeric", "data.frame", "growthObjLogIncrDeclineVar"),
		function(size,sizeNext,cov,growthObj){
			newd <- data.frame(cbind(cov,size=size),
					stringsAsFactors = FALSE)
			newd$size2 <- size^2
			newd$size3 <- size^3
			
			if (length(grep("expsize",
							names(growthObj@fit$coefficients)))>0) newd$expsize <- exp(size)
			if (length(grep("logsize",
							names(growthObj@fit$coefficients)))>0) newd$logsize=log(size)
			
			mux <- .predictMuX(grObj = growthObj, newData = newd, covPred = cov)
			sigmax2 <- growthObj@fit$sigmax2
			var.exp.coef <- growthObj@fit$var.exp.coef
			sigmax2 <- sigmax2 * exp(2 * (var.exp.coef * mux));
			
			u <- dlnorm(sizeNext - size, mux, sqrt(sigmax2), log = FALSE)  
			return(u);
		})

# =============================================================================
# =============================================================================
# growth for predicting next logincr with decline var
setMethod("growthCum", 
		c("numeric", "numeric", "data.frame", "growthObjLogIncrDeclineVar"),
		function(size, sizeNext, cov, growthObj){
			
			newd <- data.frame(cbind(cov,size = size),
					stringsAsFactors = FALSE)
			newd$size2 <- size^2
			newd$size3 <- size^3
			
			if (length(grep("expsize",
							names(growthObj@fit$coefficients)))>0) newd$expsize <- exp(size)
			if (length(grep("logsize", names(growthObj@fit$coefficients))) > 0) newd$logsize = log(size)
			
			mux <- .predictMuX(grObj = growthObj, newData=newd, covPred = cov)
			sigmax2 <- growthObj@fit$sigmax2
			var.exp.coef <- growthObj@fit$var.exp.coef
			sigmax2 <- sigmax2 * exp(2 * (var.exp.coef * mux));
			
			u <- plnorm(sizeNext - size, mux, sqrt(sigmax2), log.p =  FALSE)  
			return(u);
		})


# =============================================================================
# =============================================================================
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
			if (length(grep("expsize",
							growthObj@fit$formula))>0) newd$expsize <- exp(size)
			if (length(grep("logsize",
							growthObj@fit$formula))>0) newd$logsize=log(size)
			mux <- predict(growthObj@fit,newd,type="response")
			sigmax <-growthObj@sd
			u <- pnorm(sizeNext,mux,sigmax, log.p =FALSE)  
			return(u);
		})

# =============================================================================
# =============================================================================
# growth for predicting next incr with a polynomial or log
# using pnorm (i.e. getting cumulative at boundary points and doing difference)
setMethod("growthCum", 
		c("numeric","numeric","data.frame","growthObjIncr"),
		function(size,sizeNext,cov,growthObj){
			
			newd <- data.frame(cbind(cov,size=size),
					stringsAsFactors = FALSE)
			newd$size2 <- size^2
			newd$size3 <- size^3
			if (length(grep("expsize",
							growthObj@fit$formula))>0) newd$expsize <- exp(size)
			if (length(grep("logsize",
							growthObj@fit$formula))>0) newd$logsize=log(size)
			
			mux <- predict(growthObj@fit,newd,type="response")
			sigmax <-growthObj@sd
			u <- pnorm(sizeNext,size+mux,sigmax, log.p =FALSE)  
			return(u); 
		})

# =============================================================================
# =============================================================================
# growth for predicting next truncated incr with cumulative 
setMethod("growthCum", 
		c("numeric","numeric","data.frame","growthObjTruncIncr"),
		function(size,sizeNext,cov,growthObj){
			require(truncnorm)
			
			newd <- data.frame(cbind(cov,size=size),
					stringsAsFactors = FALSE)
			newd$size2 <- size^2
			newd$size3 <- size^3
			
			if (length(grep("expsize",
							names(growthObj@fit$coefficients)))>0) newd$expsize <- exp(size)
			if (length(grep("logsize",
							names(growthObj@fit$coefficients)))>0) newd$logsize=log(size)
			
			mux <- .predictMuX(grObj=growthObj,newData=newd,covPred=cov)
			sigmax <- sqrt(growthObj@fit$sigmax)
			
			u <- ptruncnorm(sizeNext,a=size,b=Inf,mean=size+mux,sd=sigmax)  
			return(u); 
		})

# =============================================================================
# =============================================================================
# growth for predicting next logincr with a polynomial or log
# using pnorm (i.e. getting cumulative at boundary points and doing difference)
setMethod("growthCum", 
		c("numeric","numeric","data.frame","growthObjLogIncr"),
		function(size, sizeNext, cov, growthObj){
			
			newd <- data.frame(cbind(cov,size=size),
					stringsAsFactors = FALSE)
			newd$size2 <- size^2
			newd$size3 <- size^3
			if (length(grep("expsize",
							names(growthObj@fit$coefficients)))>0) newd$expsize <- exp(size)
			if (length(grep("logsize",
							names(growthObj@fit$coefficients)))>0) newd$logsize=log(size)
			
			mux <- predict(growthObj@fit,newd,type="response")
			sigmax <-growthObj@sd
			u <- plnorm(sizeNext-size,mux,sigmax,log.p = FALSE)  
			return(u);
		})

# =============================================================================
# =============================================================================
#Simple growth methods, using  declining variance in growth
# using pnorm (i.e. getting cumulative at boundary points and doing difference)
setMethod("growthCum", 
		c("numeric","numeric","data.frame","growthObjDeclineVar"),
		function(size,sizeNext,cov,growthObj){
			
			newd <- data.frame(cbind(cov,size=size),
					stringsAsFactors = FALSE)
			newd$size2 <- size^2
			newd$size3 <- size^3
			
			if (length(grep("expsize",
							names(growthObj@fit$coefficients)))>0) newd$expsize <- exp(size)
			if (length(grep("logsize",
							names(growthObj@fit$coefficients))) > 0) newd$logsize=log(size)
			
			mux <- .predictMuX(grObj=growthObj,newData=newd,covPred=cov)
			sigmax2 <- growthObj@fit$sigmax2
			var.exp.coef<-growthObj@fit$var.exp.coef
			sigmax2<-sigmax2*exp(2*(var.exp.coef*mux));
			u <- pnorm(sizeNext,mux,sqrt(sigmax2),log.p =FALSE)  
			return(u);
		})

# =============================================================================
# =============================================================================
#Simple growth methods, using  declining variance in growth
setMethod("growth", 
		c("numeric","numeric","data.frame","growthObjDeclineVar"),
		function(size,sizeNext,cov,growthObj){
			newd <- data.frame(cbind(cov,size=size),
					stringsAsFactors = FALSE)
			newd$size2 <- size^2
			newd$size3 <- size^3
			
			if (length(grep("expsize",
							names(growthObj@fit$coefficients)))>0) newd$expsize <- exp(size)
			if (length(grep("logsize",
							names(growthObj@fit$coefficients)))>0) newd$logsize=log(size)
			
			mux <- .predictMuX(grObj=growthObj,newData=newd)
			sigmax2 <- growthObj@fit$sigmax2
			var.exp.coef<-growthObj@fit$var.exp.coef
			sigmax2<-sigmax2*exp(2*(var.exp.coef*mux));
			u <- dnorm(sizeNext,mux,sqrt(sigmax2),log=FALSE)  
			return(u);
		})

# =============================================================================
# =============================================================================
#same but with declining variance in growth on incrment
setMethod("growth", 
		c("numeric","numeric","data.frame","growthObjIncrDeclineVar"),
		function(size, sizeNext, cov, growthObj){
			newd <- data.frame(cbind(cov,size=size),
					stringsAsFactors = FALSE)
			newd$size2 <- size^2
			newd$size3 <- size^3
			
			if (length(grep("expsize",
							names(growthObj@fit$coefficients))) > 0) newd$expsize <- exp(size)
			if (length(grep("logsize",
							names(growthObj@fit$coefficients))) > 0) newd$logsize = log(size)
			
			mux <- .predictMuX(grObj=growthObj,newData=newd,covPred=cov)
			sigmax2 <- growthObj@fit$sigmax2
			var.exp.coef<-growthObj@fit$var.exp.coef
			sigmax2<-sigmax2*exp(2*(var.exp.coef*mux));
			u <- dnorm(sizeNext, size + mux, sqrt(sigmax2), log =  FALSE)  
			return(u);
		})

# =============================================================================
# =============================================================================
## Define a new growth method for Hossfeld growth (classic midpoint rule approach)
setMethod("growth", c("numeric", "numeric", "data.frame", "growthObjHossfeld"), 
		function(size, sizeNext, cov, growthObj) { 
			mux <- size+Hossfeld(size, growthObj@paras) 
			sigmax <- growthObj@sd 
			u <- dnorm(sizeNext, mux, sigmax, log =  FALSE) 
			return(u)
		}) 

# =============================================================================
# =============================================================================
## Define a new growth method for Hossfeld growth (classic midpoint rule approach)
setMethod("growthCum", c("numeric", "numeric", "data.frame", "growthObjHossfeld"), 
		function(size, sizeNext, cov, growthObj) { 
			mux <- size+Hossfeld(size, growthObj@paras) 
			sigmax <- growthObj@sd 
			u <- pnorm(sizeNext, mux, sigmax, log.p =  FALSE) 
			return(u)
		}) 

# =============================================================================
# =============================================================================
# Method combining growth and survival for doing outer (not a generic, so that don't have to
# define for all the different classes / growth and surv take care of that)
growSurv <- function(size,sizeNext,cov,growthObj,survObj){
	growth(size, sizeNext, cov, growthObj) * surv(size,cov,survObj)
}






### CLASSES AND FUNCTIONS FOR MATRICES (ENV, TMATRIX [compound or not], FMATRIX) #################
# =============================================================================
# =============================================================================
#Class for the matrix that holds the env matrix 
setClass("envMatrix",
		representation(nEnvClass = "numeric"), #number of covariate levels
		contains="matrix")

# =============================================================================
# =============================================================================
#Class for the matrix that holds the IPM
setClass("IPMmatrix",
		representation(nDiscrete = "numeric", #number of discrete classes
				nEnvClass = "numeric", #number of covariate levels
				nBigMatrix = "numeric", #the resolution of the IPM
				meshpoints = "numeric",
				env.index = "numeric",
				names.discrete = "character"),
		contains="matrix")

# =============================================================================
# =============================================================================
#Class for the matrix that holds the IPM
#setClass("DiscreteMatrix",
#		representation(nDiscrete = "numeric", #number of discrete classes
#				nEnvClass = "numeric", #number of covariate levels
#				nBigMatrix = "numeric", #the resolution of the IPM
#				meshpoints = "numeric",
#				env.index = "numeric",
#				names.discrete = "character"),
#		contains="matrix")



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
makeIPMPmatrix <- function (nEnvClass = 1, nBigMatrix = 50, minSize = -1, maxSize = 50, 
		chosenCov = data.frame(covariate = 1), growObj, survObj, 
		discreteTrans = 1, integrateType = "midpoint", correction = "none") {
	
	if (class(growObj) == "growthObjPois" | class(growObj) =="growthObjNegBin") 
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
		#remove because alden...
		#get.matrix[nBigMatrix, nBigMatrix] <- get.matrix[nBigMatrix, 
	    #			nBigMatrix] + (1 - sum(get.matrix[, nBigMatrix]))
		get.matrix <- t(t(get.matrix) * surv(size = y, cov = chosenCov, 
						survObj = survObj))
	}

	if (correction == "constant") {
		nvals <- colSums(get.matrix,na.rm=TRUE)
		loc0 <- which(nvals == 0 , arr.ind = TRUE)
		if (length(loc0) > 0) {
			print("warnings - columns that sum to 0 or that have NAs - assuming survival is along the diagonal; plot your Pmatrix to check it")
			get.matrix[,loc0] <- 0
			get.matrix[cbind(loc0, loc0)] <- surv(size = y[loc0], 
					cov = chosenCov, survObj = survObj)
		}
		nvals <- colSums(get.matrix,na.rm=TRUE)
		get.matrix <- t((t(get.matrix)/nvals) * surv(size = y, 
						cov = chosenCov, survObj = survObj))
	}

	if (correction=="discretizeExtremes"){
		tooLow <- growthCum(y,b[1], cov = chosenCov, 
						growthObj = growObj)
		tooHigh <- 1-growthCum(y,b[length(b)],cov = chosenCov, 
						growthObj = growObj)
		get.matrix[1,] <- get.matrix[1,]+tooLow
		get.matrix[nBigMatrix,] <- get.matrix[nBigMatrix,]+tooHigh
	}
	
	
	rc <- new("IPMmatrix", nDiscrete = 0, nEnvClass = 1, nBigMatrix = nBigMatrix, 
			nrow = 1 * nBigMatrix, ncol = 1 * nBigMatrix, meshpoints = y, 
			env.index = rep(1:nEnvClass, each = nBigMatrix), names.discrete = "")
	rc[, ] <- get.matrix

	if (class(discreteTrans) == "discreteTrans") {
		nDisc <- ncol(discreteTrans@meanToCont)
		moveToDiscrete <- predict(discreteTrans@moveToDiscrete, 
				data.frame(size = y, size2 = (y * y)), type = "response")
		cont.to.cont <- get.matrix * matrix(1 - moveToDiscrete, 
				nrow = nBigMatrix, ncol = nBigMatrix, byrow = TRUE)
		disc.to.disc <- discreteTrans@discreteTrans[1:nDisc, 1:nDisc]
		disc.to.cont <- matrix(0, ncol = nDisc, nrow = nBigMatrix)
		cont.to.disc <- matrix(0, nrow = nDisc, ncol = nBigMatrix)
		for (j in 1:nDisc) {
			tmp <- dnorm(y, discreteTrans@meanToCont[j], discreteTrans@sdToCont[j]) * h
			if (correction == "constant") tmp <- tmp/sum(tmp)
			tmp[which(is.na(tmp))] <- 0
			disc.to.cont[, j]  <- discreteTrans@discreteTrans["continuous", j] * tmp
			if (discreteTrans@discreteTrans[j,"continuous"]>0) {
				cont.to.disc[j, ] <- surv(y, chosenCov, survObj) * moveToDiscrete * 
					discreteTrans@discreteTrans[j,"continuous"] / sum(discreteTrans@discreteTrans[1:nDisc,"continuous"]) 
			}
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

# =============================================================================
# =============================================================================
##same function... but for discrete (no integration!)  ####
###
### NOTE ASSUMES TRANSITIONS OUT ARE DISCRETE - THIS MEANS THAT SDTOCONT
### IS THE NEGBINON ISSUE    #########
##
makeIntegerPmatrix <- function (nEnvClass = 1, 
		meshpoints=1:20,
		chosenCov = data.frame(covariate = 1), 
		growObj, survObj, 
		discreteTrans = 1) {

	y <- meshpoints
	nBigMatrix <- length(y)
	
	get.matrix <- t(outer(y, y, growSurv, cov = chosenCov, 
						growthObj = growObj, survObj = survObj))
			
	
	rc <- new("IPMmatrix", nDiscrete = 0, nEnvClass = 1, nBigMatrix = nBigMatrix, 
			nrow = 1 * nBigMatrix, ncol = 1 * nBigMatrix, meshpoints = y, 
			env.index = rep(1:nEnvClass, each = nBigMatrix), names.discrete = "")
	rc[, ] <- get.matrix
	
	if (class(discreteTrans) == "discreteTrans") {
		nDisc <- ncol(discreteTrans@meanToCont)
		moveToDiscrete <- predict(discreteTrans@moveToDiscrete, 
				data.frame(size = y, size2 = (y * y)), type = "response")
		cont.to.cont <- get.matrix * matrix(1 - moveToDiscrete, 
				nrow = nBigMatrix, ncol = nBigMatrix, byrow = TRUE)
		disc.to.disc <- discreteTrans@discreteTrans[1:nDisc, 1:nDisc]
		disc.to.cont <- matrix(0, ncol = nDisc, nrow = nBigMatrix)
		cont.to.disc <- matrix(0, nrow = nDisc, ncol = nBigMatrix)
		for (j in 1:nDisc) {
			if (discreteTrans@distToCont=="poisson") 
				tmp <- dpois(y, discreteTrans@meanToCont[j]) 
			if (discreteTrans@distToCont=="negBin") 
				tmp <- dnbinom(y, mu=discreteTrans@meanToCont[j], size=discreteTrans@thetaToCont[j]) 
		 	tmp[which(is.na(tmp))] <- 0
		 	disc.to.cont[, j] <- discreteTrans@discreteTrans["continuous", j] * tmp
			if (sum(discreteTrans@discreteTrans[j,"continuous"]>0)) {
				cont.to.disc[j, ] <- surv(y, chosenCov, survObj) * moveToDiscrete * 
					discreteTrans@discreteTrans[j,"continuous"] / sum(discreteTrans@discreteTrans[1:nDisc,"continuous"]) 
			}
		}	
		get.disc.matrix <- rbind(cbind(disc.to.disc, cont.to.disc), 
				cbind(disc.to.cont, cont.to.cont))

		rc <- new("DiscreteMatrix", nDiscrete = nDisc, nEnvClass = 1, 
				nBigMatrix = nBigMatrix, nrow = 1 * nBigMatrix + 
						nDisc, ncol = 1 * nBigMatrix + nDisc, meshpoints = y, 
				env.index = rep(1:nEnvClass, each = nBigMatrix), 
				names.discrete = rownames(discreteTrans@discreteTrans)[1:nDisc])
		rc[, ] <- get.disc.matrix
	}
	return(rc)
}

# =============================================================================
# =============================================================================
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


makeCompoundPmatrix <- function(nEnvClass = 2,
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
	if (class(discreteTrans) == "discreteTrans") nDisc <- ncol(discreteTrans@meanToCont) else nDisc <- 0
	
	
	#indexes for slotting in IPMs
	indexes <- rep(1:nEnvClass,each=(nBigMatrix+nDisc))
	
	#megamatrix
	megamatrix <- matrix(0,(nBigMatrix+nDisc)*nEnvClass,(nBigMatrix+nDisc)*nEnvClass) 
	
	#print(indexes)
	
	#loop over habitats / environments
	for (k in 1:nEnvClass) { 
		#IPM for individuals starting in env k
		get.matrix <- makeIPMPmatrix(nEnvClass = nEnvClass, 
				nBigMatrix = nBigMatrix, minSize = minSize, maxSize = maxSize, 
				chosenCov = data.frame(covariate = as.factor(k)), growObj=growObj, survObj=survObj, 
				discreteTrans = discreteTrans, 
				integrateType = integrateType, correction = correction)		
		
		# transit them
		subset <- c(1:nEnvClass)[envMatrix[,k]>0]
		for (j in subset) { 
			megamatrix[indexes==j,indexes==k] <- get.matrix*envMatrix[j,k]; 
		}
		
	}
	
	nmes<-""

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


## Extra functions for use with outer in building the IPMFmatrix ##############
# =============================================================================
# =============================================================================
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

makeIPMFmatrix <- function(fecObj,
		nEnvClass = 1,
		nBigMatrix = 50,
		minSize = -1,
		maxSize = 50,
		chosenCov = data.frame(covariate=1),
		integrateType="midpoint",
		correction="none",
		preCensus=TRUE,
		survObj=NULL,
		growObj=NULL,
		offspringObj=NULL) {
	
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
		##NOTE that the condition is necessary because you might ONLY have discrete offspring
		# in which case correction makes no sense
		if (integrateType=="midpoint"&fecObj@offspringSplitter$continuous>0) { 
			tmp <- t(outer(X=y,Y=y,.fecPreCensus,cov=chosenCov,fecObj=fecObj, offspringObj=offspringObj))*h 
		}
		##NOTE that the condition is necessary because you might ONLY have discrete offspring
		# in which case correction makes no sense
		if (integrateType=="cumul"&fecObj@offspringSplitter$continuous>0) {
			#offspring extremes (pnorm) 
			tmp.cum <- t(outer(X=y,Y=b,.offspringCum,cov=chosenCov,
							fecObj=fecObj,offspringObj=offspringObj))
			tmp <- tmp.cum[2:(nBigMatrix+1),]-tmp.cum[1:nBigMatrix,]
			#put in seed production
			tmp <- t(t(tmp)*.fecRaw(x=y,cov=chosenCov,fecObj=fecObj)[[1]])      
		}
		
		##NOTE that the condition is necessary because you might ONLY have discrete offspring
		# in which case correction makes no sense
		if (correction=="constant"&fecObj@offspringSplitter$continuous>0) {
			# in this case, column sums should equal raw fecundity
			correction.here <-.fecRaw(x=y,cov=chosenCov,fecObj=fecObj)[[1]]/colSums(tmp)
			tmp <- t(t(tmp)*correction.here)
		}
		##NOTE that the condition is necessary because you might ONLY have discrete offspring
		# in which case correction makes no sense
		if (correction=="discretizeExtremes"&fecObj@offspringSplitter$continuous>0){
			tooLow <- .offspringCum(x=y,y=b[1], cov = chosenCov, 
					fecObj = fecObj,offspringObj=offspringObj)
			tooHigh <- 1-.offspringCum(x=y,y=b[length(b)],cov = chosenCov, 
					fecObj = fecObj,offspringObj=offspringObj)
			tmp[1,] <- tmp[1,]+tooLow
			tmp[nBigMatrix,] <- tmp[nBigMatrix,]+tooHigh
		}
				
		# 2. post-census
	} else {
		# if no survObj is provided, it is assumed that all individuals survive to the time of breeding
		if (is.null(survObj)) survObj <- makeSurvObj(data.frame(size=runif(21),surv=rep(1,21)),Formula=surv~size)
		# if no growObj is provided, it is assumed that all individuals retain the same size until the time of breeding		
		if (is.null(growObj)) growObj <- makeGrowthObj(data.frame(size=seq(21),sizeNext=seq(21)),Formula=sizeNext~size)
		#print ("Warning: in the current version of IPMpack, makeIPMFmatrix still ignores the growObj you provided for your post-breeding F matrix. This will be included in a later version. Survival until breeding is already included in this version.")
		##NOTE that the condition is necessary because you might ONLY have discrete offspring
		# in which case correction makes no sense
		if (integrateType=="midpoint"&fecObj@offspringSplitter$continuous>0) { 
			tmp <- t(outer(X=y,Y=y,.fecPostCensus,
							cov=chosenCov,fecObj=fecObj, growObj=growObj,
							survObj=survObj,offspringObj=offspringObj))*h 
		}
		##NOTE that the condition is necessary because you might ONLY have discrete offspring
		# in which case correction makes no sense
		if (integrateType=="cumul"&fecObj@offspringSplitter$continuous>0) {
			#make the extreme bins offspring matrix
			tmp.cum <- t(outer(X=y,Y=b,.offspringCum,cov=chosenCov,
							fecObj=fecObj,offspringObj=offspringObj))
			tmpBabies <- tmp.cum[2:(nBigMatrix+1),]-tmp.cum[1:nBigMatrix,]
			
			#make the extreme bins growth matrix
			tmp.cum <- t(outer(X=y,Y=b,growthCum,cov=chosenCov,
							growObj=growObj,offspringObj=offspringObj))
			tmpGrowth <- tmp.cum[2:(nBigMatrix+1),]-tmp.cum[1:nBigMatrix,]
			tmpGrowth[nBigMatrix,nBigMatrix] <- tmpGrowth[nBigMatrix,nBigMatrix]+
					(1-sum(tmpGrowth[,nBigMatrix]))
			
			#put in survival and seed production
			tmp <- t(t(tmpBabies*tmpGrowth)*surv(size=y,cov=chosenCov,survObj=survObj)*.fecRaw(x=y,cov=chosenCov,fecObj=fecObj)[[1]])
			
		}
		
		##NOTE that the condition is necessary because you might ONLY have discrete offspring
		# in which case correction makes no sense
		if (correction=="constant"&fecObj@offspringSplitter$continuous>0) {
			# in this case, column sums should equal raw fecundity * survival
			correction.here <-.fecRaw(x=y,cov=chosenCov,fecObj=fecObj)[[1]]*surv(size=y,cov=chosenCov,survObj=survObj)/colSums(tmp)
			tmp <- t(t(tmp)*correction.here)
		}
		##NOTE that the condition is necessary because you might ONLY have discrete offspring
		# in which case correction makes no sense
		if (correction=="discretizeExtremes"&fecObj@offspringSplitter$continuous>0){
			tooLow <- .offspringCum(x=y,y=b[1], cov = chosenCov, 
					fecObj = fecObj,offspringObj=offspringObj)
			tooHigh <- 1-.offspringCum(x=y,y=b[length(b)],cov = chosenCov, 
					fecObj = fecObj,offspringObj=offspringObj)
			tmp[1,] <- tmp[1,]+tooLow
			tmp[nBigMatrix,] <- tmp[nBigMatrix,]+tooHigh
		}
	}
	
	get.matrix <- to.cont <- tmp
	
	#discrete classes
	nDisc <- length(fecObj@offspringSplitter)-1
	namesDiscrete <- "NA"
	if (nDisc>0) {
		namesDiscrete <- colnames(fecObj@offspringSplitter[1:nDisc])
		to.discrete <- matrix(0,nrow=nDisc,ncol=nBigMatrix)
		for (i in 1:nDisc) {
			if (length(which(fecObj@vitalRatesPerOffspringType[,namesDiscrete[i]]==1))>1) {
				to.discrete[i,] <- apply(.fecRaw(x=y,cov=chosenCov,fecObj=fecObj)[[2]][which(fecObj@vitalRatesPerOffspringType[,namesDiscrete[i]]==1),],2,prod)*unlist(fecObj@offspringSplitter[namesDiscrete[i]])
			} else {
				to.discrete[i,] <- .fecRaw(x=y,cov=chosenCov,fecObj=fecObj)[[2]][which(fecObj@vitalRatesPerOffspringType[,namesDiscrete[i]]==1),]*unlist(fecObj@offspringSplitter[namesDiscrete[i]])
			}
		}
		from.discrete <- matrix(0,ncol=nDisc,nrow=nDisc+nBigMatrix)
		if (names(fecObj@fecByDiscrete)[1]!="NA.") {
			if (sum(names(fecObj@fecByDiscrete)!=namesDiscrete)>0) stop ("Error - the names of the discrete classes as you provided for the data.frame fecByDiscrete are not 100% the same discrete class names in your data.frame offspringSplitter. They should also be in alphabetical order.")
			for (j in 1:nDisc) {
				for (i in 1:nDisc) {
					from.discrete[i,j] <- unlist(fecObj@offspringSplitter[namesDiscrete[i]]*fecObj@fecByDiscrete[namesDiscrete[j]])
				}
			}
			if (sum(fecObj@fecByDiscrete)>0&fecObj@offspringSplitter["continuous"]>0) {
				print ("WARNING - number and sizes of offspring produced by individuals in discrete classes cannot be calculated yet. The Fmatrix contains zeros instead. Only solution at this point: change the F matrix yourself afterwards.")
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

# =============================================================================
# =============================================================================
### CREATE DISCRETE F MATRIX
makeIntegerFmatrix <- function(fecObj,
		nEnvClass = 1,
		meshpoints=1:20,
		chosenCov = data.frame(covariate=1),
		preCensus=TRUE,
		survObj=NULL,
		growObj=NULL, offspringObj=NULL) {
	
	# boundary points b and mesh points y
	y <- meshpoints;
	nBigMatrix <- length(y)
	
	fecObj@fecConstants[is.na(fecObj@fecConstants)] <- 1
	
	tmp <- matrix(0,ncol=length(y),nrow=length(y))
	
	# 1. pre-census
	if (preCensus) { 
		##NOTE that the condition is necessary cos you might ONLY have discrete offspring
		if (fecObj@offspringSplitter$continuous>0) { 
			tmp <- t(outer(X=y,Y=y,.fecPreCensusInteger,cov=chosenCov,fecObj=fecObj, offspringObj= offspringObj))
		}
		##NOTE that the condition is necessary because you might ONLY have discrete offspring


		# 2. post-census
	} else {
		# if no survObj is provided, it is assumed that all individuals survive to the time of breeding
		if (is.null(survObj)) survObj <- makeSurvObj(data.frame(size=runif(21),surv=rep(1,21)),Formula=surv~size)
		# if no growObj is provided, it is assumed that all individuals retain the same size until the time of breeding		
		if (is.null(growObj)) growObj <- makeGrowthObj(data.frame(size=seq(21),sizeNext=seq(21)),Formula=sizeNext~size,Family="poisson")
		#print ("Warning: in the current version of IPMpack, makeIPMFmatrix still ignores the growObj you provided for your post-breeding F matrix. This will be included in a later version. Survival until breeding is already included in this version.")
		##NOTE that the condition is necessary because you might ONLY have discrete offspring
		# in which case correction makes no sense
		if (fecObj@offspringSplitter$continuous>0) { 
			tmp <- t(outer(X=y,Y=y,.fecPostCensusInteger,
							cov=chosenCov,fecObj=fecObj, growObj=growObj,
							survObj=survObj, offspringObj= offspringObj))
		}

		
	}
	
	get.matrix <- to.cont <- tmp
	
	#discrete classes
	nDisc <- length(fecObj@offspringSplitter)-1
	namesDiscrete <- "NA"
	if (nDisc>0) {
		namesDiscrete <- colnames(fecObj@offspringSplitter[1:nDisc])
		to.discrete <- matrix(0,nrow=nDisc,ncol=nBigMatrix)
		for (i in 1:nDisc) {
			if (length(which(fecObj@vitalRatesPerOffspringType[,namesDiscrete[i]]==1))>1) {
				to.discrete[i,] <- apply(.fecRaw(x=y,cov=chosenCov,fecObj=fecObj)[[2]][which(fecObj@vitalRatesPerOffspringType[,namesDiscrete[i]]==1),],2,prod)*unlist(fecObj@offspringSplitter[namesDiscrete[i]])
			} else {
				to.discrete[i,] <- .fecRaw(x=y,cov=chosenCov,fecObj=fecObj)[[2]][which(fecObj@vitalRatesPerOffspringType[,namesDiscrete[i]]==1),]*unlist(fecObj@offspringSplitter[namesDiscrete[i]])
			}
		}
		from.discrete <- matrix(0,ncol=nDisc,nrow=nDisc+nBigMatrix)
		if (names(fecObj@fecByDiscrete)[1]!="NA.") {
			if (sum(names(fecObj@fecByDiscrete)!=namesDiscrete)>0) stop ("Error - the names of the discrete classes as you provided for the data.frame fecByDiscrete are not 100% the same discrete class names in your data.frame offspringSplitter. They should also be in alphabetical order.")
			for (j in 1:nDisc) {
				for (i in 1:nDisc) {
					from.discrete[i,j] <- unlist(fecObj@offspringSplitter[namesDiscrete[i]]*fecObj@fecByDiscrete[namesDiscrete[j]])
				}
			}
			if (sum(fecObj@fecByDiscrete)>0&fecObj@offspringSplitter["continuous"]>0) {
				print ("WARNING - number and sizes of offspring produced by individuals in discrete classes cannot be calculated yet. The Fmatrix contains zeros instead. Only solution at this point: change the F matrix yourself afterwards.")
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

# =============================================================================
# =============================================================================
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


makeCompoundFmatrix <- function(nEnvClass = 2,
		nBigMatrix = 50,
		minSize = -1,
		maxSize = 50,
		envMatrix,
		fecObj,
		integrateType="midpoint",
		correction="none",
		preCensus=TRUE,
		survObj=NULL,
		growObj=NULL, offspringObj=NULL) {
	
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
		
		get.matrix <- makeIPMFmatrix(nEnvClass =1,
				nBigMatrix = nBigMatrix,
				minSize = minSize,
				maxSize = maxSize,
				chosenCov = data.frame(covariate=as.factor(k)),
				fecObj=fecObj,
				integrateType=integrateType,
				correction=correction,preCensus=preCensus,
				survObj=survObj,
				growObj=growObj, offspringObj= offspringObj)
		
		#print(range(get.matrix))
		#print(get.matrix[1:5,1:5])
		
		# transit them
		subset <- c(1:nEnvClass)[envMatrix[,k]>0]                 
		for (j in subset) {
			megamatrix[indexes==j,indexes==k] <- get.matrix@.Data*envMatrix[j,k]; 
		}
		
	}
	
	#warning about negative numbers should appear from makeIPMFmatrix
	
	
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



# =============================================================================
# =============================================================================
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


makeIPMCmatrix <- function(clonalObj,
		nEnvClass = 1,
		nBigMatrix = 50,
		minSize = -1,
		maxSize = 50,
		chosenCov = data.frame(covariate=1),
		integrateType="midpoint",
		correction="none", 
		preCensus=TRUE,
		survObj=NULL,
		growObj=NULL, offspringObj= NULL) {
	
	rc <- makeIPMFmatrix(fecObj=clonalObj,
			nEnvClass=nEnvClass,
			nBigMatrix=nBigMatrix,
			minSize=minSize,
			maxSize=maxSize,
			chosenCov=chosenCov,
			integrateType=integrateType,
			correction=correction, 
			preCensus=preCensus,
			survObj=survObj,
			growObj=growObj, offspringObj= offspringObj)
	
	return(rc)
}

# =============================================================================
# =============================================================================
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


makeCompoundCmatrix <- function(nEnvClass = 2,
		nBigMatrix = 50,
		minSize = -1,
		maxSize = 50,
		envMatrix,
		clonalObj,
		integrateType="midpoint",
		correction="none",
		preCensus=TRUE,
		survObj=NULL,
		growObj=NULL, offspringObj= NULL) {
	
	rc <- makeCompoundFmatrix(nEnvClass = nEnvClass,
			nBigMatrix = nBigMatrix,
			minSize = minSize,
			maxSize = maxSize,
			envMatrix = envMatrix,
			fecObj = clonalObj,
			integrateType=integrateType,
			correction=correction,
			preCensus=preCensus,
			survObj=survObj,
			growObj=growObj, offspringObj= offspringObj)
	
	return(rc) 
}

# =============================================================================
# =============================================================================
# For a single Pmatrix (!not compound and no discrete stages!), this functions makes a series
# of diagnostic plots - this is defined for growthObj,
# growthObjIncr - modification required
# if other objects used
#
# Parameters - the Pmatrix
#            - growObj - growth object used to build it
#            - survObj - survival object used to build it
#            - dff - the data from which it was built
#            - integrateType - "midpoint", or "cumul" - should
#                 correspond to what the IPM was built with
#            - do you want to implement the corrections for integration? 
# Returns - 
#

diagnosticsPmatrix <- function (Pmatrix, growObj, survObj, dff=NULL, 
		integrateType = "midpoint", 
		correction = "none", cov = data.frame(covariate = 1), 
		sizesToPlot=c(),extendSizeRange=c()) {
	print("Range of Pmatrix is ")
	print(range(c(Pmatrix)))
	if (Pmatrix@meshpoints[1] > 0) 
		new.min <- Pmatrix@meshpoints[1]/2
	else new.min <- Pmatrix@meshpoints[1] * 1.5
	new.max <- 1.5 * max(Pmatrix@meshpoints)
	
	if (length(extendSizeRange)>0 & length(extendSizeRange)!=2) print("require two values for extendSizeRange, reflecting upper and lower limits")
	if (length(extendSizeRange)>0) { new.min <- extendSizeRange[1]; new.max <- extendSizeRange[2]}
		

	#colours for 1) current; 2) bigger size range; 3) bigger no bins
	cols <- c("black","tomato","darkblue")
	ltys <- c(1,1,3)
	
	#matrix with bigger size range
	Pmatrix1 <- makeIPMPmatrix(nEnvClass = 1, nBigMatrix = length(Pmatrix@meshpoints), 
			minSize = new.min, maxSize = new.max, 
			chosenCov = cov, growObj = growObj, survObj = survObj, 
			integrateType = integrateType, correction = correction)
	
	if (sum(is.na(Pmatrix1)) > 0) {
		print("Pmatrix with extended size range returns NAs; changing these to 0, and putting survival value onto diagonal for columns that sum to zero")
		Pmatrix1[is.na(Pmatrix1)] <- 0
		bad <- which(colSums(Pmatrix1) == 0, arr.ind = TRUE)
		if (length(bad) > 0) 
			Pmatrix1[cbind(bad, bad)] <- surv(size = Pmatrix1@meshpoints[bad], 
					cov = cov, survObj = survObj)
	}
	
	
	#matrix with bigger number of bins	  
	Pmatrix2 <- makeIPMPmatrix(nEnvClass = 1, nBigMatrix = floor(length(Pmatrix@meshpoints) * 
							1.5), minSize =Pmatrix@meshpoints[1], maxSize = max(Pmatrix@meshpoints), 
			chosenCov = cov, growObj = growObj, survObj = survObj, 
			integrateType = integrateType, correction = correction)
	
	if (sum(is.na(Pmatrix2)) > 0) {
		print("Pmatrix with extended number of bins returns NAs; changing these to 0, and putting survival value onto diagonal for columns that sum to zero")
		Pmatrix2[is.na(Pmatrix2)] <- 0
		bad <- which(colSums(Pmatrix2) == 0, arr.ind = TRUE)
		if (length(bad) > 0) 
			Pmatrix2[cbind(bad, bad)] <- surv(size = Pmatrix2@meshpoints[bad], 
					cov = cov, survObj = survObj)
	}
	
	#start plots - put original Pmatrix in black
	#par(mfrow = c(3, 3), bty = "l")
	
	#save plot characters
	old.par <- par(no.readonly = TRUE) 
	
	par(mfrow = c(1, 3), bty = "l",pty="s", mar=c(5,4,4,1))
	if (!is.null(dff)) { 
	xlims <- range(c(Pmatrix@meshpoints, Pmatrix1@meshpoints,
						dff$size, dff$sizeNext), na.rm = TRUE)
		a1 <- hist(c(dff$size, dff$sizeNext), xlab = "Sizes observed", axes=FALSE,
			ylab = "Frequency", main = "", xlim = xlims, col="lightgrey", border="lightgrey",plot=TRUE)
	axis(1); axis(2)
	} else {
		a1 <- list(); a1$counts <- 1:100
		xlims <- range(c(Pmatrix@meshpoints, Pmatrix1@meshpoints))
		plot(1:100,type="n",xlab = "Sizes", ylab="", ylim=range(a1$counts), xlim=xlims)	
	} 
	
	lcs <- c(0.8,0.6,0.4)	
	lcs.x <- mean(xlims) 
	
	text(lcs.x ,max(a1$counts)*lcs[1],"Current",pos=3,col=cols[1])
	arrows(Pmatrix@meshpoints[1], max(a1$counts)*lcs[1],Pmatrix@meshpoints[length(Pmatrix@meshpoints)], max(a1$counts)*lcs[1],col=cols[1], length=0.1, code=3,lty=ltys[1])
	text(lcs.x ,max(a1$counts)*lcs[2],"Extended range",pos=3,col=cols[2])
	arrows(Pmatrix1@meshpoints[1], max(a1$counts)*lcs[2],Pmatrix1@meshpoints[length(Pmatrix1@meshpoints)], max(a1$counts)*lcs[2],col=cols[2], length=0.1, code=3,lty=ltys[1])
	text(lcs.x ,max(a1$counts)*lcs[3],"Increased no of bins",pos=3,col=cols[3])
	arrows(Pmatrix2@meshpoints[1], max(a1$counts)*lcs[3],Pmatrix2@meshpoints[length(Pmatrix2@meshpoints)], max(a1$counts)*lcs[3],col=cols[3], length=0.1, code=3,lty=ltys[1])
	
	#legend("topright", legend = "fitted range", col = "black", lty = 2, bty = "n") ##!!! change this
	title("Size range")
	
	#survival sums
	lims <- range(c(colSums(Pmatrix), surv(Pmatrix@meshpoints, cov, survObj)))
	plot(colSums(Pmatrix), surv(Pmatrix@meshpoints, cov, survObj), 
			type = "n", xlab = "Surviving in Pmatrix", ylab = "Should be surviving",col=cols[1],lty=ltys[1], 
			xlim=lims,ylim=lims)
	abline(0, 1, lwd=2,col="grey")
	points(colSums(Pmatrix), surv(Pmatrix@meshpoints, cov, survObj), type = "l", col = cols[1],lty=ltys[1])
	points(colSums(Pmatrix1), surv(Pmatrix1@meshpoints, cov, survObj), type = "l", col = cols[2],lty=ltys[2])
	points(colSums(Pmatrix2), surv(Pmatrix2@meshpoints, cov, survObj), type = "l", col = cols[3],lty=ltys[3])
	title("Survival")
	#text(lims[2],lims[2],"(0,1)",pos=1)	
	
	# mtext("Points should lie along the 0,1 line, shown in grey", side=1,outer=TRUE,line=-1)
	
	
	
	LE <- meanLifeExpect(Pmatrix)
	LE1 <- meanLifeExpect(Pmatrix1)
	LE2 <- meanLifeExpect(Pmatrix2)
	
	plot(Pmatrix@meshpoints, LE, type = "l", xlim = range(Pmatrix1@meshpoints), 
			ylim = range(c(LE, LE1)), xlab = "Sizes", ylab = "Life expectancy",col=cols[1],lty=ltys[1])
	points(Pmatrix1@meshpoints, LE1, type = "l", col = cols[2],lty=ltys[2])
	points(Pmatrix2@meshpoints, LE2, type = "l", col = cols[3],lty=ltys[3])
	legend("topleft", legend = c("Current", "Extended range", "Increased no of bins"), col = cols, lty = ltys, bty = "n")
	title("Life expectancy")
	
	#mtext("Increasing size range (red) or number of bins (blue) should not alter life expectancy estimates", side=1,outer=TRUE,line=-1)
	
	
	
	if (length(sizesToPlot)==0 & !is.null(dff)) sizesToPlot <- as.numeric(quantile(dff$size,c(0.25, 0.5,0.75),na.rm=TRUE))
	if (length(sizesToPlot)==0 & is.null(dff)) sizesToPlot <- as.numeric(quantile(Pmatrix@meshpoints,c(0.25, 0.5,0.75),na.rm=TRUE))
	
	loctest <- rep(NA,length(sizesToPlot))
	
	print("Please hit any key for the next plot")	
	scan(quiet="TRUE")
	
	par(mfrow = c(2, 3), bty = "l",pty="s")
	
	
	for (kk in c(1,3)) {
		if (kk==1) Pmat <- Pmatrix	
		if (kk==3) Pmat <- Pmatrix2	
		
		h <- diff(Pmat@meshpoints)[1]
		testSizes <- seq(min(Pmat@meshpoints), max(Pmat@meshpoints), 
				length = 5000)
		for (j in 1:3) {
			
			loctest[j] <- which(abs(Pmat@meshpoints-sizesToPlot[j])==min(abs(Pmat@meshpoints-sizesToPlot[j])),arr.ind=TRUE)[1]
			
			#print(loctest[j])
			
			ps <- surv(Pmat@meshpoints[loctest[j]], cov, survObj)
			newd <- data.frame(size = Pmat@meshpoints[loctest[j]], 
					size2 = Pmat@meshpoints[loctest[j]]^2, 
					size3 = Pmat@meshpoints[loctest[j]]^3, 
					covariate = Pmat@env.index[1])
			if (length(grep("expsize", names(growObj@fit$coefficients)))) 
				newd$expsize = log(Pmat@meshpoints[loctest[j]])
			if (length(grep("logsize", names(growObj@fit$coefficients)))) 
				newd$logsize = log(Pmat@meshpoints[loctest[j]])
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
				mux <- Pmat@meshpoints[loctest[j]] + mux
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
			
			#  range.x <- range(Pmat@meshpoints[loctest[j]] + c(-3.5 * sqrt(sigmax2), +3.5 * sqrt(sigmax2)))
			range.x <- range(mux + c(-3.5 * sqrt(sigmax2), +3.5 * sqrt(sigmax2)))
			
			plot(Pmat@meshpoints, Pmat@.Data[, loctest[j]]/h/ps, 
					type = "n", xlim = range.x, xlab = "Size next", ylab = "Kernel")
			
			#if (j == 1 & kk==1) title("Numerical resolution and growth")
			if (j == 2 & kk==1) mtext("Numerical resolution and growth",3,line=4,font=2)
			if (j ==1 & kk==1) title("Current Pmatrix")
			if (j ==1 & kk==3) title("Increased bins")
			
			
			
			for (k in 1:length(Pmat@meshpoints)) {
				points(c(Pmat@meshpoints[k]) + c(-h/2, h/2), rep(Pmat@.Data[k,loctest[j]], 2)/h/ps, type = "l",col=cols[kk])
				points(rep(Pmat@meshpoints[k] + h/2, 2), c(0, Pmat@.Data[k, loctest[j]]/h/ps), type = "l", 
						lty = 1,col=cols[kk])
				points(rep(Pmat@meshpoints[k] - h/2, 2), c(0,Pmat@.Data[k, loctest[j]]/h/ps), type = "l", 
						lty = 1,col=cols[kk])
			}
			if (length(grep("logincr", tolower(as.character(class(growObj))))) == 0 &
					length(grep("trunc", tolower(as.character(class(growObj))))) == 0) {
				points(testSizes, dnorm(testSizes, mux, sqrt(sigmax2)), type = "l", col = 2)
			} else {
				if (length(grep("trunc", tolower(as.character(class(growObj))))) > 0) {
					require(truncnorm)
					points(testSizes, dtruncnorm(testSizes, a = Pmat@meshpoints[loctest[j]], 
									b = Inf, mean = mux, sd = sqrt(sigmax2)), type = "l",col = 2)
				} else {
					points(testSizes, dlnorm(testSizes - Pmat@meshpoints[loctest[j]], 
									mux, sqrt(sigmax2)), type = "l", col = 2)
				}
			}
			
			legend("topright",legend=paste("size=",round(Pmat@meshpoints[loctest[j]],2)), col = "white", 
					lty = 1, bty = "n")
			
		}
	}
	
	#reset plot characters
	on.exit(par(old.par))
	
}



##### Functions to identify sensible numbers of bins - help file on desktop  with data-frame setup
#####

# =============================================================================
# =============================================================================
convergeIPM<-function(growObj, survObj, fecObj, nBigMatrix, minSize, maxSize, 
		discreteTrans = 1, integrateType = "midpoint", correction = "none", 
		preCensus = TRUE, tol=1e-4,
		binIncrease=5, chosenBin=1, response="lambda"){
	
	if (response=="lambda") 
		rc <- .convergeLambda(growObj=growObj, survObj=survObj, fecObj=fecObj, 
				nBigMatrix=nBigMatrix, minSize=minSize, maxSize=maxSize, 
				discreteTrans =discreteTrans, integrateType = integrateType, 
				correction = correction, preCensus = preCensus,
				tol=tol,binIncrease=binIncrease)
			
	if (response=="R0")	
		rc <- .convergeR0(growObj=growObj, survObj=survObj, fecObj=fecObj, 
				nBigMatrix=nBigMatrix, minSize=minSize, maxSize=maxSize, 
				discreteTrans =discreteTrans, integrateType = integrateType, 
				correction = correction, preCensus = preCensus,
				tol=tol,binIncrease=binIncrease)
			
	if (response=="lifeExpect")	
		rc <- .convergeLifeExpectancy(growObj=growObj, survObj=survObj, 
			nBigMatrix=nBigMatrix, minSize=minSize, maxSize=maxSize, 
			discreteTrans =discreteTrans, integrateType = integrateType, 
			correction = correction, 
			tol=tol,binIncrease=binIncrease, chosenBin=chosenBin)


	return(rc)

}
	



### FUNCTIONS FOR EXTRACTING INTEGRATED DEMOGRAPHIC MEASURES ########################
### i.e. life-expectancy and passage time

# =============================================================================
# =============================================================================
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


# =============================================================================
# =============================================================================
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




# =============================================================================
# =============================================================================
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

# =============================================================================
# =============================================================================
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

# =============================================================================
# =============================================================================
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

# =============================================================================
# =============================================================================
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
# =============================================================================
# =============================================================================
## TO DO - ADJUST TO ALLOW DISCRETE CLASSES ##

## Function to predict distribution in x time-steps given starting
## distribution and IPM (with a single covariate)
#
#
#  Parameters - startingSizes - vector of starting sizes
#             - IPM the IPM (Pmatrix if only intrested in grow surv; Pmatrix + Fmatrix otherwise)
#             - n.time.steps - number of time steps
#             - startingEnv - vector of starting env, same length as startingSizes, or length=1
#
# Returns - a list including starting numbers in each IPM size class (n.new.dist0) and
#                            final numbers in each IPM size class (n.new.dist)
#
#
predictFutureDistribution <- function(startingSizes,IPM, n.time.steps, startingEnv=1) {
	
	# turn starting sizes into the resolution of the IPM bins
	breakpoints <- c(IPM@meshpoints-(IPM@meshpoints[2]-IPM@meshpoints[1]),
			IPM@meshpoints[length(IPM@meshpoints)]+(IPM@meshpoints[2]-IPM@meshpoints[1]))
	
	# setup slightly different for coompound or non compound dists
	if (IPM@nEnvClass>1) {
		if (length(startingEnv)==1) startingEnv <- rep(startingEnv, length(startingSizes))
		compound <- TRUE
		env.index <- IPM@env.index
		n.new.dist <- rep(0,length(IPM[1,]))
		for (ev in 1:IPM@nEnvClass) { 			
			index.new.dist <- findInterval(startingSizes[startingEnv==ev],breakpoints,all.inside=TRUE)
			loc.sizes <- table(index.new.dist); 
			n.new.dist[ev==IPM@env.index][as.numeric(names(loc.sizes))] <- loc.sizes
		}
		n.new.dist0 <- n.new.dist
	} else {
		compound <- FALSE
		index.new.dist <- findInterval(startingSizes,breakpoints,all.inside=TRUE)
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

# =============================================================================
# =============================================================================
## TO DO - ADJUST TO ALLOW DISCRETE CLASSES ##

## Function to see how long it takes to get from a starting distribution to an final size
##
#  Parameters - startingSizes - vector of starting sizes (in any order)
#             - IPM the IPM (just a Pmatrix)
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







# =============================================================================
# =============================================================================
# Calculate R0 
#
# Parameters - Fmatrix, Pmatrix
#
# Returns R0
R0Calc<-function(Pmatrix, Fmatrix){
	require(MASS)
	Imatrix <- matrix(0, length(Pmatrix[1,]), length(Pmatrix[1,])); 
	diag(Imatrix) <- 1
	Nmatrix <- ginv(Imatrix - Pmatrix);
	Rmatrix <- Fmatrix %*% Nmatrix
	ave.R0 <- Re(eigen(Rmatrix)$values[1])
	return(ave.R0)
}

# =============================================================================
# =============================================================================
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
largeMatrixCalc <- function(Pmatrix, Fmatrix, tol = 1.e-8){
	require(Matrix)
	A2 <- Matrix(Pmatrix + Fmatrix);
	nt <- Matrix(1,length(Pmatrix[1,]), 1);
	nt1 <- nt; 
	
	h1 <- diff(Pmatrix@meshpoints)[1]
	
	qmax <- 1000;
	lam <- 1; 
	while(qmax > tol) {
		nt1 <- A2 %*% nt;
		qmax <- sum(abs((nt1 - lam * nt)@x));  
		lam <- sum(nt1@x); 
		nt@x <- (nt1@x) / lam; #slight cheat  
		#cat(lam,qmax,"\n");
	} 
	nt <- matrix(nt@x, length(Pmatrix[1,]), 1); 
	stableDist <- nt / (h1 * sum(nt)); #normalize so that integral=1
	lamStable <- lam; 
	
	# Check works   
	qmax <- sum(abs(lam * nt - (Pmatrix + Fmatrix) %*% nt))
	cat("Convergence: ", qmax, " should be less than ", tol, "\n")
	
	
	return(list(lam = lam, stableDist = stableDist, h1 = h1)) 
	
}

# =============================================================================
# =============================================================================
## Sensitivity of parameters - works for an IPM built out of
## growth, survival, discreteTrans, fecundity and clonality objects.
##
##
sensParams <- function (growObj, survObj, fecObj=NULL, clonalObj=NULL,
		nBigMatrix, minSize, maxSize,
		chosenCov = data.frame(covariate = 1), discreteTrans = 1,
		integrateType = "midpoint", correction = "none", 
		preCensusFec = TRUE, postCensusSurvObjFec = NULL, postCensusGrowObjFec = NULL,  
		preCensusClonal = TRUE, postCensusSurvObjClonal = NULL, postCensusGrowObjClonal = NULL,  
		delta = 1e-04, response="lambda", chosenBin=1) {
	
	if (response!="lambda" & response!="R0" & response !="lifeExpect")
		stop("response must be one of lambda or R0 or lifeExpect")
	
	nmes <- elam <- slam <- c()
	
	# get the base
	Pmatrix <- makeIPMPmatrix(nBigMatrix = nBigMatrix, minSize = minSize,
			chosenCov = chosenCov, maxSize = maxSize, growObj = growObj,
			survObj = survObj, discreteTrans = discreteTrans, integrateType = integrateType,
			correction = correction)
	if (!is.null(fecObj)) {
		Fmatrix <- makeIPMFmatrix(nBigMatrix = nBigMatrix, minSize = minSize,
				chosenCov = chosenCov, maxSize = maxSize, fecObj = fecObj,
				integrateType = integrateType, correction = correction,
				preCensus = preCensusFec, survObj = postCensusSurvObjFec, growObj = postCensusGrowObjFec)
	} else {Fmatrix <- Pmatrix*0 }
	if (!is.null(clonalObj)) {
		Cmatrix <- makeIPMCmatrix(nBigMatrix = nBigMatrix, minSize = minSize,
				chosenCov = chosenCov, maxSize = maxSize, clonalObj = clonalObj,
				integrateType = integrateType, correction = correction,
				preCensus = preCensusClonal, survObj = postCensusSurvObjClonal, 
				growObj = postCensusGrowObjClonal)
	} else {Cmatrix <- Pmatrix*0 }
	
	IPM <- Pmatrix + Fmatrix + Cmatrix
	
	if (response=="lambda") rc1 <- Re(eigen(IPM)$value[1])
	if (response=="R0") rc1 <- R0Calc(Pmatrix, Fmatrix+Cmatrix)
	if (response=="lifeExpect") rc1 <- meanLifeExpect(Pmatrix)[chosenBin]
	
	# 1. survival
	for (j in 1:length(survObj@fit$coeff)) {
		survObj@fit$coefficients[j] <- survObj@fit$coefficients[j] * (1 + delta)
		Pmatrix <- makeIPMPmatrix(nBigMatrix = nBigMatrix,
				minSize = minSize, maxSize = maxSize, growObj = growObj,
				survObj = survObj, discreteTrans = discreteTrans,
				chosenCov = chosenCov, integrateType = integrateType,
				correction = correction)
		if (!is.null(fecObj)) {
			Fmatrix <- makeIPMFmatrix(nBigMatrix = nBigMatrix,
					minSize = minSize, maxSize = maxSize, fecObj = fecObj,
					integrateType = integrateType, correction = correction,
					chosenCov = chosenCov, preCensus = preCensusFec, 
					survObj = postCensusSurvObjFec, growObj = postCensusGrowObjFec)
		} 
		if (!is.null(clonalObj)) {
			Cmatrix <- makeIPMCmatrix(nBigMatrix = nBigMatrix, minSize = minSize,
					chosenCov = chosenCov, maxSize = maxSize, clonalObj = clonalObj,
					integrateType = integrateType, correction = correction,
					preCensus = preCensusClonal, survObj = postCensusSurvObjClonal, 
					growObj = postCensusGrowObjClonal)
		} 
		
		IPM <- Pmatrix + Fmatrix + Cmatrix
		
		if (response=="lambda") rc2 <- Re(eigen(IPM)$value[1])
		if (response=="R0") rc2 <- R0Calc(Pmatrix, Fmatrix+Cmatrix)
		if (response=="lifeExpect") rc2 <- meanLifeExpect(Pmatrix)[chosenBin]
		
		survObj@fit$coefficients[j] <- survObj@fit$coefficients[j]/(1 + delta)
		
		slam <- c(slam, (rc2 - rc1)/((as.numeric(survObj@fit$coefficients[j]))* delta))
		elam <- c(elam, (rc2 - rc1)/(rc1 *delta))
		nmes <- c(nmes, as.character(paste("survival:",names(survObj@fit$coeff)[j])))
	}
	
	# 2 growth
	for (j in 1:length(growObj@fit$coeff)) {
		growObj@fit$coefficients[j] <- growObj@fit$coefficients[j] * (1 + delta)
		Pmatrix <- makeIPMPmatrix(nBigMatrix = nBigMatrix,
				minSize = minSize, maxSize = maxSize, growObj = growObj,
				chosenCov = chosenCov, survObj = survObj, discreteTrans = discreteTrans,
				integrateType = integrateType, correction = correction)
		if (!is.null(fecObj)) {
			Fmatrix <- makeIPMFmatrix(nBigMatrix = nBigMatrix,
					minSize = minSize, maxSize = maxSize, fecObj = fecObj,
					chosenCov = chosenCov, integrateType = integrateType,
					correction = correction, preCensus = preCensusFec, 
					survObj = postCensusSurvObjFec, growObj = postCensusGrowObjFec)
		}
		if (!is.null(clonalObj)) {
			Cmatrix <- makeIPMCmatrix(nBigMatrix = nBigMatrix, minSize = minSize,
					chosenCov = chosenCov, maxSize = maxSize, clonalObj = clonalObj,
					integrateType = integrateType, correction = correction,
					preCensus = preCensusClonal, survObj = postCensusSurvObjClonal, 
					growObj = postCensusGrowObjClonal)
		}
		
		IPM <- Pmatrix + Fmatrix + Cmatrix
		
		if (response=="lambda") rc2 <- Re(eigen(IPM)$value[1])
		if (response=="R0") rc2 <- R0Calc(Pmatrix, Fmatrix+Cmatrix)
		if (response=="lifeExpect") rc2 <- meanLifeExpect(Pmatrix)[chosenBin]
		
		growObj@fit$coefficients[j] <- growObj@fit$coefficients[j]/(1 + delta)
		
		slam <- c(slam, (rc2 - rc1)/(as.numeric(growObj@fit$coefficients[j]) * delta))
		elam <- c(elam, (rc2 - rc1)/(rc1 * delta))
		nmes <- c(nmes, as.character(paste("growth:",names(growObj@fit$coeff)[j])))
	}
	
	growObj@sd <- growObj@sd * (1 + delta)
	Pmatrix <- makeIPMPmatrix(nBigMatrix = nBigMatrix, minSize = minSize,
			maxSize = maxSize, growObj = growObj, survObj = survObj,
			chosenCov = chosenCov, discreteTrans = discreteTrans,
			integrateType = integrateType, correction = correction)
	if (!is.null(fecObj)) {
		Fmatrix <- makeIPMFmatrix(nBigMatrix = nBigMatrix, minSize = minSize,
				maxSize = maxSize, fecObj = fecObj, integrateType = integrateType,
				chosenCov = chosenCov, correction = correction, preCensus = preCensusFec, 
				survObj = postCensusSurvObjFec, growObj = postCensusGrowObjFec)
	}
	if (!is.null(clonalObj)) {
		Cmatrix <- makeIPMCmatrix(nBigMatrix = nBigMatrix, minSize = minSize,
				chosenCov = chosenCov, maxSize = maxSize, clonalObj = clonalObj,
				integrateType = integrateType, correction = correction,
				preCensus = preCensusClonal, survObj = postCensusSurvObjClonal, 
				growObj = postCensusGrowObjClonal)
	}
	IPM <- Pmatrix + Fmatrix + Cmatrix
	
	if (response=="lambda") rc2 <- Re(eigen(IPM)$value[1])
	if (response=="R0") rc2 <- R0Calc(Pmatrix, Fmatrix+Cmatrix)
	if (response=="lifeExpect") rc2 <- meanLifeExpect(Pmatrix)[chosenBin]
	
	growObj@sd <- growObj@sd / (1 + delta)
	
	slam <- c(slam,(rc2 - rc1)/(growObj@sd * delta))
	elam <- c(elam, (rc2 - rc1)/(rc1 * delta))
	nmes <- c(nmes, "growth: sd")
	
	# 3. DiscreteTrans
	if (class(discreteTrans)=="discreteTrans") {
		for (j in 1:(ncol(discreteTrans@discreteTrans)-1)) {
			for (i in 1:(nrow(discreteTrans@discreteTrans)-1)) {
				discreteTrans@discreteTrans[i,j]<-discreteTrans@discreteTrans[i,j] * (1 + delta)
				Pmatrix <- makeIPMPmatrix(nBigMatrix = nBigMatrix, minSize = minSize,
						maxSize = maxSize, growObj = growObj, survObj = survObj,
						chosenCov = chosenCov, discreteTrans = discreteTrans,
						integrateType = integrateType, correction = correction)
				if (!is.null(fecObj)) {
					Fmatrix <- makeIPMFmatrix(nBigMatrix = nBigMatrix, minSize = minSize,
							maxSize = maxSize, fecObj = fecObj, integrateType = integrateType,
							chosenCov = chosenCov, correction = correction, preCensus = preCensusFec, 
							survObj = postCensusSurvObjFec, growObj = postCensusGrowObjFec)
				}
				if (!is.null(clonalObj)) {
					Cmatrix <- makeIPMCmatrix(nBigMatrix = nBigMatrix, minSize = minSize,
							chosenCov = chosenCov, maxSize = maxSize, clonalObj = clonalObj,
							integrateType = integrateType, correction = correction,
							preCensus = preCensusClonal, survObj = postCensusSurvObjClonal, 
							growObj = postCensusGrowObjClonal)
				}
				IPM <- Pmatrix + Fmatrix + Cmatrix
				
				if (response=="lambda") rc2 <- Re(eigen(IPM)$value[1])
				if (response=="R0") rc2 <- R0Calc(Pmatrix, Fmatrix+Cmatrix)
				if (response=="lifeExpect") rc2 <- meanLifeExpect(Pmatrix)[chosenBin]
				
				discreteTrans@discreteTrans[i,j]<-discreteTrans@discreteTrans[i,j] / (1 + delta)
				
				slam <- c(slam,(rc2 - rc1)/(discreteTrans@discreteTrans[i,j] * delta))
				elam <- c(elam, (rc2 - rc1)/(rc1 * delta))
				nmes <- c(nmes, as.character(paste("discrete:",dimnames(discreteTrans@discreteTrans)[[2]][j],"to",dimnames(discreteTrans@discreteTrans)[[1]][i])))
			}
		}
		#if there is more than 2 discrete stages (beyond "continuous" "dead" and one discrete stage)
        #then moveToDiscrete tells you how many of surviving continuous individuals are going into 
		#discrete classes, but how they distributed also; which is the last column in discreteTrans 
		if (nrow(discreteTrans@discreteTrans)>3) {
			for (i in 1:(nrow(discreteTrans@discreteTrans)-2)) {
				discreteTrans@discreteTrans[i,"continuous"]<-discreteTrans@discreteTrans[i,"continuous"] * (1 + delta)
				Pmatrix <- makeIPMPmatrix(nBigMatrix = nBigMatrix, minSize = minSize,
						maxSize = maxSize, growObj = growObj, survObj = survObj,
						chosenCov = chosenCov, discreteTrans = discreteTrans,
						integrateType = integrateType, correction = correction)
				if (!is.null(fecObj)) {
					Fmatrix <- makeIPMFmatrix(nBigMatrix = nBigMatrix, minSize = minSize,
							maxSize = maxSize, fecObj = fecObj, integrateType = integrateType,
							chosenCov = chosenCov, correction = correction, preCensus = preCensusFec, 
							survObj = postCensusSurvObjFec, growObj = postCensusGrowObjFec)
				}
				if (!is.null(clonalObj)) {
					Cmatrix <- makeIPMCmatrix(nBigMatrix = nBigMatrix, minSize = minSize,
							chosenCov = chosenCov, maxSize = maxSize, clonalObj = clonalObj,
							integrateType = integrateType, correction = correction,
							preCensus = preCensusClonal, survObj = postCensusSurvObjClonal, 
							growObj = postCensusGrowObjClonal)
				}
				IPM <- Pmatrix + Fmatrix + Cmatrix
				
				if (response=="lambda") rc2 <- Re(eigen(IPM)$value[1])
				if (response=="R0") rc2 <- R0Calc(Pmatrix, Fmatrix+Cmatrix)
				if (response=="lifeExpect") rc2 <- meanLifeExpect(Pmatrix)[chosenBin]
				
				discreteTrans@discreteTrans[i,"continuous"]<-discreteTrans@discreteTrans[i,"continuous"] / (1 + delta)
				
				slam <- c(slam,(rc2 - rc1)/(discreteTrans@discreteTrans[i,"continuous"] * delta))
				elam <- c(elam, (rc2 - rc1)/(rc1 * delta))
				nmes <- c(nmes, as.character(paste("discrete: Continuous to",dimnames(discreteTrans@discreteTrans)[[1]][i])))
			}
		}
		
		for (j in 1:length(discreteTrans@meanToCont)) {
			discreteTrans@meanToCont[1,j]<-discreteTrans@meanToCont[1,j] * (1 + delta)
			Pmatrix <- makeIPMPmatrix(nBigMatrix = nBigMatrix, minSize = minSize,
					maxSize = maxSize, growObj = growObj, survObj = survObj,
					chosenCov = chosenCov, discreteTrans = discreteTrans,
					integrateType = integrateType, correction = correction)
			if (!is.null(fecObj)) {
				Fmatrix <- makeIPMFmatrix(nBigMatrix = nBigMatrix, minSize = minSize,
						maxSize = maxSize, fecObj = fecObj, integrateType = integrateType,
						chosenCov = chosenCov, correction = correction, preCensus = preCensusFec, 
						survObj = postCensusSurvObjFec, growObj = postCensusGrowObjFec)
			}
			if (!is.null(clonalObj)) {
				Cmatrix <- makeIPMCmatrix(nBigMatrix = nBigMatrix, minSize = minSize,
						chosenCov = chosenCov, maxSize = maxSize, clonalObj = clonalObj,
						integrateType = integrateType, correction = correction,
						preCensus = preCensusClonal, survObj = postCensusSurvObjClonal, 
						growObj = postCensusGrowObjClonal)
			}
			IPM <- Pmatrix + Fmatrix + Cmatrix
			
			if (response=="lambda") rc2 <- Re(eigen(IPM)$value[1])
			if (response=="R0") rc2 <- R0Calc(Pmatrix, Fmatrix+Cmatrix)
			if (response=="lifeExpect") rc2 <- meanLifeExpect(Pmatrix)[chosenBin]
			
			discreteTrans@meanToCont[1,j]<-discreteTrans@meanToCont[1,j] / (1 + delta)
			
			slam <- c(slam,(rc2 - rc1)/(discreteTrans@meanToCont[1,j] * delta))
			elam <- c(elam, (rc2 - rc1)/(rc1 * delta))
			nmes <- c(nmes, as.character(paste("discrete: meanToCont",dimnames(discreteTrans@meanToCont)[[2]][j])))
		}
		
		for (j in 1:length(discreteTrans@sdToCont)) {
			discreteTrans@sdToCont[1,j]<-discreteTrans@sdToCont[1,j] * (1 + delta)
			Pmatrix <- makeIPMPmatrix(nBigMatrix = nBigMatrix, minSize = minSize,
					maxSize = maxSize, growObj = growObj, survObj = survObj,
					chosenCov = chosenCov, discreteTrans = discreteTrans,
					integrateType = integrateType, correction = correction)
			if (!is.null(fecObj)) {
				Fmatrix <- makeIPMFmatrix(nBigMatrix = nBigMatrix, minSize = minSize,
						maxSize = maxSize, fecObj = fecObj, integrateType = integrateType,
						chosenCov = chosenCov, correction = correction, preCensus = preCensusFec, 
						survObj = postCensusSurvObjFec, growObj = postCensusGrowObjFec)
			}
			if (!is.null(clonalObj)) {
				Cmatrix <- makeIPMCmatrix(nBigMatrix = nBigMatrix, minSize = minSize,
						chosenCov = chosenCov, maxSize = maxSize, clonalObj = clonalObj,
						integrateType = integrateType, correction = correction,
						preCensus = preCensusClonal, survObj = postCensusSurvObjClonal, 
						growObj = postCensusGrowObjClonal)
			}
			IPM <- Pmatrix + Fmatrix + Cmatrix
			
			if (response=="lambda") rc2 <- Re(eigen(IPM)$value[1])
			if (response=="R0") rc2 <- R0Calc(Pmatrix, Fmatrix+Cmatrix)
			if (response=="lifeExpect") rc2 <- meanLifeExpect(Pmatrix)[chosenBin]
			
			discreteTrans@sdToCont[1,j]<-discreteTrans@sdToCont[1,j] / (1 + delta)
			
			slam <- c(slam,(rc2 - rc1)/(discreteTrans@sdToCont[1,j] * delta))
			elam <- c(elam, (rc2 - rc1)/(rc1 * delta))
			nmes <- c(nmes, as.character(paste("discrete: sdToCont",dimnames(discreteTrans@sdToCont)[[2]][j])))
		}
		
		for (j in 1:length(discreteTrans@moveToDiscrete$coef)) {
			discreteTrans@moveToDiscrete$coefficients[j]<-discreteTrans@moveToDiscrete$coefficients[j] * (1 + delta)
			Pmatrix <- makeIPMPmatrix(nBigMatrix = nBigMatrix, minSize = minSize,
					maxSize = maxSize, growObj = growObj, survObj = survObj,
					chosenCov = chosenCov, discreteTrans = discreteTrans,
					integrateType = integrateType, correction = correction)
			if (!is.null(fecObj)) {
				Fmatrix <- makeIPMFmatrix(nBigMatrix = nBigMatrix, minSize = minSize,
						maxSize = maxSize, fecObj = fecObj, integrateType = integrateType,
						chosenCov = chosenCov, correction = correction, preCensus = preCensusFec, 
						survObj = postCensusSurvObjFec, growObj = postCensusGrowObjFec)
			}
			if (!is.null(clonalObj)) {
				Cmatrix <- makeIPMCmatrix(nBigMatrix = nBigMatrix, minSize = minSize,
						chosenCov = chosenCov, maxSize = maxSize, clonalObj = clonalObj,
						integrateType = integrateType, correction = correction,
						preCensus = preCensusClonal, survObj = postCensusSurvObjClonal, 
						growObj = postCensusGrowObjClonal)
			}
			IPM <- Pmatrix + Fmatrix + Cmatrix
			
			if (response=="lambda") rc2 <- Re(eigen(IPM)$value[1])
			if (response=="R0") rc2 <- R0Calc(Pmatrix, Fmatrix+Cmatrix)
			if (response=="lifeExpect") rc2 <- meanLifeExpect(Pmatrix)[chosenBin]
			
			discreteTrans@moveToDiscrete$coefficients[j]<-discreteTrans@moveToDiscrete$coefficients[j] / (1 + delta)
			
			slam <- c(slam,(rc2 - rc1)/(as.numeric(discreteTrans@moveToDiscrete$coefficients[j]) * delta))
			elam <- c(elam, (rc2 - rc1)/(rc1 * delta))
			nmes <- c(nmes, as.character(paste("discrete: moveToDiscrete",names(discreteTrans@moveToDiscrete$coefficients)[j])))
		}
	}
	
	# 4. Fecundity
	if (!is.null(fecObj)) {
		for (i in 1:length(fecObj@fitFec)) {
			for (j in 1:length(fecObj@fitFec[[i]]$coefficients)) {
				fecObj@fitFec[[i]]$coefficients[j] <- fecObj@fitFec[[i]]$coefficients[j] *
						(1 + delta)
				Pmatrix <- makeIPMPmatrix(nBigMatrix = nBigMatrix,
						minSize = minSize, maxSize = maxSize, growObj = growObj,
						chosenCov = chosenCov, survObj = survObj, discreteTrans = discreteTrans,
						integrateType = integrateType, correction = correction)
				Fmatrix <- makeIPMFmatrix(nBigMatrix = nBigMatrix,
						minSize = minSize, maxSize = maxSize, fecObj = fecObj,
						chosenCov = chosenCov, integrateType = integrateType,
						correction = correction, preCensus = preCensusFec, 
						survObj = postCensusSurvObjFec, growObj = postCensusGrowObjFec)
				if (!is.null(clonalObj)) {
					Cmatrix <- makeIPMCmatrix(nBigMatrix = nBigMatrix, minSize = minSize,
							chosenCov = chosenCov, maxSize = maxSize, clonalObj = clonalObj,
							integrateType = integrateType, correction = correction,
							preCensus = preCensusClonal, survObj = postCensusSurvObjClonal, 
							growObj = postCensusGrowObjClonal)
				}
				
				IPM <- Pmatrix + Fmatrix + Cmatrix
				
				if (response=="lambda") rc2 <- Re(eigen(IPM)$value[1])
				if (response=="R0") rc2 <- R0Calc(Pmatrix, Fmatrix+Cmatrix)
				if (response=="lifeExpect") rc2 <- meanLifeExpect(Pmatrix)[chosenBin]
				
				fecObj@fitFec[[i]]$coefficients[j] <- fecObj@fitFec[[i]]$coefficients[j]/(1 + delta)
				
				slam <- c(slam,(rc2 - rc1)/((as.numeric(fecObj@fitFec[[i]]$coefficients[j]) * delta)))
				elam <- c(elam,(rc2 - rc1)/(rc1 * delta))
				nmes <- c(nmes,paste("fecundity: func", i, names(fecObj@fitFec[[i]]$coefficients)[j]))
			}
		}
		
		chs <- which(!is.na(as.numeric(fecObj@fecConstants)), arr.ind = TRUE)
		if (length(chs) > 0) {
			for (j in 1:length(chs)) {
				fecObj@fecConstants[1,chs[j]] <- fecObj@fecConstants[1,chs[j]] * (1 + delta)
				Pmatrix <- makeIPMPmatrix(nBigMatrix = nBigMatrix,
						minSize = minSize, maxSize = maxSize, growObj = growObj,
						chosenCov = chosenCov, survObj = survObj, discreteTrans = discreteTrans,
						integrateType = integrateType, correction = correction)
				Fmatrix <- makeIPMFmatrix(nBigMatrix = nBigMatrix,
						minSize = minSize, maxSize = maxSize, fecObj = fecObj,
						chosenCov = chosenCov, integrateType = integrateType,
						correction = correction, preCensus = preCensusFec, 
						survObj = postCensusSurvObjFec, growObj = postCensusGrowObjFec)
				if (!is.null(clonalObj)) {
					Cmatrix <- makeIPMCmatrix(nBigMatrix = nBigMatrix, minSize = minSize,
							chosenCov = chosenCov, maxSize = maxSize, clonalObj = clonalObj,
							integrateType = integrateType, correction = correction,
							preCensus = preCensusClonal, survObj = postCensusSurvObjClonal, 
							growObj = postCensusGrowObjClonal)
				}
				
				IPM <- Pmatrix + Fmatrix + Cmatrix
				
				if (response=="lambda") rc2 <- Re(eigen(IPM)$value[1])
				if (response=="R0") rc2 <- R0Calc(Pmatrix, Fmatrix+Cmatrix)
				if (response=="lifeExpect") rc2 <- meanLifeExpect(Pmatrix)[chosenBin]
				
				fecObj@fecConstants[1, chs[j]] <- fecObj@fecConstants[1,chs[j]]/(1 + delta)
				
				slam <- c(slam,(rc2 - rc1)/(as.numeric(fecObj@fecConstants[1,chs[j]] * delta)))
				elam <- c(elam,(rc2 - rc1)/(rc1 * delta))
				nmes <- c(nmes,paste("fecundity: constant",names(fecObj@fecConstants)[chs[j]]))
			}
		}
		
		if (max(fecObj@offspringSplitter)<1) {
			for (j in which(fecObj@offspringSplitter>0)) {
				fecObj@offspringSplitter[j] <- fecObj@offspringSplitter[j] * (1 + delta)
				Pmatrix <- makeIPMPmatrix(nBigMatrix = nBigMatrix,
						minSize = minSize, maxSize = maxSize, growObj = growObj,
						chosenCov = chosenCov, survObj = survObj, discreteTrans = discreteTrans,
						integrateType = integrateType, correction = correction)
				Fmatrix <- makeIPMFmatrix(nBigMatrix = nBigMatrix,
						minSize = minSize, maxSize = maxSize, fecObj = fecObj,
						chosenCov = chosenCov, integrateType = integrateType,
						correction = correction, preCensus = preCensusFec, 
						survObj = postCensusSurvObjFec, growObj = postCensusGrowObjFec)
				if (!is.null(clonalObj)) {
					Cmatrix <- makeIPMCmatrix(nBigMatrix = nBigMatrix, minSize = minSize,
							chosenCov = chosenCov, maxSize = maxSize, clonalObj = clonalObj,
							integrateType = integrateType, correction = correction,
							preCensus = preCensusClonal, survObj = postCensusSurvObjClonal, 
							growObj = postCensusGrowObjClonal)
				}
				
				IPM <- Pmatrix + Fmatrix + Cmatrix
				
				if (response=="lambda") rc2 <- Re(eigen(IPM)$value[1])
				if (response=="R0") rc2 <- R0Calc(Pmatrix, Fmatrix+Cmatrix)
				if (response=="lifeExpect") rc2 <- meanLifeExpect(Pmatrix)[chosenBin]
				
				fecObj@offspringSplitter[j] <- fecObj@offspringSplitter[j] / (1 + delta)
				
				slam <- c(slam,(rc2 - rc1)/(as.numeric(fecObj@offspringSplitter[j] * delta)))
				elam <- c(elam,(rc2 - rc1)/(rc1 * delta))
				nmes <- c(nmes,paste("fecundity: offspringSplitter",names(fecObj@offspringSplitter[j])))
			}
		}
		
		chs <- which(!is.na(as.numeric(fecObj@fecByDiscrete)), arr.ind = TRUE)
		if (length(chs) > 0) {
			for (j in 1:length(chs)) {
				fecObj@fecByDiscrete[1,chs[j]] <- fecObj@fecByDiscrete[1,chs[j]] * (1 + delta)
				Pmatrix <- makeIPMPmatrix(nBigMatrix = nBigMatrix,
						minSize = minSize, maxSize = maxSize, growObj = growObj,
						chosenCov = chosenCov, survObj = survObj, discreteTrans = discreteTrans,
						integrateType = integrateType, correction = correction)
				Fmatrix <- makeIPMFmatrix(nBigMatrix = nBigMatrix,
						minSize = minSize, maxSize = maxSize, fecObj = fecObj,
						chosenCov = chosenCov, integrateType = integrateType,
						correction = correction, preCensus = preCensusFec, 
						survObj = postCensusSurvObjFec, growObj = postCensusGrowObjFec)
				if (!is.null(clonalObj)) {
					Cmatrix <- makeIPMCmatrix(nBigMatrix = nBigMatrix, minSize = minSize,
							chosenCov = chosenCov, maxSize = maxSize, clonalObj = clonalObj,
							integrateType = integrateType, correction = correction,
							preCensus = preCensusClonal, survObj = postCensusSurvObjClonal, 
							growObj = postCensusGrowObjClonal)
				}
				
				IPM <- Pmatrix + Fmatrix + Cmatrix
				
				if (response=="lambda") rc2 <- Re(eigen(IPM)$value[1])
				if (response=="R0") rc2 <- R0Calc(Pmatrix, Fmatrix+Cmatrix)
				if (response=="lifeExpect") rc2 <- meanLifeExpect(Pmatrix)[chosenBin]
				
				fecObj@fecByDiscrete[1,chs[j]] <- fecObj@fecByDiscrete[1,chs[j]] / (1 + delta)
				
				slam <- c(slam,(rc2 - rc1)/(as.numeric(fecObj@fecByDiscrete[1,chs[j]] * delta)))
				elam <- c(elam,(rc2 - rc1)/(rc1 * delta))
				nmes <- c(nmes,paste("fecundity: fecByDiscrete",names(fecObj@fecByDiscrete)[chs[j]]))
			}
		}
		
		if (class(fecObj@offspringRel)=="lm") {
			for (j in 1:length(fecObj@offspringRel$coeff)) {
				fecObj@offspringRel$coefficients[j] <- fecObj@offspringRel$coefficients[j] * (1 + delta)
				Pmatrix <- makeIPMPmatrix(nBigMatrix = nBigMatrix,
						minSize = minSize, maxSize = maxSize, growObj = growObj,
						chosenCov = chosenCov, survObj = survObj, discreteTrans = discreteTrans,
						integrateType = integrateType, correction = correction)
				Fmatrix <- makeIPMFmatrix(nBigMatrix = nBigMatrix,
						minSize = minSize, maxSize = maxSize, fecObj = fecObj,
						integrateType = integrateType, correction = correction,
						chosenCov = chosenCov, preCensus = preCensusFec, 
						survObj = postCensusSurvObjFec, growObj = postCensusGrowObjFec)
				if (!is.null(clonalObj)) {
					Cmatrix <- makeIPMCmatrix(nBigMatrix = nBigMatrix, minSize = minSize,
							chosenCov = chosenCov, maxSize = maxSize, clonalObj = clonalObj,
							integrateType = integrateType, correction = correction,
							preCensus = preCensusClonal, survObj = postCensusSurvObjClonal, 
							growObj = postCensusGrowObjClonal)
				}
				
				IPM <- Pmatrix + Fmatrix + Cmatrix
				
				if (response=="lambda") rc2 <- Re(eigen(IPM)$value[1])
				if (response=="R0") rc2 <- R0Calc(Pmatrix, Fmatrix+Cmatrix)
				if (response=="lifeExpect") rc2 <- meanLifeExpect(Pmatrix)[chosenBin]
				
				fecObj@offspringRel$coefficients[j] <- fecObj@offspringRel$coefficients[j]/(1 + delta)
				
				slam <- c(slam,(rc2 - rc1)/(as.numeric(fecObj@offspringRel$coefficients[j]) *delta))
				elam <- c(elam,(rc2 - rc1)/(rc1 *delta))
				nmes <- c(nmes, as.character(paste("fecundity: offspring rel ",names(fecObj@offspringRel$coeff)[j])))
			}
			
			fecObj@sdOffspringSize <- fecObj@sdOffspringSize * (1 + delta)
			Pmatrix <- makeIPMPmatrix(nBigMatrix = nBigMatrix, minSize = minSize,
					maxSize = maxSize, growObj = growObj, chosenCov = chosenCov,
					survObj = survObj, discreteTrans = discreteTrans, integrateType = integrateType,
					correction = correction)
			Fmatrix <- makeIPMFmatrix(nBigMatrix = nBigMatrix, minSize = minSize,
					maxSize = maxSize, fecObj = fecObj, chosenCov = chosenCov,
					integrateType = integrateType, correction = correction,
					preCensus = preCensusFec, 
					survObj = postCensusSurvObjFec, growObj = postCensusGrowObjFec)
			if (!is.null(clonalObj)) {
				Cmatrix <- makeIPMCmatrix(nBigMatrix = nBigMatrix, minSize = minSize,
						chosenCov = chosenCov, maxSize = maxSize, clonalObj = clonalObj,
						integrateType = integrateType, correction = correction,
						preCensus = preCensusClonal, survObj = postCensusSurvObjClonal, 
						growObj = postCensusGrowObjClonal)
			}
			
			IPM <- Pmatrix + Fmatrix + Cmatrix
			
			if (response=="lambda") rc2 <- Re(eigen(IPM)$value[1])
			if (response=="R0") rc2 <- R0Calc(Pmatrix, Fmatrix+Cmatrix)
			if (response=="lifeExpect") rc2 <- meanLifeExpect(Pmatrix)[chosenBin]
			
			fecObj@sdOffspringSize <- fecObj@sdOffspringSize/(1 + delta)
			
			slam <- c(slam,(rc2 - rc1)/(fecObj@sdOffspringSize * delta))
			elam <- c(elam,(rc2 - rc1)/(rc1 * delta))
			nmes <- c(nmes, "fecundity: sd offspring size")
		}
	}
	
# 5. Clonality
	if (!is.null(clonalObj)) {
		for (i in 1:length(clonalObj@fitFec)) {
			for (j in 1:length(clonalObj@fitFec[[i]]$coefficients)) {
				clonalObj@fitFec[[i]]$coefficients[j] <- clonalObj@fitFec[[i]]$coefficients[j] *
						(1 + delta)
				Pmatrix <- makeIPMPmatrix(nBigMatrix = nBigMatrix,
						minSize = minSize, maxSize = maxSize, growObj = growObj,
						chosenCov = chosenCov, survObj = survObj, discreteTrans = discreteTrans,
						integrateType = integrateType, correction = correction)
				if (!is.null(fecObj)) {
					Fmatrix <- makeIPMFmatrix(nBigMatrix = nBigMatrix,
							minSize = minSize, maxSize = maxSize, fecObj = fecObj,
							chosenCov = chosenCov, integrateType = integrateType,
							correction = correction, preCensus = preCensusFec, 
							survObj = postCensusSurvObjFec, growObj = postCensusGrowObjFec)
				}
				Cmatrix <- makeIPMCmatrix(nBigMatrix = nBigMatrix, minSize = minSize,
						chosenCov = chosenCov, maxSize = maxSize, clonalObj = clonalObj,
						integrateType = integrateType, correction = correction,
						preCensus = preCensusClonal, survObj = postCensusSurvObjClonal, 
						growObj = postCensusGrowObjClonal)
				
				IPM <- Pmatrix + Fmatrix + Cmatrix
				
				if (response=="lambda") rc2 <- Re(eigen(IPM)$value[1])
				if (response=="R0") rc2 <- R0Calc(Pmatrix, Fmatrix+Cmatrix)
				if (response=="lifeExpect") rc2 <- meanLifeExpect(Pmatrix)[chosenBin]
				
				clonalObj@fitFec[[i]]$coefficients[j] <- clonalObj@fitFec[[i]]$coefficients[j]/(1 + delta)
				
				slam <- c(slam,(rc2 - rc1)/((as.numeric(clonalObj@fitFec[[i]]$coefficients[j]) * delta)))
				elam <- c(elam,(rc2 - rc1)/(rc1 * delta))
				nmes <- c(nmes,paste("clonality: func", i, names(clonalObj@fitFec[[i]]$coefficients)[j]))
			}
		}
		
		chs <- which(!is.na(as.numeric(clonalObj@fecConstants)), arr.ind = TRUE)
		if (length(chs) > 0) {
			for (j in 1:length(chs)) {
				clonalObj@fecConstants[1,chs[j]] <- clonalObj@fecConstants[1,chs[j]] * (1 + delta)
				Pmatrix <- makeIPMPmatrix(nBigMatrix = nBigMatrix,
						minSize = minSize, maxSize = maxSize, growObj = growObj,
						chosenCov = chosenCov, survObj = survObj, discreteTrans = discreteTrans,
						integrateType = integrateType, correction = correction)
				if (!is.null(fecObj)) {
					Fmatrix <- makeIPMFmatrix(nBigMatrix = nBigMatrix,
							minSize = minSize, maxSize = maxSize, fecObj = fecObj,
							chosenCov = chosenCov, integrateType = integrateType,
							correction = correction, preCensus = preCensusFec, 
							survObj = postCensusSurvObjFec, growObj = postCensusGrowObjFec)
				}
				Cmatrix <- makeIPMCmatrix(nBigMatrix = nBigMatrix, minSize = minSize,
						chosenCov = chosenCov, maxSize = maxSize, clonalObj = clonalObj,
						integrateType = integrateType, correction = correction,
						preCensus = preCensusClonal, survObj = postCensusSurvObjClonal, 
						growObj = postCensusGrowObjClonal)
				
				IPM <- Pmatrix + Fmatrix + Cmatrix
				
				if (response=="lambda") rc2 <- Re(eigen(IPM)$value[1])
				if (response=="R0") rc2 <- R0Calc(Pmatrix, Fmatrix+Cmatrix)
				if (response=="lifeExpect") rc2 <- meanLifeExpect(Pmatrix)[chosenBin]
				
				clonalObj@fecConstants[1, chs[j]] <- clonalObj@fecConstants[1,chs[j]]/(1 + delta)
				
				slam <- c(slam,(rc2 - rc1)/(as.numeric(clonalObj@fecConstants[1,chs[j]] * delta)))
				elam <- c(elam,(rc2 - rc1)/(rc1 * delta))
				nmes <- c(nmes,paste("clonality: constant",names(clonalObj@fecConstants)[chs[j]]))
			}
		}
		
		if (max(clonalObj@offspringSplitter)<1) {
			for (j in which(clonalObj@offspringSplitter>0)) {
				clonalObj@offspringSplitter[j] <- clonalObj@offspringSplitter[j] * (1 + delta)
				Pmatrix <- makeIPMPmatrix(nBigMatrix = nBigMatrix,
						minSize = minSize, maxSize = maxSize, growObj = growObj,
						chosenCov = chosenCov, survObj = survObj, discreteTrans = discreteTrans,
						integrateType = integrateType, correction = correction)
				if (!is.null(fecObj)) {
					Fmatrix <- makeIPMFmatrix(nBigMatrix = nBigMatrix,
							minSize = minSize, maxSize = maxSize, fecObj = fecObj,
							chosenCov = chosenCov, integrateType = integrateType,
							correction = correction, preCensus = preCensusFec, 
							survObj = postCensusSurvObjFec, growObj = postCensusGrowObjFec)
				}
				Cmatrix <- makeIPMCmatrix(nBigMatrix = nBigMatrix, minSize = minSize,
						chosenCov = chosenCov, maxSize = maxSize, clonalObj = clonalObj,
						integrateType = integrateType, correction = correction,
						preCensus = preCensusClonal, survObj = postCensusSurvObjClonal, 
						growObj = postCensusGrowObjClonal)
				
				IPM <- Pmatrix + Fmatrix + Cmatrix
				
				if (response=="lambda") rc2 <- Re(eigen(IPM)$value[1])
				if (response=="R0") rc2 <- R0Calc(Pmatrix, Fmatrix+Cmatrix)
				if (response=="lifeExpect") rc2 <- meanLifeExpect(Pmatrix)[chosenBin]
				
				clonalObj@offspringSplitter[j] <- clonalObj@offspringSplitter[j] / (1 + delta)
				
				slam <- c(slam,(rc2 - rc1)/(as.numeric(clonalObj@offspringSplitter[j] * delta)))
				elam <- c(elam,(rc2 - rc1)/(rc1 * delta))
				nmes <- c(nmes,paste("clonality: offspringSplitter",names(clonalObj@offspringSplitter[j])))
			}
		}
		
		chs <- which(!is.na(as.numeric(clonalObj@fecByDiscrete)), arr.ind = TRUE)
		if (length(chs) > 0) {
			for (j in 1:length(chs)) {
				clonalObj@fecByDiscrete[1,chs[j]] <- clonalObj@fecByDiscrete[1,chs[j]] * (1 + delta)
				Pmatrix <- makeIPMPmatrix(nBigMatrix = nBigMatrix,
						minSize = minSize, maxSize = maxSize, growObj = growObj,
						chosenCov = chosenCov, survObj = survObj, discreteTrans = discreteTrans,
						integrateType = integrateType, correction = correction)
				if (!is.null(fecObj)) {
					Fmatrix <- makeIPMFmatrix(nBigMatrix = nBigMatrix,
							minSize = minSize, maxSize = maxSize, fecObj = fecObj,
							chosenCov = chosenCov, integrateType = integrateType,
							correction = correction, preCensus = preCensusFec, 
							survObj = postCensusSurvObjFec, growObj = postCensusGrowObjFec)
				}
				Cmatrix <- makeIPMCmatrix(nBigMatrix = nBigMatrix, minSize = minSize,
						chosenCov = chosenCov, maxSize = maxSize, clonalObj = clonalObj,
						integrateType = integrateType, correction = correction,
						preCensus = preCensusClonal, survObj = postCensusSurvObjClonal, 
						growObj = postCensusGrowObjClonal)
				
				IPM <- Pmatrix + Fmatrix + Cmatrix
				
				if (response=="lambda") rc2 <- Re(eigen(IPM)$value[1])
				if (response=="R0") rc2 <- R0Calc(Pmatrix, Fmatrix+Cmatrix)
				if (response=="lifeExpect") rc2 <- meanLifeExpect(Pmatrix)[chosenBin]
				
				clonalObj@fecByDiscrete[1,chs[j]] <- clonalObj@fecByDiscrete[1,chs[j]] / (1 + delta)
				
				slam <- c(slam,(rc2 - rc1)/(as.numeric(clonalObj@fecByDiscrete[1,chs[j]] * delta)))
				elam <- c(elam,(rc2 - rc1)/(rc1 * delta))
				nmes <- c(nmes,paste("clonality: fecByDiscrete",names(clonalObj@fecByDiscrete)[chs[j]]))
			}
		}
		
		if (class(clonalObj@offspringRel)=="lm") {
			for (j in 1:length(clonalObj@offspringRel$coeff)) {
				clonalObj@offspringRel$coefficients[j] <- clonalObj@offspringRel$coefficients[j] * (1 + delta)
				Pmatrix <- makeIPMPmatrix(nBigMatrix = nBigMatrix,
						minSize = minSize, maxSize = maxSize, growObj = growObj,
						chosenCov = chosenCov, survObj = survObj, discreteTrans = discreteTrans,
						integrateType = integrateType, correction = correction)
				if (!is.null(fecObj)) {
					Fmatrix <- makeIPMFmatrix(nBigMatrix = nBigMatrix,
							minSize = minSize, maxSize = maxSize, fecObj = fecObj,
							chosenCov = chosenCov, integrateType = integrateType,
							correction = correction, preCensus = preCensusFec, 
							survObj = postCensusSurvObjFec, growObj = postCensusGrowObjFec)
				}
				Cmatrix <- makeIPMCmatrix(nBigMatrix = nBigMatrix, minSize = minSize,
						chosenCov = chosenCov, maxSize = maxSize, clonalObj = clonalObj,
						integrateType = integrateType, correction = correction,
						preCensus = preCensusClonal, survObj = postCensusSurvObjClonal, 
						growObj = postCensusGrowObjClonal)
				
				IPM <- Pmatrix + Fmatrix + Cmatrix
				
				if (response=="lambda") rc2 <- Re(eigen(IPM)$value[1])
				if (response=="R0") rc2 <- R0Calc(Pmatrix, Fmatrix+Cmatrix)
				if (response=="lifeExpect") rc2 <- meanLifeExpect(Pmatrix)[chosenBin]
				
				clonalObj@offspringRel$coefficients[j] <- clonalObj@offspringRel$coefficients[j]/(1 + delta)
				
				slam <- c(slam,(rc2 - rc1)/(as.numeric(clonalObj@offspringRel$coefficients[j]) *delta))
				elam <- c(elam,(rc2 - rc1)/(rc1 *delta))
				nmes <- c(nmes, as.character(paste("clonality: offspring rel ",names(clonalObj@offspringRel$coeff)[j])))
			}
			
			clonalObj@sdOffspringSize <- clonalObj@sdOffspringSize * (1 + delta)
			Pmatrix <- makeIPMPmatrix(nBigMatrix = nBigMatrix, minSize = minSize,
					maxSize = maxSize, growObj = growObj, chosenCov = chosenCov,
					survObj = survObj, discreteTrans = discreteTrans, integrateType = integrateType,
					correction = correction)
			if (!is.null(fecObj)) {
				Fmatrix <- makeIPMFmatrix(nBigMatrix = nBigMatrix,
						minSize = minSize, maxSize = maxSize, fecObj = fecObj,
						chosenCov = chosenCov, integrateType = integrateType,
						correction = correction, preCensus = preCensusFec, 
						survObj = postCensusSurvObjFec, growObj = postCensusGrowObjFec)
			}
			Cmatrix <- makeIPMCmatrix(nBigMatrix = nBigMatrix, minSize = minSize,
					chosenCov = chosenCov, maxSize = maxSize, clonalObj = clonalObj,
					integrateType = integrateType, correction = correction,
					preCensus = preCensusClonal, survObj = postCensusSurvObjClonal, 
					growObj = postCensusGrowObjClonal)
			
			IPM <- Pmatrix + Fmatrix + Cmatrix
			
			if (response=="lambda") rc2 <- Re(eigen(IPM)$value[1])
			if (response=="R0") rc2 <- R0Calc(Pmatrix, Fmatrix+Cmatrix)
			if (response=="lifeExpect") rc2 <- meanLifeExpect(Pmatrix)[chosenBin]
			
			clonalObj@sdOffspringSize <- clonalObj@sdOffspringSize/(1 + delta)
			
			slam <- c(slam,(rc2 - rc1)/(clonalObj@sdOffspringSize * delta))
			elam <- c(elam,(rc2 - rc1)/(rc1 * delta))
			nmes <- c(nmes, "clonality: sd offspring size")
		}
	}
	
	names(slam) <- nmes
	names(elam) <- nmes
	
	return(list(sens = slam, elas = elam))
}


### FUNCTIONS FOR EXTRACTING STOCH GROWTH RATE ########################
# =============================================================================
# =============================================================================
# Generic approach to get stoch rate
# by sampling list IPM
#
# Parameters - listIPMmatrix - list of IPMs corresponding to different year types
#            - nRunIn - the burnin before establishing lambda_s
#            - tMax - the total time-horizon for getting lambda_s
#
# Returns lambda_s (no density dependence)

stochGrowthRateSampleList <- function(nRunIn,tMax,listIPMmatrix=NULL,
					listPmatrix=NULL, listFmatrix=NULL, seedList = NULL,
					densDep=FALSE){
	require(MASS)
	
	if (densDep & (is.null(listPmatrix) | is.null(listFmatrix))){
		stop("Require listPmatrix & listFmatrix for densDep=TRUE")
	}
		
	if (!densDep & is.null(listIPMmatrix)) {
		stop("Require listIPMmatrix for densDep=TRUE")		
	}		
	
	if (densDep) {
		nmatrices <- length(listPmatrix)
		nBigMatrix <- length(listPmatrix[[1]][,1]) 
	} else { 
		nmatrices <- length(listIPMmatrix)
		nBigMatrix <- length(listIPMmatrix[[1]][,1])
	}
	
	nt<-rep(1,nBigMatrix)
	Rt<-rep(NA,tMax)
	
	pEst <- 1
	
	for (t in 1:tMax) {
		year.type <- sample(1:nmatrices,size=1,replace=FALSE)
		if (densDep) { 
			nseeds <- sum(listFmatrix[[year.type]]%*%nt)
			pEst <- min(seedList[min(year.type+1,nmatrices)]/nseeds,1)
			nt1 <- (listPmatrix[[year.type]]+pEst*listFmatrix[[year.type]])%*% nt
			sum.nt1 <- sum(nt1)
			Rt[t] <- log(sum.nt1)-log(sum(nt))
			nt <- nt1
		} else {
			nt1<-listIPMmatrix[[year.type]] %*% nt	
			sum.nt1<-sum(nt1)
			Rt[t]<-log(sum.nt1)
			nt<-nt1/sum.nt1
		}		
	}
		res <- mean(Rt[nRunIn:tMax],na.rm=TRUE)
	return(res)
}

# =============================================================================
# =============================================================================
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
		integrateType="midpoint",correction="none", 
		trackStruct=FALSE, plot=FALSE,...){
	require(MASS)
	
	
	nt<-rep(1,nBigMatrix)
	Rt<-rep(NA,tMax)
	fecObj@fecConstants[is.na(fecObj@fecConstants)] <- 1 
	tmp.fecObj <- fecObj
	
	#print("fec.const start")
	#print(fecObj@fecConstants)
	
	#density dep in seedling establishment 
	if (sum(nMicrosites)>0) { dd <- TRUE; seeds <- 10000 } else { dd <- FALSE; p.est <- 1}
	if (trackStruct) rc <- matrix(NA,nBigMatrix,tMax)
	
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
		
		tpS <- makeIPMPmatrix(nBigMatrix = nBigMatrix, minSize = minSize,
				maxSize = maxSize, chosenCov = covariatet,
				growObj = growthObj, survObj = survObj,
				integrateType=integrateType, correction=correction)
		
		tpF <- makeIPMFmatrix(nBigMatrix = nBigMatrix, minSize = minSize,
				maxSize = maxSize, chosenCov = covariatet,
				fecObj = tmp.fecObj,
				integrateType=integrateType, correction=correction)
		
		#total seeds for next year 
		if (dd) seeds <- sum(tpF%*%nt)
		
		IPM.here <- p.est*tpF@.Data+tpS@.Data
		nt1<-IPM.here %*% nt	
		sum.nt1<-sum(nt1)
		
		if (!trackStruct){
			if (!dd) { 
				Rt[t]<-log(sum.nt1)
				nt<-nt1/sum.nt1
			} else {
				Rt[t]<-log(sum.nt1)-log(sum(nt))
				nt<-nt1	
			}} else {
			nt<-nt1	
			rc[,t] <- nt1
		}
		
		
		
	}
	
	if (trackStruct & plot){
		.plotResultsStochStruct(tVals=1:tMax, meshpoints=tpS@meshpoints,
				st=rc,covTest=covariate, nRunIn = nRunIn,  ...) 	
	}
	
	if (!trackStruct) {res <- mean(Rt[nRunIn:tMax],na.rm=TRUE); return(list(Rt=res))}
	if (trackStruct) return(list(rc=rc,IPM.here=tpS))
	
}






## FUNCTIONS FOR SENSITIVITY AND ELASTICITY ###################################

# =============================================================================
# =============================================================================
#parameters - an IPM (with full survival and fecundity complement)
# returns - the sensitivity of every transition for pop growth
sens<-function(A) {
	w<-Re(eigen(A)$vectors[,1]); 
	v<-Re(eigen(t(A))$vectors[,1]);
	vw<-sum(v*w);
	s<-outer(v,w)
	return(s/vw); 
}   

# =============================================================================
# =============================================================================
#parameters - an IPM (with full survival and fecundity complement)
# returns - the elasticity of every transition for pop growth
elas<-function(A) {
	s<-sens(A)
	lam<-Re(eigen(A)$values[1]);
	return((s*A)/lam);
}



## FUNCTIONS TO BE DEPRECATED #################################################
# =============================================================================
# =============================================================================
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
	
	print("warning: deprecated! use makeIPMFmatrix instead")
	
 rc <- makeIPMFmatrix(fecObj=fecObj,
		 nEnvClass = nEnvClass,
		 nBigMatrix = nBigMatrix,
		 minSize = minSize,
		 maxSize = maxSize,
		 chosenCov = chosenCov,
		 integrateType=integrateType,
		 correction=correction,
		 preCensus=preCensus,
		 survObj=survObj,
		 growObj=growObj)
	 
	 return(rc)

}

# =============================================================================
# =============================================================================
createIPMCmatrix <- function(clonalObj,
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
	
	print("warning: deprecated! use makeIPMCmatrix instead")
	
	rc <- makeIPMFmatrix(fecObj=clonalObj,
			nEnvClass = nEnvClass,
			nBigMatrix = nBigMatrix,
			minSize = minSize,
			maxSize = maxSize,
			chosenCov = chosenCov,
			integrateType=integrateType,
			correction=correction,
			preCensus=preCensus,
			survObj=survObj,
			growObj=growObj)
		
		return(rc)
	
}

# =============================================================================
# =============================================================================
createIPMPmatrix <- function (nEnvClass = 1, 
		nBigMatrix = 50, minSize = -1, maxSize = 50, 
		chosenCov = data.frame(covariate = 1), growObj, survObj, 
		discreteTrans = 1, integrateType = "midpoint", correction = "none") {
	
	print("warning: deprecated! use makeIPMPmatrix instead")
	
	rc <- makeIPMPmatrix(nEnvClass = nEnvClass,
			nBigMatrix = nBigMatrix,
			minSize = minSize,
			maxSize = maxSize,
			chosenCov = chosenCov,
			growObj=growObj, survObj=survObj,
			discreteTrans =discreteTrans,
			integrateType=integrateType,
			correction=correction)
	
	return(rc)
}


# =============================================================================
# =============================================================================
createIntegerPmatrix <- function (nEnvClass = 1, 
		meshpoints=1:20,
		chosenCov = data.frame(covariate = 1), 
		growObj, survObj, 
		discreteTrans = 1){
	
	print("warning: deprecated! use makeIntegerPmatrix instead")
	
	rc <- makeIntegerPmatrix(nEnvClass = nEnvClass,meshpoints=meshpoints,
			chosenCov = chosenCov,growObj=growObj, survObj=survObj,discreteTrans =discreteTrans)
	
	return(rc)
	
}

createIntegerFmatrix <- function(fecObj,
		nEnvClass = 1,
		meshpoints=1:20,
		chosenCov = data.frame(covariate=1),
		preCensus=TRUE,
		survObj=NULL,
		growObj=NULL){
	
	print("warning: deprecated! use makeIntegerFmatrix instead")
	
	rc <- makeIntegerFmatrix(fecObj=fecObj,
			nEnvClass = nEnvClass,meshpoints=meshpoints,
			chosenCov = chosenCov,preCensus=preCensus, 
			growObj=growObj, survObj=survObj)
	
	
	return(rc)
	
}

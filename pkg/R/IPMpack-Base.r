

##can you do this to get rid of teh warnings? check 
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
		representation(fit="lm"))

setClass("growthObjMultiCov",
		representation(fit="lm"))

# Create a generic growth object with normal errors on increment
setClass("growthObj.incr",
		representation(fit="lm"))

setClass("growthObjMultiCov.incr",
		representation(fit="lm"))

# Create a generic growth object with truncated normal errors on increment
setClass("growthObj.truncincr",
		representation(fit=c("numeric")))

# Create a generic growth object with log normal errors on increment
setClass("growthObj.logincr",
		representation(fit="lm"))

setClass("growthObjMultiCov.logincr",
		representation(fit="lm"))

# Create a generic growth object with declining errors 
setClass("growthObj.declinevar",
		representation(fit="gls"))

setClass("growthObjMultiCov.declinevar",
		representation(fit="gls"))

# Create a generic growth object with declining errors for increment
setClass("growthObj.incr.declinevar",
		representation(fit="gls"))

setClass("growthObjMultiCov.incr.declinevar",
		representation(fit="gls"))

# Create a generic growth object with declining errors for logincrement
setClass("growthObj.logincr.declinevar",
		representation(fit="gls"))

setClass("growthObjMultiCov.logincr.declinevar",
		representation(fit="gls"))


# Create a generic growth object containing the Hossfeld parameters 
setClass("growthObj.Hossfeld",
		representation(paras="numeric",
				sd="numeric", 
				logLik="numeric"))

## SURVIVAL OBJECTS ##
# Create a generic survival object
setClass("survObj",
		representation(fit="glm"))

setClass("survObjMultiCov",
		representation(fit="glm"))

# Create a generic survival object with a log covariate, and 5 years
setClass("survObjLog.multiyear",
		representation(fit="glm"))

# Create a generic survival object for use where over-dispersion
# modeled, using Diggles approximate correction for the transform
setClass("survObjOverDisp",
		representation(fit="glm"))



## FECUNDITY OBJECTS ##
# Create a generic fecundity object
setClass("fecObj",
		representation(fit.fec = "list",
				fec.constants = "numeric",
				OffspringSplitter = "data.frame",
				mean.offspring.size = "numeric",
				var.offspring.size = "numeric",
				fec.by.discrete = "data.frame",
				Transform = "character")
)

# Create a generic fecundity object for multiple covariates
setClass("fecObjMultiCov",
		representation(fit.fec = "list",
				fec.constants = "numeric",
				OffspringSplitter = "data.frame",
				mean.offspring.size = "numeric",
				var.offspring.size = "numeric",
				fec.by.discrete = "data.frame",
				Transform = "character")
)



## DISCRETE TRANSITION MATRIX OBJECTS ##
# Create a generic discrete transition matrix object
setClass("discreteTrans",
		representation(nclasses="numeric",
				discrete.trans="matrix",
				discrete.surv="matrix",
				mean.to.cont="matrix",
				sd.to.cont="matrix",
				distrib.to.discrete="matrix",
				surv.to.discrete="glm"))





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
		c("numeric","numeric","survObj"),
		function(size,cov,survObj){
			newd <- data.frame(size=size,size2=size^2,size3=size^3,
					covariate=as.factor(rep(cov,length(size)))) 
			if (length(grep("logsize",
							survObj@fit$formula))>0) newd$logsize=log(size)
			if (length(grep("logsize2",
							survObj@fit$formula))>0) newd$logsize2=(log(size))^2
			
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
		c("numeric","numeric","survObjOverDisp"),
		function(size,cov,survObj){
			newd <- data.frame(size=size,size2=size^2,size3=size^3,
					covariate=as.factor(rep(cov,length(size)))) 
			if (length(grep("logsize",
							survObj@fit$formula))>0) newd$logsize=log(size)
			if (length(grep("logsize2",
							survObj@fit$formula))>0) newd$logsize2=(log(size))^2
			
			u <- predict(survObj@fit,newd,type="link")
			c2 <- ((16 * sqrt(3))/(15 * pi))^2  #from MCMCglmm course notes, search for c2
			u <- logit(u/sqrt(1+c2)) 
			return(u);
		})


#Method to obtain probability of survival
# logistic regression on size with potentially many
# continuous covariates
#
#Parameters -
#   size = current size (vector)
#   cov = a vector of potentially many covariates
#   survObj = a survival object, containig e.g. a glm fitted to data
#
#Returns -
#  survival probability for given sizes and covariate level
setMethod("surv", 
		c("numeric","data.frame","survObjMultiCov"),
		function(size,cov,survObj){
			newd <- cov
			newd[2:length(size),] <- rep(as.numeric(cov[1,]), each=(length(size)-1))
			newd$size <- size
			newd$size2 <- size^2
			newd$size3 <- size^3
			
			#print(survObj)
			
			if (length(grep("logsize",
							survObj@fit$formula))>0) newd$logsize=log(size)
			if (length(grep("logsize2",
							survObj@fit$formula))>0) newd$logsize2=(log(size))^2
			u <- predict(survObj@fit,newd,type="response")
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
		c("numeric","numeric","numeric","growthObj"),
		function(size,sizeNext,cov,growthObj){
			newd <- data.frame(size=size,size2=size^2,size3=size^3,
					covariate=as.factor(rep(cov,length(size))))
			#print(head(newd))
			if (length(grep("logsize",
							growthObj@fit$formula))>0) { newd$logsize <- log(size)}
			
			mux <- predict(growthObj@fit,newd,type="response")
			sigmax <- summary(growthObj@fit)$sigma
			u <- dnorm(sizeNext,mux,sigmax,log=FALSE)  
			return(u);
		})

# Same for many covariates
setMethod("growth", 
		c("numeric","numeric","data.frame","growthObjMultiCov"),
		function(size,sizeNext,cov,growthObj){
			newd <- cov
			newd[2:length(size),] <- rep(as.numeric(cov[1,]), each=(length(size)-1))
			newd$size <- size
			newd$size2 <- size^2
			newd$size3 <- size^3
			
			if (length(grep("logsize",
							growthObj@fit$formula))>0) newd$logsize=log(size)
			if (length(grep("logsize2",
							growthObj@fit$formula))>0) newd$logsize2=(log(size))^2
			
			mux <- predict(growthObj@fit,newd,type="response")
			sigmax <- summary(growthObj@fit)$sigma
			u <- dnorm(sizeNext,mux,sigmax,log=FALSE)  
			return(u);
		})


# growth for predicting next incr 
setMethod("growth", 
		c("numeric","numeric","numeric","growthObj.incr"),
		function(size,sizeNext,cov,growthObj){
			newd <- data.frame(size=size,size2=size^2,size3=size^3,
					covariate=as.factor(rep(cov,length(size))))
			
			if (length(grep("logsize",
							growthObj@fit$formula))>0) newd$logsize=log(size)
			
			mux <- predict(growthObj@fit,newd,type="response")
			
			#print(mux)
			
			sigmax <- summary(growthObj@fit)$sigma
			u <- dnorm(sizeNext,size+mux,sigmax,log=FALSE)  
			return(u); 
		})



# growth for predicting next truncated incr 
setMethod("growth", 
		c("numeric","numeric","numeric","growthObj.truncincr"),
		function(size,sizeNext,cov,growthObj){
			require(truncnorm)
			newd <- data.frame(size=size,size2=size^2,size3=size^3,
					logsize=log(size),logsize2=(log(size))^2,
					covariate=as.factor(rep(cov,length(size))))
			
			
			m1 <-match(names(growthObj@fit),colnames(newd)); m1 <- m1[!is.na(m1)]
			mux <- colSums(growthObj@fit[1]+t(newd[,m1])*growthObj@fit[2:(length(growthObj@fit)-1)])
			#print(range(mux))
			sigmax <- exp(growthObj@fit["logSigma"])
			u <- dtruncnorm(sizeNext,a=size,b=Inf,mean=size+mux,sd=sigmax)  
			return(u); 
		})


# Same for many covariates on increment
setMethod("growth", 
		c("numeric","numeric","numeric","growthObjMultiCov.incr"),
		function(size,sizeNext,cov,growthObj){
			newd <- cov
			newd[2:length(size),] <- rep(as.numeric(cov[1,]), each=(length(size)-1))
			newd$size <- size
			newd$size2 <- size^2
			newd$size3 <- size^3
			
			if (length(grep("logsize",
							growthObj@fit$formula))>0) newd$logsize=log(size)
			if (length(grep("logsize2",
							growthObj@fit$formula))>0) newd$logsize2=(log(size))^2
			
			mux <- predict(growthObj@fit,newd,type="response")
			sigmax <- summary(growthObj@fit)$sigma
			u <- dnorm(sizeNext,size+mux,sigmax,log=FALSE)  
			return(u);
		})


# growth for predicting next logincr with a polynomial or log
setMethod("growth", 
		c("numeric","numeric","numeric","growthObj.logincr"),
		function(size,sizeNext,cov,growthObj){
			newd <- data.frame(size=size,size2=size^2,size3=size^3,
					covariate=as.factor(rep(cov,length(size))))
			
			if (length(grep("logsize",
							growthObj@fit$formula))>0) newd$logsize=log(size)
			
			mux <- predict(growthObj@fit,newd,type="response")
			sigmax <- summary(growthObj@fit)$sigma
			u <- dlnorm(sizeNext-size,mux,sigmax,log=FALSE)  
			return(u);
		})


# Same for many covariates on log increment
setMethod("growth", 
		c("numeric","numeric","numeric","growthObjMultiCov.logincr"),
		function(size,sizeNext,cov,growthObj){
			newd <- cov
			newd[2:length(size),] <- rep(as.numeric(cov[1,]), each=(length(size)-1))
			newd$size <- size
			newd$size2 <- size^2
			newd$size3 <- size^3
			
			if (length(grep("logsize",
							growthObj@fit$formula))>0) newd$logsize=log(size)
			if (length(grep("logsize2",
							growthObj@fit$formula))>0) newd$logsize2=(log(size))^2
			
			mux <- predict(growthObj@fit,newd,type="response")
			sigmax <- summary(growthObj@fit)$sigma
			u <- dlnorm(sizeNext,size+mux,sigmax,log=FALSE)  
			return(u);
		})






# Slightly alternative approach, using cumulative probs ##

# growth for predicting next size with a polynomial or log
# using pnorm (i.e. getting cumulative at boundary points and doing difference)
setMethod("growthCum", 
		c("numeric","numeric","numeric","growthObj"),
		function(size,sizeNext,cov,growthObj){
			newd <- data.frame(size=size,size2=size^2,size3=size^3,
					covariate=as.factor(rep(cov,length(size))))
			if (length(grep("logsize",
							growthObj@fit$formula))>0) newd$logsize=log(size)
			mux <- predict(growthObj@fit,newd,type="response")
			sigmax <- summary(growthObj@fit)$sigma
			u <- pnorm(sizeNext,mux,sigmax,log=FALSE)  
			return(u);
		})

# growth for predicting next incr with a polynomial or log
# using pnorm (i.e. getting cumulative at boundary points and doing difference)
setMethod("growthCum", 
		c("numeric","numeric","numeric","growthObj.incr"),
		function(size,sizeNext,cov,growthObj){
			newd <- data.frame(size=size,size2=size^2,size3=size^3,#logsize=log(size),
					covariate=as.factor(rep(cov,length(size))))
			if (length(grep("logsize",
							growthObj@fit$formula))>0) newd$logsize=log(size)
			
			mux <- predict(growthObj@fit,newd,type="response")
			sigmax <- summary(growthObj@fit)$sigma
			u <- pnorm(sizeNext,size+mux,sigmax,log=FALSE)  
			return(u); 
		})



# growth for predicting next truncated incr with cumulative 
setMethod("growthCum", 
		c("numeric","numeric","numeric","growthObj.truncincr"),
		function(size,sizeNext,cov,growthObj){
			require(truncnorm)
			newd <- data.frame(size=size,size2=size^2,size3=size^3,
					logsize=log(size),logsize2=(log(size))^2,
					covariate=as.factor(rep(cov,length(size))))
			m1 <-match(names(growthObj@fit),colnames(newd)); m1 <- m1[!is.na(m1)]
			mux <- colSums(growthObj@fit[1]+t(newd[,m1])*growthObj@fit[2:(length(growthObj@fit)-1)])
			#print(range(mux))
			sigmax <- exp(growthObj@fit["logSigma"])
			u <- ptruncnorm(sizeNext,a=size,b=Inf,mean=size+mux,sd=sigmax)  
			return(u); 
		})


# growth for predicting next logincr with a polynomial or log
# using pnorm (i.e. getting cumulative at boundary points and doing difference)
setMethod("growthCum", 
		c("numeric","numeric","numeric","growthObj.logincr"),
		function(size,sizeNext,cov,growthObj){
			newd <- data.frame(size=size,size2=size^2,size3=size^3,#logsize=log(size),
					covariate=as.factor(rep(cov,length(size))))
			if (length(grep("logsize",
							growthObj@fit$formula))>0) newd$logsize=log(size)
			
			mux <- predict(growthObj@fit,newd,type="response")
			sigmax <- summary(growthObj@fit)$sigma
			u <- plnorm(sizeNext-size,mux,sigmax,log=FALSE)  
			return(u);
		})


#Simple growth methods, using  declining variance in growth
# using pnorm (i.e. getting cumulative at boundary points and doing difference)
setMethod("growthCum", 
		c("numeric","numeric","numeric","growthObj.declinevar"),
		function(size,sizeNext,cov,growthObj){
			newd <- data.frame(size=size,size2=size^2,covariate=as.factor(rep(cov,length(size))))
			mux <- predict(growthObj@fit,newd,type="response")
			sigmax2 <- summary(growthObj@fit)$sigma^2
			var.exp.coef<-as.numeric(growthObj@fit$modelStruct$varStruct[1])
			sigmax2<-sigmax2*exp(2*(var.exp.coef*mux));
			u <- pnorm(sizeNext,mux,sqrt(sigmax2),log=FALSE)  
			return(u);
		})




#Simple growth methods, using  declining variance in growth
setMethod("growth", 
		c("numeric","numeric","numeric","growthObj.declinevar"),
		function(size,sizeNext,cov,growthObj){
			newd <- data.frame(size=size,size2=size^2,size3=size^3,covariate=as.factor(rep(cov,length(size))))
			mux <- predict(growthObj@fit,newd,type="response")
			sigmax2 <- summary(growthObj@fit)$sigma^2
			var.exp.coef<-as.numeric(growthObj@fit$modelStruct$varStruct[1])
			sigmax2<-sigmax2*exp(2*(var.exp.coef*mux));
			u <- dnorm(sizeNext,mux,sqrt(sigmax2),log=FALSE)  
			return(u);
		})




#same but with declining variance in growth on incrment
setMethod("growth", 
		c("numeric","numeric","numeric","growthObj.incr.declinevar"),
		function(size,sizeNext,cov,growthObj){
			newd <- data.frame(size=size,size2=size^2,size3=size^3,covariate=as.factor(rep(cov,length(size))))
			mux <- predict(growthObj@fit,newd,type="response")
			sigmax2 <- summary(growthObj@fit)$sigma^2
			var.exp.coef<-as.numeric(growthObj@fit$modelStruct$varStruct[1])
			sigmax2<-sigmax2*exp(2*(var.exp.coef*mux));
			u <- dnorm(sizeNext,size+mux,sqrt(sigmax2),log=FALSE)  
			return(u);
		})


## Define a new growth method for Hossfeld growth (classic midpoint rule approach)
setMethod("growth", c("numeric", "numeric", "numeric", "growthObj.Hossfeld"), 
		function(size, sizeNext, cov, growthObj) { 
			mux <- size+Hossfeld(size, growthObj@paras) 
			sigmax <- growthObj@sd 
			u <- dnorm(sizeNext, mux, sigmax, log = F) 
			return(u)
		}) 




# Method combining growth and survival for doing outer (not a generic, so that don't have to
# define for all the different classes / growth and surv take care of that)
growSurv <- function(size,sizeNext,cov,growthObj,survObj){
	growth(size,sizeNext,cov,growthObj)*surv(size,cov,survObj)
}






### CLASSES AND FUNCTIONS FOR MATRICES (ENV, TMATRIX [compound or not], FMATRIX) #################

#Class for the matrix that holds the env matrix 
setClass("envMatrix",
		representation(nEnvClass = "numeric"), #number of covariate levels
		contains="matrix")


#Class for the matrix that holds the IPM
setClass("IPM.matrix",
		representation(n.discrete = "numeric", #number of discrete classes
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
#   chosen.cov - current level of covariate - can be vector for continuous
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
create.IPM.Tmatrix <- function(nEnvClass = 1,
		nBigMatrix = 50,
		minSize = -1,
		maxSize = 50,
		chosen.cov = 1,
		growObj,
		survObj,
		discreteTrans =1,
		integrateType = "midpoint",
		correction="none") {
	# boundary points b and mesh points y
	b<-minSize+c(0:nBigMatrix)*(maxSize-minSize)/nBigMatrix;
	y<-0.5*(b[1:nBigMatrix]+b[2:(nBigMatrix+1)]);
	
	# step size for mid point rule, see equations 4 and 5
	h<-y[2]-y[1]
	
	# fill in matrix                          
	if (integrateType=="midpoint") { 
		get.matrix <- 
				t(outer(y,y,growSurv,cov=chosen.cov,
								growthObj=growObj,survObj=survObj))*h  
		
	}
	if (integrateType=="cumul") {
		get.matrix.cum <- 
				t(outer(y,b,growthCum,cov=chosen.cov,
								growthObj=growObj))
		get.matrix <- get.matrix.cum[2:(nBigMatrix+1),]-get.matrix.cum[1:nBigMatrix,]
		#fix last size
		get.matrix[nBigMatrix,nBigMatrix] <- get.matrix[nBigMatrix,nBigMatrix]+
				(1-sum(get.matrix[,nBigMatrix]))
		#put in survival
		get.matrix <- t(t(get.matrix)*surv(size=y,cov=chosen.cov,survObj=survObj))
	}
	
	#fix any integration issues reducing survival by dividing by col sums and multiply by survival
	if (correction=="constant") { 
		nvals <- colSums(get.matrix); nvals[nvals==0] <- 1
		get.matrix <- t((t(get.matrix)/nvals)*surv(size=y,cov=chosen.cov,survObj=survObj))    
	}
	
	rc <- new("IPM.matrix",
			n.discrete =0,
			nEnvClass = 1, 
			nBigMatrix = nBigMatrix,
			nrow = 1*nBigMatrix,
			ncol =1*nBigMatrix,
			meshpoints = y,
			env.index = rep(1:nEnvClass,each=nBigMatrix),
			names.discrete="")
	rc[,] <-get.matrix
	
	# In case of discrete classes, take the IPM constructed above and add discrete classes defined in discreteTrans
	if (class(discreteTrans)=="discreteTrans") {
		
		ndisc <- ncol(discreteTrans@discrete.surv)
		surv.to.discrete <- predict(discreteTrans@surv.to.discrete,data.frame(size=y,size2=(y*y)),type="response")
		cont.to.cont <- get.matrix*matrix(1-surv.to.discrete,nrow=nBigMatrix,ncol=nBigMatrix,byrow=T)
		disc.to.disc <- discreteTrans@discrete.trans[1:ndisc,1:ndisc]*matrix(c(discreteTrans@discrete.surv),nrow=ndisc,ncol=ndisc,byrow=T)
		disc.to.cont <- matrix(NA,ncol=ndisc,nrow=nBigMatrix)
		cont.to.disc <- matrix(NA,nrow=ndisc,ncol=nBigMatrix)
		
		for (j in 1:ndisc) {
			tmp<-dnorm(y,discreteTrans@mean.to.cont[j],discreteTrans@sd.to.cont[j])*h
			if (correction=="constant") tmp<-tmp/sum(tmp) 
			disc.to.cont[,j] <- discreteTrans@discrete.surv[,j]*discreteTrans@discrete.trans["continuous",j]*tmp
			cont.to.disc[j,] <- discreteTrans@distrib.to.discrete[j,]*surv(y,chosen.cov,survObj)*surv.to.discrete
		}
		
		get.disc.matrix <- rbind(cbind(disc.to.disc,cont.to.disc),
				cbind(disc.to.cont,cont.to.cont))
		
	
		rc <- new("IPM.matrix",
				n.discrete = ndisc,
				nEnvClass = 1, 
				nBigMatrix = nBigMatrix,
				nrow = 1*nBigMatrix+ndisc,
				ncol =1*nBigMatrix+ndisc,
				meshpoints = y,
				env.index = rep(1:nEnvClass,each=nBigMatrix),
				names.discrete=rownames(discreteTrans@discrete.trans)[1:ndisc])
		
		rc[,] <-get.disc.matrix   
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


create.compound.Tmatrix <- function(nEnvClass = 2,
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
	if (class(discreteTrans)=="discreteTrans") ndisc <- ncol(discreteTrans@discrete.surv) else ndisc <- 0
	
	#indexes for slotting in IPMs
	indexes <- rep(1:nEnvClass,each=(nBigMatrix+ndisc))
	
	#megamatrix
	megamatrix <- matrix(0,(nBigMatrix+ndisc)*nEnvClass,(nBigMatrix+ndisc)*nEnvClass) 
	
	#print(indexes)
	
	#loop over habitats / environments
	for (k in 1:nEnvClass) { 
		#IPM for individuals starting in env k
		
		if (integrateType=="midpoint") { 
			get.matrix <- (maxSize-minSize)*
					t(outer(y,y,growSurv,cov=as.factor(k),
									growthObj=growObj,survObj=survObj))/nBigMatrix  
		}
		if (integrateType=="cumul") {
			get.matrix.cum <- 
					t(outer(y,b,growthCum,cov=as.factor(k),
									growthObj=growObj))
			get.matrix <- get.matrix.cum[2:(nBigMatrix+1),]-get.matrix.cum[1:nBigMatrix,]
			get.matrix <- t(t(get.matrix)*surv(size=y,cov=k,survObj=survObj))
			
		}
		
		#fix any integration issues reducing survival by dividing by col sums and multiply by survival
		if (correction=="constant") { 
			nvals <- colSums(get.matrix); nvals[nvals==0] <- 1
			get.matrix <- t((t(get.matrix)/nvals)*surv(size=y,cov=as.factor(k),survObj=survObj))    
		}
		
		#names of discrete classes default
		nmes <- ""
		
		
		# In case of discrete classes, take the IPM constructed above and add discrete classes defined in discreteTrans
		if (class(discreteTrans)=="discreteTrans") {
			nmes <- rownames(discreteTrans@discrete.trans)
			surv.to.discrete <- predict(discreteTrans@surv.to.discrete,data.frame(size=y,size2=(y*y)),type="response")
			cont.to.cont <- get.matrix*matrix(1-surv.to.discrete,nrow=nBigMatrix,ncol=nBigMatrix,byrow=T)
			disc.to.disc <- discreteTrans@discrete.trans[1:ndisc,1:ndisc]*matrix(c(discreteTrans@discrete.surv),
					nrow=ndisc,ncol=ndisc,byrow=T)
			disc.to.cont <- matrix(NA,ncol=ndisc,nrow=nBigMatrix)
			cont.to.disc <- matrix(NA,nrow=ndisc,ncol=nBigMatrix)
			
			#print(discreteTrans@mean.to.cont)
			#print(discreteTrans@sd.to.cont)
			
			for (j in 1:ndisc) {
				tmp<-dnorm(y,discreteTrans@mean.to.cont[j],discreteTrans@sd.to.cont[j])*h
				if (correction=="constant") tmp<-tmp/sum(tmp)
				
				#print(tmp)
				
				disc.to.cont[,j] <- discreteTrans@discrete.surv[,j]*discreteTrans@discrete.trans["continuous",j]*tmp
				cont.to.disc[j,] <- discreteTrans@distrib.to.discrete[j,]*surv(y,cov=k,survObj)*surv.to.discrete
			}
			
			get.matrix <- rbind(cbind(disc.to.disc,cont.to.disc),cbind(disc.to.cont,cont.to.cont))
		}
		
		# transit them
		subset <- c(1:nEnvClass)[envMatrix[,k]>0]
		for (j in subset) { 
			megamatrix[indexes==j,indexes==k] <- get.matrix*envMatrix[j,k]; 
		}
		
	}
	
	rc <- new("IPM.matrix",
			nEnvClass = nEnvClass, 
			nBigMatrix = nBigMatrix,
			nrow = nEnvClass*(nBigMatrix+ndisc),
			ncol =nEnvClass*(nBigMatrix+ndisc),
			meshpoints = y,
			env.index = rep(1:nEnvClass,each=nBigMatrix,
					names.discrete=nmes))
	
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
#   chosen.cov - current level of covariate
#   fecObj - a fecundity object
#   integrateType - options include "midpoint" "cumul" 
#   etc...
#Returns -
#  an IPM object


create.IPM.Fmatrix <- function(fecObj,
		nEnvClass = 1,
		nBigMatrix = 50,
		minSize = -1,
		maxSize = 50,
		chosen.cov = 1,
		integrateType="midpoint",
		correction="none") {
	
	# boundary points b and mesh points y
	b<-minSize+c(0:nBigMatrix)*(maxSize-minSize)/nBigMatrix;
	y<-0.5*(b[1:nBigMatrix]+b[2:(nBigMatrix+1)]);
	
	# step size for mid point rule, see equations 4 and 5
	h<-y[2]-y[1]
	#size<-y
	newd <- data.frame(size=y,size2=y^2,size3=y^3)
	if (length(as.numeric(chosen.cov))==1) newd$covariate <- as.factor(rep(chosen.cov,length(y)))
	#print(head(newd))
	for (i in 1:length(fecObj@fit.fec)) if (length(grep("logsize",fecObj@fit.fec[[i]]$formula))>0) newd$logsize <- log(y)
	
	fecObj@fec.constants[is.na(fecObj@fec.constants)] <- 1
	
	fecValues <- matrix(c(rep(1,length(fecObj@fit.fec)),fecObj@fec.constants),ncol=nBigMatrix,nrow=length(fecObj@fit.fec)+length(fecObj@fec.constants))
	
	for (i in 1:length(fecObj@fit.fec)) fecValues[i,]<-predict(fecObj@fit.fec[[i]],newd,type="response")
	
	#Transforms
	if (length(grep("log",fecObj@Transform))>0) for (i in grep("log",fecObj@Transform)) fecValues[i,]<-exp(fecValues[i,])
	if (length(grep("sqrt",fecObj@Transform))>0) for (i in grep("sqrt",fecObj@Transform)) fecValues[i,]<-(fecValues[i,])^2
	if (length(grep("-1",fecObj@Transform))>0) for (i in grep("-1",fecObj@Transform)) fecValues[i,]<-fecValues[i,]+1
	fecValues[!is.finite(fecValues)] <- exp(200)
	prodFecValues<-apply(fecValues,2,prod)
	
	#Kids normal dist
	tmp<-dnorm(y,fecObj@mean.offspring.size,sqrt(fecObj@var.offspring.size))*h
	if (integrateType=="cumul") { 
		tmp1 <- dnorm(b,fecObj@mean.offspring.size,sqrt(fecObj@var.offspring.size))
		tmp <- tmp1[2:(nBigMatrix+1)]-tmp1[1:nBigMatrix]
	}
	if (correction=="constant") tmp<-tmp/sum(tmp)
	to.cont<-tmp%*%t(as.numeric(fecObj@OffspringSplitter["continuous"])*prodFecValues)
	get.matrix <- to.cont
	ndisc <- length(fecObj@OffspringSplitter)-1
	
	namesDiscrete <- "NA"
	if (ndisc>0) {
		namesDiscrete <- colnames(fecObj@OffspringSplitter[1:ndisc])
		
		to.discrete <- as.numeric(fecObj@OffspringSplitter)[1:ndisc]%*%t(prodFecValues)
		
		from.discrete <- matrix(0,ncol=ndisc,nrow=ndisc+nBigMatrix)
		if (names(fecObj@fec.by.discrete)[1]!="NA.") {
			if (sum(names(fecObj@fec.by.discrete)!=namesDiscrete)>0) stop ("Error - the names of the discrete classes as you provided for the data.frame fec.by.discrete are not 100% the same discrete class names in your data.frame OffspringSplitter. They should also be in alphabetical order.")
			from.discrete <- c(as.numeric(fecObj@OffspringSplitter)[1:ndisc],as.numeric(fecObj@OffspringSplitter)[ndisc+1]*tmp)%*%as.matrix(fecObj@fec.by.discrete)
		}
		get.matrix <- cbind(from.discrete,rbind(to.discrete,to.cont))
	}
	
	#warning about negative numbers
	if (min(get.matrix)<0) print("Warning: fertility values < 0 exist in matrix, consider transforms.") 
	
	rc <- new("IPM.matrix",
			n.discrete = ndisc,
			nEnvClass = 1, 
			nBigMatrix = nBigMatrix,
			nrow = 1*nBigMatrix+ndisc,
			ncol =1*nBigMatrix+ndisc,
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


create.compound.Fmatrix <- function(nEnvClass = 2,
		nBigMatrix = 50,
		minSize = -1,
		maxSize = 50,
		envMatrix,
		fecObj,
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
	if (ncol(fecObj@OffspringSplitter)>1) ndisc <- ncol(fecObj@OffspringSplitter)-1 else ndisc <- 0
	
	#indexes for slotting in IPMs
	indexes <- rep(1:nEnvClass,each=(nBigMatrix+ndisc))
	#print(indexes)
	
	
	#megamatrix
	megamatrix <- matrix(0,(nBigMatrix+ndisc)*nEnvClass,(nBigMatrix+ndisc)*nEnvClass) 
	
	#loop over habitats / environments
	for (k in 1:nEnvClass) {
		
		get.matrix <- create.IPM.Fmatrix(nEnvClass =1,
				nBigMatrix = nBigMatrix,
				minSize = minSize,
				maxSize = maxSize,
				chosen.cov = k,
				fecObj=fecObj,
				integrateType=integrateType,
				correction=correction)
		
		#print(range(get.matrix))
		#print(get.matrix[1:5,1:5])
		
		# transit them
		subset <- c(1:nEnvClass)[envMatrix[,k]>0]                 
		for (j in subset) {
			megamatrix[indexes==j,indexes==k] <- get.matrix@.Data*envMatrix[j,k]; 
		}
		
	}
	
	#warning about negative numbers
	if (min(megamatrix)<0) print("Warning: fertility values < 0 exist in matrix, consider transforms.") 
	
	
	rc <- new("IPM.matrix",
			nEnvClass = nEnvClass, 
			nBigMatrix = nBigMatrix,
			nrow = nEnvClass*(nBigMatrix+ndisc),
			ncol =nEnvClass*(nBigMatrix+ndisc),
			meshpoints = y,
			env.index = rep(1:nEnvClass,each=nBigMatrix),
			names.discrete=get.matrix@names.discrete)
	
	
	rc[,] <- megamatrix
	
	return(rc) 
}




# For a single Tmatrix (!not compound and no discrete stages!), this functions makes a series
# of diagnostic plots - this is defined for growthObj,
# growthObj.incr - modification required
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
diagnosticsTmatrix <- function(Tmatrix,growObj,survObj, dff, integrateType="midpoint", correction="none") {
	
	#Print the range of the Tmatrix (should be on 0-1)
	print("Range of Tmatrix is "); print(range(c(Tmatrix)))
	
	#Create a new Tmatrix with 0.5*min and 1.5*max and 1.5*meshpoints
	if (Tmatrix@meshpoints[1]>0)
		new.min <- Tmatrix@meshpoints[1]/2 else new.min <- Tmatrix@meshpoints[1]*1.5
	Tmatrix1 <- create.IPM.Tmatrix(nEnvClass = 1,
			nBigMatrix = floor(length(Tmatrix@meshpoints)*1.5),
			minSize = new.min, maxSize = 1.5*max(Tmatrix@meshpoints),
			chosen.cov = 1, growObj=growObj,survObj=survObj,
			integrateType=integrateType, correction=correction)
	
	#Is the size range sufficient? 
	par(mfrow=c(2,3),bty="l")
	hist(c(dff$size,dff$sizeNext),xlab="Sizes observed",
			ylab="Frequency",main="", xlim=range(c(Tmatrix@meshpoints,dff$size,dff$sizeNext),na.rm=TRUE))
	abline(v=c(Tmatrix@meshpoints[1],Tmatrix@meshpoints[length(Tmatrix@meshpoints)]),lty=2,col=2)
	legend("topright",legend="fitted range", col="red",lty=2,bty="n")
	title("Size range")
	
	#Are there losses due to survival? 
	plot(colSums(Tmatrix),surv(Tmatrix@meshpoints,1,survObj),
			type="l",xlab="Surviving in Tmatrix", ylab="Should be surviving")
	points(colSums(Tmatrix1),surv(Tmatrix1@meshpoints,1,survObj), type="l",col=2)
	abline(0,1,lty=2)
	title("Survival")
	
	#Do resolution and size range affect results for Life Expect? 
	LE <- MeanLifeExpect(Tmatrix)
	LE1 <- MeanLifeExpect(Tmatrix1)
	plot(Tmatrix@meshpoints,LE,type="l",
			xlim=range(Tmatrix1@meshpoints),ylim=range(c(LE,LE1)),xlab="Sizes", ylab="Life expectancy")
	points(Tmatrix1@meshpoints,LE1,type="l",col=2)
	legend("topleft", legend=c("current", "extended"),col=1:2,lty=1,bty="n")
	title("Extend size range + resolution")
	
	
	#Does the resolution adequately reflect normal distribution of growth? (Zuidema issue)
	loctest <- floor(quantile(1:Tmatrix@nBigMatrix,c(0.25,0.5,0.75)))
	h <- diff(Tmatrix@meshpoints)[1]
	testSizes <- seq(min(Tmatrix@meshpoints),max(Tmatrix@meshpoints),length=5000)
	
	#growth plots for midpoint integration
	if (integrateType=="midpoint") {
		
		for (j in 1:3) {
			#prob survive
			ps <- surv(Tmatrix@meshpoints[loctest[j]],Tmatrix@env.index[1],survObj)
			#plot template
			plot(Tmatrix@meshpoints,Tmatrix@.Data[,loctest[j]]/h/ps, type="n",
					xlim=range(Tmatrix@meshpoints[loctest[j]]+
									c(-3.5*summary(growObj@fit)$sigma,+3.5*summary(growObj@fit)$sigma)),
					xlab="Size next", ylab="pdf")
			if (j==1) title("Numerical resolution and growth")
			for (k in 1:length(Tmatrix@meshpoints)) {
				points(c(Tmatrix@meshpoints[k])+c(-h/2,h/2),
						rep(Tmatrix@.Data[k,loctest[j]],2)/h/ps,type="l")
				points(rep(Tmatrix@meshpoints[k]+h/2,2),c(0,Tmatrix@.Data[k,loctest[j]]/h/ps),type="l",lty=1)
				points(rep(Tmatrix@meshpoints[k]-h/2,2),c(0,Tmatrix@.Data[k,loctest[j]]/h/ps),type="l",lty=1)
			}
			newd <- data.frame(size=Tmatrix@meshpoints[loctest[j]],
					size2=Tmatrix@meshpoints[loctest[j]]^2,
					covariate=Tmatrix@env.index[1])
			
			if(length(grep("logsize",growObj@fit$formula)))
				newd$logsize=log(Tmatrix@meshpoints[loctest[j]])
			
			if (length(growObj@fit$model$covariate)>0)
				if (is.factor(growObj@fit$model$covariate))
					newd$covariate <- as.factor(newd$covariate)
			
			#predict mean
			mux <- predict(growObj@fit,newd,type="response")
			#add to size if it is a incr object
			if (class(growObj)=="growthObj.incr") mux <- Tmatrix@meshpoints[loctest[j]]+mux
			#define variance if it is not a declinevar object
			if (class(growObj)!="growthObj.declinevar" &
					class(growObj)!="growthObj.incr.declinevar" &
					class(growObj)!="growthObj.logincr.declinevar")
				sigmax2 <- summary(growObj@fit)$sigma^2 else (print("undefined growth variance class"))
			#plot using mean 
			if (class(growObj)!="growthObj.logincr"){
				points(testSizes,dnorm(testSizes,mux,sqrt(sigmax2)),type="l",col=2)
			} else {
				points(testSizes,dlnorm(testSizes-Tmatrix@meshpoints[loctest[j]],mux,
								sqrt(sigmax2)),type="l",col=2)
			}
			if (j==1) legend("topright", legend=c("Small"),col="white",lty=1,bty="n")
			if (j==2) legend("topright", legend=c("Medium"),col="white",lty=1,bty="n")
			if (j==3) legend("topright", legend=c("Large"),col="white",lty=1,bty="n")
			
		}
	}
	
	#growth plots for cumul integration
	if (integrateType=="cumul")  {
		if (class(growObj)!="growthObj.logincr") rval <- 3.5 else rval <- 0.5
		for (j in 1:3) {
			#prob survive
			ps <- surv(Tmatrix@meshpoints[loctest[j]],Tmatrix@env.index[1],survObj)
			#plot template
			plot(Tmatrix@meshpoints,Tmatrix@.Data[,loctest[j]]/ps, type="n",
					xlim=range(Tmatrix@meshpoints[loctest[j]]+
									c(-rval*summary(growObj@fit)$sigma,+rval*summary(growObj@fit)$sigma)),
					xlab="Size next", ylab="pdf")
			if (j==1) title("Numerical resolution and growth")
			for (k in 1:length(Tmatrix@meshpoints)) {
				points(c(Tmatrix@meshpoints[k])+c(-h/2,h/2),
						rep(Tmatrix@.Data[k,loctest[j]],2)/ps,type="l")
				points(rep(Tmatrix@meshpoints[k]+h/2,2),c(0,Tmatrix@.Data[k,loctest[j]]/ps),type="l",lty=1)
				points(rep(Tmatrix@meshpoints[k]-h/2,2),c(0,Tmatrix@.Data[k,loctest[j]]/ps),type="l",lty=1)
			}
			newd <- data.frame(size=Tmatrix@meshpoints[loctest[j]],
					size2=Tmatrix@meshpoints[loctest[j]]^2,
					covariate=Tmatrix@env.index[1])
			
			if(length(grep("logsize",growObj@fit$formula)))
				newd$logsize=log(Tmatrix@meshpoints[loctest[j]])
			
			
			if (length(growObj@fit$model$covariate)>0)
				if (is.factor(growObj@fit$model$covariate))
					newd$covariate <- as.factor(newd$covariate)
			
			mux <- predict(growObj@fit,newd,type="response")
			if (class(growObj)=="growthObj.incr") mux <- Tmatrix@meshpoints[loctest[j]]+mux
			sigmax2 <- summary(growObj@fit)$sigma^2
			if (class(growObj)!="growthObj.logincr"){
				points(testSizes,dnorm(testSizes,mux,sqrt(sigmax2))*h,type="l",col=2)
			} else {
				points(testSizes,dlnorm(testSizes-Tmatrix@meshpoints[loctest[j]],mux,sqrt(sigmax2))*h,type="l",col=2)
			}
			if (j==1) legend("topright", legend=c("Small"),col="white",lty=1,bty="n")
			if (j==2) legend("topright", legend=c("Medium"),col="white",lty=1,bty="n")
			if (j==3) legend("topright", legend=c("Large"),col="white",lty=1,bty="n")
			
		}
		
	}
	
}







### FUNCTIONS FOR EXTRACTING INTEGRATED DEMOGRAPHIC MEASURES ########################
### i.e. life-expectancy and passage time


#Generic for mean life expectancy
#parameters - an IPM
# returns - the life expectancy for every starting size. 
setGeneric("MeanLifeExpect",
		function(IPM.matrix) standardGeneric("MeanLifeExpect"))

setMethod("MeanLifeExpect",
		c("IPM.matrix"),
		function(IPM.matrix){
			require(MASS)
			nBigMatrix <- length(IPM.matrix@.Data[1,]) #this nBigMatrix actually contains discrete, env, etc
			#tmp <-  ginv(diag(IPM.matrix@nEnvClass*nBigMatrix)-IPM.matrix)
			tmp <-  ginv(diag(nBigMatrix)-IPM.matrix)
			lifeExpect <- colSums(tmp)
			return(lifeExpect)
		})



#Generic for variance life expectancy (p119 Caswell)
#parameters - an IPM
# returns - the variance in life expectancy for every starting size. 
setGeneric("VarLifeExpect",
		function(IPM.matrix) standardGeneric("VarLifeExpect"))

setMethod("VarLifeExpect",
		c("IPM.matrix"),
		function(IPM.matrix){
			require(MASS)
			nBigMatrix <- length(IPM.matrix@.Data[1,])
			#tmp <-  ginv(diag(IPM.matrix@nEnvClass*nBigMatrix)-IPM.matrix)
			tmp <-  ginv(diag(nBigMatrix)-IPM.matrix)
			#varlifeExpect <- (2*diag(tmp)-diag(length(IPM.matrix[,1])))%*%tmp-(tmp*tmp)
			#varLifeExpect <- colSums(varlifeExpect)
			varLifeExpect <- colSums(2*(tmp%*%tmp)-tmp)-colSums(tmp)*colSums(tmp)                  
			return(varLifeExpect)
		})





#Generic for survivorship
#parameters - IPM.matrix - an IPM
#           - size1 - a size at age 1
#           - maxAge - a maxAge
# returns - a list including the survivorship up to the max age,
#                      this broken down by stage,
#                       and mortality over age 

## WON'T WORK WITH DISCRETE STAGES AS IS!!
setGeneric("Survivorship",
		function(IPM.matrix, size1, maxAge) standardGeneric("Survivorship"))

setMethod("Survivorship",
		c("IPM.matrix","numeric","numeric"),
		function(IPM.matrix, size1, maxAge=300){
			nBigMatrix <- length(IPM.matrix@.Data[1,])
			#n <- IPM.matrix@nEnvClass*nBigMatrix
			n <- nBigMatrix
			A1 <- tmp <-  IPM.matrix
			stage.agesurv <- matrix(NA,n,maxAge)
			surv.curv <- rep (NA,maxAge)
			
			#identify the starting size you want to track
			loc <- which(abs(size1-IPM.matrix@meshpoints)==min(abs(size1-IPM.matrix@meshpoints)),arr.ind=T)
			popvec <- matrix(0,n,1)
			popvec[loc,1] <- 1
			
			for (a in 1:maxAge) {
				surv.curv[a]<-sum(A1[,loc]); 
				stage.agesurv[c(1:n),a]<-A1[,]%*%popvec
				A1<-A1%*%tmp
			}
			
			mortality <- -log(surv.curv[2:length(surv.curv)]/surv.curv[1:(length(surv.curv)-1)])
			
			return(list(surv.curv=surv.curv,stage.agesurv=stage.agesurv, mortality = mortality))
		})





#Generic for first passage time 
#parameters - an IPM
#           - a size for which passage time is required            
# returns - the passage time to this size from each of the sizes in the IPM 
setGeneric("PassageTime",
		function(chosen.size,IPM.matrix) standardGeneric("PassageTime"))

setMethod("PassageTime",
		c("numeric","IPM.matrix"),
		function(chosen.size,IPM.matrix){
			require(MASS)
			
			loc <- which(abs(chosen.size-IPM.matrix@meshpoints) ==
							min(abs(chosen.size - IPM.matrix@meshpoints)),arr.ind=TRUE)[1]
			matrix.dim <- length(IPM.matrix[1,])
			
			Tprime <- IPM.matrix
			Tprime[,loc] <- 0
			
			Mprime <- 1-colSums(IPM.matrix)
			Mprime[loc]<-0
			Mprime <- rbind(Mprime,rep(0,matrix.dim))
			Mprime[2,loc] <- 1
			
			Bprime <- Mprime%*% ginv(diag(matrix.dim)-Tprime)
			#print(round(Bprime[2,],2))
			#print(sum(Bprime[2,]<1e-6))
			Bprime[2,][Bprime[2,]==0] <- 1
			
			diagBprime <- diag(Bprime[2,])
			#plot(IPM.matrix@meshpoints,diag(diagBprime),type="l",log="y")
			#abline(v=chosen.size)
			Tc <- diagBprime%*%Tprime%*%ginv(diagBprime)
			eta1 <- ginv(diag(matrix.dim)-Tc)             
			
			time.to.absorb <- colSums(eta1)
			time.to.absorb[loc:length(time.to.absorb)] <- 0
			return(time.to.absorb)
		})




#Generic for first variance first passage time (not sure!!!)
#parameters - an IPM
#           - a size for which passage time is required            
# returns - the variance passage time to this size from each of the sizes in the IPM 
setGeneric("varPassageTime",
		function(chosen.size,IPM.matrix) standardGeneric("varPassageTime"))

setMethod("varPassageTime",
		c("numeric","IPM.matrix"),
		function(chosen.size,IPM.matrix){
			require(MASS)
			
			loc <- which(abs(chosen.size-IPM.matrix@meshpoints)==min(abs(chosen.size-IPM.matrix@meshpoints)),arr.ind=TRUE)
			matrix.dim <- length(IPM.matrix[1,])
			
			Tprime <- IPM.matrix
			Tprime[,loc] <- 0
			
			Mprime <- 1-colSums(IPM.matrix)
			Mprime[loc]<-0
			Mprime <- rbind(Mprime,rep(0,matrix.dim))
			Mprime[2,loc] <- 1
			
			Bprime <- Mprime%*% solve(diag(matrix.dim)-Tprime)
			
			Tc <- diag(Bprime[2,])%*%Tprime%*%ginv(diag(Bprime[2,]))
			eta1 <- ginv(diag(matrix.dim)-Tc)             
			
			vartimeAbsorb <- colSums(2*(eta1%*%eta1)-eta1)-colSums(eta1)*colSums(eta1)                  
			
			return(vartimeAbsorb)
		})




#Generic for life expectancy in a modeled stoch env  -NEEDS DOUBLE-CHECKING
#parameters - a compound IPM of dim nEnvClass*nBigMatrix
#           - an environmental matrix
# returns - the life expectancy for each of the sizes in the IPM (columns)
#           for each of the starting env states
setGeneric("LifeExpect",
		function(IPM.matrix,envMatrix) standardGeneric("LifeExpect"))

setMethod("LifeExpect",
		c("IPM.matrix", "envMatrix"),
		function(IPM.matrix,envMatrix){
			require(MASS)
			
			matrix.dim <- length(IPM.matrix[1,])
			nstages <- IPM.matrix@nBigMatrix
			nstates <- IPM.matrix@nEnvClass
			
			pis <- Re(eigen(envMatrix)$vector[,1])
			pis <- pis/(sum(pis))
			
			#print(pis)
			
			#ckron <- kronecker(envMatrix,diag(nstages))  #doesnt work??
			m <- IPM.matrix  #%*%ckron  #eq 29 in carols paper
			
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
				qatildas[indext,,i] <- IPM.matrix[indext,indext]/envMatrix[i,i]
				#print( IPM.matrix[cbind(indext,indext)]/envMatrix[i,i] )
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
		})






#Generic for first passage time in a stochastic environment
#parameters - a compound IPM (megamatrix)
#           - an environmental transition matrix 
#           - a size for which passage time is required            
# returns - the passage time to this size from each of the sizes in the IPM 
setGeneric("StochPassageTime",
		function(chosen.size,IPM.matrix,envMatrix) standardGeneric("StochPassageTime"))

setMethod("StochPassageTime",
		c("numeric","IPM.matrix","envMatrix"),
		function(chosen.size,IPM.matrix,envMatrix){
			require(MASS)
			#get the right index for the size you want
			loc <- which(abs(chosen.size-
											IPM.matrix@meshpoints)==min(abs(chosen.size-
													IPM.matrix@meshpoints)),arr.ind=TRUE)
			#expand out to find that size in every env
			locs.all <- loc*c(1:IPM.matrix@nEnvClass)
			matrix.dim <- length(IPM.matrix[1,])
			
			Tprime <- IPM.matrix
			Tprime[,locs.all] <- 0
			
			dhat <- 1-colSums(IPM.matrix)
			dhat[locs.all]<-0
			dhat <- rbind(dhat,rep(0,matrix.dim))
			dhat[2,locs.all] <- 1
			
			bhat <- dhat%*% solve(diag(matrix.dim)-Tprime)
			
			Mc <- diag(bhat[2,])%*%Tprime%*%solve(diag(bhat[2,]))
			eta1 <- solve(diag(matrix.dim)-Mc)             
			
			time.to.absorb <-colSums(eta1)
			
			return(time.to.absorb)
		})



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
	
	par(mfrow=c(2,2),bty="l")
	plot(IPM@meshpoints,n.new.dist0[env.index==1],type="l",xlab="size", ylab="n", ylim=range(c(n.new.dist0,n.new.dist)))
	points(IPM@meshpoints,n.new.dist[env.index==1],type="l",col=2)
	abline(v=IPM@meshpoints[cutoff],lty=3)
	
	plot(survivorship[1:t], xlab="Time", ylab="Survivorship", type="l")
	
	if (time.reach>5) { 
		image(1:time.reach,IPM@meshpoints,t(log(ts.dist)),xlab="Time steps", ylab="Number in each size class")
		contour(1:time.reach,IPM@meshpoints,t(log(ts.dist[,1:time.reach])),add=TRUE)
	}
	print(paste("Time to reach:",time.reach))
	
	return(list(ts.dist=ts.dist, time.reach=time.reach, survivorship=survivorship))     
}








# Calculate R0 
#
# Parameters - Fmatrix, Tmatrix
#
# Returns R0
R0.calc<-function(Tmatrix,Fmatrix){
	require(MASS)
	Imatrix <- length(Tmatrix[1,])
	Nmatrix <- ginv(Imatrix-Tmatrix);
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
largeMatrixCalc <- function(Tmatrix,Fmatrix, tol=1.e-8){
	require(Matrix)
	A2 <- Matrix(Tmatrix+Fmatrix);
	nt <- Matrix(1,length(Tmatrix[1,]),1);
	nt1 <- nt; 
	
	h1 <- diff(Tmatrix@meshpoints)[1]
	
	qmax <- 1000;
	lam <- 1; 
	while(qmax>tol) {
		nt1 <- A2%*%nt;
		qmax <- sum(abs((nt1-lam*nt)@x));  
		lam <- sum(nt1@x); 
		nt@x <- (nt1@x)/lam; #slight cheat  
		#cat(lam,qmax,"\n");
	} 
	nt <- matrix(nt@x,length(Tmatrix[1,]),1); 
	stableDist <- nt/(h1*sum(nt)); #normalize so that integral=1
	lamStable <- lam; 
	
	# Check works   
	qmax <- sum(abs(lam*nt-(Tmatrix+Fmatrix)%*%nt)); 
	cat("Convergence: ",qmax," should be less than ",tol,"\n"); 
	
	return(list(lam=lam,stableDist=stableDist,h1=h1)) 
	
}





## Sensitivity of parameters - works for an IPM built out of
## growth, survival, and fecundity objects with a fitted regression
## as one of their slots. Note that fecObj and growObj need to reflect
## the 
##
##

sensParams <- function(growObj,survObj,fecObj,
		nBigMatrix,minSize,maxSize,
		integrateType="midpoint",
		correction="none") { 
	
	delta<-0.0001;
	
	#nfec objects
	nfec <- 0
	fec.coeff.names <- c()
	for (i in 1:length(fecObj@fit.fec)){
		nfec <- nfec + length(fecObj@fit.fec[[i]]$coeff)	
	    fec.coeff.names <- c(fec.coeff.names,
				paste("reprod",i,names(fecObj@fit.fec[[i]]$coeff)))
	
	}
	
	# define numbers of parameters
	npar <- length(growObj@fit$coeff)+1+
			length(survObj@fit$coeff)+
			(sum(!is.na(fecObj@fec.constants)))+nfec	
	#print(npar)
	
	# create a named vector to hold them 
	elam <- rep(0,npar);
	nmes <- c(paste("grow",names(growObj@fit$coeff)), "sd growth",paste("surv",names(survObj@fit$coeff)),
			paste("reprod constant",which(!is.na(fecObj@fec.constants), arr.ind=TRUE)),
			fec.coeff.names)
	
	names(elam) <- nmes[1:npar]
	
	slam <- elam
	
	# build the IPM and get the lamdba value
	Tmatrix <- create.IPM.Tmatrix(nBigMatrix = nBigMatrix, minSize=minSize,maxSize=maxSize,
			growObj=growObj,survObj=survObj, integrateType=integrateType, correction=correction)
	Fmatrix <- create.IPM.Fmatrix(nBigMatrix = nBigMatrix, minSize=minSize,maxSize=maxSize,
			fecObj=fecObj,integrateType=integrateType, correction=correction)
	IPM <- Tmatrix+Fmatrix
	lambda1 <- Re(eigen(IPM)$value[1]);#print(lambda1)
	
	# change the growth parameters
	for (param.test in 1:length(growObj@fit$coeff)){
		growObj@fit$coefficients[param.test] <- growObj@fit$coefficients[param.test]*(1+delta);
		Tmatrix <- create.IPM.Tmatrix(nBigMatrix = nBigMatrix, minSize=minSize,maxSize=maxSize,
				growObj=growObj,survObj=survObj, integrateType=integrateType, correction=correction)
		Fmatrix <- create.IPM.Fmatrix(nBigMatrix = nBigMatrix, minSize=minSize,maxSize=maxSize,
				fecObj=fecObj, integrateType=integrateType, correction=correction)
		IPM <- Tmatrix+Fmatrix
		lambda2 <- Re(eigen(IPM)$value[1]); #print(lambda2)
		growObj@fit$coefficients[param.test] <- growObj@fit$coefficients[param.test]/(1+delta);
		slam[param.test]<-(lambda2-lambda1)/(growObj@fit$coefficients[param.test]*delta);
		elam[param.test]<-(lambda2-lambda1)/(log(1+delta));
	}
	# change the variance in growth
	param.test <- param.test+1
	resids <- growObj@fit$residuals
	growObj@fit$residuals <- rnorm(length(growObj@fit$residuals),0,sd(growObj@fit$residuals)*(1+delta))
	Tmatrix <- create.IPM.Tmatrix(nBigMatrix = nBigMatrix, minSize=minSize,maxSize=maxSize,
			growObj=growObj,survObj=survObj, integrateType=integrateType, correction=correction)
	Fmatrix <- create.IPM.Fmatrix(nBigMatrix = nBigMatrix, minSize=minSize,maxSize=maxSize,
			fecObj=fecObj, integrateType=integrateType, correction=correction)
	IPM <- Tmatrix + Fmatrix
	lambda2 <- Re(eigen(IPM)$value[1]); #print(lambda2)
	growObj@fit$residuals <- resids
	slam[param.test]<-(lambda2-lambda1)/(sd(growObj@fit$residuals)*delta);
	elam[param.test]<-(lambda2-lambda1)/(log(1+delta));
	
	# change the survival parameters
	count <- param.test
	for (param.test in 1:length(survObj@fit$coeff)){
		survObj@fit$coefficients[param.test] <- survObj@fit$coefficients[param.test]*(1+delta);
		Tmatrix <- create.IPM.Tmatrix(nBigMatrix = nBigMatrix, minSize=minSize,maxSize=maxSize,
				growObj=growObj,survObj=survObj, integrateType=integrateType, correction=correction)
		Fmatrix <- create.IPM.Fmatrix(nBigMatrix = nBigMatrix, minSize=minSize,maxSize=maxSize,
				fecObj=fecObj,integrateType=integrateType, correction=correction)
		IPM <- Tmatrix+Fmatrix
		lambda2 <- Re(eigen(IPM)$value[1]);
		survObj@fit$coefficients[param.test] <- survObj@fit$coefficients[param.test]/(1+delta);
		slam[param.test+count]<-(lambda2-lambda1)/(survObj@fit$coefficients[param.test]*delta);
		elam[param.test+count]<-(lambda2-lambda1)/(log(1+delta));
	}
	
	#change the constant fecundity objects
	chs <- which(!is.na(fecObj@fec.constants), arr.ind=TRUE)
	if (length(chs)>0) { 
		count <- count + param.test;
		for (param.test in 1:length(chs)) {
			fecObj@fec.constants[chs[param.test]] <- fecObj@fec.constants[chs[param.test]]*(1+delta);
			Tmatrix <- create.IPM.Tmatrix(nBigMatrix = nBigMatrix, minSize=minSize,maxSize=maxSize,
					growObj=growObj,survObj=survObj, integrateType=integrateType, correction=correction)
			Fmatrix <- create.IPM.Fmatrix(nBigMatrix = nBigMatrix, minSize=minSize,maxSize=maxSize,
					fecObj=fecObj,integrateType=integrateType, correction=correction)
			IPM <- Tmatrix+Fmatrix
			lambda2 <- Re(eigen(IPM)$value[1]);
			fecObj@fec.constants[chs[param.test]] <- fecObj@fec.constants[chs[param.test]]/(1+delta);
			slam[param.test+count]<-(lambda2-lambda1)/(fecObj@fec.constants[param.test]*delta);
			elam[param.test+count]<-(lambda2-lambda1)/(log(1+delta));
			
		}
	}
	
	
	# change the reprod prob parameters in sequence
	count <- count + param.test;
	for (i in 1:length(fecObj@fit.fec)){
	for (param.test in 1:length(fecObj@fit.fec[[i]]$coeff)){
		fecObj@fit.fec[[i]]$coefficients[param.test] <- fecObj@fit.fec[[i]]$coefficients[param.test]*(1+delta);
		Tmatrix <- create.IPM.Tmatrix(nBigMatrix = nBigMatrix, minSize=minSize,maxSize=maxSize,
				growObj=growObj,survObj=survObj, integrateType=integrateType, correction=correction)
		Fmatrix <- create.IPM.Fmatrix(nBigMatrix = nBigMatrix, minSize=minSize,maxSize=maxSize,
				fecObj=fecObj,integrateType=integrateType, correction=correction)
		IPM <- Tmatrix+Fmatrix
		lambda2 <- Re(eigen(IPM)$value[1]);
		fecObj@fit.fec[[i]]$coefficients[param.test] <- fecObj@fit.fec[[i]]$coefficients[param.test]/(1+delta);
		slam[param.test+count]<-(lambda2-lambda1)/(fecObj@fit.fec[[i]]$coefficients[param.test]*delta);
		elam[param.test+count]<-(lambda2-lambda1)/(log(1+delta));		
	}
	count <- count + param.test;
	}	
	
	return(list(slam=slam,elam=elam))
	
}


### FUNCTIONS FOR EXTRACTING STOCH GROWTH RATE ########################

# Generic approach to get stoch rate
# by sampling list IPM
#
# Parameters - list.IPM.matrix - list of IPMs corresponding to different year types
#            - n.runin - the burnin before establishing lambda_s
#            - Tmax - the total time-horizon for getting lambda_s
#
# Returns lambda_s (no density dependence)

setGeneric("StochGrowthRateSampleList",
		function(list.IPM.matrix, n.runin, Tmax) standardGeneric("StochGrowthRateSampleList"))

setMethod("StochGrowthRateSampleList",
		c("list","numeric","numeric"),
		function(list.IPM.matrix,n.runin,Tmax){
			require(MASS)
			
			nmatrices <- length(list.IPM.matrix)
			
			nt<-rep(1,length(list.IPM.matrix[[1]][,1]))
			Rt<-rep(NA,Tmax)
			
			for (t in 1:Tmax) {
				year.type <- sample(1:nmatrices,size=1,replace=FALSE)
				nt1<-list.IPM.matrix[[year.type]] %*% nt	
				sum.nt1<-sum(nt1)
				Rt[t]<-log(sum.nt1)
				nt<-nt1/sum.nt1
				
			}
			
			res <- mean(Rt[n.runin:Tmax],na.rm=TRUE)
			return(res)
		})


# Approach to get stoch rate
# with time-varying covariates
#
# Parameters - covariate - the covariates (temperature, etc)
#            - n.runin - the burnin before establishing lambda_s
#            - Tmax - the total time-horizon for getting lambda_s
#
# Returns lambda_s 


StochGrowthRateManyCov <- function(covariate,n.runin,Tmax,
		growthObj,survObj,fecObj,
		nBigMatrix,minSize,maxSize, n.microsites,
		integrateType="midpoint",correction="none"){
	require(MASS)
	
	
	nt<-rep(1,nBigMatrix)
	Rt<-rep(NA,Tmax)
	fecObj@fec.constants[is.na(fecObj@fec.constants)] <- 1 
	tmp.fecObj <- fecObj
	
	#print("fec.const start")
	#print(fecObj@fec.constants)
	
	#density dep? 
	if (sum(n.microsites)>0) { dd <- TRUE } else { dd <- FALSE}
	
	
	for (t in 1:Tmax) {
		if (dd) tmp.fecObj@fec.constants <- c(fecObj@fec.constants, 
					min(n.microsites[min(t,length(n.microsites))]/nt[1],1))
		
		
		tpS <- create.IPM.Tmatrix(nBigMatrix = nBigMatrix, minSize = minSize,
				maxSize = maxSize, chosen.cov = covariate[t,],
				growObj = growthObj, survObj = survObj,
				integrateType=integrateType, correction=correction)
		
		tpF <- create.IPM.Fmatrix(nBigMatrix = nBigMatrix, minSize = minSize,
				maxSize = maxSize, #chosen.cov = covariate[t,],
				fecObj = tmp.fecObj,
				integrateType=integrateType, correction=correction)
		
		#print(range(tpF))
		#print(tmp.fecObj)	
		
		IPM.here <- tpF@.Data+tpS@.Data
		nt1<-IPM.here %*% nt	
		sum.nt1<-sum(nt1)
		Rt[t]<-log(sum.nt1)
		nt<-nt1/sum.nt1
		
	}
	
	res <- mean(Rt[n.runin:Tmax],na.rm=TRUE)
	return(res)
}


# Approach to get track pop struct
# with time-varying covariates; assuming density
# dep in seedling establishment (i.e limited no microsites)
#
# Parameters - covariate - the covariates (temperature, etc)
#            - n.runin - the burnin before establishing lambda_s
#            - Tmax - the total time-horizon for getting lambda_s
#
# Returns matrix with time as columns, and pop struct as rows


TrackPopStructManyCov<-function(covariate,n.runin,Tmax,
		growthObj,survObj,fecObj,
		nBigMatrix,minSize,maxSize,
		n.microsites,integrateType="midpoint",correction="none"){
	require(MASS)
	
	nt <- rep(1,nBigMatrix)
	rc <- matrix(NA,nBigMatrix,Tmax)
	fecObj@fec.constants[is.na(fecObj@fec.constants)] <- 1 
	tmp.fecObj <- fecObj
	#density dep? 
	if (sum(n.microsites)>0) { dd <- TRUE } else { dd <- FALSE}
	
	
	for (t in 1:Tmax) {
		if (dd) tmp.fecObj@fec.constants <- c(fecObj@fec.constants, 
					min(n.microsites[min(t,length(n.microsites))]/nt[1],1))
		
		#print(tmp.fecObj@fec.constants)
		
		tpS <- create.IPM.Tmatrix(nBigMatrix = nBigMatrix, minSize = minSize,
				maxSize = maxSize, chosen.cov = covariate[t,],
				growObj = growthObj, survObj = survObj,
				integrateType=integrateType, correction=correction)
		tpF <- create.IPM.Fmatrix(nBigMatrix = nBigMatrix, minSize = minSize,
				maxSize = maxSize, #chosen.cov = covariate[t,],
				fecObj = tmp.fecObj,
				integrateType=integrateType, correction=correction)
		
		IPM.here <- tpF+tpS
		nt1<-IPM.here %*% nt	
		rc[,t] <- nt1
		nt<-nt1
		
		#print("Smat")
		#print(range(tpS))
		#print("Fmat")
		#print(range(tpF))
		
		
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







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
                                 explanatoryVariables="size+size2",
                                 responseType="sizeNext",
                                 regType="constantVar"){


    if (responseType=="incr" & length(dataf$incr)==0) {
        print("building incr as sizeNext-size")
        dataf$incr <- dataf$sizeNext-dataf$size
    }

    if (responseType=="logincr" & length(dataf$logincr)==0) {
        print("building logincr as log(sizeNext-size) - pre-build if this is not appropriate")
        dataf$logincr <- log(dataf$sizeNext-dataf$size)
    }

    Formula<-paste(responseType,'~',explanatoryVariables,sep='')
    
    #create appropriate size based covariates
    dataf$size2 <- dataf$size^2
    dataf$size3 <- dataf$size^3
    if (length(grep("logsize",Formula))>0) dataf$logsize <- log(dataf$size)

    #setup for discrete covariates if data suggests may be implemented by the
    #presence of "covariate" and "covariatenext"
    if (length(grep("covariate",Formula))>0 &
        length(grep("covariatenext",colnames(dataf)))>0) {
            dataf$covariate <- as.factor(dataf$covariate)
             #convert to 1:n for indexing later
            levels(dataf$covariate) <- 1:length(unique(dataf$covariate))
    }
    
    #eval fit
    if (regType=="constantVar")  {
        fit <-lm(Formula, data=dataf)
    } else { 
        if (regType=="declineVar"){
            fit<-gls(formula(Formula),
                     na.action=na.omit,weight=varExp(form=~fitted(.)),data=dataf)
        }
    }

    #make the objects
    #with sizeNext as response
    if (responseType=="sizeNext") { 

        if (class(fit)=="lm") { 
            gr1 <- new("growthObj")
            gr1@fit <- fit
        } else {
            if (class(fit)=="gls") { 
                gr1 <- new("growthObj.declinevar")
                gr1@fit <- fit
            } else {
                print("unknown formula;
                 please use lm or gls for declining variance models")
            }
        }
    } else {
        if (responseType=="incr") { 

            if (class(fit)=="lm") { 
                gr1 <- new("growthObj.incr")
                gr1@fit <- fit
            } else {
                if (class(fit)=="gls") { 
                    gr1 <- new("growthObj.incr.declinevar")
                    gr1@fit <- fit
                } else {
                    print("unknown formula;
                 please use lm or gls for declining variance models")
                }
            }
  
        } else {
            if (responseType=="logincr") {

                 if (class(fit)=="lm") { 
                    gr1 <- new("growthObj.logincr")
                    gr1@fit <- fit
                } else {
                    if (class(fit)=="gls") { 
                        gr1 <- new("growthObj.logincr.declinevar")
                        gr1@fit <- fit
                    } else {
                        print("unknown formula;
                 please use lm or gls for declining variance models")
                    }
                }
                
            }
        }}
    
    
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
					explanatoryVariables="size+size2+covariate1",
					responseType="sizeNext",
					regType="constantVar"){
    
				
	if (responseType=="incr" & length(dataf$incr)==0) {
		print("building incr as sizeNext-size")
		dataf$incr <- dataf$sizeNext-dataf$size
	}
				
	if (responseType=="logincr" & length(dataf$logincr)==0) {
		print("building logincr as log(sizeNext-size) - pre-build if this is not appropriate")
		dataf$logincr <- log(dataf$sizeNext-dataf$size)
	}
				
				
    Formula<-paste(responseType,'~',explanatoryVariables,sep='')
    
    #create appropriate size based covariates
    dataf$size2 <- dataf$size^2
    dataf$size3 <- dataf$size^3
    if (length(grep("logsize",Formula))>0) dataf$logsize <- log(dataf$size)
    
    #eval fit
    if (regType=="constantVar")  {
        fit <-lm(Formula, data=dataf)
    } else { 
        if (regType=="declineVar"){
            fit<-gls(formula(Formula),
                     na.action=na.omit,weight=varExp(form=~fitted(.)),data=dataf)
    }}

    #make the objects
    #with sizeNext as response
    if (responseType=="sizeNext") { 

        if (class(fit)=="lm") { 
            gr1 <- new("growthObjMultiCov")
            gr1@fit <- fit
        } else {
            if (class(fit)=="gls") { 
                gr1 <- new("growthObjMultiCov.declinevar")
                gr1@fit <- fit
            } else {
                print("unknown formula;
                 please use lm or gls for declining variance models")
            }
        }
    } else {
        if (responseType=="incr") { 

            if (class(fit)=="lm") { 
                gr1 <- new("growthObjMultiCov.incr")
                gr1@fit <- fit
            } else {
                if (class(fit)=="gls") { 
                    gr1 <- new("growthObjMultiCov.incr.declinevar")
                    gr1@fit <- fit
                } else {
                    print("unknown formula;
                 please use lm or gls for declining variance models")
                }
            }
  
        } else {
            if (responseType=="logincr") {

                 if (class(fit)=="lm") { 
                    gr1 <- new("growthObjMultiCov.logincr")
                    gr1@fit <- fit
                } else {
                    if (class(fit)=="gls") { 
                        gr1 <- new("growthObjMultiCov.logincr.declinevar")
                        gr1@fit <- fit
                    } else {
                        print("unknown formula;
                 please use lm or gls for declining variance models")
                    }
                }
                
            }
        }}
    
    
    return(gr1)
    
}




## Function to create a new Hossfeld growth object
#
# Parameters - dataf - a dataframe
#
# Returns - a Hossfeld growth object

#no covariate, and one polynom, linear regression
makeGrowthObjHossfeld <- function(dataf) {  
    if (length(dataf$incr)==0) dataf$incr <- dataf$sizeNext-dataf$size
    dataf$incr[dataf$incr<0] <- 0
    tmp <- optim(c(1, 1, 1), wrapHossfeld, dataf = dataf, method = "Nelder-Mead")
    print(tmp$convergence)
    gr1 <- new("growthObj.Hossfeld")
    gr1@paras <- tmp$par
    resids <- Hossfeld(dataf$size, tmp$par) - dataf$incr 
    gr1@sd <- sd(resids, na.rm = T)
    return(gr1)
}


# Function to create a truncated increment growth model 
#
# Paramteres - dataf - a dataframe
#
# Retrurns - 
#
#no covariate, and one polynom, linear regression on increment
makeGrowthObjIncrTrunc <- function(dataf) {
    require(censReg)
    dataf$size2 <- dataf$size^2
    if (length(dataf$incr)==0) dataf$incr <- dataf$sizeNext-dataf$size
    fit <- censReg(incr~logsize+logsize2,data=dataf, left=0)
    gr1 <- new("growthObj.truncincr")
    gr1@fit <- fit$estimate
    return(gr1)
}




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
    #presence of "covariate" and "covariatenext"
    if (length(grep("covariate",formula))>0 &
        length(grep("covariatenext",colnames(dataf)))>0) {
        dataf$covariate <- as.factor(dataf$covariate)
        #convert to 1:n for indexing later
        levels(dataf$covariate) <- 1:length(unique(dataf$covariate))
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
                              fec.constants=as.numeric(NA),
                              explanatoryVariables="size",
                              Family="gaussian",
                              Transform="none",
                              mean.offspring.size=NA,
                              var.offspring.size=NA,
                              offspring.splitter=data.frame(continuous=1),
                              fec.by.discrete=matrix(NA,nrow=0,ncol=0)){

	#order stage names from discrete to continuous
	stages <- names(tapply(c(levels(dataf$stage),levels(dataf$stageNext)),c(levels(dataf$stage),levels(dataf$stageNext)),length))
	stages <- stages[stages!="dead"] 
	stages <- c(stages[stages!="continuous"],"continuous") 
	if ((sum(names(offspring.splitter)%in%stage.names)/length(offspring.splitter))<1) {
		stop("Error - the variable names in your offspring.splitter data.frame are not all part of the levels of stage or stageNext in your data file. Please fix this by adjusting your offspring.splitter entry to include the correct variable names, e.g. offspring.splitter=data.frame(continuous=.7,seed.age.1=.3)")
	}
	dummy<-rep(0,length(stages));names(dummy)<-stages;dummy<-as.data.frame(t(as.matrix(dummy)))
	for (i in names(offspring.splitter)) dummy[i]<-offspring.splitter[i]
	offspring.splitter <- dummy
  
	##warnings
    if (length(dataf$stage)==0) {
        print("Warning - no column named stage - assuming all continuous")
        dataf$stageNext <- dataf$stage <- rep("continuous", length(dataf[,1]))
        dataf$stage[is.na(dataf$size)] <- NA
        dataf$stageNext[is.na(dataf$sizeNext)] <- "dead"
    }
	
    if (ncol(offspring.splitter)>1 & (ncol(offspring.splitter)-1)!=ncol(fec.by.discrete)) {
        print("Warning - offspring splitter indicates more than just continuous stages. No fecundity by the discrete stages supplied in fec.by.discrete; assumed that is 0")
        fec.by.discrete <- matrix(0,col(offspring.splitter)-1,col(offspring.splitter)-1)
    }
    
	if (sum(offspring.splitter)!=1) {
		print("Warning - offspring splitter does not sum to 1. It is now rescaled to sum to 1.")
		offspring.splitter <- offspring.splitter / sum(offspring.splitter) 
	}
	
	if ("covariate"%in%strsplit(explanatoryVariables,"[+-\\*]")[[1]]) dataf$covariate <- as.factor(dataf$covariate)
	if ("covariatenext"%in%strsplit(explanatoryVariables,"[+-\\*]")[[1]]) dataf$covariatenext <- as.factor(dataf$covariatenext)
    
	f1 <- new("fecObj")
    dataf$size2 <- dataf$size^2
    if (length(grep("logsize",explanatoryVariables))>0) dataf$logsize <- log(dataf$size)
    
    fecnames <- names(dataf)[grep("fec",names(dataf))]
    if (length(fecnames)>length(explanatoryVariables)) {
        misE <- (length(explanatoryVariables)+1):length(fecnames)
        print(c("number in explanatoryVariables not the same as the number of fecundity columns in the data file, using default of `size' for missing ones which are:",fecnames[misE],". (which might be exactly what you want)"))
        explanatoryVariables <- c(explanatoryVariables,rep("size",length(fecnames)-length(explanatoryVariables)))
    }
    if (length(fecnames)>length(Family)) {
        misE <- (length(Family)+1):length(fecnames)
        print(c("number of families not the same as the number of fecundity columns in the data file, using default of `gaussian' for missing ones which are:",fecnames[misE],". (which might be exactly what you want)"))
        Family <- c(Family,rep("gaussian",length(fecnames)-length(Family)))
    }
    if (length(fecnames)>length(Transform)) {
        misE <- (length(Transform)+1):length(fecnames)
        print(c("number of transforms not the same as the number of fecundity columns in the data file, using default of `none' for missing ones which are:",fecnames[misE],". (which might be exactly what you want)"))
        Transform <- c(Transform,rep("none",length(fecnames)-length(Transform)))
    }
    	
    for (i in 1:length(fecnames)) {
        if (Transform[i]=="log") dataf[,fecnames[i]] <- log(dataf[,fecnames[i]])
        if (Transform[i]=="sqrt") dataf[,fecnames[i]] <- sqrt(dataf[,fecnames[i]])
        dataf[!is.finite(dataf[,fecnames[i]]),fecnames[i]] <- NA
        fit <- glm(paste(fecnames[i],'~',explanatoryVariables[i],sep=''),family=Family[i],data=dataf)
        if (i==1) f1@fit.fec1 <- fit
        if (i==2) f1@fit.fec2 <- fit
        if (i==3) f1@fit.fec3 <- fit
        if (i==4) f1@fit.fec4 <- fit
        if (i==5) f1@fit.fec5 <- fit
        if (i==6) f1@fit.fec6 <- fit
        if (i==7) f1@fit.fec7 <- fit
        if (i==8) f1@fit.fec8 <- fit
        if (i==9) f1@fit.fec9 <- fit
    }
    
	if (is.na(mean.offspring.size)) {
        offspringdata<-subset(dataf,is.na(dataf$stage)&dataf$stageNext=="continuous")
        mean.offspring.size <- mean(offspringdata$sizeNext)
        var.offspring.size <- var(offspringdata$sizeNext) }
    f1@fec.constants <- fec.constants
    f1@mean.offspring.size <- mean.offspring.size
    f1@var.offspring.size <- var.offspring.size
    f1@offspring.splitter <- offspring.splitter 
    f1@fec.by.discrete <- fec.by.discrete
    f1@Transform <- Transform
    return(f1)
}
		 


## NO different from the above yet, except in what it produces

makeFecObjManyCov <- function(dataf,
		fec.constants=as.numeric(NA),
		explanatoryVariables="size",
		Family="gaussian",
		Transform="none",
		mean.offspring.size=NA,
		var.offspring.size=NA,
		offspring.splitter=data.frame(continuous=1),
		fec.by.discrete=matrix(NA,nrow=0,ncol=0)){
	
	##warnings
	if (length(dataf$stage)==0) {
		print("Warning - no column named stage - assuming all continuous")
		dataf$stageNext <- dataf$stage <- rep("continuous", length(dataf[,1]))
		dataf$stage[is.na(dataf$size)] <- NA
		dataf$stageNext[is.na(dataf$sizeNext)] <- "dead"
	}
	
	if(ncol(offspring.splitter)>1 & (ncol(offspring.splitter)-1)!=ncol(fec.by.discrete)) {
		print("Warning - offspring splitter indicates more than just continuous stages. No fecundity by the discrete stages supplied in fec.by.discrete; assumed that is 0")
		fec.by.discrete <- matrix(0,col(offspring.splitter)-1,col(offspring.splitter)-1)
	}
	
	if(sum(offspring.splitter)!=1) {
		print("Warning - offspring splitter does not sum to 1. It is now rescaled to sum to 1.")
		
	}
	
	if (length(grep("covariate",explanatoryVariables))>0) {
		dataf$covariate <- as.factor(dataf$covariate)
		dataf$covariatenext <- as.factor(dataf$covariatenext)
		levels(dataf$covariate) <- 1:length(unique(dataf$covariate))
	}
	
	f1 <- new("fecObjMultiCov")
	dataf$size2 <- dataf$size^2
	if (length(grep("logsize",explanatoryVariables))>0) dataf$logsize <- log(dataf$size)
	
	fecnames <- names(dataf)[grep("fec",names(dataf))]
	if (length(fecnames)>length(explanatoryVariables)) {
		misE <- length(explanatoryVariables):length(fecnames)
		print(c("number in explanatoryVariables not the same as the number of fecundity columns in the data file, using default of `size' for missing ones which are:",fecnames[misE]))
		explanatoryVariables <- c(explanatoryVariables,rep("size",length(fecnames)-length(explanatoryVariables)))
	}
	if (length(fecnames)>length(Family)) {
		misE <- length(Family):length(fecnames)
		print(c("number of families not the same as the number of fecundity columns in the data file, using default of `gaussian' for missing ones which are:",fecnames[misE]))
		Family <- c(Family,rep("gaussian",length(fecnames)-length(Family)))
	}
	if (length(fecnames)>length(Transform)) {
		misE <- length(Transform):length(fecnames)
		print(c("number of transforms not the same as the number of fecundity columns in the data file, using default of `none' for missing ones which are:",fecnames[misE]))
		Transform <- c(Transform,rep("none",length(fecnames)-length(Transform)))
	}
	
	for (i in 1:length(fecnames)) {
		
		if (Transform[i]=="log") dataf[,fecnames[i]] <- log(dataf[,fecnames[i]])
		if (Transform[i]=="sqrt") dataf[,fecnames[i]] <- sqrt(dataf[,fecnames[i]])
		dataf[!is.finite(dataf[,fecnames[i]]),fecnames[i]] <- NA
		
		#print(range(dataf[,fecnames[i]]))
		#print(range(dataf[,"size"], na.rm=TRUE))
		#print(range(dataf[,"covariate"]))
		
		fit <- glm(paste(fecnames[i],'~',explanatoryVariables[i],sep=''),family=Family[i],data=dataf)
		if (i==1) f1@fit.fec1 <- fit
		if (i==2) f1@fit.fec2 <- fit
		if (i==3) f1@fit.fec3 <- fit
		if (i==4) f1@fit.fec4 <- fit
		if (i==5) f1@fit.fec5 <- fit
		if (i==6) f1@fit.fec6 <- fit
		if (i==7) f1@fit.fec7 <- fit
		if (i==8) f1@fit.fec8 <- fit
		if (i==9) f1@fit.fec9 <- fit
	}
	
	if (is.na(mean.offspring.size)) {
		offspringdata<-subset(dataf,is.na(dataf$stage)&dataf$stageNext=="continuous")
		mean.offspring.size <- mean(offspringdata$sizeNext)
		var.offspring.size <- var(offspringdata$sizeNext) }
	f1@fec.constants <- fec.constants
	f1@mean.offspring.size <- mean.offspring.size
	f1@var.offspring.size <- var.offspring.size
	f1@offspring.splitter <- offspring.splitter 
	f1@fec.by.discrete <- fec.by.discrete
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
# Returns - an object of class DiscreteTrans
#
makeDiscreteTrans <- function(dataf) {

	#order stage names from discrete to continuous
	stages <- names(tapply(c(levels(dataf$stage),levels(dataf$stageNext)),c(levels(dataf$stage),levels(dataf$stageNext)),length))
	stages <- stages[stages!="dead"] 
	stages <- c(stages[stages!="continuous"],"continuous") 
    #define the number of classes
    nclasses <- length(stages)
    #define matrices to hold the transition between all classes
    discrete.trans <- matrix(NA,nrow=nclasses,ncol=nclasses, dimnames=list(stages,stages))
    #define matrix to hold sd and mean of re-entry into continous + matrix of  survival for all discrete stages
    sd.to.cont <- mean.to.cont <- discrete.surv <- matrix(NA,nrow=1,ncol=nclasses-1,dimnames=list(1,stages[1:length(stages)-1]))
    # define matrix to hold transitions from the continuous to the discrete
    distrib.to.discrete <- matrix(NA,ncol=1,nrow=nclasses-1,dimnames=list(stages[1:length(stages)-1],"continuous"))

    #loop over discrete stages and fill 
    for (j in stages[1:(length(stages)-1)]) {
      for (i in stages) discrete.trans[i,j] <- sum(dataf[dataf$stage==j & dataf$stageNext==i,]$number,na.rm=T)
      discrete.surv[,j] <- sum(discrete.trans[,j],na.rm=T)/sum(dataf[dataf$stage==j,]$number,na.rm=T)
      discrete.trans[,j] <- discrete.trans[,j]/sum(discrete.trans[,j],na.rm=T)
      mean.to.cont[,j] <- mean(dataf[dataf$stage==j&dataf$stageNext==i,]$sizeNext,na.rm=T)
      sd.to.cont[,j] <- sd(dataf[dataf$stage==j&dataf$stageNext==i,]$sizeNext,na.rm=T)
      }

    for (i in stages[1:(length(stages)-1)])
        distrib.to.discrete[i,] <- sum(dataf[dataf$stage=="continuous"&dataf$stageNext==i,]$number,na.rm=T)
    distrib.to.discrete <- distrib.to.discrete/sum(distrib.to.discrete,na.rm=T)

    
    subdata <- subset(dataf,dataf$stage=="continuous"&dataf$surv==1)
    subdata$cont.to.discrete <- 1
    subdata$cont.to.discrete[subdata$stageNext=="continuous"] <- 0
    subdata$size2 <- subdata$size^2
    surv.to.discrete <- glm(cont.to.discrete~size+size2,family=binomial,data=subdata)

    #define new object
    disTrans <- new("DiscreteTrans")
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

DeathDataAugment <- function (dataf, size.thresh, prop.dead) { 
    
    n.now <- sum(dataf$size>size.thresh)
    n.new.dead <- ceiling(prop.dead*n.now/(1-prop.dead))
    new.size <- rnorm(n.new.dead,size.thresh,sd(dataf$size[dataf$size>size.thresh]))

    datanew <- data.frame(size =new.size, sizeNext=rep(NA,n.new.dead), surv=rep(0,n.new.dead), 
                        covariate = rep(0,n.new.dead), covariatenext = rep(0,n.new.dead),
                          fec = rep(NA,n.new.dead), age = rep(NA,n.new.dead)) 

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
                               explanatoryVariables="size+size2+covariate",
                               responseType="sizeNext",
                               meanB=rep(0,3),varB=rep(1e10),nitt=50000) {
    
    require(MCMCglmm)
    
    if (responseType=="incr" & length(dataf$incr)==0) {
        print("building incr as sizeNext-size")
        dataf$incr <- dataf$sizeNext-dataf$size
    }

    if (responseType=="logincr" & length(dataf$logincr)==0) {
        print("building logincr as log(sizeNext-size) - pre-build if this is not appropriate")
        dataf$logincr <- log(dataf$sizeNext-dataf$size)
    }
    
    dataf$size2 <- dataf$size^2
    dataf$size3 <- dataf$size^3
	if (length(dataf$covariate)>0){
		dataf$covariate <- as.factor(dataf$covariate)
    	levels(dataf$covariate) <- 1:length(unique(dataf$covariate))
	}
    if (length(grep("logsize",explanatoryVariables))>0) dataf$logsize <- log(dataf$size)
    
	#get rid of NAs
	dataf <- dataf[!is.na(dataf$size) & !is.na(dataf$sizeNext),]
	
    #fit growth model
    Formula<-as.formula(paste(responseType,'~',explanatoryVariables,sep=''))
    fit<-MCMCglmm(Formula, data=dataf, verbose=FALSE, nitt=nitt)
    dummy.fit <- lm(Formula,data=dataf)

    #create list of growth models reflecing posterior
    gr <- list()
    for (k in 1:length(fit$Sol[,1])) {
        dummy.fit$coefficients <- fit$Sol[k,]
        dummy.fit$residuals <- rnorm(length(dummy.fit$residuals),0,sqrt(fit$VCV[k,1]))
        if (responseType=="sizeNext") gr[[k]] <-  new("growthObj")
        if (responseType=="incr") gr[[k]] <-  new("growthObj.incr")
        if (responseType=="logincr") gr[[k]] <-  new("growthObj.logincr")
        gr[[k]]@fit <- dummy.fit       
    }

    return(gr)

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
                                 meanB=rep(0,2),varB=rep(1e10), nitt=50000) {

    require(MCMCglmm)
    #build appropriate size based covariates
    dataf$size2 <- dataf$size^2
    dataf$size3 <- dataf$size^3
    dataf$covariate <- as.factor(dataf$covariate)
    levels(dataf$covariate) <- 1:length(unique(dataf$covariate))
    if (length(grep("logsize",explanatoryVariables))>0) dataf$logsize <- log(dataf$size)

	#get rid of NAs
	dataf <- dataf[!is.na(dataf$size) & !is.na(dataf$sizeNext),]
	
	
    #build formula
    Formula<-paste('surv','~',explanatoryVariables,sep='')
    Formula <- as.formula(Formula)
    
    #fit survival model 
    #avoid fitting a prior unless you need to (takes longer with)
    if (sum(meanB!=0)>1 | sum(varB!=1e10)>1) { 
        fit<-MCMCglmm(Formula, data=dataf[!is.na(dataf$surv),],
                      verbose=FALSE,family="categorical",
                      prior=list(B=list(mu=meanB,V=diag(length(meanB))*varB)), nitt=nitt)
    } else { 
        fit<-MCMCglmm(Formula, data=dataf[!is.na(dataf$surv),], verbose=FALSE,family="categorical")        
    }
    dummy.fit <- glm(Formula, data=dataf,family=binomial)
    
    #create list survival models reflecting posterior
    sv <- list()
    for (k in 1:length(fit$Sol[,1])) {
        dummy.fit$coefficients <- fit$Sol[k,]
        sv[[k]] <-  new("survObj")
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
                            fec.constants=as.numeric(NA),
                            explanatoryVariables="size",
                            Family="gaussian",
                            Transform="none",
                            mean.offspring.size=NA,
                            var.offspring.size=NA,
                            offspring.splitter=data.frame(continuous=1),
                            fec.by.discrete=matrix(NA,nrow=0,ncol=0),nitt=50000) {

    require(MCMCglmm)

    ##warnings
    if (length(dataf$stage)==0) {
        print("Warning - no column named stage - assuming all continuous")
        dataf$stageNext <- dataf$stage <- rep("continuous", length(dataf[,1]))
        dataf$stage[is.na(dataf$size)] <- NA
        dataf$stageNext[is.na(dataf$sizeNext)] <- "dead"
    }
    
    if(ncol(offspring.splitter)>1 & (ncol(offspring.splitter)-1)!=ncol(fec.by.discrete)) {
        print("Warning - offspring splitter indicates more than just continuous stages. No fecundity by the discrete stages supplied in fec.by.discrete; assumed that is 0")
        fec.by.discrete <- matrix(0,col(offspring.splitter)-1,col(offspring.splitter)-1)
    }
    
    if (length(grep("covariate",explanatoryVariables))>0) {
        dataf$covariate <- as.factor(dataf$covariate)
        dataf$covariatenext <- as.factor(dataf$covariatenext)
        levels(dataf$covariate) <- 1:length(unique(dataf$covariate))
    }
    
    dataf$size2 <- dataf$size^2
    if (length(grep("logsize",explanatoryVariables))>0) dataf$logsize <- log(dataf$size)
    fecnames <- names(dataf)[grep("fec",names(dataf))]
    if (length(fecnames)>length(explanatoryVariables)) {
        print("number in explanatoryVariables not the same as the number of fecundity columns in the data file, taking first for all")
        explanatoryVariables <- rep(explanatoryVariables[1],length(fecnames))
    }
    if (length(fecnames)>length(Family)) {
        print("number of families not the same as the number of fecundity columns in the data file, taking first for all")
        Family <- rep(Family[1],length(fecnames))
    }
    if (length(fecnames)>length(Transform)) {
        print("number of transforms not the same as the number of fecundity columns in the data file, taking first for all")
        Transform <- rep(Transform[1],length(fecnames))
    }

    if (is.na(mean.offspring.size)) {
        offspringdata<-subset(dataf,is.na(dataf$stage)&dataf$stageNext=="continuous")
        mean.offspring.size <- mean(offspringdata$sizeNext)
        var.offspring.size <- var(offspringdata$sizeNext)
    }
 
	#get rid of NAs
	dataf <- dataf[!is.na(dataf$size) & !is.na(dataf$sizeNext),]
	
    fit <- dummy.fit <- list()
    
    for (i in 1:length(fecnames)) {

        if (Transform[i]=="log") dataf[,fecnames[i]] <- log(dataf[,fecnames[i]])
        if (Transform[i]=="sqrt") dataf[,fecnames[i]] <- sqrt(dataf[,fecnames[i]])
        dataf[!is.finite(dataf[,fecnames[i]]),fecnames[i]] <- NA

        Formula <- paste(fecnames[i],'~',explanatoryVariables[i],sep='')
        Formula <- as.formula(Formula)
         
        fit[[i]] <- MCMCglmm(Formula,
                        data=dataf[!is.na(dataf[,fecnames[i]]),], verbose=FALSE, nitt=nitt,family=Family[i])

        dummy.fit[[i]] <- glm(Formula,
                             data=dataf[!is.na(dataf[,fecnames[i]]),],family=Family[i])
        
    }

        
    #create list of growth models reflecing posterior
    fv <- list()
    for (k in 1:length(fit[[1]]$Sol[,1])) {
        fv[[k]] <-  new("fecObj")

        for (i in 1:length(fecnames)) { 
            dummy.fit[[i]]$coefficients <- fit[[i]]$Sol[k,]
            dummy.fit[[i]]$residuals <- rnorm(length(dummy.fit[[i]]$residuals),0,sqrt(fit[[i]]$VCV[k,1]))
            if (i==1) fv[[k]]@fit.fec1 <- dummy.fit[[i]]
            if (i==2) fv[[k]]@fit.fec2 <- dummy.fit[[i]]
            if (i==3) fv[[k]]@fit.fec3 <- dummy.fit[[i]]
            if (i==4) fv[[k]]@fit.fec4 <- dummy.fit[[i]]
            if (i==5) fv[[k]]@fit.fec5 <- dummy.fit[[i]]
            if (i==6) fv[[k]]@fit.fec6 <- dummy.fit[[i]]
            if (i==7) fv[[k]]@fit.fec7 <- dummy.fit[[i]]
            if (i==8) fv[[k]]@fit.fec8 <- dummy.fit[[i]]
            if (i==9) fv[[k]]@fit.fec9 <- dummy.fit[[i]]
            }
            
        fv[[k]]@fec.constants <- fec.constants
        fv[[k]]@mean.offspring.size <- mean.offspring.size
        fv[[k]]@var.offspring.size <- var.offspring.size
        fv[[k]]@offspring.splitter <- offspring.splitter
        fv[[k]]@fec.by.discrete <- fec.by.discrete
        fv[[k]]@Transform <- Transform 
    }

    return(fv)
    
}


# Function to take a list of growth and survival objects and make a list of Tmatrices
#
# Parameters - growObjList - a list of growth objects
#            - survObjList - a list of survival objects
#            - n.big.matrix - the number of bins
#            - minsize - the minimum size
#            - maxsize - the maximum size
#            - cov - is a discrete covariate considered
#            - env.mat - enviromental matrix for transition between
# 
# Returns    - a list of Tmatrices
makeListTmatrix <- function(growObjList,survObjList,
                            n.big.matrix,minsize,maxsize, cov=FALSE, env.mat=NULL) {

    if (length(growObjList)>length(survObjList)) { 
        survObjList <- sample(survObjList,size=length(growObjList),replace=TRUE)
    } else { 
        if (length(growObjList)<length(survObjList))  
            growObjList <- sample(growObjList,size=length(survObjList),replace=TRUE)
    }
    
    nsamp <- length(growObjList)
    
    Tmatrixlist <- list()
    for ( k in 1:length(growObjList)) { 
        if (!cov) {
            Tmatrixlist[[k]] <- create.IPM.Tmatrix(n.big.matrix = n.big.matrix, minsize = minsize, 
                                                   maxsize = maxsize, growObj = growObjList[[k]],
                                                   survObj = survObjList[[k]]) 
        } else {
            Tmatrixlist[[k]] <- create.compound.Tmatrix(n.env.class = length(env.mat[1,]),
                                                        n.big.matrix = n.big.matrix, minsize = minsize, 
                                                        maxsize = maxsize, envMatrix=env.mat,
                                                        growObj = growObjList[[k]],
                                                        survObj = survObjList[[k]])    
        }
    }
        
    return(Tmatrixlist)
}

# Function to take a list of growth and survival objects and make a list of Fmatrices

makeListFmatrix <- function(growObjList,survObjList,fecObjList,
                            n.big.matrix,minsize,maxsize, cov=FALSE, env.mat=NULL) {

    nsamp <- max(length(growObjList),length(survObjList),length(fecObjList))
    if (length(survObjList)<nsamp)  
        survObjList <- sample(survObjList,size=nsamp,replace=TRUE)
    if (length(fecObjList)<nsamp)  
        fecObjList <- sample(fecObjList,size=nsamp,replace=TRUE)
    if (length(growObjList)<nsamp)  
        growObjList <- sample(growObjList,size=nsamp,replace=TRUE)
  
    Fmatrixlist <- list()
    for ( k in 1:nsamp) {
        if (!cov) { 
            Fmatrixlist[[k]] <- create.IPM.Fmatrix(n.big.matrix = n.big.matrix, minsize = minsize, 
                                                   maxsize = maxsize, 
                                                   fecObj=fecObjList[[k]])
        } else {
            Fmatrixlist[[k]] <- create.compound.Fmatrix(n.env.class = length(env.mat[1,]),
                                                        n.big.matrix = n.big.matrix, minsize = minsize, 
                                                        maxsize = maxsize, envMatrix=env.mat,
                                                        fecObj=fecObjList[[k]])
        }

         Fmatrixlist[[k]] <-  Fmatrixlist[[k]]
    }
    return(Fmatrixlist)
}



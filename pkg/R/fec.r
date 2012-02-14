
# Create a generic fertility object #in base
setClass("fecObj",
         representation(fit.fec1 = "glm",
                        fit.fec2 = "glm",
                        fit.fec3 = "glm",
                        fit.fec4 = "glm",
                        fit.fec5 = "glm",
                        fit.fec6 = "glm",
                        fit.fec7 = "glm",
                        fit.fec8 = "glm",
                        fit.fec9 = "glm",
                        fec.constants = "numeric",
                        offspring.splitter = "data.frame",
						mean.offspring.size = "numeric",
                        var.offspring.size = "numeric",
						fec.by.discrete = "matrix",
						Transform = "character")
         )


makeFecObjGeneral <- function(dataf,
                              fec.constants=1,
                              explanatoryVariables="size+size2",
                              Family="gaussian",
                              Transform="none",
                              mean.offspring.size=NA,
                              var.offspring.size=NA,
                              offspring.splitter=data.frame(continuous=1),
                              fec.by.discrete=matrix(NA,nrow=0,ncol=0)){

    ##warnings
    if (length(dataf$stage)==0) {
        print("Warning - no column named stage - assuming all continuous")
        dataf$stagenext <- dataf$stage <- rep("continuous", length(dataf[,1]))
        dataf$stage[is.na(dataf$size)] <- NA
        dataf$stagenext[is.na(dataf$sizenext)] <- "dead"
    }
    
    if(ncol(offspring.splitter)>1 & (ncol(offspring.splitter)-1)!=ncol(fec.by.discrete)) {
        print("Warning - offspring splitter indicates more than just continuous stages. No fertility by the discrete stages supplied in fec.by.discrete; assumed that is 0")
        fec.by.discrete <- matrix(0,col(offspring.splitter)-1,col(offspring.splitter)-1)
    }

   if (length(grep("covariate",explanatoryVariables))>0) {
       dataf$covariate <- as.factor(dataf$covariate)
       dataf$covariatenext <- as.factor(dataf$covariatenext)
       levels(dataf$covariate) <- 1:length(unique(dataf$covariate))
   }
    
    f1 <- new("fecObj")
    dataf$size2 <- dataf$size^2
   if (length(grep("logsize",explanatoryVariables))>0) dataf$logsize <- log(dataf$size)

    fecnames <- names(dataf)[grep("fec",names(dataf))]
    if (length(fecnames)>length(explanatoryVariables)) {
        print("number in explanatoryVariables not the same as the number of fecundity columns in the data file, taking first for all")
        explanatoryVariables <- rep(explanatoryVariable[1],length(fecnames)) }
    if (length(fecnames)>length(Family)) {
        print("number in families not the same as the number of fecundity columns in the data file, taking first for all")
        family <- rep(Family[1],length(fecnames)) }
    
    for (i in 1:length(fecnames)) {
        if (Transform[i]=="log") dataf[,fecnames[i]] <- log(dataf[,fecnames[i]])
        if (Transform[i]=="sqrt") dataf[,fecnames[i]] <- sqrt(dataf[,fecnames[i]])
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
        offspringdata<-subset(dataf,is.na(dataf$stage)&dataf$stagenext=="continuous")
        mean.offspring.size <- mean(offspringdata$sizenext)
        var.offspring.size <- var(offspringdata$sizenext) }
    f1@fec.constants <- fec.constants
    f1@mean.offspring.size <- mean.offspring.size
    f1@var.offspring.size <- var.offspring.size
    f1@offspring.splitter <- offspring.splitter / max(offspring.splitter) 
    f1@fec.by.discrete <- fec.by.discrete
    f1@Transform <- Transform
    return(f1)
}
		 

		 
create.IPM.Fmatrix <- function(n.env.class = 1,
                               n.big.matrix = 50,
                               minsize = -1,
                               maxsize = 50,
                               chosen.cov = 1,
                               fecObj,
                               integrate.type="midpoint",
                               correction="none") {
    
    # boundary points b and mesh points y
    b<-minsize+c(0:n.big.matrix)*(maxsize-minsize)/n.big.matrix;
    y<-0.5*(b[1:n.big.matrix]+b[2:(n.big.matrix+1)]);
    
    # step size for mid point rule, see equations 4 and 5
    h<-y[2]-y[1]
    #size<-y
    newd <- data.frame(size=y,size2=y^2,size3=y^3,
                       covariate=as.factor(rep(chosen.cov,length(y))))
    if (length(grep("logsize",fecObj@fit.fec1$formula))>0&length(grep("logsize",fecObj@fit.fec2$formula))>0&
        length(grep("logsize",fecObj@fit.fec3$formula))>0
        &length(grep("logsize",fecObj@fit.fec4$formula))>0&length(grep("logsize",fecObj@fit.fec5$formula))>0&length(grep("logsize",fecObj@fit.fec6$formula))>0
        &length(grep("logsize",fecObj@fit.fec7$formula))>0&length(grep("logsize",fecObj@fit.fec8$formula))>0
        &length(grep("logsize",fecObj@fit.fec9$formula))>0) { newd$logsize <- log(y)}
    
    fec.values <- matrix(c(rep(1,9),fecObj@fec.constants),ncol=n.big.matrix,nrow=(9+length(fecObj@fec.constants)))
    if (is.null(fecObj@fit.fec1)==FALSE) fec.values[1,]<-predict(fecObj@fit.fec1,newd,type="response")
    if (is.null(fecObj@fit.fec2)==FALSE) fec.values[2,]<-predict(fecObj@fit.fec2,newd,type="response")
    if (is.null(fecObj@fit.fec3)==FALSE) fec.values[3,]<-predict(fecObj@fit.fec3,newd,type="response")
    if (is.null(fecObj@fit.fec4)==FALSE) fec.values[4,]<-predict(fecObj@fit.fec4,newd,type="response")
    if (is.null(fecObj@fit.fec5)==FALSE) fec.values[5,]<-predict(fecObj@fit.fec5,newd,type="response")
    if (is.null(fecObj@fit.fec6)==FALSE) fec.values[6,]<-predict(fecObj@fit.fec6,newd,type="response")
    if (is.null(fecObj@fit.fec7)==FALSE) fec.values[7,]<-predict(fecObj@fit.fec7,newd,type="response")
    if (is.null(fecObj@fit.fec8)==FALSE) fec.values[8,]<-predict(fecObj@fit.fec8,newd,type="response")
    if (is.null(fecObj@fit.fec9)==FALSE) fec.values[9,]<-predict(fecObj@fit.fec9,newd,type="response")
    if (length(grep("log",fecObj@Transform))>0) for (i in grep("log",fecObj@Transform)) fec.values[i,]<-exp(fec.values[i,])
    if (length(grep("sqrt",fecObj@Transform))>0) for (i in grep("sqrt",fecObj@Transform)) fec.values[i,]<-(fec.values[i,])^2
    prod.fec.values<-apply(fec.values,2,prod)
    
    tmp<-dnorm(y,fecObj@mean.offspring.size,sqrt(fecObj@var.offspring.size))*h
    if (correction=="constant") tmp<-tmp/sum(tmp)
    to.cont<-tmp%*%t(as.numeric(fecObj@offspring.splitter["continuous"])*prod.fec.values)
    get.matrix <- to.cont
    ndisc <- length(fecObj@offspring.splitter)-1
	
    if (ndisc>0) {
        to.discrete <- as.numeric(fecObj@offspring.splitter)[1:ndisc]%*%t(prod.fec.values)
        from.discrete <- matrix(0,ncol=ndisc,nrow=ndisc+n.big.matrix)
        if (length(fecObj@fec.by.discrete)>0)
            from.discrete <- c(as.numeric(fecObj@offspring.splitter)[1:ndisc],
                               as.numeric(fecObj@offspring.splitter)[ndisc+1]*tmp)%*%t(fecObj@fec.by.discrete)
        

        get.matrix <- cbind(from.discrete,rbind(to.discrete,to.cont)) }
    
    rc <- new("IPM.matrix",
              n.discrete = ndisc,
              n.env.class = 1, 
              n.big.matrix = n.big.matrix,
              nrow = 1*n.big.matrix+ndisc,
              ncol =1*n.big.matrix+ndisc,
              meshpoints = y,
              env.index = rep(1:n.env.class,each=n.big.matrix))
    rc[,] <-get.matrix   
    
    return(rc)
}


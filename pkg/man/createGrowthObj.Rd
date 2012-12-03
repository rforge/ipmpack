\name{createGrowthObj}
\alias{createGrowthObj}
\alias{createSurvObj}
\alias{createFecObj}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Function to create growth or survival objects with user-defined parameters 
}
\description{
supplied with a Formula, coefficients, and for growth a standard
deviation, and for fertility, a range of other necessary arguments (for many of which defaults may be sufficient), creates a growth, survival or fertility object 
}
\usage{

createGrowthObj(Formula=sizeNext~size, coeff=c(1,1), sd=1)
createSurvObj(Formula=surv~size, coeff=c(1,1))
createFecObj(Formula=list(fec1~size,fec2~size+size2), 
							coeff=list(c(1,1),c(1,1,1)),
							Family = c("gaussian","binomial"),
							Transform = c("log","none"),
							meanOffspringSize = NA, sdOffspringSize = NA, 
							offspringSplitter = data.frame(continuous = 1), 
							vitalRatesPerOffspringType = data.frame(NA), 
							fecByDiscrete = data.frame(NA), 
							offspringSizeExplanatoryVariables = "1",
							fecConstants = data.frame(NA)){ 


}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{Formula}{an object of class formula for createGrowthObj or createSurvObj, or a list of formulas for createFecObj}
    \item{coeff}{a numeric vector of length the number of coefficients
      supplied by Formula, or a list for createFecObj}
    \item{sd}{a numeric vector of length 1, for createGrowthObj}
 
  \item{fecConstants}{for createFecObj, a list containing the value by which each of the fecundity rates will be multiplied \emph{in the order defined by the order in \code{Formula}}.  This data frame adjusts the probability of establishment of seeds or other stages in  sexual reproduction that are not explicitly incorporated  for each parent (e.g., 25\% of seeds across all individuals germinate).  The default is NA if no constants are used (equivalent to multiplying by 1).}
     \item{Family}{for createFecObj, a character vector containing the names of the families to be used for the glms, e.g., binomial, poisson, etc. Again, these must appear  \emph{in the order defined by \code{Formula}} }
  \item{Transform}{for createFecObj, a character vector containing the names of the
    transforms to be used for the response variables, e.g., log, sqrt, -1, etc. Again, these must appear  \emph{in the order defined by \code{Formula}} }
  \item{meanOffspringSize}{for createFecObj, numeric vector, defining mean offspring size. Defaults to NA, in which case the function will use to data to assess the mean offspring size according to the relationship defined in offspringSizeExplanatoryVariables (which either simply fits a mean, or may fit more complex relationships linking maternal size to offspring size).}
  \item{sdOffspringSize}{for createFecObj, numeric vector, defining standard deviation of offspring size. Defaults to NA, in which case the function will use the data to assess the standard deviation of offspring size; as described for meanOffspringSize}
  \item{offspringSplitter}{for createFecObj, dataframe with values defining the number of
    offspring going into the indicated offspring category; will be re-scaled to
    sum to 1 within the function. This argument needs to be entered as a data.frame, and the names in the data.frame need to precisely match the used stage names in the data file.}
  \item{vitalRatesPerOffspringType}{for createFecObj, dataframe defining which fecundity rates (both functions and constants) apply to which offspring category. This only needs to be specified when some fecundity rates do not apply to all offspring categories. The offspring categories in the column names of this dataframe should match those in the offspringSplitter exactly. The row names of the dataframe should match the fecundity column names in the data file and the supplied fecundity constants, in that order. In the dataframe, a '1' indicates that a fecundity rate applies to an offspring category, while a '0' indicates an omission. For instance, establishment and seedling survival rates may be applicable to seedlings, but not to seeds that go into a seedbank (depending on the life cycle and definition of vital rates).}
  \item{fecByDiscrete}{for createFecObj, data.frame defining number of offspring produced by each
    discrete class ; defaults to 0. If specified, ALL discrete classes MUST appear in alphabetical order, 
    so NO "continuous". e.g. fecByDiscrete=data.frame(dormant=0,seedAge1=4.2,seedOld=0)}    
  \item{offspringSizeExplanatoryVariables}{for createFecObj, a character defining the relationship defining offspring size; the default is "1", indicating simply fitting a mean and a variance; alternatives would including defining offspring size as a function of maternal size (i.e., offspringSizeExplanatoryVariables="size") or more complex polynomials of size (i.e., offspringSizeExplanatoryVariables="size+size2"). The corresponding relationship is fitted to the data contained in dataf, taking as the response variable the column "sizeNext" in dataf for rows where the column "offspringNext" is equal to "sexual" and the column "stageNext" is equal to "continuous". }  
    
    
}
  
\details{
  These functions create growth and survival objects where direct
  fitting is not desired; note that currently they are not capable of
  managing Formula that encompass factors. 
}

\value{
  an object of class growthObj / survObj
}
\author{
C. Jessica E. Metcalf, Sean M. McMahon, Roberto Salguero-Gomez & Eelke Jongejans.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{makeSurvObj}}, \code{\link{makeGrowthObj}},\code{\link{makeFecObj}}, \code{\link{coerceGrowthObj}},\code{\link{coerceSurvObj}}

}
\examples{
dff<-generateData()

#for growth
gr1 <- createGrowthObj(Formula=sizeNext~size,coeff=c(0.1,0.9),sd=1)

#halve the slope
gr2 <- createGrowthObj(Formula=sizeNext~size,coeff=c(0.1,0.45),sd=1)


#for survival
sv1 <- createSurvObj(Formula=surv~size,coeff=c(-10,3))


#for fertility 
fv1<-createFecObj(Formula=list(fec1~size,fec1~size),
				coeff=list(c(0.1,0.9),c(-2,10)),
				Transform=c("log","none"),
				Family=c("gaussian","binomial"),
				meanOffspringSize=1,sdOffspringSize=1)


}

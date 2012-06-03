\name{makePostFecObjs}
\alias{makePostFecObjs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Make list of fecundity objects corresponding to posterior of Bayes fit
}
\description{
   Allows a series of different glms to be fit all on the way to
   fecundity, e.g., probability of flowering, number of flower heads
   produced, etc. using a Bayesian framework via the MCMCglmm
   package. Needs work on priors.  

}
\usage{
makePostFecObjs(dataf,
		fecConstants=data.frame(NA),fecNames=NA,
		explanatoryVariables="size",
		Family="gaussian",
		Transform="none",
		meanOffspringSize=NA,
		sdOffspringSize=NA,
		offspringSplitter=data.frame(continuous=1),
		offspringTypeRates=data.frame(NA),
		fecByDiscrete=data.frame(NA),
		offspringSizeExplanatoryVariables="1",
		burnin=3000,nitt=50000) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{dataf}{a dataframe with columns "size", "sizeNext", "stage",
      "stageNext", "fec1", "fec2" or other names containing fec; up to
      nine are currently possible (and note that ANY name containing
      "fec" will be included and these MUST NOT be transformed via log,
      etc, since the fecundity object produced will assume that it has total reproductive output.)}
 \item{fecNames}{a character vector containing the names of the columns of the dataframe from which fecundity data will be analysed (one column at a time). If the fecNames argument is not specified (as by default, NA), fecNames will be filled by the names of all column names that contain "fec".} 
   \item{fecConstants}{a dataframe containing the value by which
    each of the fecundity transforms will be multplied \emph{in the
      order defined by the alphabetical order of column names
      containing "fec"}; these might capture for example the probability
    of establishment of seeds; default is NA if no constants are used.}
  \item{explanatoryVariables}{a character vector containing formulas describing the desired
    explanatory variables (interactions, etc) in classical R
    style, i.e. separated by `+', `*', `:'. Possible covariates include `size', 'size2' (size^2), `size3'
    (size^3),`logsize' (log(size)), and `covariate'. Again, these must appear  \emph{in the
      order defined by the alphabetical order of column names
      containing "fec"}}
  \item{Family}{a character vector containing the names of the
    families to be used for the glms, e.g., binomial, poisson, etc. Again, these must appear  \emph{in the
      order defined by the alphabetical order of column names
      containing "fec"} }
  \item{Transform}{a character vector containing the names of the
    transforms to be used for the response variables, e.g., log, sqrt,
    etc. Again, these must appear  \emph{in the
      order defined by the alphabetical order of column names
      containing "fec"} }
  \item{meanOffspringSize}{numeric vector, defining mean offspring size. Defaults to NA, in which case the
    function will use to data to assess the mean offspring size}
  \item{sdOffspringSize}{numeric vector, defining sd of offspring size. Defaults to NA, in which case the
    function will use to data to assess the mean offspring size}
  \item{offspringSplitter}{dataframe, with column headings for EVERY
    category, so always at minimum containing a column called "continuous" but also
    including every discrete category that will enter into the Fmatrix (and, if
    compatibility is required, the Tmatrix), with values defining the number of
    offspring going into each offspring category; will be re-scaled to
    sum to 1 within the function. THESE MUST appear in alphabetical
    order.}
  \item{offspringTypeRates}{dataframe defining which fecundity rates (both functions and constants) apply to which offspring category. This only needs to be 
    specified when some fecundity rates do not apply to all offspring categories. The offspring categories in the column names of this dataframe should match those 
    in the offspringSplitter exactly. The row names of the dataframe should match the fecundity column names in the data file and the supplied fecundity constants, 
    in that order. In the dataframe, a '1' indicates that a fecundity rate applies to an offspring category, while a '0' indicates an omission. For instance, establishment
    and seedling survival rates may be applicable to seedlings, but not to seeds that go into a seedbank (depending on the life cycle and definition of vital rates).}
  \item{fecByDiscrete}{data.frame defining number of offspring produced by each
    discrete class ; defaults to 0. If specified, ALL discrete classes MUST appear in alphabetical order, so NO "continuous". e.g. fecByDiscrete=data.frame(dormant=0,seedAge1=4.2,seedOld=0)}     
   \item{offspringSizeExplanatoryVariables}{a character defining the relationship defining offspring size; the default is "1", indicating simply fitting a mean and a variance; alternatives would including defining offspring size as a function of maternal size (i.e., offspringSizeExplanatoryVariables="size") or more complex polynomials of size (i.e., offspringSizeExplanatoryVariables="size+size2"). The corresponding relationship is fitted to the data contained in dataf, taking as the response variable the column "sizeNext" in dataf for rows where the column "offspringNext" is equal to "sexual" and the column "stageNext" is equal to "continuous". }  
  \item{burnin}{
    numeric defining burn in for MCMCglmm
}
  \item{nitt}{
    number of iterations
    }
}
\details{
NOTE: priors not implemented; also uncertainty / variation in mean seedling size,
variance in seedling size is not implemented. Requires MCMCglmm.
Function does not return output of Bayes analysis which might be useful for assessment; in development. 

}
\value{
  list of fecundity objects
}
\author{
C. Jessica E. Metcalf, Sean M. McMahon, Roberto Salguero-Gomez & Eelke Jongejans
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{makeListIPMs}}, \code{\link{makePostGrowthObjs}}, \code{\link{makeFecObj}}
}
\examples{
dff <- generateData()
fvlist <- makePostFecObjs(dataf=dff,
	explanatoryVariables="size+size2+covariate",
                               burnin=100,nitt = 500)
                               
fvlist[[1]]                               
                               

}

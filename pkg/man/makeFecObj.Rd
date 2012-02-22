\name{makeFecObj}
\alias{makeFecObj}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Function to build fecundity objects 
}
\description{
  Allows a series of different glms to be fit all on the way to fecundity, e.g., probability of
  flowering, number of flower heads produced, etc. 
}
\usage{
makeFecObj(dataf, fec.constants=as.numeric(NA),
	explanatoryVariables="size", Family="gaussian", 
	Transform="none", fecNames=NA, mean.offspring.size=NA, 
	var.offspring.size=NA, offspring.splitter=data.frame(continuous=1), 
	fec.by.discrete=data.frame(NA))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{dataf}{a dataframe with columns "size", "sizeNext", "stage", "stageNext", and any number of columns with fecundity data. Unless the fecundity columns are specified in the argument fecNames, fecundity data will be extracted from columns "fec1", "fec2" or any other names containing fec; (note that ANY name containing "fec" will be included). Fecundity data MUST NOT be transformed via log, etc, since the fecundity object produced will assume that it has total reproductive output.)}
  \item{fec.constants}{a numeric vector containing the value by which each of the product of the fecundity rates will be multiplied \emph{in the order defined by the alphabetical order of column names containing "fec"}; these might capture for example the probability of establishment of seeds or other steps in the sexual reproductive pathway that are not measured for each parent; default is NA if no constants are used.}
  \item{explanatoryVariables}{a character vector containing formulas describing the desired
    explanatory variables (interactions, etc) in classical R style, i.e. separated by `+', `*', `:'. Possible covariates include `size', 'size2' (size^2), `size3'  (size^3),`logsize' (log(size)), and `covariate'. Again, these must appear  \emph{in the order defined by the alphabetical order of column names containing "fec" or in the order of column names specified in fecNames}}
  \item{Family}{a character vector containing the names of the families to be used for the glms, e.g., binomial, poisson, etc. Again, these must appear  \emph{in the
      order defined by the alphabetical order of column names 
      containing "fec" or in the order of column names specified in fecNames} }
  \item{Transform}{a character vector containing the names of the
    transforms to be used for the response variables, e.g., log, sqrt, -1, etc. Again, these must appear  \emph{in the order defined by the alphabetical order of column names
      containing "fec" or in the order of column names specified in fecNames} 
  \item{fecNames}{a character vector containing the names of the columns of the dataframe from which fecundity data will be analysed (one column at a time). If the fecNames argument is not specified (as by defaul, NA), fecNames will be filled by the names of all column names that contain "fec".} }
  \item{mean.offspring.size}{numeric vector, defining mean offspring size. Defaults to NA, in which case the function will use to data to assess the mean offspring size}
  \item{var.offspring.size}{numeric vector, defining mean offspring size. Defaults to NA, in which case the function will use to data to assess the mean offspring size}
  \item{offspring.splitter}{dataframe with values defining the number of
    offspring going into the indicated offspring category; will be re-scaled to
    sum to 1 within the function. This argument needs to be entered as a data.frame, and the names in the data.frame need to precisely match the used stage names in the data file.}
  \item{fec.by.discrete}{data.frame defining number of offspring produced by each
    discrete class ; defaults to 0. If specified, ALL discrete classes MUST appear in alphabetical order, 
    so NO "continuous". e.g. fec.by.discrete=data.frame(dormant=0,seed.age.1=4.2,seed.old=0)}     
}
  
\details{
  This function fits a suite of subfunctions of fecundity towards creating
  a Fecundity transition  projection model; e.g., the probability of flowering as a function of size, the number of seeds produced as a function of size, etc. Users can
  define the functional form of each relationship, as well as the
  distribution and any transforms. There is also a possibility of
  defining reproduction from discrete sizes into each of
  the subfunction outcomes; defined in the matrix fec.by.discrete. NOTE
  THAT by default all column headings containing "fec" will be considered as
  contributing to the fecundity product, unless otherwise specified in the argument fecNames. 
}

\value{
  an object of class fecObj
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
C. Jessica E. Metcalf, Sean McMahon, Roberto Salguero-Gomez & Eelke Jongejans
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{

}
\examples{
dff<-generateDataDiscrete()
fv1 <- makeFecObj(dff,offspring.splitter=data.frame(continuous=0.2,seed.age.1=0.8),
	fec.by.discrete=data.frame(dormant=0,seed.age.1=4.2,seed.old=0))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

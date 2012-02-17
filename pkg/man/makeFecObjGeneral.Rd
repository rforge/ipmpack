\name{makeFecObjGeneral}
\alias{makeFecObjGeneral}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Function to build fertility objects 
}
\description{
  Allows a series of different glms to be fit all on the way to fertility, e.g., probability of
  flowering, number of flower heads produced, etc. 
}
\usage{
makeFecObjGeneral(dataf, fec.constants,explanatoryVariables,Family,
Transform,mean.offspring.size, var.offspring.size, offspring.splitter,
fec.by.discrete)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{dataf}{a dataframe with columns "size", "sizenext", "stage",
      "stagenext", "fec1", "fec2" or other names containing fec; up to
      nine are currently possible (and note that ANY name containing
      "fec" will be included and these MUST NOT be transformed via log,
      etc, since the fertility object produced will assume that it has total reproductive output.)}
  \item{fec.constants}{a numeric vector containing the value by which
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
  \item{mean.offspring.size}{numeric vector, defining mean offspring size. Defaults to NA, in which case the
    function will use to data to assess the mean offspring size}
  \item{var.offspring.size}{numeric vector, defining mean offspring size. Defaults to NA, in which case the
    function will use to data to assess the mean offspring size}
  \item{offspring.splitter}{dataframe, with column headings for EVERY
    category, so always at minimum containing a column called "continuous" but also
    including every discrete category that will enter into the Fmatrix (and, if
    compatibility is required, the Tmatrix), with values defining the number of
    offspring going into each offspring category; will be re-scaled to
    sum to 1 within the function. THESE MUST appear in alphabetical
    order.}
  \item{fec.by.discrete}{matrix defining number produced by each
    discrete class (rows) into each of the "fec" types (columns);
    defaults to 0}    
}
  
\details{
  This function fits a suite of subfunctions of fertility towards creating
  a Fecundity transition  projection model; e.g., the probability of flowering as a function of size,
  the number of seeds produced as a function of size, etc. Users can
  define the functional form of each relationship, as well as the
  distribution and any transforms. There is also a possibility of
  defining reproduction from discrete sizes into each of
  the subfunction outcomes; defined in the matrix fec.by.discrete. NOTE
  THAT all column headings containing "fec" will be considered as
  contributing to the fecundity product; re-name data-frame
  accordingly. 
}

\value{
  an object of class fecObj
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{

}
\examples{
dff<-generateDataDiscrete()
fv1 <- makeFecObjGeneral(dff)
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

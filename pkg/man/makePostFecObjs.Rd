\name{makePostFecObjs}
\alias{makePostFecObjs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Make list of fecundity objects corresponding to posterior of Bayes fit
}
\description{
   Allows a series of different glms to be fit all on the way to
   fertility, e.g., probability of flowering, number of flower heads
   produced, etc. using a Bayesian framework via the MCMCglmm
   package. Needs work on priors.  

}
\usage{
makePostFecObjs(dataf,fec.constants,explanatoryVariables,Family,Transform,mean.offspring.size,var.offspring.size,offspring.splitter,fec.by.discrete,nitt)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{dataf}{a dataframe with columns "size", "sizeNext", "stage",
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
  \item{nitt}{
    number of iterations
    }
}
\details{
NOTE: priors not implemented; also uncertainty / variation in mean seedling size,
variance in seedling size is not implemented. Requires MCMCglmm.
}
\value{
  list of fertility objects
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
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
dff<-generateData()
fvlist <- makePostFecObjs(dataf=dff,explanatoryVariables="size+size2+covariate",
                               nitt = 5000)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

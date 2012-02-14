\name{makeDiscreteTrans}
\alias{makeDiscreteTrans}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Functions to build matrix for transitions between discrete stages
  and between discrete and continuous stages}
\description{
  Function that takes a dataframe reflecting a population that includes continuous and discrete
  classes of individuals, and builds a transition matrix from this for movement
  between discrete and continuous stages (providing just a single value
  for the continuous stages; the detail of movement between continuous
  stages (e.g., sizes) is generated elsewhere). This object can then be used as an
  argument in the function to "create.IPM.Tmatrix" to build a T matrix
  that contains both discrete and continuous strategies. 

  The data-frame must contain columns  "stage" and "stagenext"
  with values of the names of the discrete classes and the name
  "continuous" where appropriate, in the current and
  subsequent time-step. Continuous categories are also defined by their
  measurements, contained in "size" and "sizenext" (necessary for
  defining the mean and variance of the size of individuals emerging
  from discrete classes)
}
\usage{
makeDiscreteTrans(dataf)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataf}{a dataframe with columns "stage", "stagenext" (both
    character vectors containing either names of discrete stages or
  "continuous" or "dead"), "size",
    "sizenext", and "surv" (continuous and booean variables respectively)}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  an object of class "DiscreteTrans"
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
\code{\link{create.IPM.Tmatrix}}

}
\examples{
dff <- generateDataDiscrete()
discTrans <- makeDiscreteTrans(dff)
create.IPM.Tmatrix(n.big.matrix = 10,growObj=makeGrowthObjGeneral(dff),survObj=makeSurvObjGeneral(dff),DiscreteTrans=discTrans)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

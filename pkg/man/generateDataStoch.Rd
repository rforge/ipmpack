\name{generateDataStoch}
\alias{generateDataStoch}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Function to generate data with stochastic covariates
}
\description{
  Simulates growth, survival, fertility to create a data-frame of the
  form required by the functions and methods used in the package
}
\usage{
generateDataStoch()
}
%- maybe also 'usage' for other objects documented here.
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  a dataframe with headings size (continuous variable, indicating current size),
  sizenext (continuous variable, indicating size in the next time step), surv
  (boolean, indicating whether individual survived or not to the next time-step), covariate1
  covariate2,  covariate3, ..., fec (continuous variable, indicating fecundity).
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
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
dff<-generateDataStoch()
head(dff)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

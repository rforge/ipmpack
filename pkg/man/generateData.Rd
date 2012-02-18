\name{generateData}
\alias{generateData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Function to generate random data in the form used by IPMpack
}
\description{
  Simulates growth, survival and fertility to create a data-frame of the form required by the functions and methods used in IPMpack. Note that the data created includes only one continuous covariate structuring the population; user can include additional covariates if necessary.
}
\usage{
generateData()
}
%- maybe also 'usage' for other objects documented here.
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  a dataframe with headings size (continuous variable, indicating current size),
  sizeNext (continuous variable, indicating size in the next time step), surv
  (boolean, indicating whether individual survived or not to the next time-step), covariate
  (discrete covariate), covariatenext (discrete covariate in the next
  time-step),  fec (continuous variable, indicating fecundity).
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
\code{\link{generateDataStoch}},\code{\link{generateDataDiscrete}}
}
\examples{
dff<-generateData()
head(dff)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

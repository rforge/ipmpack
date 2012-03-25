\name{generateData}
\alias{generateData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Generates random data in the form used by IPMpack.
}
\description{
  Simulates growth, survival and fecundity to create a dataframe of the form required by the functions and methods used in IPMpack. Demographic stage data is only continuous.
}
\usage{
generateData(nSamp=1000)
}
\arguments{
	\item{nSamp}{
		number of samples desired in the base population, defaults to 1000
	}
}
%- maybe also 'usage' for other objects documented here.
\details{
Note that the random data created includes only one discrete covariate structuring the population (e.g., habitat type [shaded, open], age [0, 1, 2, 3 years old...], time since fire [0, 1, 2, 3 years...], etc); user can include into the data frame additional discrete covariates if necessary. To include continuous covariates (e.g., precipitation, temperature, light incidence, hurricane recurrence interval, etc) in randomized data, refer to generateDataStoch(). For randomized datasets with complex life cycles that include discrete stages (e.g., seedbank), refer to generateDataDiscrete().
}
\value{
  A dataframe with headings:
  	- "size": continuous variable, indicating current size.
  	- "sizeNext" continuous variable, indicating size in the next time step.
  	- "surv": boolean, indicating whether individual survived or not to the next time step.
  	- "covariate": discrete covariate.
  	- "covariateNext": discrete covariate in the next
  time step.
  	- "fec": continuous variable, indicating fecundity.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
C. Jessica E. Metcalf, Sean McMahon, Roberto Salguero-Gomez & Eelke Jongejans.
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{generateDataStoch}},\code{\link{generateDataDiscrete}}
}
\examples{
dff <- generateData(nSamp=2000)
head(dff)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

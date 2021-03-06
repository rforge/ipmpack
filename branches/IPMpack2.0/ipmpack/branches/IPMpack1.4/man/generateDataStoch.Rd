\name{generateDataStoch}
\alias{generateDataStoch}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Generates data with stochastic covariates.
}
\description{
  Simulates growth, survival, fecundity to create a dataframe of the
  form required by the functions and methods used in IPMpack.
}
\usage{
generateDataStoch()
}
%- maybe also 'usage' for other objects documented here.
\value{
  A dataframe with headings:
  	- "size": continuous variable, indicating current size.
  	- "sizeNext": continuous variable, indicating size in the next time step.
  	- "surv": boolean, indicating whether individual survived or not to the next time step.
  	- "covariate1", "covariate2", "covariate3", ...: discrete or continuous covariates.
  	- "fec": continuous variable, indicating fecundity, that is, per-capita contributions via sexual reproduction.
}
\references{
Caswell. 2001. Matrix population models: construction, analysis, and interpretation. 2nd ed. Sinauer. p377-502.

Rees & Ellner. 2009. Integral projection models for populations in temporally varying environments. The American Naturalist 79, p575-594.

}
\author{
C. Jessica E. Metcalf, Sean M. McMahon, Roberto Salguero-Gomez & Eelke Jongejans
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{generateData}},\code{\link{generateDataDiscrete}}
}
\examples{
dff <- generateDataStoch()
head(dff)
}

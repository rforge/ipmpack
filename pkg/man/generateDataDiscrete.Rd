\name{generateDataDiscrete}
\alias{generateDataDiscrete}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Function to generate data including discrete stages (here seed stages)
}
\description{
  Simulates growth, survival, fertility to create a data-frame of the
  form required by the functions and methods used in the package,
  including three possible discrete stages, "dormant", "seed.age.1", "seed.old"
}
\usage{
generateDataDiscrete()
}
%- maybe also 'usage' for other objects documented here.
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  a dataframe with headings "size" (continuous variable, indicating current size),
  "sizenext" (continuous variable, indicating size in the next time step), "surv"
  (boolean, indicating whether individual survived or not to the next
  time-step, whatever stage the individual ends up in), "fec" (continuous
  variable, indicating fecundity), "stage" (character vector, containing names of the discrete stages in
  that time-step, or "continuous"), "stagenext" (character vector, containing names of the discrete stages in
  the following time-step, or "continuous"), "number" corresponds to the number of
  individuals moving between stages - so "number" is 1 for all movements
  out of the "continous" stage, and is >1 for all movements out of
  discrete stages - this is just so that you don't need to have an
  individual line for every movement between discrete stages. 
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
Note that the data-frame generated contains no covariates (either discrete or continuous) 
and so may be inappropriate for some functions. %%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{generateDataStoch}},\code{\link{generateData}}
}
\examples{
dff<-generateDataDiscrete()
head(dff)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

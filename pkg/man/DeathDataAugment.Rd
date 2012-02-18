\name{DeathDataAugment}
\alias{DeathDataAugment}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to augment data on mortality (where rare)
}
\description{
If mortality is rare (e.g., in large trees) this method allows the
proportion dead above a chosen threshold to be adjusted by extending
the data-set (e.g., using a prior).
}
\usage{
DeathDataAugment(dataf, size.thresh, prop.dead)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataf}{a dataframe with columns `size',
    `sizenext', and `surv' (`size' is size at t, `sizenext' is size at t+1,
    and `surv' is 0s and 1s, indicating if the individual
    survived or not); facultatively, dataf may include `covariate'
    and `covariatel' for discrete covariates, indicating values at t,
    and at t+1, respectively; these must take values of sequential integers, starting at
    `1'; alternatively, dataf may contain columns 
    `covariate1', `covariate2', etc., for continuous covariates. }
  \item{size.thresh}{
    the size threshold above which prop.dead is required to have died
  }
  \item{prop.dead}{
    the required proportion for sizes> size.thresh expected to have died
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
a data-frame with the same headings as dataf but with extra rows to
capture extra mortality events. 
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
%dff<-generateData()
%dffnew<-DeathDataAugment(dff,size.thresh=10,prop.dead=0.5)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

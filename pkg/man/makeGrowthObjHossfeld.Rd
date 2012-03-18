\name{makegrowthObjHossfeld}
\alias{makegrowthObjHossfeld}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to make a Hossfeld Growth Object
}
\description{
Takes a data-frame with at minimum columns size and sizeNext; and fits a Hossfeld type growth function
to increment - if length(dataf$incr) is zero, it will calculate dataf$increment as the difference between 
size and sizeNext; otherwise it will take the column provided
}
\usage{
makegrowthObjHossfeld(dataf)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataf}{a dataframe with columns `size' and
    `sizeNext'(`size' is size at t, `sizeNext' is size at t+1);
	}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns a growth object of class growthObjHossfeld for which growth methods exist. 
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
dff <- generateData()
gr1 <- makegrowthObjHossfeld(dff)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

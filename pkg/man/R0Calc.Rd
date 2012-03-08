\name{R0Calc}
\alias{R0Calc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculates net reproductive rate (R0) from an IPM.
}
\description{
Estimates lifetime reproductive success from a full IPM, including survival, growth and fecundity. 
}
\usage{
R0Calc(Tmatrix, Fmatrix)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Tmatrix}{
    a matrix (not necessarily of class IPMmatrix).
}
  \item{Fmatrix}{
    a matrix (not necessarily of class IPMmatrix).
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  numeric
}
\references{
Caswell, 2001. Matrix population models: analysis, construction and
interpretation. 2nd ed. Sinauer. p126.
}
\author{
C. Jessica E. Metcalf, Sean McMahon, Roberto Salguero-Gomez & Eelke Jongejans.
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{createIPMTmatrix}},\code{\link{createIPMFmatrix}}  
}
\examples{
dff <- generateData()
Tmatrix <- createIPMTmatrix(minSize = min(dff$size, na.rm = TRUE),
maxSize = max(dff$size, na.rm = TRUE), growObj = makeGrowthObj(dff), 
survObj = makeSurvObj(dff), correction="constant")
Fmatrix <- createIPMFmatrix(minSize = min(dff$size, na.rm = TRUE),
maxSize = max(dff$size, na.rm = TRUE), fecObj = makeFecObj(dff), correction="constant")
R0Calc(Tmatrix, Fmatrix)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

\name{R0.calc}
\alias{R0.calc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculating R0 from an IPM
}
\description{
Function to estimate lifetime reproductive success from a full IPM
(including growth, survival and fertility transitions). 
}
\usage{
R0.calc(Tmatrix, Fmatrix)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Tmatrix}{
    a matrix (not necessarily of class IPM.matrix)
}
  \item{Fmatrix}{
    a matrix (not necessarily of class IPM.matrix)
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
interpretation. Sinauer. p126
}
\author{
C. Jessica E. Metcalf, Sean McMahon, Roberto Salguero-Gómez & Eelke Jongejans
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{create.IPM.Tmatrix}},\code{\link{create.IPM.Fmatrix}}  
}
\examples{
dff<-generateData()
Tmatrix <-create.IPM.Tmatrix(minsize=min(dff$size, na.rm = TRUE),
maxsize=max(dff$size, na.rm = TRUE),growObj=makeGrowthObjGeneral(dff), survObj=makeSurvObjGeneral(dff))
Fmatrix <-create.IPM.Fmatrix(minsize=min(dff$size, na.rm = TRUE),
maxsize=max(dff$size, na.rm = TRUE),fecObj=makeFecObjGeneral(dff))
R0.calc(Tmatrix, Fmatrix)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

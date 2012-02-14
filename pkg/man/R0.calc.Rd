\name{R0.calc}
\alias{R0.calc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculating R0 from an IPM
}
\description{
Methods to estimate lifetime reproductive success from a full IPM
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
dff<-generateData()
Tmatrix <-create.IPM.Tmatrix(minsize=min(dff$size),
maxsize=max(dff$size),growObj=makeGrowthObj(dff), survObj=makeSurvObj(dff))
Fmatrix <-create.IPM.Fmatrix(minsize=min(dff$size),
maxsize=max(dff$size),growObj=makeGrowthObj(dff), survObj=makeSurvObj(dff), fecObj=makeFecObj(dff))
R0.calc(Tmatrix, Fmatrix)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

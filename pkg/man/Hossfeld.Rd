\name{Hossfeld}
\alias{Hossfeld}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Hossfeld function defining growth 
}
\description{
Functional form describing growth according to a Hossfeld function
}
\usage{
Hossfeld(size, par)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{size}{
    vector of sizes
}
  \item{par}{
    vector of length 3
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Zuidema et al. 2010. Integral Projection Models for trees: a new parametrization method and a validation of model output. 
Journal of Ecology, 98, 345-355
}
\author{
C. Jessica E. Metcalf, Sean McMahon, Roberto Salguero-Gómez & Eelke Jongejans
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
Hossfeld(1:10,rep(1,3))

plot(1:10,Hossfeld(1:10,rep(1,3)), type="l",
ylab="predicted increment as a function of size")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

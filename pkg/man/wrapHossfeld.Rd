\name{wrapHossfeld}
\alias{wrapHossfeld}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fitting Hossfeld growth function
}
\description{
Function to define deviance of a Hossfeld function for use with optim to
generate a Hossfeld growth object
}
\usage{
wrapHossfeld(par, dataf)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{par}{
    vector of length three
}
  \item{dataf}{
    dataframe with columns size and incr
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
C. Jessica E. Metcalf, Sean McMahon, Roberto Salguero-Gomez & Eelke Jongejans
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
 \code{\link{help}} {Hossfeld}
 }
\examples{

#simulate data and create a column for growth increment
dff<-generateData()
dff$incr <- dff$sizeNext-dff$size

#current sum of squares
wrapHossfeld(c(1,1,1), dff)

#use optim to get best parameters values [not run]
#tmp <- optim(c(1,1,1),wrapHossfeld,dataf=dff,method="Nelder-Mead")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

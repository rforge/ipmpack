\name{growthObj.logincr-class}
\Rdversion{1.1}
\docType{class}
\alias{growthObj.logincr-class}
\alias{growth,numeric,numeric,numeric,growthObj.logincr-method}
\alias{growthCum,numeric,numeric,numeric,growthObj.logincr-method}

\title{Class \code{"growthObj.logincr"}}
\description{
%%  ~~ A concise (1-5 lines) description of what the class is. ~~
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("growthObj.logincr", ...)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{fit}:}{Object of class \code{"lm"} ~~ }
  }
}
\section{Methods}{
  \describe{
    \item{growth}{\code{signature(size = "numeric", sizeNext = "numeric", cov = "numeric", growthObj = "growthObj.logincr")}: ... }
    \item{growthCum}{\code{signature(size = "numeric", sizeNext = "numeric", cov = "numeric", growthObj = "growthObj.logincr")}: ... }
	 }
}
\references{
%%  ~~put references to the literature/web site here~~
}
\author{
C. Jessica E. Metcalf, Sean McMahon, Roberto Salguero-Gomez & Eelke Jongejans.
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%%  ~~objects to See Also as \code{\link{~~fun~~}}, ~~~
%%  ~~or \code{\linkS4class{CLASSNAME}} for links to other classes ~~~
}
\examples{
showClass("growthObj.logincr")
}
\keyword{classes}

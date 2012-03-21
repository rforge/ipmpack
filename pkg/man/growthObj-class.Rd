\name{growthObj-class}
\Rdversion{1.1}
\docType{class}
\alias{growthObj-class}
%\alias{growth,numeric,numeric,numeric,growthObj-method}
%\alias{growthIntegrate,numeric,numeric,numeric,numeric,growthObj-method}

\title{Class "growthObj"}
\description{
A class object description
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("growthObj", ...)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{fit}:}{Object of class \code{"lm"} ~~ }
  }
}
\section{Methods}{
  \describe{
    \item{growth}{\code{signature(size = "numeric", sizeNext = "numeric", cov = "numeric", growthObj = "growthObj")}: ... }
    \item{growthCum}{\code{signature(size = "numeric", sizeNext = "numeric",  cov = "numeric", growthObj = "growthObj")}: ... }
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
showClass("growthObj")
}
\keyword{classes}

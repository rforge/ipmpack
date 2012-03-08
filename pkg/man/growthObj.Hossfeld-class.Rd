\name{growthObj.Hossfeld-class}
\Rdversion{1.1}
\docType{class}
\alias{growthObj.Hossfeld-class}
\alias{growth,numeric,numeric,numeric,growthObj.Hossfeld-method}
\alias{growthIntegrate,numeric,numeric,numeric,numeric,growthObj.Hossfeld-method}

\title{Class "growthObj.Hossfeld"}
\description{
%%  ~~ A concise (1-5 lines) description of what the class is. ~~
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("growthObj.Hossfeld", ...)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{logLik}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{paras}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{sd}:}{Object of class \code{"numeric"} ~~ }
  }
}
\section{Methods}{
  \describe{
    \item{growth}{\code{signature(size = "numeric", sizeNext = "numeric", cov = "numeric", growthObj = "growthObj.Hossfeld")}: ... }
    \item{growthIntegrate}{\code{signature(size = "numeric", sizeNext = "numeric", h = "numeric", cov = "numeric", growthObj = "growthObj.Hossfeld")}: ... }
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
showClass("growthObj.Hossfeld")
}
\keyword{classes}

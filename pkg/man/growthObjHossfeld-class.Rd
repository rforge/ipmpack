\name{growthObjHossfeld-class}
\Rdversion{1.1}
\docType{class}
\alias{growthObjHossfeld-class}
%\alias{growth,numeric,numeric,numeric,growthObjHossfeld-method}
%\alias{growthIntegrate,numeric,numeric,numeric,numeric,growthObjHossfeld-method}

\title{Class "growthObjHossfeld"}
\description{
A class object description
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("growthObjHossfeld", ...)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{logLik}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{paras}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{sd}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{hessian}:}{Object of class \code{"matrix"} ~~ }
  }
}
\section{Methods}{
  \describe{
    \item{growth}{\code{signature(size = "numeric", sizeNext = "numeric", cov = "numeric", growthObj = "growthObjHossfeld")}: ... }
    \item{growthIntegrate}{\code{signature(size = "numeric", sizeNext = "numeric", h = "numeric", cov = "numeric", growthObj = "growthObjHossfeld")}: ... }
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
showClass("growthObjHossfeld")
}
\keyword{classes}

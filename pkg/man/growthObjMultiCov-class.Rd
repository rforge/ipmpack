\name{growthObjMultiCov-class}
\Rdversion{1.1}
\docType{class}
\alias{growthObjMultiCov-class}
%\alias{growth,numeric,numeric,data.frame,growthObjMultiCov-method}

\title{Class \code{"growthObjMultiCov"}}
\description{
A class object description
%%  ~~ A concise (1-5 lines) description of what the class is. ~~
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("growthObjMultiCov", ...)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{fit}:}{Object of class \code{"lm"} ~~ }
  }
}
\section{Methods}{
  \describe{
    \item{growth}{\code{signature(size = "numeric", sizeNext = "numeric", cov = "data.frame", growthObj = "growthObjMultiCov")}: ... }
	 }
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
showClass("growthObjMultiCov")
}
\keyword{classes}

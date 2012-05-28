\name{growthObjMultiCovIncr-class}
\Rdversion{1.1}
\docType{class}
\alias{growthObjMultiCovIncr-class}
%\alias{growth,numeric,numeric,data.frame,growthObjMultiCovIncr-method}

\title{Class \code{"growthObjMultiCovIncr"}}
\description{
A class object description
%%  ~~ A concise (1-5 lines) description of what the class is. ~~
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("growthObjMultiCovIncr", ...)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{fit}:}{Object of class \code{"lm"} ~~ }
     \item{\code{sd}:}{Object of class \code{"numeric"} ~~ }
  }
}
\section{Methods}{
  \describe{
    \item{growth}{\code{signature(size = "numeric", sizeNext = "numeric", cov = "data.frame", growthObj = "growthObjMultiCovIncr")}: ... }
	 }
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
showClass("growthObjMultiCovIncr")
}
\keyword{classes}

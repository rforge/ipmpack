\name{survObjMultiCov-class}
\Rdversion{1.1}
\docType{class}
\alias{survObjMultiCov-class}
\alias{surv,numeric,data.frame,survObjMultiCov-method}

\title{Class \code{"survObjMultiCov"}}
\description{
A class object description
%%  ~~ A concise (1-5 lines) description of what the class is. ~~
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("survObjMultiCov", ...)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{fit}:}{Object of class \code{"glm"} ~~ }
  }
}
\section{Methods}{
  \describe{
    \item{surv}{\code{signature(size = "numeric", cov = "data.frame", survObj = "survObjMultiCov")}: ... }
	 }
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
showClass("survObjMultiCov")
}
\keyword{classes}

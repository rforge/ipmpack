\name{IPMmatrix-class}
\Rdversion{1.1}
\docType{class}
%\alias{IPMmatrix-class}
%\alias{lifeExpect,IPMmatrix,envMatrix-method}
%\alias{meanLifeExpect,IPMmatrix-method}
%\alias{passageTime,numeric,IPMmatrix-method}
%\alias{stochpassageTime,numeric,IPMmatrix,envMatrix-method}
%\alias{survivorship,IPMmatrix-method}
%\alias{varLifeExpect,IPMmatrix-method}
%\alias{varpassageTime,numeric,IPMmatrix-method}

\title{Class "IPMmatrix"}
\description{
A class object description
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("IPMmatrix", data, nrow, ncol, byrow, dimnames, ...)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{.Data}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{nEnvClass}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{nBigMatrix}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{meshpoints}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{env.index}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{nDiscrete}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{names.discrete}:}{Object of class \code{"character"} ~~ }
  }
}
\section{Extends}{
Class \code{"\linkS4class{matrix}"}, from data part.
Class \code{"\linkS4class{array}"}, by class "matrix", distance 2.
Class \code{"\linkS4class{structure}"}, by class "matrix", distance 3.
Class \code{"\linkS4class{vector}"}, by class "matrix", distance 4, with explicit coerce.
}
\section{Methods}{
  \describe{
    \item{lifeExpect}{\code{signature(IPMmatrix = "IPMmatrix", envMatrix = "envMatrix")}: ... }
    \item{meanLifeExpect}{\code{signature(IPMmatrix = "IPMmatrix")}: ... }
    \item{passageTime}{\code{signature(chosenSize = "numeric", IPMmatrix = "IPMmatrix")}: ... }
    \item{stochpassageTime}{\code{signature(chosenSize = "numeric", IPMmatrix = "IPMmatrix", envMatrix = "envMatrix")}: ... }
    \item{survivorship}{\code{signature(IPMmatrix = "IPMmatrix")}: ... }
    \item{varLifeExpect}{\code{signature(IPMmatrix = "IPMmatrix")}: ... }
    \item{varpassageTime}{\code{signature(chosenSize = "numeric", IPMmatrix = "IPMmatrix")}: ... }
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
showClass("IPMmatrix")
}
\keyword{classes}

\name{envMatrix-class}
\Rdversion{1.1}
\docType{class}
\alias{envMatrix-class}
\alias{lifeExpect,IPMmatrix,envMatrix-method}
\alias{stochPassageTime,numeric,IPMmatrix,envMatrix-method}

\title{Class "envMatrix"}
\description{
%%  ~~ A concise (1-5 lines) description of what the class is. ~~
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("envMatrix", data, nrow, ncol, byrow, dimnames, ...)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{.Data}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{nEnvClass}:}{Object of class \code{"numeric"} ~~ }
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
    \item{stochPassageTime}{\code{signature(chosenSize = "numeric", IPMmatrix = "IPMmatrix", envMatrix = "envMatrix")}: ... }
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
showClass("envMatrix")
}
\keyword{classes}

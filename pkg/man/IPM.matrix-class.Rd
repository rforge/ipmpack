\name{IPM.matrix-class}
\Rdversion{1.1}
\docType{class}
\alias{IPM.matrix-class}
\alias{LifeExpect,IPM.matrix,env.matrix-method}
\alias{MeanLifeExpect,IPM.matrix-method}
\alias{PassageTime,numeric,IPM.matrix-method}
\alias{StochPassageTime,numeric,IPM.matrix,env.matrix-method}
\alias{Survivorship,IPM.matrix-method}
\alias{VarLifeExpect,IPM.matrix-method}
\alias{varPassageTime,numeric,IPM.matrix-method}

\title{Class "IPM.matrix"}
\description{
%%  ~~ A concise (1-5 lines) description of what the class is. ~~
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("IPM.matrix", data, nrow, ncol, byrow, dimnames, ...)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{.Data}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{n.env.class}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{n.big.matrix}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{meshpoints}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{env.index}:}{Object of class \code{"numeric"} ~~ }
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
    \item{LifeExpect}{\code{signature(IPM.matrix = "IPM.matrix", env.matrix = "env.matrix")}: ... }
    \item{MeanLifeExpect}{\code{signature(IPM.matrix = "IPM.matrix")}: ... }
    \item{PassageTime}{\code{signature(chosen.size = "numeric", IPM.matrix = "IPM.matrix")}: ... }
    \item{StochPassageTime}{\code{signature(chosen.size = "numeric", IPM.matrix = "IPM.matrix", env.matrix = "env.matrix")}: ... }
    \item{Survivorship}{\code{signature(IPM.matrix = "IPM.matrix")}: ... }
    \item{VarLifeExpect}{\code{signature(IPM.matrix = "IPM.matrix")}: ... }
    \item{varPassageTime}{\code{signature(chosen.size = "numeric", IPM.matrix = "IPM.matrix")}: ... }
	 }
}
\references{
%%  ~~put references to the literature/web site here~~
}
\author{
%%  ~~who you are~~
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
showClass("IPM.matrix")
}
\keyword{classes}

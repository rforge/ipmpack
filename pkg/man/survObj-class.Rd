\name{survObj-class}
\Rdversion{1.1}
\docType{class}
\alias{survObj-class}
\alias{surv,numeric,numeric,survObj-method}

\title{Class "survObj"}
\description{
A class object description
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("survObj", ...)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{fit}:}{Object of class \code{"glm"} ~~ }
  }
}
\section{Methods}{
  \describe{
    \item{surv}{\code{signature(size = "numeric", cov = "numeric", survObj = "survObj")}: ... }
	 }
}
\references{
%%  ~~put references to the literature/web site here~~
}
\author{
C. Jessica E. Metcalf, Sean M. McMahon, Roberto Salguero-Gomez & Eelke Jongejans
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
showClass("survObj")
}
\keyword{classes}

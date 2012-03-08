\name{logit}
\alias{logit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Implements a logistic transform.
}
\description{
Provided a vector or numeric changes it into a vector on 0-1 using the logit transform.
}
\usage{
logit(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{x}{
vector of numbers for which the transform is required.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
C. Jessica E. Metcalf, Sean McMahon, Roberto Salguero-Gomez & Eelke Jongejans.
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
x <- rnorm(100)
plot(sort(x), logit(sort(x)), type = "l", xlab = "State x", ylab = 
"Transformed state x", col = "dark gray")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

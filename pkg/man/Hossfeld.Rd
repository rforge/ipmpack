\name{Hossfeld}
\alias{Hossfeld}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Hossfeld function defining growth 
}
\description{
Functional form describing growth according to a Hossfeld function
}
\usage{
Hossfeld(size, par)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{size}{
    vector of sizes
}
  \item{par}{
    vector of length 3
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
Zuidema, Jongejans, Chien, During & Schieving. 2010. Integral projection models for trees: a new parameterization method and a validation of model output. 
Journal of Ecology 98 p345-355

Rivas, Gonzalez, Gonzalez & von Gadow. 2004. Compatible height and site index models for five pine species in El Salto, Durango (Mexico). Forest Ecology and Management 201 p145-160

}
\author{
C. Jessica E. Metcalf, Sean McMahon, Roberto Salguero-Gomez & Eelke Jongejans
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
Hossfeld(1:10,rep(1,3))

plot(1:10,Hossfeld(1:10,rep(1,3)), type="l",ylab="Predicted increment t to t+1",xlab="Size in t")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

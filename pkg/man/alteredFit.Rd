\name{alteredFit}
\alias{alteredFit}
	\title{
		Allows the insertion of new coefficients and new variances into the \code{fit} slot of growth objects.
		}
\description{
	Function replaces a growth object fit with new coefficients and new residuals that are read by \code{predict} when making the IPM.
	}
\usage{
	alteredFit(dummyFit, 
		newCoef = dummyFit$coefficients, 
		desiredSd = 1) 
}
%- maybe also 'usage' for other objects documented here.
		\arguments{
			\item{dummyFit}{
				object of class \code{lm} that has coefficients and/or residuals to be replaced. 
			}
			\item{newCoef}{
				vector, the value of the new coefficients. Defaults to old coefficients.
			}
			\item{desiredSd}{
				numeric, the standard deviation wanted for the output object. Defaults to 1.
			}
}
\details{
	Note that the longer the vector of dummyFit$residuals, the more accurate the desired sd will be.  To make a highly accurate sd, merely use a grObjOld with a high sample size, e.g., by increasing nsamp in generateData.
}
\value{
	\code{dummyFit}, an object of class \code{lm} with new residuals and, potentially, new coefficients.
}
\references{
}
\author{
	C. Jessica E. Metcalf, Sean McMahon, Roberto Salguero-Gomez & Eelke Jongejans.
}
\note{
	%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~
\examples{
# Data with size and sizeNext
	dff <- generateData(nsamp=100000)
	gr1 <- makeGrowthObj(dff)
	print(sd(gr1@fit$residuals))
	print(coef(gr1@fit))
	newFit <- alteredFit(dummyFit = gr1@fit, newCoef = c(0.24, 0.55), desiredSd = 4)
	print(sd(newFit$residuals))
	print(coef(newFit))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
		% R documentation directory.
		\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
		
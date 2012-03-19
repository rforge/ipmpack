\name{alteredFit}
\alias{alteredFit}

		\title{
		Allows the insertion of new coefficients and new variances into the \code{fit} slot of growth objects.
		}
\description{
	Function replaces a growth object fit with new coefficients and new residuals that are read by \code{predict} when making the IPM.
	}
\usage{
	alteredFit <- function(grObjOld, 
		newCoef = groObjOld@fit$coefficients, 
		desiredSd = 1) 
}
%- maybe also 'usage' for other objects documented here.
		\arguments{
			\item{grObjOld}{
				object of class \code{grObj} that has a \code{fit} slot of class \code{lm} that will have coefficients and/or residuals replaced. 
			}
			\item{newCoef}{
				vector, the value of the new coefficients. Defaults to old coefficients.
			}
			\item{desiredSd}{
				numeric, the standard deviation wanted for the output object. Defaults to 1.
			}
}
\details{
	Note that the longer the vector of dummy residuals, the more accurate the desired sd will be.  To make a highly accurate sd, merely use a grObjOld with a high sample size. This cannot yet be adjusted in generateData().
}
\value{
	\code{grObjNew}, an object of class \code{grObj} with a \code{fit} slot of class \code{lm} with identical values as those of the \code{groObjOld}, save for new residuals and, potentially, new coefficients.
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
	dff <- generateData()
	gr1 <- makeGrowthObj(dff)
	print(sigma.hat(gr1@fit))
	print(coef(gr1@fit))
	gr1New <- alteredFit(grObjOld = gr1, newCoef = c(0.24, 0.55), desiredSd = 4)
	print(sigma.hat(gr1New@fit))
	print(coef(gr1New@fit))
	
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
		% R documentation directory.
		\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
		
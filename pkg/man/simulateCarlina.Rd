\name{simulateCarlina}
\alias{simulateCarlina}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Generates random data in the form used by IPMpack based on the population dynamics of Carlina vulgaris
}
\description{
  Simulates growth, survival and fecundity and density dependent seedling establishment to create a dataframe of the form required by the functions and methods used in IPMpack. Demographic stage data is only continuous.
}
\usage{
simulateCarlina(nSamp=2000,nYrs=1000,nSampleYrs=15,
		m0=-1.37,ms=0.59,
		b0=-12.05,bs=3.64,
		A=-1,B=2,
		ag=1.13,bg=0.74,sig=sqrt(0.095),
		mean.kids=3.0,sd.kids=0.52,
		meanYear=c(0,0,0),
		matVarYear=matrix(c(1.34,0.1,0,0.1,0.04,0,0,0,0.01),3,3))
}
\arguments{
	\item{nSamp}{
	number of samples desired in the base population, defaults to 2000
	}
	\item{nYrs}{
	number of years in the simulation, defaults to 1000
	}
	\item{nSampleYrs}{
	number of years sampled, defaults to 15
	}
	\item{m0}{
intercept survival
	}
	\item{ms}{
	slope survival
	}
\item{b0}{
	intercept flowering
	}
\item{bs}{
	slope flowering
	}
\item{A}{
	intercept reproductive allometry seed production
	}
\item{B}{
	slope reproductive allometry seed production
	}
\item{ag}{
	intercept growth
	}
\item{bg}{
	slope growth
	}
\item{sig}{
	variance growth
	}
\item{mean.kids}{
mean kid size
	}
\item{sd.kids}{
variance kid size
	}
\item{meanYear}{
mean year effects
	}
\item{matVarYear}{
covariance covariance in eyar effects
	}



}
%- maybe also 'usage' for other objects documented here.
\value{
 
A list including: 
  dataf:
  A dataframe with headings:
  	- "size": continuous variable, indicating current size.
  	- "sizeNext" continuous variable, indicating size in the next time step.
  	- "surv": boolean, indicating whether individual survived or not to the next time step.
  	- "covariate": discrete covariate.
  	- "covariateNext": discrete covariate in the next
  time step.
  	- "fec": continuous variable, indicating fecundity.
	- nSeedlings: number seedlings corresponding to that year
	- m.year: intercept of mortality for that year
	- cg.year: intercept of growth for that year
	- b.year: intercept of offspring size for that year
	- offspringNext: where the row corresponds to offspring, this takes the value offspringNexxt
	- year: year of the sample
	
 listpar: 
 	- a list of all the other parameters	
		
}
\author{
C. Jessica E. Metcalf, Sean M. McMahon, Roberto Salguero-Gomez & Eelke Jongejans.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{generateDataStoch}},\code{\link{generateDataDiscrete}}
}
\examples{
dff <- simulateCarlina(nSamp=1000)
head(dff$dataf)

}

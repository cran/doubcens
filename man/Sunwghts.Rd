
\name{Sunwghts}
\alias{Sunwghts}
\title{Estimate Survivor Function using Doubly Interval-Censored Left-Truncated Failure Time Data}
\usage{Sunwghts(Ei, Ri, Li, Ui, Bi1, Bi2)}
\arguments{
\item{Ei}{The left endpoint of the censoring interval for the initial event time}
\item{Ri}{The right endpoint of the censoring interval for the initial event time}
\item{Li}{The left endpoint of the censoring interval for the terminating event time}
\item{Ui}{The right endpoint of the censoring interval for the terminating event time}
\item{Bi1}{The left endpoint of the truncation interval for the terminating event time}
\item{Bi2}{The right endpoint of the truncation interval for the terminating event time}
}
\description{Estimates the discrete Survivor Function from doubly interval-censored left-truncated failure time data using the algorithm of Sun.}
\details{
Set Ei = Ri if the initial event is observed. Set Li = Ui if the 
terminating event is observed. Set Ei = -INF if the initial event is left-censored. Set Ri = INF if the terminating event is right-censored.
}
\value{
Sunwghts returns a list containing the following components
\item{uj}{A vector of mass points for survival lengths}
\item{fnew}{A vector of estimated probabilities for uj}
\item{counter}{Number of iterations required for convergence}
}
\references{Sun, J. (1995). Empirical Estimation of a Distribution Function with Truncated and Doubly Interval-Censored Data and Its Applications to AIDS Studies. Biometrics 51 (3): 1096-1104.}

\examples{
test <- Sunwghts(c(1,2,1), c(4,4,4), c(5,8,9), c(5,9,10), 
		c(4.4,4.5,8), c(Inf, Inf, Inf))

}



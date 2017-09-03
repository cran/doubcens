
\name{DGLwghts}
\alias{DGLwghts}
\title{Estimate Survivor Function using Doubly Interval-Censored Failure Time Data}
\usage{DGLwghts(X_L, X_R, Z_L, Z_R)}
\arguments{
\item{X_L}{The left endpoint of the censoring interval for the initial event time}
\item{X_R}{The right endpoint of the censoring interval for the initial event time}
\item{Z_L}{The left endpoint of the censoring interval for the terminating event time}
\item{Z_R}{The right endpoint of the censoring interval for the terminating event time}
}
\description{Estimates the discrete Survivor Function from doubly interval-censored failure time data using the algorithm of De Gruttola and Lagakos.}
\details{
Set X_L = X_R if the initial event is observed. Set Z_L = Z_R if the 
terminating event is observed. Set X_L = -INF if the initial event is left-censored. Set Z_R = INF if the terminating event is right-censored.
}
\value{
DGLwghts returns a list containing the following components
\item{x_val}{A vector of mass points for initial event}
\item{w_new}{A vector of estimated probabilities for x_val}
\item{t_val}{A vector of mass points for terminating event}
\item{f_new}{A vector of estimated probabilities for t_val}
\item{counter}{Number of iterations required for convergence}
}


\references{De Gruttola, V. and Lagakos, S. (1989). Analysis of Doubly-Censored Survival Data, with Applications to AIDS. Biometrics 45 (1): 1-11.}

\examples{
test <- DGLwghts(c(1,1,1), c(1,2,1), c(1,2,3), c(Inf, Inf, Inf))

}

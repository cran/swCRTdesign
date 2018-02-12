\name{swSim}
\alias{swSim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulating individual-level data for specified study design of Stepped Wedge Cluster Randomized Trial (SW CRT)
}
\description{
\code{swSim} returns individual-level data set of a SW CRT study design for the specified number of clusters per wave, fractional treatment effect at each time after crossing over from control, time (standard SW CRT time computed; extending trial beyond standard time needs to be specified by the user), family (and link function), number of individuals per cluster per wave, mean in control arm, mean in treatment arm, time effect, pooled standard deviation for both arms, standard deviation of random intercepts, standard deviation of random slopes, correlation between random intercepts and random slopes, time point label, specified random seed, and option to simulate data with time on treatment lag. 
}
\usage{
swSim(design, family, 
n, mu0, mu1, time.effect, sigma, tau, eta, rho = 0, time.lab = NULL, seed = NULL,
retTimeOnTx = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{design}{
list: A stepped wedge design object, typically from swDsn, that includes at least the following components: 
##             swDsn, clusters, n.clusters, total.time
}
  \item{family}{
character: Used in typical way. However, only Gaussian, Binomial, and Poisson families accepted. Also, only identity, logit, and log links accepted, but any combination of the three families and three link functions (i.e., do not have to specify canonical link). ***NOTE: It is the users responsibility to make sure specified parameters (mu0, mu1, time.effect, sigma, tau, eta, rho) are ALL on SAME scale as specified link function.***
}
  \item{n}{
integer (scalar, vector, or matrix): Number of observations: (scalar) for all clusters and all time points; (vector) for each cluster at all time points; and (matrix) for each cluster at each time point.
}
  \item{mu0}{
numeric (scalar): Mean in control arm/group.
}
  \item{mu1}{
numeric (scalar): Mean in treatment arm/group.
}
  \item{time.effect}{
integer (scalar or vector): Time effect at each time point (added to mean at each time).
}
  \item{sigma}{
numeric (scalar): Pooled treatment and control arm standard deviation. Ignored if family != Gaussian
}
  \item{tau}{
numeric (scalar): Standard deviation of random intercepts.
}
  \item{eta}{
numeric (scalar): Standard deviation of random slopes.
}
  \item{rho}{
numeric (scalar): Correlation between random intercepts and random slopes. The default value is 0 (i.e., uncorrelated)
}
  \item{time.lab}{
character (vector): Labeling for time points when output is display; choice of labeling does not affect results.
}
  \item{seed}{
integer (scalar): An initial random seed.
}
  \item{retTimeOnTx}{
logical: If \code{TRUE}, simulates data with time on treatment lag. If \code{FALSE}, simulates data ignoring time on treatment lag. The default is \code{FALSE}.
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
numeric (data frame): Returns the following (individual-level) variables corresponding to the specified SW CRT design:
%% *****SPECIFY WHAT THESE ARE!!!!!**********
	\item{$response.var}{ 
		numeric (vector):  Response variable based on specified SW CRT design of interest (including family and link function) for each observation in the data frame/set.
	}
	\item{$tx.var}{  
		numeric (vector): Predictor of interest. (Fractional) treatment effect corresponding to 0=control, 1=treatment, and value between 0 and 1 corresponding to treatment arm with fractional treatment effect (for each observation in the data frame/set).
	}
	\item{$timeOnTx.var}{  
		numeric (vector): Predictor of interest when interested in time on treatment lag effect. (Fractional) treatment effect corresponding to 0=control, 1=treatment, and value between 0 and 1 corresponding to treatment arm with fractional treatment effect (for each observation in the data frame/set).
	}
	\item{$time.var}{  
		numeric (vector): Time point id for each observation in the data frame/set.
	}
	\item{$cluster.var}{  
		numeric (vector): Grouping variable. Cluster id for each observation in the data frame/set.
	}
}
\references{
Hussey MA, Hughes JP. Design and analysis of stepped wedge cluster randomized trials. Contemporary Clinical Trials 2007;28:182-191.
}
\author{
James P Hughes and Navneet R Hakhu
}
%\note{
%%  ~~further notes~~
%}
%
%% ~Make other sections like Warning with \section{Warning }{....} ~
%
%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
library(swCRTdesign)
# Example 1 [ n = scalar; can be vector (for different n for each cluster, 
# n=rep(120,22)) or matrix (different n for each cluster at each time point, 
# n=matrix(120,22,5)) ]

# generate SW data (fractional treatment effect)
design <- swDsn(clusters=c(6,6,6,4), tx.effect=c(0.8,0.9,1.0), extra.time=0, all.ctl.time0=TRUE)
swGenData.nScalar <- swSim( design,
family=binomial(link="logit"), n=120, 
mu0=log(0.1/0.9), mu1=log(0.9) + log(0.1/0.9),
time.effect=c(0,.2,.3,.4,.5)*0, sigma=1, tau=0.2, eta=0, 
rho=0, time.lab=seq(0,12,3), seed=5, retTimeOnTx=FALSE)

# summarize SW data by wave
swSummary(response.var, tx.var, time.var, cluster.var, swGenData.nScalar, 
type="mean", digits=3)$response.wave

swSummary(response.var, tx.var, time.var, cluster.var, swGenData.nScalar, 
type="mean", digits=3)$swDsn

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ design }
\keyword{ stepped wedge }
\keyword{ cluster randomized trial }% __ONLY ONE__ keyword per line

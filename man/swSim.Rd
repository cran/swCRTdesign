\name{swSim}
\alias{swSim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulating individual-level data for specified study design of Stepped Wedge Cluster Randomized Trial (SW CRT)
}
\description{
\code{swSim} returns individual-level data set of a SW CRT study design for the specified number of clusters per wave, treatment effect (possibly varying according to the exposure time - time after crossing over from control), time effect, family (and link function), number of individuals per cluster per time period, mean in control arm, mean in treatment arm(s), standard deviation (if applicable), standard deviation of random intercepts, standard deviation of random treatment effects, correlation between random intercepts and random treatment effects, standard deviation of random time effects, and, for closed cohorts, standard deviation of individual random effects. Alternatively, for a Gaussian family, standard deviations of random effects may be specified using ICC, CAC and, for closed cohorts, IAC; see \code{\link{swPwr}} details. An option to add time point labels is also included.  
}
\usage{
swSim(design, family, log.gaussian = FALSE, n, mu0, mu1, time.effect, 
sigma, tau=0, eta=0, rho=0, gamma=0, zeta=0, icc=0, cac=1, iac=0, 
time.lab = NULL, retTimeOnTx=TRUE, silent = FALSE, nocheck=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{design}{
list: A stepped wedge design object, typically from \code{swDsn}, that includes at least the following components: swDsn, swDsn.unique.clusters, clusters, n.clusters, total.time, nTxLev, TxLev. Fractional treatment effects specified in \code{swDsn} are used in the simulation.
}
  \item{family}{
character: Used in typical way. However, only Gaussian, Binomial, and Poisson families accepted. Also, only identity, logit, and log links accepted.  Logit link is only available for Binomial family, and log link is only available for Binomial and Poisson.  Currently, 'Binomial' implies Bernoulli. Default links are identity for Gaussian, logit for Binomial, log for Posson. ***NOTE: It is the users responsibility to make sure specified parameters (mu0, mu1, time.effect, tau, eta, rho, gamma) are ALL on the SAME scale as the specified link function; see example.***
}
  \item{log.gaussian}{
character: When \code{TRUE} with a Gaussian family, simulates data whose log follows a Gaussian distribution; all parameters (mu0, mu1, time.effect, variance parameters) refer to the log scale.  Default is \code{FALSE}.
}
  \item{n}{
integer (scalar, vector, or matrix): Number of observations: (scalar) for all clusters and all time points; (vector) for each cluster at all time points; and (matrix) for each cluster at each time point, where rows correspond to clusters and columns correspond to time.  \code{n} can also be used to specify a design with transition periods were no data is collected; see \code{swPwr}.
}
  \item{mu0}{
numeric (scalar): Mean outcome in the control group on the appropriate scale.
}
  \item{mu1}{
numeric (scalar or vector): Mean outcome in the treatment group on the appropriate scale. If design$nTxLev > 1 (number of treatment levels greater than 1) then \code{mu1} should have length equal to \code{design$nTxLev}; if \code{design$nTxLev = 1} then \code{mu1} should be either a scalar (constant treatment effect) or a vector with length equal to the maximum number of exposure times (time-varying treatment effect). Do not specify time-varying or multiple treatment effects if the design includes fractional treatment effects.
}
  \item{time.effect}{
integer (scalar or vector): Time effect at each time point on the appropriate scale (added to mean at each time).
}
  \item{sigma}{
numeric (scalar): Pooled treatment and control arm standard deviation on the appropriate scale. Ignored if family != Gaussian.
}
  \item{tau}{
numeric (scalar): Standard deviation of random intercepts on the appropriate scale. Default is 0.
}
  \item{eta}{
numeric (scalar): Standard deviation of random treatment effects on the appropriate scale. Default is 0.
}
  \item{rho}{
numeric (scalar): Correlation between random intercepts and random treatment effects on the appropriate scale. Default is 0.
}
  \item{gamma}{
numeric (scalar): Standard deviation of random time effects on the appropriate scale. Default is 0.
}
  \item{zeta}{
numeric (scalar): Standard deviation of individual effects for closed cohort sampling on the appropriate scale. Default is 0 which implies cross-sectional sampling. Values greater than 0 imply closed cohort sampling.
}
  \item{icc}{
numeric (scalar): Within-period intra-cluster correlation on the appropriate scale.  Can be used with CAC instead of tau, eta, rho, and gamma when the outcome is Gaussian. Default is 0.
}
  \item{cac}{
numeric (scalar): Cluster auto-correlation on the appropriate scale.  Can be used with ICC instead of tau, eta, rho, and gamma when the outcome is Gaussian. Default is 1.
}
  \item{iac}{
numeric (scalar): Individual auto-correlation for closed cohort sampling on the appropriate scale. Can be used with ICC and CAC instead of tau, eta, rho, gamma, and zeta when the outcome is gaussian. Default is 0 which implies cross-sectional sampling. Values greater than 0 imply closed cohort sampling.
}
  \item{time.lab}{
character (vector): Labeling for time points when output is display; choice of labeling does not affect results.
}
  \item{retTimeOnTx}{
  This argument has been deprecated. TimeOnTx is always returned.
  }
  \item{silent}{
logical: if TRUE, hides reminder about order of entries in \code{n} when \code{n} is not a scalar.  Default value is \code{FALSE}.
}
  \item{nocheck}{
  logical: if TRUE, various checks on the validity of the arguments are not done. For internal use only.
  }

}
\details{
%%  ~~ If necessary, more details than the description above ~~
When simulating from a Gaussian distribution and if eta (and rho) are 0, instead of using tau, eta, rho, gamma and zeta, the icc, cac and iac can be used to define the variability of the random effects.  In this model,
\deqn{ ICC = \frac{\tau^2+\gamma^2}{\tau^2+\gamma^2+\zeta^2+\sigma^2} }
\deqn{ CAC = \frac{\tau^2}{\tau^2+\gamma^2} }
\deqn{ IAC = \frac{\zeta^2}{\zeta^2+\sigma^2}}
Choose one parameterization or the other, do not mix parameterizations.
}
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
		numeric (vector): Predictor of interest when interested in time on treatment lag effect. Total time spent on treatment for each observation in the data frame/set, with 0=control, 1=first time period on treatment, 2=second time period on treatment, etc.
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

Voldal EC, Hakhu NR, Xia, F, Heagerty PJ, Hughes JP. swCRTdesign: An R Package for Stepped Wedge Trial Design and Analysis. Computer Methods and Programs in Biomedicine 2020;196:105514.
}
\author{
James P Hughes, Navneet R Hakhu, and Emily C Voldal
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

# generate cross-sectional SW data 
design <- swDsn(clusters=c(6,6,6,4))
set.seed(5)
swGenData.nScalar <- swSim( design, family=binomial(link="logit"), n=120,
mu0=log(0.1/0.9), mu1=log(0.9) + log(0.1/0.9), time.effect=0, tau=0.2, eta=0,
rho=0, gamma=0, time.lab=seq(0,12,3))

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

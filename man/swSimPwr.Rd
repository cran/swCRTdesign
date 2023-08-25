\name{swSimPwr}
\alias{swSimPwr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulating power for Stepped Wedge Cluster Randomized Trials (SW CRT)
}
\description{
\code{swSimPwr} uses a simulation approach to compute (two-sided) power for the treatment effect(s) from an immediate treatment effect (IT) model or an exposure time indicator (ETI) model (Kenny et al, 2022) for the specified SW CRT design. A cross-sectional or closed cohort sampling scheme can be specified. The outcome of interest can be Gaussian, binomial or Poisson. Individual-level data are generated using \code{swSim} for the specified number of clusters per wave, treatment effect (possibly varying according to the exposure time - time after crossing over from control), time effect, family (and link function), number of individuals per cluster per time period, mean in control arm, mean in treatment arm(s), standard deviation (if applicable), standard deviation of random intercepts, standard deviation of random treatment effects, correlation between random intercepts and random treatment effects, standard deviation of random time effects, and, for closed cohorts, standard deviation of individual random effects. For a Gaussian family, standard deviations of random effects may be specified using ICC, CAC and, for closed cohorts, IAC; see \code{\link{swPwr}} details. These data are then analyzed by either \code{lmer} or \code{glmer} and the proportion of p-values less than or equal to the specified alpha are computed over the nsim simulations.
}
\usage{
swSimPwr(design, family, n, mu0, mu1, time.effect=0, H = NULL, 
sigma, tau=0, eta=0, rho=0, gamma=0,zeta=0, icc=0, cac=1, iac=0, 
alpha = 0.05, nsim = 500, retDATA = FALSE, silent = FALSE, counter=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{design}{
list: A stepped wedge design object, typically from \code{swDsn}, that includes at least the following components: swDsn, swDsn.unique.clusters, clusters, n.clusters, total.time, nTxLev, TxLev. Fractional treatment effects specified in \code{swDsn} are used in the simulation.
}
  \item{family}{
character: Used in typical way. Only Gaussian, Binomial, and Poisson families accepted. Also, only identity, logit, and log links accepted.  Logit link is only available for Binomial family, and log link is only available for Binomial and Poisson.  Currently, 'Binomial' implies Bernoulli. Default links are identity for Gaussian, logit for Binomial, log for Posson. ***NOTE: It is the users responsibility to make sure specified parameters (mu0, mu1, time.effect, tau, eta, rho, gamma, zeta) are ALL on the SAME scale as the specified link function; see example.***
}
  \item{n}{
integer (scalar, vector, or matrix): Number of observations: (scalar) for all clusters and all time points; (vector) for each cluster at all time points; and (matrix) for each cluster at each time point, where rows correspond to clusters and columns correspond to time.  \code{n} can also be used to specify a design with transition periods were no data is collected (a transition period is a period, typically, immediately after treatment is introduced were no data are collected). Simply define \code{n} as a matrix with a sample size of 0 during every transition period. This is equivalent to putting an NA in design$swDsn.
}
  \item{mu0}{
numeric (scalar): Mean outcome in the control group on the appropriate scale.
}
  \item{mu1}{
numeric (scalar or vector): Mean outcome in the treatment group on the appropriate scale. If design$nTxLev > 1 (number of treatment levels greater than 1) then mu1 should have length equal to design$nTxLev; if design$nTxLev = 1 then mu1 should be either a scalar (constant treatment effect) or a vector with length equal to the maximum number of exposure times (time-varying treatment effect). Do not specify time-varying or multiple treatment effects if the design includes fractional treatment effects.
}
  \item{time.effect}{
integer (scalar or vector): Time effect at each time point on the appropriate scale (added to mean at each time). Default is 0.
}
  \item{H}{
numeric (vector): If NULL, then swPwr assumes an immediate, constant treatment effect (IT) model. If not NULL, then an exposure time indicator (ETI) model is assumed and H should be a vector as long as the longest exposure time (typically, number of time periods minus one). H specifies the desired linear combination of exposure time treatment effects. For example, in a stepped wedge trial with 5 time periods and four exposure times, H = rep(.25,4) gives the average treatment effect over the four exposure times; H = c(0,0,.5,.5) ignores the first two periods after the intervention is introduced and averages the remaining periods. Typically, the sum of H is 1.0; if not, it is renormalized to sum to 1.0. H can only be specified when there is a single intervention level (i.e. design$swDsn includes only NA,0,1; see cautions in help for \code{swDsn} about including NA periods when H is used).
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
  \item{alpha}{
numeric (scalar): Desired alpha level for a two-sided test.
}
  \item{nsim}{
numeric (scalar): Number of simulations.
}
  \item{retDATA}{
logical: if FALSE, only power is returned; if TRUE, then estimates, standard errors, p-values and convergence information from all simulations are also returned.
  }
  \item{silent}{
logical: if TRUE, hides most warning messages.  Default value is \code{FALSE}.
  }
  \item{counter}{
logical: if TRUE, show a counter for the simulations. Default value is \code{TRUE}
  }
}

\details{
\code{swSimPwr} uses a simulation approach to compute the two-sided statistical power for testing the hypotheses \eqn{Ho: \mu_1 - \mu_0 = 0} if H is NULL or \eqn{Ho: \sum H*(\mu_1 - \mu_0) = 0} if H is non-NULL. swSim is called to generate \code{nsim} simulations and then \code{lmer} or \code{glmer} is used to analyze each simulated dataset. The simulated power is equal to the proportion of p-values from \code{lmer}/\code{glmer} that are less than or equal to alpha. For small sample sizes it is recommended that the user run a simulation under the null to check the accuracy of the simulated type I error rate. Optionally, a counter is printed as each simulation is completed.
}

\value{
numeric (vector): if \code{retDATA = FALSE}, \code{swSimPwr} returns the power for the two-sided design$nTxLev hypothesis test(s) treatment effect(s) against a null of 0.0.

numeric (list): if \code{retDATA = TRUE}, \code{swSimPwr} returns the following objects of a list.
	\item{power}{
	 Power for the two-sided design$nTxLev hypothesis test(s) of treatment effect(s) against a null of 0.0.
	}
	\item{nsim}{
	 Number of simulation power is based on
	}
	\item{est}{
		matrix: nsim x design$nTxLev matrix of estimated treatment effects. If H is specified then this is the linear combination of the exposure time estimates.
	}
	\item{se}{
		matrix: nsim x design$nTxLev matrix of standard errros of the estimated treatment effects. If H is specified then this is the standard error of the linear combination of the exposure time estimates.
	}
	\item{sdcor}{
	  matrix: matrix with nsim rows and number of columns equal to the number of variance components in the model (expressed on the random effects scale) giving the estimated variance components (expressed as standard deviations) from each fit. Note that eta is expressed as a correlation.
	}
	\item{df}{
		matrix: nsim x design$nTxLev matrix of degrees of freedom for testing the hypotheses on the treatment effect(s). Only returned if lmer is called for estimation; otherwise NULL.
	}
	\item{pval}{
		matrix: nsim x design$nTxLev matrix of p-values for testing the hypotheses on the treatment effect(s). 
	}
	\item{code}{
	  vector: vector of length \code{nsim} giving either 0 (convergence) or the value of rslt@optinfo$conv$lme4$code for each simulation 
	}
	\item{msg}{
	vector: vector of length \code{nsim} giving all messages from rslt@optinfo$conv$lme4$messages for each simulation
	}
}
\references{
Kenny A, Voldal E, Xia F, Heagerty PJ, Hughes JP. Analysis of stepped wedge cluster randomized trials in the presence of a time-varying treatment effect. Statistics in Medicine, in press, 2022.
}
\examples{
logit = function(x){log(x/(1 - x))}
stdy1 <- swDsn(c(6,6,6,6),swBlk=matrix(c(0,1,1,1,1,0,0,1,1,1,0,0,0,1,1,0,0,0,0,1),4,5,byrow=TRUE))
# analytic power by swPwr
swPwr(stdy1, distn="binomial",n=120, mu0=0.45, mu1=0.5, 
      tau=0.1, alpha=0.05, retDATA=FALSE)
# 
# equivalent analysis on logit scale with transformation of variables
swGlmPwr(stdy1,distn="binomial",n=120,fixed.intercept=logit(0.45),
         fixed.treatment.effect = logit(0.5)-logit(0.45), fixed.time.effect=0.0,
         tau=.1/(.475*.525))
#
# Same problem by simulation using three different approaches
#
#swSimPwr(stdy1,family=binomial(link="logit"),n=120,
#         mu0=logit(.45),mu1=logit(.5),time.effect=0.0,
#         tau=.1/(.475*.525),nsim=500)  
#
#swSimPwr(stdy1,family="gaussian",n=120,
#         mu0=.45,mu1=.5,time.effect=0.0,
#         sigma=sqrt(.475*.525),tau=.1,nsim=500) 
#
#swSimPwr(stdy1,family=binomial(link="identity"),n=120,
#         mu0=.45,mu1=0.5,time.effect=0.0,
#         tau=.1,nsim=500)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.

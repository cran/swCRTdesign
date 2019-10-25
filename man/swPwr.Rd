\name{swPwr}
\alias{swPwr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Power of Stepped Wedge Cluster Randomized Trial (SW CRT)
}
\description{
\code{swPwr} returns (two-sided) power of treatment effect (\eqn{\theta}) for the specified SW CRT design using a linear models weighted least squares (WLS) approach. The response/outcome of interest can be binomial or Gaussian distributed and is assumed to come from a model with random intercepts, random treatment effects, and random cluster-specific time effects. Variance components can be specified using either tau, eta, rho, and gamma, or icc and cac (see details).  If a random intercepts only model is used (i.e., eta and gamma are 0 and n is constant over clusters and time), then the power calculation is comparable to the closed-form formula of [Hussey and Hughes, 2007].
}
\usage{
swPwr(design, distn, n, mu0, mu1, sigma, tau, eta, rho, gamma,
icc, cac, alpha=0.05, retDATA=FALSE, silent=FALSE)
%
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{design}{
list: A stepped wedge design object, typically from swDsn, that includes at least the following components:
##             swDsn, swDsn.unique.clusters, clusters, n.clusters, total.time
}
  \item{distn}{
character: Distribution assumed (\code{gaussian} or \code{binomial}).  Currently, 'Binomial' implies Bernoulli.
}
  \item{n}{
integer (scalar, vector, or matrix): Number of observations: (scalar) for all clusters and all time points; (vector) for each cluster at all time points; and (matrix) for each cluster at each time point, where rows correspond to clusters and columns correspond to time.  \code{n} can also be used to specify a design with transition periods (e.g. in the first time period that each sequence receives treatment, no observations are collected from that sequence).  Simply define \code{n} as a matrix with a sample size of 0 during every transition period.
}
  \item{mu0}{
numeric (scalar): Mean outcome in the control group.
}
  \item{mu1}{
numeric (scalar): Mean outcome in the treatment group. Note: Treatment effect is the difference in means \eqn{\theta = \mu_1 - \mu_0}.
}
  \item{sigma}{
numeric (scalar): Standard deviation when assuming Gaussian distribution (\code{distn=gaussian}). For binomial distribution \eqn{\sigma^2} is automatically set to \eqn{\bar{\mu}(1-\bar{\mu})} where \eqn{\bar{\mu} = (\mu_1 + \mu_0)/2}
}
  \item{tau}{
numeric (scalar): Standard deviation of random intercepts.
}
  \item{eta}{
numeric (scalar): Standard deviation of random treatment effects.
}
  \item{rho}{
numeric (scalar): Correlation between random intercepts and random treatment effects.
}
  \item{gamma}{
numeric (scalar): Standard deviation of random time effects.
}
  \item{icc}{
numeric (scalar): Within-period intra-cluster correlation.  Can be used with CAC instead of tau, eta, rho, and gamma; see details.
}
  \item{cac}{
numeric (scalar): Cluster auto-correlation.  Can be used with ICC instead of tau, eta, rho, and gamma; see details.
}
  \item{alpha}{
numeric (scalar): Two-sided statistical significance level.
}
  \item{retDATA}{
logical: if TRUE, all stored (input, intermediate, and output) values of \code{swPwr} are returned. Default value is \code{FALSE}.
}
  \item{silent}{
logical: if TRUE, hides a warning about differences in argument order between version 3.0 and prior versions. When \code{n} is not a scalar, also hides reminder about order of entries in \code{n}.  Default value is \code{FALSE}.
}
}
\details{
%%We assume the outcome variable of interest can be modeled using (possibly correlated) random intercepts and random treatment effects
%%\deqn{ Y_{ij} = \mu_0 + \beta_j + \theta X_{ij} + a_i + r_i X_{ij} + e_{ij} }
%%where
%%\deqn{ a_i \sim N(0, \tau^2) }
%%\deqn{ r_i \sim N(0, \eta^2) }
%%\deqn{ Cov(a_i, r_i) = \rho \tau \eta }
%%\deqn{ e_{ij} \sim N(0, \sigma^2 / K) }
%%for \eqn{K} observations in cluster \eqn{i=1,\ldots,I}  and  at time \eqn{j=1,\ldots,J}.
%%
The two-sided statistical power of treatment effect \eqn{(\theta = \mu_1 - \mu_0)} is
\deqn{ Pwr(\theta) = \Phi( Z - z_{1 - \alpha /2} ) + \Phi( -Z - z_{1 - \alpha /2} ) }
where
\deqn{ Z = \frac{ |\theta| }{ \sqrt{Var(\hat{\theta}_{WLS})} } }
and \eqn{\Phi} is the cumulative distribution function of the standard normal \eqn{N(0,1)} distribution.

When eta (and rho) are 0, instead of using tau, eta, rho, and gamma, the icc and cac can be used to define the variability of the random intercepts and time effects.  In this model,
\deqn{ ICC = \frac{\tau^2+\gamma^2}{\tau^2+\gamma^2+\sigma^2} }
\deqn{ CAC = \frac{\tau^2}{\tau^2+\gamma^2} }
}
\value{
numeric (matrix): \code{swPwr} returns the power of treatment effect (\eqn{\theta}), where the variance of treatment effect is computed by WLS.  %%%% (if \eqn{\eta = 0}, then power computed using closed-form formula for variance is also given) if \code{retDATA = FALSE}.
\item{}{

}
numeric (list): \code{swPwr} returns all specified and computed items as objects of a list if \code{retDATA = TRUE}.
%% *****SPECIFY WHAT THESE ARE!!!!!**********
	\item{design}{
		list: The stepped wedge design object as input.
	}
	\item{distn}{
		character: Distribution assumed (\code{gaussian} or \code{binomial}).
	}
	\item{n}{
		integer (scalar, vector, or matrix): Number of observations: (scalar) for all clusters and all time points; (vector) for each cluster at all time points; and (matrix) for each cluster at each time point, where rows correspond to clusters and columns correspond to time.
	}
	\item{mu0}{
		numeric (scalar): Mean outcome in the control group.
	}
	\item{mu1}{
		numeric (scalar): Mean outcome in intervention group. Note: treatment effect is difference in means \eqn{\theta = \mu_1 - \mu_0}.
	}
	\item{sigma}{
		numeric (scalar): Standard deviation input. For binomial distribution, sigma = NA
	}
	\item{tau}{
		numeric (scalar): Standard deviation of random intercepts.
	}
	\item{eta}{
		numeric (scalar): Standard deviation of random treatment effects.
	}
	\item{rho}{
		numeric (scalar): Correlation between random intercepts and random treatment effects.
	}
	\item{gamma}{
numeric (scalar): Standard deviation of random time effects.
  }
  \item{icc}{
numeric (scalar): Within-period intra-cluster correlation.  Can be used with CAC instead of tau, eta, rho, and gamma; see details.
  }
  \item{cac}{
numeric (scalar): Cluster auto-correlation.  Can be used with ICC instead of tau, eta, rho, and gamma; see details.
  }
	\item{alpha}{
		numeric (scalar): Statistical significance level.
	}
	\item{Xmat}{
		numeric (matrix): Design matrix for this SW CRT design.
	}
	\item{Wmat}{
		numeric (matrix): Covariance matrix for this SW CRT design.
	}
	\item{var.theta.WLS}{
		numeric (scalar): Variance estimate of \eqn{\theta} using weighted least squares (WLS) for this SW CRT design.
	}
	\item{pwrWLS}{
		numeric (scalar): Power of treatment effect (\eqn{\theta}) using weighted least squares (WLS) for this SW CRT design.
	}
%	\item{U}{
%		numeric (scalar): (Returned if \eqn{\eta \neq 0}--no random slopes. See closed-form formula from Hughes et al. (2003).)
%	}
%	\item{W}{
%		numeric (scalar): (Returned if \eqn{\eta \neq 0}--no random slopes. See closed-form formula from Hughes et al. (2003).)
%	}
%	\item{V}{
%		numeric (scalar): (Returned if \eqn{\eta \neq 0}--no random slopes. See closed-form formula from Hughes et al. (2003).)
%	}
%	\item{var.theta.CLOSED}{
%		numeric (scalar): Variance estimate of \eqn{\theta} using closed-form formula from Hughes et al. (2003) for this SW CRT design. Returned if \eqn{\eta \neq 0}--no random slopes.
%	}
	\item{pwrCLOSED}{
		numeric (scalar): Power of treatment effect (\eqn{\theta}) using closed-form formula from Hughes, et al. (2003) for this SW CRT design. Returned if eta and gamma are 0 and n is constant over clusters and time.
	}
}
\references{

Hussey MA, Hughes JP. Design and analysis of stepped wedge cluster randomized trials. Contemporary Clinical Trials 2007;28:182-191.
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
# Example 1 (Random Intercepts Only, standard Stepped Wedge (SW) design)
swPwr.Ex1.RIO.std <- swPwr(swDsn(c(6,6,6,6)), distn="binomial",
n=120, mu0=0.05, mu1=0.035, tau=0.01, eta=0, rho=0, gamma=0, alpha=0.05, retDATA=FALSE)
swPwr.Ex1.RIO.std

# Example 2 (Random Intercepts Only, extended SW design)
swPwr.Ex1.RIO.extend <- swPwr(swDsn(c(6,6,6,6), extra.time=3), distn="binomial",
n=120, mu0=0.05, mu1=0.035, tau=0.01, eta=0, rho=0, gamma=0,
alpha=0.05, retDATA=FALSE)
swPwr.Ex1.RIO.extend

# Example 3 (Independent Random Intercepts and Treatment effects, standard SW design)
swPwr.Ex1.IRIS <- swPwr(swDsn(c(6,6,6,6)), distn="binomial",
n=120, mu0=0.05, mu1=0.035, tau=0.01, eta=0.0045, rho=0, gamma=0,
alpha=0.05, retDATA=FALSE)
swPwr.Ex1.IRIS

# Example 4 (Correlated Random Intercepts and Slopes, standard SW design)
swPwr.Ex1.CRIS <- swPwr(swDsn(c(6,6,6,6)), distn="binomial",
n=120, mu0=0.05, mu1=0.035, tau=0.01, eta=0.0045, rho=0.4, gamma=0,
alpha=0.05, retDATA=FALSE)
swPwr.Ex1.CRIS

# Example 5 (Random time effect and correlated Random Intercepts and Slopes, standard SW design)
swPwr.Ex1.RTCRIS <- swPwr(swDsn(c(6,6,6,6)), distn="binomial",
n=120, mu0=0.05, mu1=0.035, tau=0.01, eta=0.0045, rho=0.4, gamma = 0.1,
alpha=0.05, retDATA=FALSE)
swPwr.Ex1.RTCRIS

#Example 6 (Sample size varying by cluster)
sample.size.vector <- c(35219,53535,63785,456132,128670,96673,
           51454,156667,127440,68615,56502,17719,
           75931,58655,52874,75936)
swPwr.Ex1.vector <- swPwr(swDsn(c(4,3,5,4)), distn="gaussian",
n=sample.size.vector, mu0=2.66, mu1=2.15,
sigma=sqrt(1/2.66), tau=0.31, eta=0.2, rho=0, gamma = 0.15,
alpha=0.05, retDATA=FALSE)
swPwr.Ex1.vector

#Example 7 (Sample size varying by cluster and time)
sample.size.matrix <- matrix(c(26, 493,  64,  45,  48, 231, 117,  17,  49,  36,  19,  77, 67, 590,
261, 212,  67, 318, 132,  58,  44,  57,  59,  78, 115, 532, 176, 199,  73, 293, 129,  79,  51,
62, 109,  94, 174, 785, 133,  79, 120, 305, 224,  99,  83,  79, 122, 122, 94, 961,  90, 131, 166,
352, 316,  59,  54, 131, 101, 133),nrow=12,ncol=5, byrow=FALSE)
swPwr.Ex1.matrix <- swPwr(swDsn(c(3,3,3,3)), distn="binomial",
n=sample.size.matrix, mu0=0.08, mu1=0.06, tau=0.017, eta=0.006, rho=-0.5, gamma = 0,
alpha=0.05, retDATA=FALSE)
swPwr.Ex1.matrix

#Example 8 (Using ICC and CAC)
swPwr.Ex1.icccac <- swPwr(swDsn(c(6,6,6,6)), distn="gaussian",
n=120, mu0=0.05, mu1=0.035, sigma=0.1, icc=0.02, cac=0.125, alpha=0.05, retDATA=FALSE)
swPwr.Ex1.icccac
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ power }
\keyword{ stepped wedge }% __ONLY ONE__ keyword per line
\keyword{ cluster randomized trial }

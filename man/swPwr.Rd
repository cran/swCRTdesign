\name{swPwr}
\alias{swPwr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Power of Stepped Wedge Cluster Randomized Trial (SW CRT)
}
\description{
\code{swPwr} returns (two-sided) power of treatment effect (\eqn{\theta}) for the specified SW CRT design using a linear models weighted least squares (WLS) approach. The response/outcome of interest can be binomial or Gaussian distributed and is assumed to come from a model with random intercepts and random treatment effects. If a random intercepts only model is used (i.e., \eqn{\eta=0}), then the power calculation is comparable to the closed-form formula of [Hussey and Hughes, 2007].  
}
\usage{
swPwr(design, distn, n, mu0, mu1, tau, eta, rho=0, sigma, alpha=0.05, retDATA=FALSE)
%
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{design}{
list: A stepped wedge design object, typically from swDsn, that includes at least the following components: 
##             swDsn, swDsn.unique.clusters, clusters, n.clusters, total.time
}
  \item{distn}{
character: Distribution assumed (\code{gaussian} or \code{binomial}).
}
  \item{n}{
integer (scalar): Number of observations for each cluster at each time point (same for all clusters at all time points).
}
  \item{mu0}{
numeric (scalar): Mean (outcome) in the control group.
}
  \item{mu1}{
numeric (scalar): Mean (outcome) in the treatment group. Note: Treatment effect is the difference in means \eqn{\theta = \mu_1 - \mu_0}.
}
  \item{tau}{
numeric (scalar): Standard deviation of random intercepts.
}
  \item{eta}{
numeric (scalar): Standard deviation of random treatment effects.
}
  \item{rho}{
numeric (scalar): Correlation between random intercepts and random treatment effects. Default value is \code{0}.
}
  \item{sigma}{
numeric (scalar): Standard deviation when assuming Gaussian distribution (\code{distn=gaussian}). For binomial distribution \eqn{\sigma^2} is automatically set to \eqn{\bar{\mu}(1-\bar{\mu})} where \eqn{\bar{\mu} = (\mu_1 + \mu_0)/2}
}
  \item{alpha}{
numeric (scalar): Two-sided statistical significance level.
}
  \item{retDATA}{
logical: if TRUE, all stored (input, intermediate, and output) values of \code{swPwr} are returned. Default value is \code{FALSE}.
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
		integer (scalar): Number of observations for each cluster at each time point (same for all clusters at all time points).
	}
	\item{mu0}{
		numeric (scalar): Mean (outcome) in the control group.
	}
	\item{mu1}{
		numeric (scalar): Mean (outcome) in intervention group. Note: treatment effect is difference in proportions \eqn{\theta = \mu_1 - \mu_0}.
	}
	\item{tau}{
		numeric (scalar): Standard deviation of random intercepts.
	}
	\item{eta}{
		numeric (scalar): Standard deviation of random slopes.
	}
	\item{rho}{
		numeric (scalar): Correlation between random intercepts and random slopes.
	}
	\item{sigma}{
		numeric (scalar): Standard deviation input. For binomial distribution, sigma = NA
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
		numeric (scalar): Returned if \eqn{\eta \neq 0}. Power of treatment effect (\eqn{\theta}) using closed-form formula from Hughes, et al. (2003) for this SW CRT design. Returned if \eqn{\eta \neq 0}--no random slopes. 
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
# Example 1 (Random Intercepts Only, standard Stepped Wedge (SW) design)
swPwr.Ex1.RIO.std <- swPwr(swDsn(c(6,6,6,6)), distn="binomial", 
n=120, mu0=0.05, mu1=0.035, tau=0.01, eta=0, rho=0, alpha=0.05, retDATA=FALSE)
swPwr.Ex1.RIO.std

# Example 2 (Random Intercepts Only, extended SW design)
swPwr.Ex1.RIO.extend <- swPwr(swDsn(c(6,6,6,6), extra.time=3), distn="binomial", 
n=120, mu0=0.05, mu1=0.035, tau=0.01, eta=0, rho=0, 
alpha=0.05, retDATA=FALSE)
swPwr.Ex1.RIO.extend

# Example 3 (Independent Random Intercepts and Treatment effects, standard SW design)
swPwr.Ex1.IRIS <- swPwr(swDsn(c(6,6,6,6)), distn="binomial", 
n=120, mu0=0.05, mu1=0.035, tau=0.01, eta=0.0045, rho=0,  
alpha=0.05, retDATA=FALSE)
swPwr.Ex1.IRIS


# Example 4 (Correlated Random Intercepts and Slopes, standard SW design)
swPwr.Ex1.CRIS <- swPwr(swDsn(c(6,6,6,6)), distn="binomial", 
n=120, mu0=0.05, mu1=0.035, tau=0.01, eta=0.0045, rho=0.4,  
alpha=0.05, retDATA=FALSE)
swPwr.Ex1.CRIS

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ power }
\keyword{ stepped wedge }% __ONLY ONE__ keyword per line
\keyword{ cluster randomized trial }
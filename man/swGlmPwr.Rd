\name{swGlmPwr}
\alias{swGlmPwr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Power of Stepped Wedge Cluster Randomized Trial with Discrete Outcomes
}
\description{
\code{swGlmPwr} returns (two-sided) power of the treatment effect for the specified SW CRT design in the context of generalized linear models by adopting the Laplace approximation detailed in Breslow and Clayton (1993) to obtain the covariance matrix of the estimated parameters. The response/outcome of interest can be binomial or Poisson distributed. 
The outcome is assumed to come from a model with fixed treatment effect (using an immediate treatment (IT) or exposure time indicator (ETI) model - see Kenny et al (2022)), fixed time effect, random intercepts, random treatment effects, and random cluster-specific time effects. The coefficients for fixed effects can be specified using \code{fixed.intercept}, \code{fixed.treatment.effect}, and \code{fixed.time.effect}. Variance components can be specified using \code{tau}, \code{eta}, \code{rho}, and \code{gamma}. 
}
\usage{
swGlmPwr(design, distn, n, fixed.intercept,
fixed.treatment.effect, fixed.time.effect, H = NULL,
tau = 0, eta = 0, rho = 0, gamma = 0, alpha=0.05, retDATA = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{design}{
list: A stepped wedge design object, typically from swDsn, that includes at least the   following components: swDsn, swDsn.unique.clusters, clusters, n.clusters, total.time
}
\item{distn}{
character: Distribution assumed (binomial or Poisson). "binomial" implies binomial outcomes and "poisson" implies Poisson outcome.
}
\item{n}{
integer (scalar, vector, or matrix): Number of observations: (scalar) for all clusters and all time points; (vector) for each cluster at all time points; and (matrix) for each cluster at each time point, where rows correspond to clusters and columns correspond to time. \code{n} can also be used to specify a design with transition periods (e.g. in the first time period that each sequence receives treatment, no observations are collected from that sequence). Simply define \code{n} as a matrix with a sample size of 0 during every transition period.
}
\item{fixed.intercept}{
numeric (scalar): Intercept for the fixed effect on canonical scales (logit for binomial outcomes and log for Poisson outcomes). It is the mean outcome under the control condition in the first time point transformed to the canonical scales.
}
\item{fixed.treatment.effect}{
numeric (scalar, vector): If H = NULL then an IT model is assumed and and fixed.treatment.effect is the scalar coefficient for the treatment in the fixed effect model on canonical scales (logit for binomial outcomes and log for Poisson outcomes). If H is non-NULL then an ETI model is assumed and fixed.treatment.effect is a vector as long as the longest treatment effect lag (typically, number of time periods minus one) giving the coefficient for the treatment effect on the canonical scale.
}
\item{fixed.time.effect}{
numeric(scalar, vector): Coefficients for the time (as dummy variables) in the fixed effect model on canonical scales (logit for binomial outcomes and log for Poisson outcomes). The first time point is always used as reference.
Specify a common time effect for all time points after the first (scalar) or differnt time effects for each time point after the first (vector of length (total time-1)).
}
  \item{H}{
numeric (vector): If NULL, then swGlmPwr assumes an immediate, constant treatment effect (IT) model. If not NULL, then an exposure time indicator (ETI) model is assumed and H should be a vector as long as the longest treatment effect lag (typically, number of time periods minus one). H specifies the desired linear combination of fixed.treatment.effect. For example, in a stepped wedge trial with 5 time periods and four exposure times, H = rep(.25,4) gives the average treatment effect over the four exposure times; H = c(0,0,.5,.5) ignores the first two periods after the intervention is introduced and averages the remaining periods. Typically, the sum of H is 1.0.
}
\item{tau}{
numeric (scalar): Standard deviation of random intercepts on canonical scales (logit for binomial outcomes and log for Poisson outcomes).
}
\item{eta}{ 
numeric (scalar): Standard deviation of random treatment effects on canonical scales (logit for binomial outcomes and log for Poisson outcomes).
}
\item{rho}{ 
numeric (scalar): Correlation between random intercepts and random treatment effects on canonical scales (logit for binomial outcomes and log for Poisson outcomes).
}
\item{gamma}{
numeric (scalar): Standard deviation of random time effects on canonical scales (logit for binomial outcomes and log for Poisson outcomes).
}
\item{alpha}{
numeric (scalar): Statistical significance level. Default is 0.05.
}
\item{retDATA}{
logical: if \code{TRUE}, all stored (input, intermediate, and output) values of \code{swGlmPwr} are returned. Default value is \code{FALSE}.
}
}
\details{
The two-sided statistical power of treatment effect \eqn{\theta} (equal to H\%*\%fixed.treatment.effect if H is non-NULL) is
\deqn{Pwr(\theta) = \Phi( \frac{Z - z_{1 - \alpha /2} \sqrt{V_0(\hat{\theta})}}{\sqrt{V_\alpha(\hat{\theta})}}) + 1 - \Phi( \frac{Z+ z_{1 - \alpha /2} \sqrt{V_0(\hat{\theta})}}{\sqrt{V_\alpha(\hat{\theta})}})},
where \eqn{\Phi} is the cumulative distribution function of the standard normal distribution. 

The variance of  \eqn{\hat{\theta}} under the null is denoted as \eqn{V_0(\hat{\theta})}, and the variance of \eqn{\hat{\theta}} under the alternative is denoted as \eqn{V_\alpha(\hat{\theta})}). Both variances are approximated by simplifying the Laplace approximation that marginalizes the random effects in the generalized linear mixed models. For more details, see Xia et al. (2020). 

When the outcome is Gaussian, the method adopted by \code{swGlmPwr} coincides with that of \code{swPwr}, so power calculation for Gaussian outcomes is not included in \code{swGlmPwr} to avoid repetition. When the outcome is binomial, \code{swGlmPwr} performs power calculation on the natural scale (logit), while \code{swPwr} performs power calculation on the linear scale. 

}
\value{
numeric (scalar): \code{swGlmPwr} returns the power of treatment effect if retDATA = FALSE.

numeric (list): \code{swGlmPwr} returns all specified and computed items as objects of a list if retDATA = TRUE.  

\item{$design}{
list: A stepped wedge design object, typically from swDsn, that includes at least the   following components: swDsn, swDsn.unique.clusters, clusters, n.clusters, total.time
}
\item{$distn}{
character: Distribution assumed (binomial or Poisson). "binomial" implies binomial outcomes and "poisson" implies Poisson outcome.  
}
\item{$n}{
integer (scalar, vector, or matrix): Number of observations: (scalar) for all clusters and all time points; (vector) for each cluster at all time points; and (matrix) for each cluster at each time point, where rows correspond to clusters and columns correspond to time. \code{n} can also be used to specify a design with transition periods (e.g. in the first time period that each sequence receives treatment, no observations are collected from that sequence). Simply define \code{n} as a matrix with a sample size of 0 during every transition period.
}
\item{$fixed.intercept}{
numeric (scalar): Intercept for the fixed effect on canonical scales (logit for binomial outcomes and log for Poisson outcomes).
}

\item{$fixed.treatment.effect}{
numeric (scalar): Coefficient for the treatment in the fixed effect model on canonical scales (logit for binomial outcomes and log for Poisson outcomes).
}
\item{$fixed.time.effect}{
numeric(scalar, vector): Coefficients for the time (as dummy variables) in the fixed effect model on canonical scales (logit for binomial outcomes and log for Poisson outcomes). The first time point is always used as reference.
 A common time effect for all time points after the first (scalar) or differnt time effects for each time point after the first (vector of length (total time-1)).
}
\item{$tau}{
numeric (scalar): Standard deviation of random intercepts on canonical scales (logit for binomial outcomes and log for Poisson outcomes).
}
\item{$eta}{ 
numeric (scalar): Standard deviation of random treatment effects on canonical scales (logit for binomial outcomes and log for Poisson outcomes).
}
\item{$rho}{ 
numeric (scalar): Correlation between random intercepts and random treatment effects on canonical scales (logit for binomial outcomes and log for Poisson outcomes).
}
\item{$gamma}{
numeric (scalar): Standard deviation of random time effects on canonical scales (logit for binomial outcomes and log for Poisson outcomes).
}
\item{$alpha}{
numeric (scalar): Statistical significance level. Default is 0.05.
}
\item{$var.theta.null}{
numeric (scalar): Variance estimate of the estimated treatment effect under the null for this SW CRT design.
}
\item{$var.theta.alt}{
numeric (scalar): Variance estimate of the estimated treatment effect under the alternative for this SW CRT design.
}
\item{$pwrGLM}{
numeric (scalar): Power of treatment effect using a simplified Laplace approximation.
}
}
\references{
Breslow NE and Clayton DG (1993). Approximate inference in generalized linear mixed models. Journal of the American Statistical Association, 88(421):9-25.

Kenny A, Voldal E, Xia F, Heagerty PJ, Hughes JP. Analysis of stepped wedge cluster randomized trials in the presence of a time-varying treatment effect. Statistics in Medicine, in press, 2022.

Xia F, Hughes JP, Voldal EC, Heagerty PJ. Power and sample size calculation for stepped-wedge designs with discrete outcomes. Trials. 2021 Dec;22(1):598.

}
\author{
Fan Xia, James P Hughes, and Emily C Voldal
}
%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%%}
\examples{
##test-case large clusters
library(swCRTdesign)
#specify large cluster sizes
size = c(35219,53535,63785,456132,128670,96673,
51454,156667,127440,68615,56502,17719,75931,58655,52874,75936)
#calculate power
swGlmPwr(design=swDsn(c(4,3,5,4)),distn="binomial",n=size,
fixed.intercept=log(28.62/(2*100000)),fixed.time.effect = 1,fixed.treatment.effect = log(.6),
tau=.31,eta=abs(0.4*log(.6)),rho=0,gamma=.15,alpha=.05,retDATA = FALSE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ design }
\keyword{ power }
\keyword{ discrete outcomes }
\keyword{ stepped wedge }
\keyword{ cluster randomized trial }% __ONLY ONE__ keyword per line

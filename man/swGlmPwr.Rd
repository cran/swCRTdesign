\name{swGlmPwr}
\alias{swGlmPwr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Power of Stepped Wedge Cluster Randomized Trial with Discrete Outcomes
}
\description{
\code{swGlmPwr} returns (two-sided) power of the treatment effect for the specified SW CRT design in the context of generalized linear models by adopting the Laplace approximation detailed in Breslow and Clayton (1993) to obtain the covariance matrix of the estimated parameters. The response/outcome of interest can be binomial (logit link only) or Poisson (log link only) distributed. A cross-sectional or closed cohort sampling scheme can be specified.
The outcome is assumed to come from a model with fixed treatment effect(s) (using an immediate treatment (IT) or exposure time indicator (ETI) model - see Kenny et al (2022)), fixed time effect, random intercepts, random treatment effects, random cluster-specific time effects  and, in the case of closed cohort sampling, an individual random effect. The coefficients for fixed effects can be specified using \code{fixed.intercept}, \code{fixed.treatment.effect}, and \code{fixed.time.effect}. Variance components can be specified using \code{tau}, \code{eta}, \code{rho}, \code{gamma} and \code{zeta}. 
}
\usage{
swGlmPwr(design, distn, n, fixed.intercept, fixed.treatment.effect, 
fixed.time.effect, H = NULL, tau = 0, eta = 0, rho = 0, gamma = 0, zeta = 0, 
alpha=0.05, retDATA = FALSE, silent=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{design}{
list: A stepped wedge design object, typically from \code{swDsn}, that includes at least the   following components: swDsn, swDsn.unique.clusters, clusters, n.clusters, total.time, nTxLev, TxLev. Fractional treatment effects specified in \code{swDsn} are used in the calculation.
}
\item{distn}{
character: Distribution assumed (binomial or Poisson). "binomial" implies binomial outcomes with logit link and "poisson" implies Poisson outcome with log link. ***NOTE: It is the users responsibility to make sure specified parameters (fixed.intercept, fixed.treatment effect, fixed.time effect, tau, eta, rho, gamma, zeta) are ALL on SAME scale as specified link function; see example.***
}
\item{n}{
integer (scalar, vector, or matrix): Number of observations: (scalar) for all clusters and all time points that are not NA in design$swDsn; (vector) for each cluster at all time points that are not NA in design$swDsn; and (matrix) for each cluster at each time point, where rows correspond to clusters and columns correspond to time. \code{n} can also be used to specify a design with transition periods (a transition period is a period, typically, immediately after treatment is introduced were no data are collected).  Simply define \code{n} as a matrix with a sample size of 0 during every transition period. This is equivalent to putting an NA in design$swDsn.
}
\item{fixed.intercept}{
numeric (scalar): Intercept for the fixed effect on canonical scales (logit for binomial outcomes and log for Poisson outcomes). It is the mean outcome under the control condition in the first time point transformed to the canonical scales.
}
\item{fixed.treatment.effect}{
numeric (scalar, vector): Gives the coefficients for the treatment effect on the canonical scale (logit for binomial outcomes and log for Poisson outcomes). 

If H = NULL then an IT model is assumed and \code{fixed.treatment.effect} gives the coefficient(s) for the \code{design$nTxLev} treatment effect(s). If \code{tx.effect.frac} is specified in \code{design} then the assumed treatment effect at exposure time k is \code{tx.effect.frac[k]*fixed.treatment.effect}. 

If H is non-NULL then an ETI model is assumed and fixed.treatment.effect is a vector as long as the longest treatment effect exposure time; in a classic stepped wedge design, the longest exposure time is the number of time periods minus one. 
}
\item{fixed.time.effect}{
numeric(scalar, vector): Coefficients for the time (as dummy variables) in the fixed effect model on canonical scales (logit for binomial outcomes and log for Poisson outcomes). The first time point is always used as reference.
Specify a common time effect for all time points after the first (scalar) or differnt time effects for each time point after the first (vector of length (total time-1)).
}
  \item{H}{
numeric (vector): If NULL, then swGlmPwr assumes an immediate, constant treatment effect (IT) model. If not NULL, then an exposure time indicator (ETI) model is assumed and H should be a vector as long as the longest exposure time (typically, number of time periods minus one). H specifies the desired linear combination of fixed.treatment.effect. For example, in a stepped wedge trial with 5 time periods and four exposure times, H = rep(.25,4) gives the average treatment effect over the four exposure times; H = c(0,0,.5,.5) ignores the first two periods after the intervention is introduced and averages the remaining periods. Typically, the sum of H is 1.0; if not, it is renormalized to sum to 1.0. H can only be specified when there is a single intervention level (i.e. design$swDsn includes only NA,0,1; see cautions in help for \code{swDsn} about including NA periods when H is used). Do not use an ETI model (non-null \code{H}) if you specified a fractional treatment effect in \code{swDsn}.
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
\item{zeta}{
numeric (scalar): Standard deviation of individual effects for closed cohort sampling on canonical scales (logit for binomial outcomes and log for Poisson outcomes). Default is 0 which implies cross-sectional sampling. Values greater than 0 imply closed cohort sampling.
}
\item{alpha}{
numeric (scalar): Statistical significance level. Default is 0.05.
}
\item{retDATA}{
logical: if \code{TRUE}, all stored (input, intermediate, and output) values of \code{swGlmPwr} are returned. Default value is \code{FALSE}.
}
\item{silent}{
logical: if TRUE, hides most warning messages.  Default value is \code{FALSE}.
}

}
\details{
The two-sided statistical power of treatment effect \eqn{\theta} (equal to H\%*\%fixed.treatment.effect if H is non-NULL) is
\deqn{Pwr(\theta) = \Phi( \frac{Z - z_{1 - \alpha /2} \sqrt{V_0(\hat{\theta})}}{\sqrt{V_\alpha(\hat{\theta})}}) + 1 - \Phi( \frac{Z+ z_{1 - \alpha /2} \sqrt{V_0(\hat{\theta})}}{\sqrt{V_\alpha(\hat{\theta})}})},
where \eqn{\Phi} is the cumulative distribution function of the standard normal distribution. 

The variance of  \eqn{\hat{\theta}} under the null is denoted as \eqn{V_0(\hat{\theta})}, and the variance of \eqn{\hat{\theta}} under the alternative is denoted as \eqn{V_\alpha(\hat{\theta})}). Both variances are approximated by simplifying the Laplace approximation that marginalizes the random effects in the generalized linear mixed models. For more details, see Xia et al. (2020). 

When the outcome is Gaussian, the method adopted by \code{swGlmPwr} coincides with that of \code{swPwr}, so power calculation for Gaussian outcomes is not included in \code{swGlmPwr} to avoid repetition. When the outcome is binomial, \code{swGlmPwr} performs power calculation on the natural scale (logit), while \code{swPwr} performs power calculation on the linear scale. 

The value of zeta defines the samling scheme. When zeta = 0, cross-sectional sampling is assumed; if zeta > 0 then closed cohort sampling is assumed.
}
\value{
numeric (scalar): \code{swGlmPwr} returns the power of treatment effect if retDATA = FALSE.

numeric (list): \code{swGlmPwr} returns all specified and computed items as objects of a list if retDATA = TRUE.  

\item{$design}{
list: A stepped wedge design object, typically from swDsn, that includes at least the   following components: swDsn, swDsn.unique.clusters, clusters, n.clusters, total.time, nTxLev
}
\item{$distn}{
character: Distribution assumed (binomial or Poisson). "binomial" implies binomial outcomes and "poisson" implies Poisson outcome.  
}
\item{$n}{
integer (scalar, vector, or matrix): Number of observations: (scalar) for all clusters and all time points; (vector) for each cluster at all time points; and (matrix) for each cluster at each time point, where rows correspond to clusters and columns correspond to time. 
}
\item{$fixed.intercept}{
numeric (scalar): Intercept for the fixed effect on canonical scales (logit for binomial outcomes and log for Poisson outcomes).
}

\item{$fixed.treatment.effect}{
numeric (scalar): Coefficient(s) for the treatment(s) in the fixed effect model on canonical scales (logit for binomial outcomes and log for Poisson outcomes).
}
\item{$fixed.time.effect}{
numeric(scalar, vector): Coefficients for the time (as dummy variables) in the fixed effect model on canonical scales (logit for binomial outcomes and log for Poisson outcomes). The first time point is always used as reference.
 A common time effect for all time points after the first (scalar) or differnt time effects for each time point after the first (vector of length (total time-1)).
}
\item{H}{
  numeric (vector): H specifies the desired linear combination of exposure time treatment effects for a ETI model-based estimate.
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
\item{$zeta}{
numeric (scalar):  Standard deviation of individual random effects for closed cohort sampling on canonical scales (logit for binomial outcomes and log for Poisson outcomes).
}
\item{$alpha}{
numeric (scalar): Statistical significance level. Default is 0.05.
}
\item{$var.theta.null}{
numeric (martix): Variance-covariance martix of the estimated treatment effect(s) under the null for this SW CRT design.
}
\item{$var.theta.alt}{
numeric (marix): Variance-covariance matrix of the estimated treatment effect(s) under the alternative for this SW CRT design.
}
\item{$pwrGLM}{
numeric (scalar): Power of treatment effect using a simplified Laplace approximation.
}
}
\references{
Breslow NE and Clayton DG (1993). Approximate inference in generalized linear mixed models. Journal of the American Statistical Association, 88(421):9-25.

Kenny A, Voldal E, Xia F, Heagerty PJ, Hughes JP. Analysis of stepped wedge cluster randomized trials in the presence of a time-varying treatment effect. Statistics in Medicine, in press, 2022.

Voldal EC, Hakhu NR, Xia, F, Heagerty PJ, Hughes JP. swCRTdesign: An R Package for Stepped Wedge Trial Design and Analysis. Computer Methods and Programs in Biomedicine 2020;196:105514.

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
logit <- function(x){log(x/(1 - x))}
#specify large cluster sizes
size = c(35219,53535,63785,456132,128670,96673,
51454,156667,127440,68615,56502,17719,75931,58655,52874,75936)

# Cross-sectional example
swGlmPwr(design=swDsn(c(4,3,5,4)),distn="binomial",n=size,
        fixed.intercept=log(28.62/(2*100000)),fixed.time.effect = 1,
        fixed.treatment.effect = log(.6),
        tau=.31,eta=abs(0.4*log(.6)),rho=0,gamma=.15,alpha=.05,retDATA = FALSE)

# Closed cohort example, comparing average of intervention period exposure 
# lags 3 and 4 to control period
swGlmPwr(design=swDsn(c(5,5,5,5,5),extra.ctl.time=3,extra.trt.time=5),
         distn="binomial",n=20,
         fixed.intercept=logit(0.40),
         fixed.treatment.effect=c(0,0,rep(logit(0.60)-logit(0.40),8)),
         fixed.time.effect=0.08,
         H=c(0,0,1,1,0,0,0,0,0,0),
         tau=sqrt(.1316),gamma=sqrt(.1974),eta=0,zeta=sqrt(2.5))

# Example wih periods with no data and multiple treatment levels
stdy <- swDsn(c(6,6,6,6),swBlk=matrix(c(0,1,2,2,2,2,
                                         NA,0,1,2,2,2,
                                         NA,NA,0,1,2,2,
                                         NA,NA,NA,0,1,2),4,6,byrow=TRUE))
swGlmPwr(stdy, distn="binomial",n=10000, 
         fixed.intercept=log(28.62/(2*100000)),fixed.time.effect = 1,
         fixed.treatment.effect = c(log(.6),log(.5)),
         tau=.31,eta=0,rho=0,gamma=0,alpha=.05,retDATA = FALSE)

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

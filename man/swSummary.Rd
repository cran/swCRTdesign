\name{swSummary}
\alias{swSummary}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Summary of Response/Outcome for Stepped Wedge Cluster Randomized Trial (SW CRT)
}
\description{
\code{swSummary} returns the mean, sum, and number of non-missing values for the response/outcome variable of interest for each cluster at each time point from a SW CRT.
}
\usage{
swSummary(response.var, tx.var, time.var, cluster.var, data, 
type="mean", digits=16, fcn.Call=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{response.var}{
numeric(vector): Response (Outcome) variable.
}
  \item{tx.var}{
numeric(vector): Treatment (Predictor of Interest) variable. Typically, 0=placebo, 1=intervention, and values greater than 1 correspond to other treatment levels. 
}
  \item{time.var}{
integer(vector): Time (points) variable, corresponding to the time points when data were collected during the SW CRT.
}
  \item{cluster.var}{
integer(vector): Cluster (identification) variable, corresponding to the cluster where an individual is from.
}
  \item{data}{
An optional data frame containing (at least) the response, treatment (tx), time, and cluster variables.
}
  \item{type}{
character (scalar): Specify which summary measure is of interest from \code{"mean"}, \code{"sum"}, and \code{"n"}. (Note: The default returns "mean" as the summary measure of response. Note that all summary measures may not be scientifically relevant in every situation.)
}
  \item{digits}{
integer (scalar): Number of places right of the decimal. The default value is 16.
}
  \item{fcn.Call}{
logical: Only \code{TRUE} when calling \code{swSummary} from within \code{swPlot}. The default is \code{FALSE}.
}
}
\details{
Returns a list containing a matrix of dimension \code{length(unique(data$cluster))} by \code{length(} \code{unique(data$time))} that summarizes \code{data$response} for specified \code{type}. Either the mean, sum, or the number of non-missing \code{data$response} values may be requested using \code{type}. \code{dimnames} of the matrix correspond to the unique values of cluster and time. Note that the stepping pattern in the data may be obtained by specifying the treatment variable name as the response and \code{type = "mean"}.
}
\value{
numeric (list): \code{swSummary} returns a list containing the following
	\item{type}{
One of user-specified options \code{"mean"}, \code{"sum"}, or \code{"n"}.
}
	\item{swDsn}{
The SW design.		
}
	\item{swDsn.unique.clusters}{
The unique clusters (i.e., rows) SW design.		
}
	\item{n.waves}{
Number of waves. Note that this is usually, but not necessarily, the same as data$n.waves. See example 2.
}
	\item{clusters}{
Clusters per wave.		
}
	\item{n.clusters}{
Total number of clusters.		
}
	\item{time.at.each.wave}{
Time at each wave.		
}
	\item{total.time}{
Total time points.		
}
	\item{response.cluster}{
numeric (matrix): Response variable summarized according to \code{type} for all clusters, with dimension 
\code{length(data$cluster)} by \code{length(unique(data$time))}.
}
	\item{response.wave}{
numeric (matrix): Response variable summarized according to \code{type}, for all waves (all clusters of a particular wave are combined), with dimension \code{length(} \code{unique(data$cluster))} by \code{length(unique(data$time))}.
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
# Example 1 
# Generate binary response with 5 clusters, 4 time points
data.Ex1 <- swSim(swDsn(c(2,2,1)),family=binomial(link="identity"), n=120,
                     mu0=0.25, mu1=0.30,time.effect=0, tau=0.05, gamma=.01)

## Example 1 (type="mean", by cluster and by wave)
swSummary.Ex1.mean <- swSummary(response.var, tx.var, time.var, cluster.var, 
data=data.Ex1, type="mean", digits=3)

swSummary.Ex1.mean$response.cluster

swSummary.Ex1.mean$response.wave

# Example 1 (type="sum", by cluster and by wave)
swSummary.Ex1.sum <- swSummary(response.var, tx.var, time.var, cluster.var, 
data=data.Ex1, type="sum")

swSummary.Ex1.sum$response.cluster

swSummary.Ex1.sum$response.wave

## Example 1 (type="n", by cluster and by wave)
swSummary.Ex1.n <- swSummary(response.var, tx.var, time.var, cluster.var, 
data=data.Ex1, type="n")

swSummary.Ex1.n$response.cluster

swSummary.Ex1.n$response.wave

## Example 2
design2 <- swDsn(clusters=c(6,6,6,6)) 
# design2$n.waves says 4 waves
nmat2 = rbind(matrix(rep(c(120,120,120,120,120),6),6,5,byrow=TRUE),
              matrix(rep(c(120,120,120,120,120),3),3,5,byrow=TRUE),
              matrix(rep(c(0,120,120,120,120),3),3,5,byrow=TRUE),
              matrix(rep(c(0,0,120,120,120),6),6,5,byrow=TRUE),
              matrix(rep(c(0,0,0,120,120),6),6,5,byrow=TRUE))
swGenData2 <- swSim( design2,family=binomial(link="logit"), n=nmat2,
                     mu0=log(0.1/0.9), mu1=log(0.9) + log(0.1/0.9),
                     time.effect=0, tau=0.2, time.lab=seq(0,12,3))
swSummary(response.var, tx.var, time.var, cluster.var, swGenData2,
          type="mean", digits=3)
# This summary create 5 waves
#
# The following is equivalent to above
design2b <- swDsn(c(6,3,3,6,6),
                   swBlk=matrix(c(0,1,1,1,1,
                                  0,0,1,1,1,
                                  NA,0,1,1,1,
                                  NA,NA,0,1,1,
                                  NA,NA,NA,0,1),5,5,byrow=TRUE)) 
# design2b$n.waves says 5 waves
swGenData2b <- swSim(design2b,family=binomial(link="logit"), n=120,
                      mu0=log(0.1/0.9), mu1=log(0.9) + log(0.1/0.9),
                      time.effect=0, tau=0.2, time.lab=seq(0,12,3))
swSummary(response.var, tx.var, time.var, cluster.var, swGenData2b,
          type="mean", digits=3)
# This summary create 5 waves          
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ summary }
\keyword{ stepped wedge }
\keyword{ cluster randomized trial }% __ONLY ONE__ keyword per line

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
numeric(vector): Treatment (Predictor of Interest) variable. Typically, 0=placebo, 1=intervention, values between 0 and 1 correspond to fractional treatment/intervention effects, and values greater than 1 correspond to other treatment options. 
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
Number of waves.
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
# Example 1 (Generating SW CRT data)
# (binary response with 1 missing value, 5 clusters, 4 time points)

n.Ex1 <- 120
p0.Ex1 <- 0.05
clusters.Ex1 <- c(2,2,1)
dsn.Ex1 <- swDsn( clusters=clusters.Ex1 )
time.Ex1 <- c(1:dsn.Ex1$total.time)*4 - 4
response.Ex1 <- rbinom(n.Ex1 * dsn.Ex1$n.clusters * dsn.Ex1$total.time, 1, p0.Ex1)
response.Ex1[1] <- NA
tx.Ex1 <- as.vector( apply( dsn.Ex1$swDsn, 1, function(z){rep(z, n.Ex1)}) )
time.Ex1 <- rep( time.Ex1, n.Ex1 * dsn.Ex1$n.clusters )
cluster.Ex1 <- rep( 1:dsn.Ex1$n.clusters, each=n.Ex1 * dsn.Ex1$total.time )
data.Ex1 <- data.frame(response.Ex1, tx.Ex1, time.Ex1, cluster.Ex1)

## Example 1 (type="mean", by cluster and by wave)
swSummary.Ex1.mean <- swSummary(response.Ex1, tx.Ex1, time.Ex1, cluster.Ex1, 
data=data.Ex1, type="mean", digits=3)

swSummary.Ex1.mean$response.cluster

swSummary.Ex1.mean$response.wave

# Example 1 (type="sum", by cluster and by wave)
swSummary.Ex1.sum <- swSummary(response.Ex1, tx.Ex1, time.Ex1, cluster.Ex1, 
data=data.Ex1, type="sum")

swSummary.Ex1.sum$response.cluster

swSummary.Ex1.sum$response.wave

## Example 1 (type="n", by cluster and by wave)
swSummary.Ex1.n <- swSummary(response.Ex1, tx.Ex1, time.Ex1, cluster.Ex1, 
data=data.Ex1, type="n")

swSummary.Ex1.n$response.cluster

swSummary.Ex1.n$response.wave

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ summary }
\keyword{ stepped wedge }
\keyword{ cluster randomized trial }% __ONLY ONE__ keyword per line

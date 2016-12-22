\name{swPlot}
\alias{swPlot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot of Mean Response/Outcome for Stepped Wedge Cluster Randomized Trial (SW CRT)
}
\description{
\code{swPlot} returns plot of the mean response versus time based on waves and/or clusters from a SW CRT.
}
\usage{
swPlot(response.var, tx.var, time.var, cluster.var, data, choose.mfrow=NULL, 
by.wave=TRUE, combined.plot=TRUE, choose.xlab="Time", choose.main=NULL, 
choose.pch=NULL, choose.cex=1, choose.tx.col=NULL, choose.ncol=2, 
choose.tx.pos="topright", choose.legend.pos=locator(1))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{response.var}{
numeric(vector): Response (Outcome) variable.
}
  \item{tx.var}{
numeric(vector): Treatment (Predictor of Interest) variable. Typically, 0=control, 1=intervention, values between 0 and 1 correspond to fractional treatment/intervention effects, and values greater than 1 correspond to other treatment options. 
}
  \item{time.var}{
integer(vector): Time (points) variable, corresponding to the time points when data were collected during the SW CRT.
}
  \item{cluster.var}{
integer(vector): Cluster (identification) variable, corresponding to the cluster where an observation is from.
}
  \item{data}{
An optional data frame containing (at least) the response, treatment (tx), time, and cluster variables.
}
  \item{choose.mfrow}{
numeric (vector): Choose \code{mfrow} for plot. If \code{NULL}, \code{mfrow} is automatically assigned based on the plot created. The default is \code{NULL}.
}
  \item{by.wave}{
logical: If \code{TRUE}, plot mean response according to each wave. If \code{FALSE}, plot mean response according to each cluster. The default is \code{TRUE}.
}
  \item{combined.plot}{
logical: If \code{TRUE}, plot mean response on same plot (what is plotted depends on \code{by.wave}). If \code{FALSE}, plot mean response on separate plots for each wave (what is plotted depends on \code{by.wave}). The default is \code{TRUE}.
}
  \item{choose.xlab}{
Choose \code{xlab} for plot. The default is "Time".
}
  \item{choose.main}{
Choose \code{main} for plot. If \code{NULL}, \code{main} is chosen for the user; which is highly recommended. The default is \code{NULL}.
}
  \item{choose.pch}{
Choose \code{pch} for plot. If \code{NULL}, \code{pch} are chosen for the user; which is highly recommended. The default is \code{NULL}.
}
  \item{choose.cex}{
Choose \code{cex} for \code{choose.pch}. Standard \code{cex} option in \code{points()} applies. The default is \code{1}.
}
  \item{choose.tx.col}{
Choose colors for different treatment options (at least two: 0=placebo, 1=intervention). If \code{NULL}, colors are chosen for the unique treatment options in the data supplied. The default is \code{NULL}.
}
  \item{choose.ncol}{
Choose number of columns for non-treatment legend. Standard \code{ncol} option in \code{legend()} applies. The default is \code{2}.
}
  \item{choose.tx.pos}{
Choose where to place treatment colors "legend". Standard \code{legend()} positions apply. The default is "topright".
}
  \item{choose.legend.pos}{
Choose where to place the non-treatment legend. Standard \code{legend()} positions apply. The default is \code{locator(1)}, user specifies location with mouse/trackpad. Specify \code{NULL} to remove legend.
}
}
\details{
Returns a plot of the mean response versus time with a combination of by wave (\code{TRUE} / \code{FALSE}) and combined plot (\code{TRUE} / \code{FALSE}) from a SW CRT.
}
%\value{
%
%}
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

# Example 1 (Mean Response vs Time, by.wave=TRUE, combined.plot=TRUE)
swPlot(response.Ex1, tx.Ex1, time.Ex1, cluster.Ex1, data.Ex1, by.wave=TRUE, 
combined.plot=TRUE, choose.tx.pos="bottomright", choose.legend.pos="bottom")

# Example 2 (Mean Response vs Time, by.wave=TRUE, combined.plot=FALSE)
swPlot(response.Ex1, tx.Ex1, time.Ex1, cluster.Ex1, data.Ex1, by.wave=TRUE, 
combined.plot=FALSE, choose.tx.pos="bottomright", choose.legend.pos="bottom")

# Example 3 (Mean Response vs Time, by.wave=FALSE, combined.plot=TRUE)
swPlot(response.Ex1, tx.Ex1, time.Ex1, cluster.Ex1, data.Ex1, by.wave=FALSE, 
combined.plot=TRUE, choose.tx.pos="bottomright", choose.legend.pos="bottom")

# Example 4 (Mean Response vs Time, by.wave=FALSE, combined.plot=FALSE)
swPlot(response.Ex1, tx.Ex1, time.Ex1, cluster.Ex1, data.Ex1, by.wave=FALSE, 
combined.plot=FALSE, choose.tx.pos="bottomright", choose.legend.pos="bottom")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ plot }
\keyword{ stepped wedge }
\keyword{ cluster randomized trial }% __ONLY ONE__ keyword per line

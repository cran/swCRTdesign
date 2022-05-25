\name{swCRTdesign-package}
\alias{swCRTdesign-package}
\alias{swCRTdesign}
\docType{package}
\title{
Stepped Wedge Cluster Randomized Trial (SW CRT) Design
}
\description{
This package includes functions for the design and analysis of stepped wedge cluster randomized trials according to a repeated cross-sectional sampling scheme. For additional guidance, see (Voldal EC, Hakhu NR, Xia F, Heagerty PJ, Hughes JP. swCRTdesign: An R package for stepped wedge trial design and analysis. Computer Methods and Programs in Biomedicine 2020;196:105514. <doi:10.1016/j.cmpb.2020.105514>). Four primary functions - \code{swPwr}, \code{swSim}, \code{swSummary}, and \code{swPlot} - and two support functions - \code{blkDiag}, \code{swDsn} - are included. The \code{blkDiag} function creates a block diagonal matrix from a specified array or list of block-matrices. The \code{swDsn} function creates a stepped wedge (SW) design object based on specified information on clusters, time points, and the two arms of the cluster randomized trial (CRT). The \code{swPwr} function computes the (two-sided) power of treatment effect (\eqn{\theta}) for the specified SW CRT design via weighted least squares (WLS), where the response/outcome of interest is assumed to come from a mixed effects model with random time effects and (possibly correlated) random intercepts and random treatment effects. The random time effects apply to all time points, and time is treated as categorical.  The \code{swSim} function generates individual-level data consisting of response, treatment, time, and cluster variables based on a specified SW CRT design. The \code{swSummary} function computes the mean, sum, or number of non-missing response values for clusters separately or aggregated by wave at each time point from stepped wedge data that includes, at least, response, treatment, time, and cluster variables. The \code{swPlot} function plots mean response as a combined or separate plot, for waves and clusters.  Some features of the package are also available as a shiny app, available online (https://swcrtdesign.shinyapps.io/stepped_wedge_power_calculation/) or to download and run locally (https://github.com/swCRTdesign/Stepped-wedge-power-calculation).
}
% If assuming the response of interest comes from a random intercepts only model, then \code{swPwr} returns (two-sided) power that is equivalent to the closed-form formula of Hughes et al., 2003.
\details{
\tabular{ll}{
Package: \tab swCRTdesign\cr
Type: \tab Package\cr
Version: \tab 3.2\cr
Date: \tab 2022-05-24\cr
License: \tab GPL (>= 2)\cr
}
%~~ An overview of how to use the package, including the most important functions ~~
}
\author{
James P Hughes, Navneet R Hakhu, Emily C Voldal, and Fan Xia

Maintainer: James P Hughes <jphughes@uw.edu>
}
\references{

Hussey MA, Hughes JP. Design and analysis of stepped wedge cluster randomized trials. Contemporary Clinical Trials 2007;28:182-191.

Voldal EC, Hakhu NR, Xia F, Heagerty PJ, Hughes JP. swCRTdesign: An R package for stepped wedge trial design and analysis. Computer Methods and Programs in Biomedicine 2020;196:105514.
}
%~~ Optionally other standard keywords, one per line, from file KEYWORDS in the R documentation directory ~~
\keyword{ stepped wedge }
\keyword{ cluster randomized trial }
\keyword{ repeated cross-sectional sampling }
\keyword{ design }
\keyword{ power }
\keyword{ simulated data set }
\keyword{ summary }
\keyword{ plot }
\keyword{ block diagonal }
%\seealso{
%~~ Optional links to other man pages, e.g. ~~
%~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
%}
%\examples{
%~~ simple examples of the most important functions ~~
%}

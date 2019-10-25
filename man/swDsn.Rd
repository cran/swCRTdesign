\name{swDsn}
\alias{swDsn}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Study design of Stepped Wedge Cluster Randomized Trial (SW CRT)
}
\description{
\code{swDsn} returns a SW CRT study design object based on a repeated cross-sectional sampling scheme. All clusters that start the intervention at a given time point are collectively referred to as a wave or sequence. There many be a variable number of clusters in each wave. By default, all clusters are assumed to start in the control condition. Fractional treatment effect may be specified for each time after the intervention is introduced. Additional observations may be added to the end of the trial after the intervention has been introduced in all clusters. For incorporating transition periods where no data is collected, see \code{swPwr}.  \code{swDsn} is used by other functions in this package.
}
\usage{
swDsn(clusters, tx.effect.frac = 1, extra.time = 0, all.ctl.time0 = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{clusters}{
integer (vector): Number of clusters for each wave (e.g. c(6,6,6,6) specifies four waves with 6 clusters in each wave). A value of 0 in the vector means that no clusters introduce the intervention at a given time (see examples).
}
  \item{tx.effect.frac}{
numeric (scalar or vector): Fractional treatment effect upon crossing over from control. Note that this is not the treatment effect!  If a scalar with value of 1, the standard SW CRT treatment effect will be presumed. If a scalar with a fractional value between 0 and 1, then only the first time point upon crossing over from control will have fractional treatment effect; the remaining time points in SW CRT design will have value of 1. If a vector of fractional treatment effect is specified, each element of the vector corresponds to the (fractional) treatment effect upon crossing over from control; if length of vector less than total number of time points after crossing over, the remaining time points will have treatment effect value of 1; if length of vector greater than total number of time points after crossing over, not all elements of vector will be used. The default value is (scalar) 1.
}
  \item{extra.time}{
integer (scalar): Number of additional time steps beyond the standard SW CRT design (standard + extended times corresponds to total time). The default value is 0.
}
  \item{all.ctl.time0}{
logical: If \code{TRUE}, all clusters receive control at the first time point. If \code{FALSE}, clusters in the first wave (i.e., the first element of \code{clusters}) receive intervention at the first time. The default is \code{TRUE}.
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
numeric (list): Returns the following user-specified and computed objects
%% *****SPECIFY WHAT THESE ARE!!!!!**********
	\item{swDsn}{
		numeric (matrix):  schematic representation of the specified SW CRT design. Number of clusters is equal to the number of rows of the matrix and number of time intervals is equal to the number of columns of the matrix. \code{swDsn[i,j]} gives the intervention status for cluster \code{i} at time \code{j} and has possible values 0 (control), 1 (intervention) or a fractional value as specified by \code{tx.effect.frac}.
	}
	\item{swDsn.unique.clusters}{
		numeric (matrix): Truncated SW CRT design of interest, with one row for each wave.
	}
	\item{n.waves}{
		numeric (scalar): Number of waves for the SW CRT design of interest.
	}
	\item{clusters}{
		numeric (vector): Number of clusters changing from control to intervention at each wave for the SW CRT design of interest.
	}
	\item{n.clusters}{
		numeric (scalar): Total number of clusters for the SW CRT design of interest.
	}
	\item{tx.effect.frac}{
		numeric (scalar or vector): Fractional treatment effect for time points upon crossing over from control of SW CRT design of interest.
	}
	\item{total.time}{
		numeric (scalar): Total number of time points for the SW CRT design of interest.
	}
	\item{extra.time}{
		numeric (matrix): Number of time points added on to the standard SW CRT time points for the user-specified values (i.e., extending the trial).
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
# Example 1 (Equal clusters per wave, standard SW design)
swDsn.Ex1.std <- swDsn( clusters=c(3,3,3) )
swDsn.Ex1.std$swDsn

# Example 2 (Equal clusters per wave, extended SW design)
swDsn.Ex1.extend <- swDsn( clusters=c(3,3,3), extra.time=2 )
swDsn.Ex1.extend$swDsn

# Example 3 (Equal clusters per wave, not all ctl at time 0, "standard" for time SW design)
swDsn.Ex1.std.noAllctl <- swDsn( clusters=c(3,3,3), all.ctl.time0=FALSE )
swDsn.Ex1.std.noAllctl$swDsn

# Example 4 (Equal clusters per wave, not all ctl at time 0, extended SW design)
swDsn.Ex1.extend.noAllctl <- swDsn( clusters=c(3,3,3), extra.time=2, all.ctl.time0=FALSE )
swDsn.Ex1.extend.noAllctl$swDsn

# Example 5 (Unequal clusters per wave, standard SW design)
swDsn.Ex1.std.unequal <- swDsn( clusters=c(3,0,2) )
swDsn.Ex1.std.unequal$swDsn

# Example 6 (Unequal clusters per wave, extended SW design)
swDsn.Ex1.extend.unequal <- swDsn( clusters=c(3,0,2), extra.time=2 )
swDsn.Ex1.extend.unequal$swDsn

# Example 7 (Unequal clusters per wave, extended SW design)
swDsn.Ex1.extend.unequal.varyTxEffect <- swDsn( clusters=c(3,0,2), tx.effect.frac=c(0.8,0.9,1.0),
extra.time=2 )
swDsn.Ex1.extend.unequal.varyTxEffect$swDsn

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ design }
\keyword{ stepped wedge }
\keyword{ cluster randomized trial }% __ONLY ONE__ keyword per line

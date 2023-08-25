\name{swDsn}
\alias{swDsn}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Study design of Stepped Wedge Cluster Randomized Trial (SW CRT)
}
\description{
\code{swDsn} returns a SW CRT study design object. 

All calls to \code{swDsn} require a \code{clusters} argument, a vector that specifies the number of clusters in each wave (all clusters that start the intervention at a given time point are collectively referred to as a wave or sequence). Used by itself, \code{clusters} specifies a classic stepped wedge design with number of waves equal to the length of \code{clusters} and number of time points equal to number of waves plus one. Further design modifications may be specified using the \code{extra.ctl.time}, \code{extra.trt.time}, and \code{all.ctl.time0} arguments. Fractional treatment effects may be specified for each time after the intervention is introduced using the \code{tx.effect.frac} argument.

Alternatively, the \code{swBlk} argument can be used to specify a completely arbitrary stepped wedge design with multiple treatment levels and/or periods with no data collection (alternatively, periods with no data collection can  be defined in the sample size arguments of \code{swPwr}, \code{swGlmPwr} and \code{swSimPwr}).

One can choose to use the \code{extra.ctl.time}, \code{extra.trt.time}, and \code{all.ctl.time0} arguments or the \code{swBlk} argument, but not both. Also, \code{tx.effect.frac} may not be used with the \code{swBlk} argument.

\code{swDsn} is used by other functions in this package. 
}
\usage{
swDsn(clusters, tx.effect.frac=1, extra.ctl.time=0, extra.trt.time=0, 
all.ctl.time0=TRUE, swBlk, extra.time=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{clusters}{
integer (vector): Number of clusters for each wave (e.g. c(6,6,6,6) specifies four waves with 6 clusters in each wave). A value of 0 in the vector means that no clusters introduce the intervention at a given time (see examples). Note that 0 clusters in a wave is not allowed when using the swBlk specification.
}
  \item{tx.effect.frac}{
numeric (scalar or vector): Relative treatment effect upon crossing over from control. Note that this is not the treatment effect (the treatment effect is specified in the call to other functions such as \code{swPwr} or \code{swSim})!  If \code{tx.effect.frac} is set to a scalar with value of 1, the standard SW CRT treatment effect will be presumed. If code{tx.effect.frac} is set to a scalar with a fractional value between 0 and 1, then only the first time point upon crossing over from control will have fractional relative treatment effect; the remaining exposure times will have a relative treatment effect of 1. If a vector of fractional treatment effect is specified, each element of the vector corresponds to the (fractional) treatment effect upon crossing over from control; if the length of the vector is less than total number of exposure times, the remaining time points will have a relative treatment effect value of 1; if the length of the vector is greater than total number of time points after crossing over, not all elements of \code{tx.effect.frac} will be used. The default value is (scalar) 1.
}
  \item{extra.ctl.time}{
integer (scalar): Number of additional time steps during the control period beyond the standard SW CRT design. The default value is 0.
}
  \item{extra.trt.time}{
integer (scalar): Number of additional time steps during the treatment period beyond the standard SW CRT design. The default value is 0.
}
  \item{all.ctl.time0}{
logical: If \code{TRUE}, all clusters receive control at the first time point. If \code{FALSE}, clusters in the first wave (i.e., the first element of \code{clusters}) receive intervention at the first time. The default is \code{TRUE}.
}
  \item{swBlk}{
integer (matrix): matrix with number of rows equal to the number of waves and number of columns equal to the total number of time periods of the study. Each row of \code{swBlk} corresponds to a wave and describes the introduction and levels of the intervention. Each element of \code{swBlk} is either NA, meaning no data will be collected in that cluster-time period, or an integer. Integers correspond to levels of the intervention where 0 is the control condition and values >0 correspond to various levels of the intervention. A classic SW design includes only 2 values (0 and 1). More generally, one can specify integers corresponding to different levels of the intervention (e.g., a learning period, full intervention effect period and followup period could be coded as 1, 2 and 3, respectively). Alternatively, level 2 might denote the introduction of a second intervention that replaces the first, or is added to the first. A separate parameter is assigned for each intervention level so power can be computed for comparisons between any level and the control, or between different levels. 

Note that the ETI model of Kenny et al (2022) implies a different intervention level for each exposure lag time. However, if you plan to use an ETI-based estimator of the treatment effect, then you should specify a standard 2-value design (0/1) in \code{swDsn} and use the H argument in \code{swPwr} or \code{swGlmPwr} to determine power for the estimator of interest. Also, note that NA periods that are interspersed with intervention periods (e.g. 0, 1, NA, 1) may cause unpredictable results for power calculations for ETI-based estimators since the lag for periods after the NA is ill-defined; NA periods that occur after the last intervention lag period of interest (e.g. 0, 1, 1, NA) are okay, however.}	
  \item{extra.time}{
  integer (scalar): deprecated; if used, replaces extra.trt.time.
  }
}
\details{
The \code{clusters} argument is required. After that, a design may be specified using the \code{tx.effect.frac}, \code{extra.ctrl.time}, \code{extra.trt.time}, \code{all.ctl.time0} arguments or the \code{swBlk} argument:

swDsn(clusters, tx.effect.frac=1, extra.ctl.time=0, extra.trt.time=0, all.ctl.time0=TRUE)
 
 or 
 
swDsn(clusters, swBlk)
}
\value{
numeric (list): Returns the following user-specified and computed objects
%% *****SPECIFY WHAT THESE ARE!!!!!**********
	\item{swDsn}{
		numeric (matrix):  schematic representation of the specified SW CRT design. Number of clusters is equal to the number of rows of the matrix and total number of time intervals is equal to the number of columns of the matrix. \code{swDsn[i,j]} gives the intervention status for cluster \code{i} at time \code{j} and has possible values NA (no data), 0 (control), or >=1 (intervention levels).
	}
	\item{swDsn.unique.clusters}{
		numeric (matrix): Truncated SW CRT design of interest, with one row for each wave.
	}
	\item{n.waves}{
		numeric (scalar): Number of waves for the SW CRT design of interest.
	}
	\item{clusters}{
		numeric (vector): Number of clusters for each wave for the SW CRT design of interest.
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
	\item{nTxLev}{
		numeric (integer): Number of treatment levels. 
	}
	\item{TxLev}{
		numeric (integer): Vector of unique treatment levels. 
	}
}
\references{
Hussey MA, Hughes JP. Design and analysis of stepped wedge cluster randomized trials. Contemporary Clinical Trials 2007;28:182-191.

Voldal EC, Hakhu NR, Xia, F, Heagerty PJ, Hughes JP. swCRTdesign: An R Package for Stepped Wedge Trial Design and Analysis. Computer Methods and Programs in Biomedicine 2020;196:105514.
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
# Example 1 (Equal clusters per wave, standard SW design)
swDsn.Ex1.std <- swDsn( clusters=c(3,3,3) )
swDsn.Ex1.std$swDsn

# Example 2 (Equal clusters per wave, extended SW design)
swDsn.Ex1.extend <- swDsn( clusters=c(3,3,3), extra.trt.time=2 )
swDsn.Ex1.extend$swDsn

# Example 3 (Equal clusters per wave, not all ctl at time 0, "standard" for time SW design)
swDsn.Ex1.std.noAllctl <- swDsn( clusters=c(3,3,3), all.ctl.time0=FALSE )
swDsn.Ex1.std.noAllctl$swDsn

# Example 4 (Equal clusters per wave, not all ctl at time 0, extended SW design)
swDsn.Ex1.extend.noAllctl <- swDsn( clusters=c(3,3,3), extra.trt.time=2, all.ctl.time0=FALSE )
swDsn.Ex1.extend.noAllctl$swDsn

# Example 5 (Unequal clusters per wave, standard SW design)
swDsn.Ex1.std.unequal <- swDsn( clusters=c(3,0,2) )
swDsn.Ex1.std.unequal$swDsn

# Example 6 (Unequal clusters per wave, extended SW design)
swDsn.Ex1.extend.unequal <- swDsn( clusters=c(3,0,2), extra.trt.time=2 )
swDsn.Ex1.extend.unequal$swDsn

# Example 7 (Unequal clusters per wave, extended SW design)
swDsn.Ex1.extend.unequal.varyTxEffect <- swDsn( clusters=c(3,0,2), tx.effect.frac=c(0.8,0.9,1.0),
extra.trt.time=2 )
swDsn.Ex1.extend.unequal.varyTxEffect$swDsn

#Example 8 (Equal clusters per wave, extra control time and extra treatment time)
swDsn.Ex1.extendctrl.extendtrt <- swDsn(c(5,5,5,5,5),extra.ctl.time=3,extra.trt.time=5)
swDsn.Ex1.extendctrl.extendtrt$swDsn

# Examples with swBlk
# Example 9 (Equal clusters per wave, standard SW design)
mat9 = matrix(c(0,1,1,1,1,0,0,1,1,1,0,0,0,1,1,0,0,0,0,1),4,5,byrow=TRUE)
swDsn.Ex9 <- swDsn(c(6,6,6,6),swBlk=mat9)
swDsn.Ex9$swDsn

# Example 10 (Unequal clusters per wave, periods with no data collection)
mat10 = matrix(c(0,1,1,1,1,NA,0,1,1,1,NA,NA,0,1,1,NA,NA,NA,0,1),4,5,byrow=TRUE)
swDsn.Ex10 <- swDsn(c(5,6,6,5),swBlk=mat10)
swDsn.Ex10$swDsn

# Example 11 (Unequal clusters per wave, periods with no data collection, 
# multiple intervention levels)
mat11 = matrix(c(0,1,2,2,2,2,NA,0,1,2,2,2,NA,NA,0,1,2,2,NA,NA,NA,0,1,2),4,6,byrow=TRUE)
swDsn.Ex11 <- swDsn(c(5,6,6,5),swBlk=mat11)
swDsn.Ex11$swDsn

# Example 12 (Dogleg design)
swDsn.Ex12 <- swDsn(c(4,4,4), swBlk=matrix(c(1,NA,0,1,NA,0),3,2,byrow=TRUE))
swDsn.Ex12$swDsn
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ design }
\keyword{ stepped wedge }
\keyword{ cluster randomized trial }% __ONLY ONE__ keyword per line

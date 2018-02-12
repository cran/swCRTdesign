\name{blkDiag}
\alias{blkDiag}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Block Diagonal Matrix Generator
}
\description{
\code{blkDiag} returns block diagonal matrix based on specified square blocks (either as an array or a list).
}
\usage{
blkDiag(z)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{z}{
	numeric (array or list): User-specified matrices to be combined into one block diagonal matrix.
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
numeric (matrix): \code{blkDiag} gives a block diagonal matrix.
}
%\references{
%% ~put references to the literature/web site here ~
%}
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
# Example 1 (input: array)
blkDiag.Ex1.array <- blkDiag( z=array(1:12, c(2,2,3)) )
blkDiag.Ex1.array

# Example 2 (input: list)
blkDiag.Ex2.list <- blkDiag( z=list(matrix(1,2,2), diag(2,3), matrix(3,4,4)) )
blkDiag.Ex2.list

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ block diagonal matrix }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

swSim <-
function(design, family, n, mu0, mu1, time.effect, sigma, tau, eta, rho=0, time.lab=NULL, seed=NULL, retTimeOnTx=FALSE)
{ 
##
## clusters
## extra.time
## design = stepped wedge design, as from swDsn, including, at least, swDsn, clusters, n.clusters, total.time
## distn = distribution family 
## link = scale on which you have linear predictor
## n = number of observations for each cluster at each time point (scalar, vector--each element is for each cluster, or matrix)
## mu0 = mean in control group
## mu1 = mean in intervention/experimental group
## time.effect = vector--length (total time points - 1); specify all time effects after first time point
## sigma
## tau
## eta
## rho
## time.lab
## seed
##
####
## taken from beginning of glmer() in lme4 package
mc <- match.call()
if( is.character(family) ) {
family <- get(family, mode="function", envir=parent.frame(2))
}
if( is.function(family) ){
family <- family()
}
if( !is.list(family) || is.null(family$family) ){
stop(gettextf("family '%s' not recognized", deparse(substitute(family)), domain="R-swCRTdesign"))
}
distn <- family$family
link <- family$link
####
p0 <- mu0
p1 <- mu1
#**#
## treatment effect [SCALAR]
theta <- p1 - p0
## standard deviation of random intercepts [SCALAR]
CV.p0 <- tau / p0
## standard deviation of random slopes [SCALAR]
CV.theta <- eta / abs(theta)
## SW CRT design [MATRIX]
X.ij <- design$swDsn
##********************************************************************************************************************
if( length(time.effect)==1 ) {
time.effectVec <- rep( time.effect, design$total.time )
}else if( length(time.effect) == design$total.time ) {
time.effectVec <- time.effect
}else {
warning("Invalid time.effects length. Either specify a scalar fixed time effect (i.e., the same fixed time effect at each time point), or specify a vector of fixed time effects for each time point. It is best to ignore any results which follow as a result of this warning message, and correctly assign value(s) to the 'time.effect' function argument of swSim().")
}
## time effect [MATRIX] 
## (for all clusters, all time points, but only 1 observation per (i,j)-pair)
beta.ij <- matrix( rep(time.effectVec, design$n.clusters), design$n.clusters, design$total.time, byrow=TRUE )
##********************************************************************************************************************
##********************************************************************************************************************
##********************************************************************************************************************
## treatment effect [MATRIX] 
## (for all clusters, all time points, but only 1 observation per (i,j)-pair)
thetaX.ij <- X.ij * theta
#
timeOnTx.ij <- swDsn(clusters=design$clusters, tx.effect=(1:design$total.time))$swDsn
## covariance matrix for random intercepts and random slopes [MATRIX]
sigMat <- matrix( c(tau^2, rho*tau*eta, rho*tau*eta, eta^2), 2, 2 )
## set random seed
set.seed(seed)
if(tau!=0 & eta!=0) {
## generate 2 standard normals for all clusters [VECTOR]
z <- rnorm( 2 * design$n.clusters )
## standard normals [MATRIX] (row=(a.i, r.i), col=(each cluster))
zMat <- matrix(z, nrow=2, byrow=FALSE)
## linear transformation of 
## randomly generated univariate standard normals with Cholesky decomposition
## to get randomly generated bivariate normals [MATRIX]
ar.i <- chol( sigMat ) %*% zMat
## random normally generated random intercepts [VECTOR]
a.i <- rep( ar.i[1,], design$total.time )
## random normally generated random slopes [VECTOR]
r.i <- rep( ar.i[2,], design$total.time )
}else if(tau!=0 & eta==0) {
a.i <- rep( rnorm( design$n.clusters, 0, tau ), design$total.time )
r.i <- 0
}else if(tau==0 & eta!=0) {
a.i <- 0
r.i <- rep( rnorm( design$n.clusters, 0, eta ), design$total.time )
}else if(tau==0 & eta==0) {
a.i <- 0
r.i <- 0
}
## random intercepts [MATRIX]
## (for all clusters, all time points, but only 1 observation per (i,j)-pair)
a.ij <- matrix(a.i, nrow=design$n.clusters, ncol=design$total.time, byrow=FALSE )
## random slopes [MATRIX]
## (for all clusters, all time points, but only 1 observation per (i,j)-pair)
rX.ij <- X.ij * matrix(r.i, nrow=design$n.clusters, ncol=design$total.time, byrow=FALSE )
## link(mu.ijk) [MATRIX]
## (for all clusters, all time points, but only 1 observation per (i,j)-pair)
link.mu.ij <- mu0 + beta.ij + thetaX.ij + a.ij + rX.ij
##
## expit [FUNCTION]
expit <- function(x){ 1 / (1 + exp(-x)) }
## logit [FUNCTION] (not needed)
## logit <- function(x){ log(1 / (1/x - 1)) }
##
## time [MATRIX]
## (for all clusters, all time points, but only 1 observation per (i,j)-pair)
##
## time.between provides the spacing between time points
##
if( is.null(time.lab) ){
time.lab <- 1:design$total.time
time.ij <- matrix( rep(c(time.lab),design$n.clusters), design$n.clusters, design$total.time, byrow=TRUE )
}else if( length(time.lab)==design$total.time ) {
time.ij <- matrix( rep(c(time.lab),design$n.clusters), design$n.clusters, design$total.time, byrow=TRUE )
}else {
warning("The length of the specified 'time.lab' vector is not equal to the total time points determined based on the SW design from specifying 'cluters', 'extra.time', and 'all.ctl.time0'. Ignore any results that follow and re-specify the vector 'time.lab' if you want to specify the label at each time point; specify 'NULL' to have the default labeling 1, 2, ... to the total number of time points.")
}
## cluster [MATRIX]
## (for all clusters, all time points, but only 1 observation per (i,j)-pair)
cluster.ij <- matrix( rep(c(1:design$n.clusters),each=design$total.time), design$n.clusters, design$total.time, byrow=TRUE )
##**********************************************************************************************************************
## creating response, tx, time, and cluster variables
## for specified observations 'n'
##   - scalar (same 'n' for each cluster AND each time point)
##   + vector (each element of vector is unique 'n[k]' for all time points within each cluster)
##   + matrix (each element of matrix is unique 'n[k,l]' for each cluster and each time points)
if( is.vector(n) & length(n)==1 ){
##
## link(mu.ijk) [VECTOR]
## (for all clusters, all time points, AND all observations)
tmplink.mu.ijk <- rep( as.vector( t(link.mu.ij) ), each=n )
##
tmpTx <- rep( as.vector( t(X.ij) ), each=n )
tmpTime <- rep( as.vector( t(time.ij) ), each=n )
#
#
tmpTimeOnTx <- rep( as.vector( t(timeOnTx.ij) ), each=n )
#
#
tmpCluster <- rep( as.vector( t(cluster.ij) ), each=n )
##
## generate Y.ijk from specified distn [VECTOR]
if(distn=="gaussian") {
## use inverse-link function [VECTOR]
if(link=="identity") {
mu.ijk <- tmplink.mu.ijk
}else if(link=="logit") {
mu.ijk <- expit( tmplink.mu.ijk )
## count number of mu.ijk < 0 & mu.ijk > 1
mu.ijk.lessthan0 <- length( mu.ijk[mu.ijk < 0] )
mu.ijk.greaterthan1 <- length( mu.ijk[mu.ijk > 1] )
## change "< 0" to 0 and "> 1" to 1
mu.ijk[mu.ijk < 0] <- 0
mu.ijk[mu.ijk > 1] <- 1
if( mu.ijk.lessthan0 > 0 | mu.ijk.greaterthan1 > 0 ) warning( paste("When using a logit link function with family=gaussian, the probability of outcome for person k in cluster i at time j, denoted by mu.ijk in swSim(), cannot be less than 0 or greater than 1. Out of the ", length(mu.ijk), " probabilities computed, ", mu.ijk.lessthan0, " (", mu.ijk.lessthan0/length(mu.ijk)*100, "%) were less than 0 [ACTION TAKEN BY swSim: forced these to be 0] and ", mu.ijk.greaterthan1, " (", mu.ijk.greaterthan1/length(mu.ijk)*100, "%) were greater than 1 [ACTION TAKEN BY swSim: forced these to be 1].", sep="") )
}else if(link=="log"){
mu.ijk <- exp( tmplink.mu.ijk )
## count number of mu.ijk < 0
mu.ijk.lessthan0 <- length( mu.ijk[mu.ijk < 0] )
## change "< 0" to 0
mu.ijk[mu.ijk < 0] <- 0
if( mu.ijk.lessthan0 > 0 ) warning( paste("When using a log link function with family=gaussian, the mean of outcome for person k in cluster i at time j, denoted by mu.ijk in swSim(), cannot be less than 0. Out of the ", length(mu.ijk), " means computed, ", mu.ijk.lessthan0, " (", mu.ijk.lessthan0/length(mu.ijk)*100, "%) were less than 0 [ACTION TAKEN BY swSim: forced these to be 0].", sep="") )
}
## generating response variable (Y.ijk)
Y.ijk <- rnorm( n *design$n.clusters*design$total.time, mu.ijk, sd=sigma )
##
##
}else if(distn=="binomial") {
## use inverse-link function [VECTOR]
if(link=="identity") {
mu.ijk <- tmplink.mu.ijk
## count number of mu.ijk < 0 & mu.ijk > 1
mu.ijk.lessthan0 <- length( mu.ijk[mu.ijk < 0] )
mu.ijk.greaterthan1 <- length( mu.ijk[mu.ijk > 1] )
## change "< 0" to 0 and "> 1" to 1
mu.ijk[mu.ijk < 0] <- 0
mu.ijk[mu.ijk > 1] <- 1
if( mu.ijk.lessthan0 > 0 | mu.ijk.greaterthan1 > 0 ) warning( paste("When using an identity link function with family=binomial, the probability of outcome for person k in cluster i at time j, denoted by mu.ijk in swSim(), cannot be less than 0 or greater than 1. Out of the ", length(mu.ijk), " probabilities computed, ", mu.ijk.lessthan0, " (", mu.ijk.lessthan0/length(mu.ijk)*100, "%) were less than 0 [ACTION TAKEN BY swSim: forced these to be 0] and ", mu.ijk.greaterthan1, " (", mu.ijk.greaterthan1/length(mu.ijk)*100, "%) were greater than 1 [ACTION TAKEN BY swSim: forced these to be 1].", sep="") )
}else if(link=="logit") {
mu.ijk <- expit( tmplink.mu.ijk )
}else if(link=="log"){
mu.ijk <- exp( tmplink.mu.ijk )
#**#
## count number of mu.ijk > 0
mu.ijk.greaterthan1 <- length( mu.ijk[mu.ijk < 0] )
## change "> 1" to 1
mu.ijk[mu.ijk > 1] <- 1
if( mu.ijk.greaterthan1 > 1 ) warning( paste("When using a log link function with family=binomial, the mean of outcome for person k in cluster i at time j, denoted by mu.ijk in swSim(), cannot be greater than 1. Out of the ", length(mu.ijk), " means computed, ", mu.ijk.greaterthan1, " (", mu.ijk.greaterthan1/length(mu.ijk)*100, "%) were greater than 1 [ACTION TAKEN BY swSim: forced these to be 1].", sep="") )
#**#
}
## generating response variable (Y.ijk)
Y.ijk <- rbinom( n *design$n.clusters*design$total.time, 1, mu.ijk )
##
##
}else if(distn=="poisson") {
## use inverse-link function [VECTOR]
if(link=="identity") {
mu.ijk <- tmplink.mu.ijk
## count number of mu.ijk < 0 & mu.ijk > 1
mu.ijk.lessthan0 <- length( mu.ijk[mu.ijk < 0] )
mu.ijk.greaterthan1 <- length( mu.ijk[mu.ijk > 1] )
## change "< 0" to 0 and "> 1" to 1
mu.ijk[mu.ijk < 0] <- 0
mu.ijk[mu.ijk > 1] <- 1
if( mu.ijk.lessthan0 > 0 | mu.ijk.greaterthan1 > 0 ) warning( paste("When using an identity link function with family=poisson, the probability of outcome for person k in cluster i at time j, denoted by mu.ijk in swSim(), cannot be less than 0 or greater than 1. Out of the ", length(mu.ijk), " probabilities computed, ", mu.ijk.lessthan0, " (", mu.ijk.lessthan0/length(mu.ijk)*100, "%) were less than 0 [ACTION TAKEN BY swSim: forced these to be 0] and ", mu.ijk.greaterthan1, " (", mu.ijk.greaterthan1/length(mu.ijk)*100, "%) were greater than 1 [ACTION TAKEN BY swSim: forced these to be 1].", sep="") )
}else if(link=="logit") {
mu.ijk <- expit( tmplink.mu.ijk )
## count number of mu.ijk < 0 & mu.ijk > 1
mu.ijk.lessthan0 <- length( mu.ijk[mu.ijk < 0] )
mu.ijk.greaterthan1 <- length( mu.ijk[mu.ijk > 1] )
## change "< 0" to 0 and "> 1" to 1
mu.ijk[mu.ijk < 0] <- 0
mu.ijk[mu.ijk > 1] <- 1
if( mu.ijk.lessthan0 > 0 | mu.ijk.greaterthan1 > 0 ) warning( paste("When using a logit link function with family=poisson, the probability of outcome for person k in cluster i at time j, denoted by mu.ijk in swSim(), cannot be less than 0 or greater than 1. Out of the ", length(mu.ijk), " probabilities computed, ", mu.ijk.lessthan0, " (", mu.ijk.lessthan0/length(mu.ijk)*100, "%) were less than 0 [ACTION TAKEN BY swSim: forced these to be 0] and ", mu.ijk.greaterthan1, " (", mu.ijk.greaterthan1/length(mu.ijk)*100, "%) were greater than 1 [ACTION TAKEN BY swSim: forced these to be 1].", sep="") )
}else if(link=="log"){
mu.ijk <- exp( tmplink.mu.ijk )
}
## generating response variable (Y.ijk)
Y.ijk <- rpois( n *design$n.clusters*design$total.time, mu.ijk )
##
##
}
}else if( (is.vector(n) & length(n) > 1) | is.matrix(n) ) { 
## create nMat [MATRIX]
if( (is.vector(n) & length(n) > 1)  ) {
nMat <- matrix( rep(n, each=design$total.time), ncol=design$total.time, byrow=TRUE )  ## nrow=design$n.clusters
}else if( is.matrix(n) ){
nMat <- n
}
## convert nMat to nMat2Vec [VECTOR]
nMat2Vec <- as.vector( t(nMat) )
##
tmplink.mu.ijk <- NULL
tmpTx <- NULL
tmpTime <- NULL
#
tmpTimeOnTx <- NULL
#
#
tmpCluster <- NULL
for( k in 1:length(nMat2Vec) ){
## link(mu.ijk) [VECTOR]## (for all clusters, all time points, AND all observations)
tmplink.mu.ijk <- c( tmplink.mu.ijk, rep( as.vector( t(link.mu.ij) )[k], each= nMat2Vec[k] ) )
##
tmpTx <- c( tmpTx, rep( as.vector( t(X.ij) )[k], each=nMat2Vec[k] ) )
tmpTime <- c( tmpTime, rep( as.vector( t(time.ij) )[k], each=nMat2Vec[k] ) )
#
#
tmpTimeOnTx <- c( tmpTimeOnTx, rep( as.vector( t(timeOnTx.ij) )[k], each=nMat2Vec[k] ) )
#
#
tmpCluster <- c( tmpCluster, rep( as.vector( t(cluster.ij) )[k], each=nMat2Vec[k] ) )
}
##
## generate Y.ijk from specified distn [VECTOR]
if(distn=="gaussian") {
## use inverse-link function [VECTOR]
if(link=="identity") {
mu.ijk <- tmplink.mu.ijk
}else if(link=="logit") {
mu.ijk <- expit( tmplink.mu.ijk )
## count number of mu.ijk < 0 & mu.ijk > 1
mu.ijk.lessthan0 <- length( mu.ijk[mu.ijk < 0] )
mu.ijk.greaterthan1 <- length( mu.ijk[mu.ijk > 1] )
## change "< 0" to 0 and "> 1" to 1
mu.ijk[mu.ijk < 0] <- 0
mu.ijk[mu.ijk > 1] <- 1
if( mu.ijk.lessthan0 > 0 | mu.ijk.greaterthan1 > 0 ) warning( paste("When using a logit link function with family=gaussian, the probability of outcome for person k in cluster i at time j, denoted by mu.ijk in swSim(), cannot be less than 0 or greater than 1. Out of the ", length(mu.ijk), " probabilities computed, ", mu.ijk.lessthan0, " (", mu.ijk.lessthan0/length(mu.ijk)*100, "%) were less than 0 [ACTION TAKEN BY swSim: forced these to be 0] and ", mu.ijk.greaterthan1, " (", mu.ijk.greaterthan1/length(mu.ijk)*100, "%) were greater than 1 [ACTION TAKEN BY swSim: forced these to be 1].", sep="") )
}else if(link=="log"){
mu.ijk <- exp( tmplink.mu.ijk )
## count number of mu.ijk < 0
mu.ijk.lessthan0 <- length( mu.ijk[mu.ijk < 0] )
## change "< 0" to 0
mu.ijk[mu.ijk < 0] <- 0
if( mu.ijk.lessthan0 > 0 ) warning( paste("When using a log link function with family=gaussian, the mean of outcome for person k in cluster i at time j, denoted by mu.ijk in swSim(), cannot be less than 0. Out of the ", length(mu.ijk), " means computed, ", mu.ijk.lessthan0, " (", mu.ijk.lessthan0/length(mu.ijk)*100, "%) were less than 0 [ACTION TAKEN BY swSim: forced these to be 0].", sep="") )
}
##
Y.ijk <- rnorm( n *design$n.clusters*design$total.time, mu.ijk, sd=sigma )
##
##
}else if(distn=="binomial") {
## use inverse-link function [VECTOR]
if(link=="identity") {
mu.ijk <- tmplink.mu.ijk
## count number of mu.ijk < 0 & mu.ijk > 1
mu.ijk.lessthan0 <- length( mu.ijk[mu.ijk < 0] )
mu.ijk.greaterthan1 <- length( mu.ijk[mu.ijk > 1] )
## change "< 0" to 0 and "> 1" to 1
mu.ijk[mu.ijk < 0] <- 0
mu.ijk[mu.ijk > 1] <- 1
if( mu.ijk.lessthan0 > 0 | mu.ijk.greaterthan1 > 0 ) warning( paste("When using an identity link function with family=binomial, the probability of outcome for person k in cluster i at time j, denoted by mu.ijk in swSim(), cannot be less than 0 or greater than 1. Out of the ", length(mu.ijk), " probabilities computed, ", mu.ijk.lessthan0, " (", mu.ijk.lessthan0/length(mu.ijk)*100, "%) were less than 0 [ACTION TAKEN BY swSim: forced these to be 0] and ", mu.ijk.greaterthan1, " (", mu.ijk.greaterthan1/length(mu.ijk)*100, "%) were greater than 1 [ACTION TAKEN BY swSim: forced these to be 1].", sep="") )
}else if(link=="logit") {
mu.ijk <- expit( tmplink.mu.ijk )
}else if(link=="log"){
mu.ijk <- exp( tmplink.mu.ijk )
#**#
## count number of mu.ijk > 0
mu.ijk.greaterthan1 <- length( mu.ijk[mu.ijk < 0] )
## change "> 1" to 1
mu.ijk[mu.ijk > 1] <- 1
if( mu.ijk.greaterthan1 > 1 ) warning( paste("When using a log link function with family=binomial, the mean of outcome for person k in cluster i at time j, denoted by mu.ijk in swSim(), cannot be greater than 1. Out of the ", length(mu.ijk), " means computed, ", mu.ijk.greaterthan1, " (", mu.ijk.greaterthan1/length(mu.ijk)*100, "%) were greater than 1 [ACTION TAKEN BY swSim: forced these to be 1].", sep="") )
#**#
}
## generating response variable (Y.ijk)
Y.ijk <- rbinom( n *design$n.clusters*design$total.time, 1, mu.ijk )
##
##
}else if(distn=="poisson") {
## use inverse-link function [VECTOR]
if(link=="identity") {
mu.ijk <- tmplink.mu.ijk
## count number of mu.ijk < 0 & mu.ijk > 1
mu.ijk.lessthan0 <- length( mu.ijk[mu.ijk < 0] )
mu.ijk.greaterthan1 <- length( mu.ijk[mu.ijk > 1] )
## change "< 0" to 0 and "> 1" to 1
mu.ijk[mu.ijk < 0] <- 0
mu.ijk[mu.ijk > 1] <- 1
if( mu.ijk.lessthan0 > 0 | mu.ijk.greaterthan1 > 0 ) warning( paste("When using an identity link function with family=poisson, the probability of outcome for person k in cluster i at time j, denoted by mu.ijk in swSim(), cannot be less than 0 or greater than 1. Out of the ", length(mu.ijk), " probabilities computed, ", mu.ijk.lessthan0, " (", mu.ijk.lessthan0/length(mu.ijk)*100, "%) were less than 0 [ACTION TAKEN BY swSim: forced these to be 0] and ", mu.ijk.greaterthan1, " (", mu.ijk.greaterthan1/length(mu.ijk)*100, "%) were greater than 1 [ACTION TAKEN BY swSim: forced these to be 1].", sep="") )
}else if(link=="logit") {
mu.ijk <- expit( tmplink.mu.ijk )
## count number of mu.ijk < 0 & mu.ijk > 1
mu.ijk.lessthan0 <- length( mu.ijk[mu.ijk < 0] )
mu.ijk.greaterthan1 <- length( mu.ijk[mu.ijk > 1] )
## change "< 0" to 0 and "> 1" to 1
mu.ijk[mu.ijk < 0] <- 0
mu.ijk[mu.ijk > 1] <- 1
if( mu.ijk.lessthan0 > 0 | mu.ijk.greaterthan1 > 0 ) warning( paste("When using a logit link function with family=poisson, the probability of outcome for person k in cluster i at time j, denoted by mu.ijk in swSim(), cannot be less than 0 or greater than 1. Out of the ", length(mu.ijk), " probabilities computed, ", mu.ijk.lessthan0, " (", mu.ijk.lessthan0/length(mu.ijk)*100, "%) were less than 0 [ACTION TAKEN BY swSim: forced these to be 0] and ", mu.ijk.greaterthan1, " (", mu.ijk.greaterthan1/length(mu.ijk)*100, "%) were greater than 1 [ACTION TAKEN BY swSim: forced these to be 1].", sep="") )
}else if(link=="log"){
mu.ijk <- exp( tmplink.mu.ijk )
}
## generating response variable (Y.ijk)
Y.ijk <- rpois( n *design$n.clusters*design$total.time, mu.ijk )
##
##
}
##
}
##**********************************************************************************************************************
## response variable [VECTOR] (for each cluster, each time point, and each observation)
response.var <- Y.ijk
## treatment variable [VECTOR] (for each cluster, each time point, and each observation)
tx.var <- tmpTx
#
#

## time-on-treatment variable [VECTOR] (for each cluster, each time point, and each observation)
#**#
tx.var <- ifelse( tmpTx > 0, 1, 0 )

#**#
timeOnTx.var <- tmpTimeOnTx
#
## time variable [VECTOR] (for each cluster, each time point, and each observation)
time.var <- tmpTime
## cluster variable [VECTOR] (for each cluster, each time point, and each observation)
cluster.var <- tmpCluster
## swData [DATAFRAME] (result of function)
##
##
if( !retTimeOnTx ) {
swData <- data.frame(response.var, tx.var, time.var, cluster.var)
}else if( retTimeOnTx ) {
swData <- data.frame(response.var, tx.var, timeOnTx.var, time.var, cluster.var)
}
##
## returning swData [DATAFRAME]
swData
}

swSim <- function (design, family, log.gaussian=FALSE, n, mu0, mu1, time.effect, sigma, tau=0,
                         eta=0, rho=0, gamma=0, zeta=0, icc=0, cac=1, iac=0, time.lab = NULL, 
                         retTimeOnTx=TRUE, silent = FALSE, nocheck=FALSE)
{
  #
  #Last update: 4/28/23, v. 4.0, Jim Hughes
  # Note: timeOnTx always returned (argument retTimeOnTx removed)
  #In this file, 'random slopes' and 'random treatment effects' are used interchangeably
  ##########
  # Helper functions
  ##########
  replaceNAwithZero = function(x){
    ifelse(is.na(x),0,x)
  }
  expit <- function(x) {
    1/(1 + exp(-x))
  }
  ##########
  #Warnings, unpack family, translate ICC/CAC
  ##########
  if (!missing(retTimeOnTx)) warning("argument retTimeOnTx always set to TRUE")
  ## taken from beginning of glmer() in lme4 package
  mc <- match.call()
  if (is.character(family)) {
    family <- get(family, mode = "function", envir = parent.frame(2))
  }
  if (is.function(family)) {
    family <- family()
  }
  if (!is.list(family) || is.null(family$family)) {
    stop(gettextf("family '%s' not recognized", deparse(substitute(family)),
                  domain = "R-swCRTdesign"))
  }
  distn <- family$family
  if (!nocheck){
  #Basic input checks
  if(! all(n%%1 == 0)){
    stop("n (either scalar, vector, or matrix) must consist only of integers.")
  }
  if (design$nTxLev==1){
    maxET = max(as.vector(apply(replaceNAwithZero(design$swDsn),1,cumsum)))
    if (!(length(mu1)==1 | length(mu1)==maxET)){
      stop("Length of mu1 must be 1 or maximum number of exposure time periods")  
    }
    } else {
    if(design$nTxLev!=length(mu1)){
      stop("Length of mu1 must correspond to number of treatment levels specified in swDsn1")
    }
    }
  if (distn == 'gaussian' & missing(sigma)) stop("If distribution is gaussian, a value for sigma must be entered")
  #Checks to make sure people are using random effects OR ICC/CAC.
    param.icc.any <- !missing(icc) | !missing(cac) | !missing(iac)
    param.re.any <- !missing(tau) | !missing(eta) | !missing(rho) | !missing(gamma) | !missing(zeta)
    if(param.re.any == TRUE & param.icc.any == TRUE){
      stop("The two parameterizations (random effects and ICC/CAC/IAC) are mutually exclusive.  Either enter values for ICC, CAC (and IAC, if cohort design), or tau, eta, gamma, rho (and zeta, if cohort design).")
    }
    #If using ICC/CAC/IAC, translate to random effects
  if (param.icc.any == TRUE){
    #Check range restrictions
    if(icc < 0 | icc > 1 | cac < 0 | cac > 1 | iac < 0 | iac > 1){
      stop("ICC, CAC and IAC must be between 0 and 1.")
    }
    if(icc == 1){
      stop("There are multiple combinations of random effects that can make the ICC be 1; if you believe this is a realistic scenario, use the random effect parameterization.")
    }
    if(iac == 1){
      stop("There are multiple combinations of random effects that can make the IAC be 1; if you believe this is a realistic scenario, use the random effect parameterization.")
    }
    #Assume eta=0 and rho=0
    eta <- 0
    rho <- 0
    if (distn == 'gaussian'){
      sigmasq.temp <- sigma^2
      if(sigma == 0){
        stop("When sigma is 0, the random effects parameterization must be used.")
      }
    }
    if (distn == 'binomial'){
      stop("For a Binomial outcome, please use the random effects parameterization instead of ICC/CAC.")
    }
    if (distn == 'poisson'){
      stop("For a Poisson outcome, please use the random effects parameterization instead of ICC/CAC.")
    }
    if(cac == 1){
      zeta <- sqrt(sigmasq.temp*iac/(1-iac))
      gamma <- 0
      tau <- sqrt((sigmasq.temp+zeta^2)*icc/(1-icc))
    }else{
      zeta <- sqrt(sigmasq.temp*iac/(1-iac)) 
      gamma <- sqrt(icc*(sigmasq.temp+zeta^2)*(1-cac)/(1-icc))
      tau <- sqrt(gamma^2*cac/(1-cac))
    }
  }
  #Basic restrictions on newly defined variance components
  if(rho < -1 | rho > 1){
    stop("Rho must be a numeral between -1 and 1.")
  }
  if(tau < 0 | eta < 0 | gamma < 0 | zeta < 0){
    stop("Tau, eta, gamma and zeta must be non-negative numerals.")
  }
  if((tau == 0 | eta == 0) & rho != 0){
    stop("If either tau or eta is zero, rho must be zero, since fixed effects cannot be correlated.")
  }
  if(distn == 'gaussian' && sigma == 0 & gamma == 0 & zeta == 0){
    stop("For a non-deterministic Gaussian outcome, at least one of sigma, zeta or gamma needs to be non-zero.")
  }
  if (distn != "gaussian" & missing(sigma) == FALSE & silent==FALSE){
    warning("sigma is only used when family is gaussian")
  }
  if (distn != "gaussian" & log.gaussian == TRUE){
    stop("The log.gaussian=TRUE option must be paired with a Gaussian family.")
  }
  if (length(tau) > 1 | length(eta) > 1 | length(zeta) > 1) {
    stop("Function cannot simulate stepped wedge design data for tau-vector, zeta-vector or eta-vector; tau, zeta and eta must be scalars.")
  }
  if (any(rowSums(design$swDsn,na.rm=TRUE) == 0)) {
    warning("For the specified total number of clusters (I), total number of time periods (J), and number of cluster repetitions (I.rep), the specified stepped wedge design has at least one cluster which does not crossover from control(0) to treatment(1) arm.")
  }
  if (is.vector(n) & length(n) > 1 & length(n) != design$n.clusters){
    stop("The number of clusters in 'design' (design$n.clusters) and 'n' (length(n)) do not match.")
  }
  if (is.matrix(n) && ((nrow(n) != design$n.clusters)|(ncol(n) != design$total.time))){
    stop("The number of clusters and/or time steps in 'design' (design$n.clusters and design$total.time) and 'n' (number of rows and columns, respectively) do not match.")
  }
  if (length(n)>1 & silent == FALSE){
    warning("When sample sizes are not uniform, power depends on order of clusters (see documentation).")
  }
  if (zeta > 0){
    if (silent==FALSE) warning("Closed cohort design assumed")
    if (is.matrix(n)) {
      for (i in 1:nrow(n)){
        if (var(n[i,n[i,]>0])>0) stop("For closed cohort design (zeta > 0) sample size must be constant across time for each cluster (exception: sample size can be 0 to denote time periods which will not be included in the analysis)")
      }
    }
  } else {
    if (silent==FALSE) warning("Cross-sectional design assumed")
  }
  }
  ##########
  ##########
  link <- family$link
  p0 <- mu0
  p1 <- mu1
  theta <- p1 - p0## treatment effect 
#  CV.p0 <- tau/p0## standard deviation of random intercepts for log/logit scale[SCALAR]
#  CV.theta <- eta/abs(theta)## standard deviation of random slopes  for log/logit scale
  X.ij <- design$swDsn## SW CRT design [MATRIX]
  if (length(time.effect) == 1) {
    time.effectVec <- rep(time.effect, design$total.time)
  }
  else if (length(time.effect) == design$total.time) {
    time.effectVec <- time.effect
  } else {
    stop("Invalid time.effects length. Either specify a scalar fixed time effect (i.e., the same fixed time effect at each time point), or specify a vector of fixed time effects for each time point. It is best to ignore any results which follow as a result of this warning message, and correctly assign value(s) to the 'time.effect' function argument of swSim().")
  }
  ## time effect [MATRIX]
  ## (for all clusters, all time points, but only 1 observation per (i,j)-pair)
  beta.ij <- matrix(rep(time.effectVec, design$n.clusters),design$n.clusters, design$total.time, byrow = TRUE)
  #Make matrix of time on treatment (Note: time on any treatment, regardless of fractional treatment effect)
  timeOnTx.ij.pre <- replaceNAwithZero(design$swDsn)
  timeOnTx.ij <- t(apply(timeOnTx.ij.pre>0, 1, cumsum))
  ## treatment effect [MATRIX]
  if (design$nTxLev > 1){
# multiple treatments or multiple levels
    thetaX.ij = matrix(0,nrow(X.ij),ncol(X.ij))
    TxLev = design$TxLev
    for (k in 1:design$nTxLev){
      thetaX.ij = thetaX.ij + (X.ij==TxLev[k])*theta[k]
    }
  } else {
    if (length(theta)==1){
# immediate, constant treatment effect; length(theta) is 1
# But note that X.ij may incorporate fractional treatment effects
      thetaX.ij <- X.ij * theta 
    } else {
# time-varying treatment effect; length(theta) is maxET
      thetaX.ij = matrix(NA,nrow(X.ij),ncol(X.ij))
      for (i in 1:nrow(X.ij)){
        for (j in 1:ncol(X.ij)){
        thetaX.ij[i,j] <- ifelse(X.ij[i,j]==1,theta[timeOnTx.ij[i,j]],0)
      }}
    }
  }
  ## covariance matrix for random intercepts and random treatment effects [MATRIX] (can do gamma separately because assumed no correlation)
  sigMat <- matrix(c(tau^2, rho * tau * eta, rho * tau * eta, eta^2), 2, 2)
  ##########
  #Generate random effects at the coarsest level
  ##########
  if (tau != 0 & eta != 0) {
    ## generate 2 standard normals for all clusters [VECTOR]
    z <- rnorm(2 * design$n.clusters)
    ## standard normals [MATRIX] (row=(a.i, r.i), col=(each cluster))
    zMat <- matrix(z, nrow = 2, byrow = FALSE)
    ## linear transformation of
    ## randomly generated univariate standard normals with Cholesky decomposition
    ## to get randomly generated bivariate normals [MATRIX]
    ar.i <- chol(sigMat) %*% zMat
    ## random normally generated random intercepts [VECTOR]
    a.i <- rep(ar.i[1, ], design$total.time)
    ## random normally generated random slopes [VECTOR]
    r.i <- rep(ar.i[2, ], design$total.time)
  }
  else if (tau != 0 & eta == 0) {
    a.i <- rep(rnorm(design$n.clusters, 0, tau), design$total.time)
    r.i <- 0
  }
  else if (tau == 0 & eta != 0) {
    a.i <- 0
    r.i <- rep(rnorm(design$n.clusters, 0, eta), design$total.time)
  }
  else if (tau == 0 & eta == 0) {
    a.i <- 0
    r.i <- 0
  }
  g.i <- 0
  if (gamma != 0){#generate one g.i per time point per cluster (random time effects)
    g.i <- rnorm(design$n.clusters*design$total.time,0,gamma)
  }
  ## random intercepts [MATRIX]
  ## (for all clusters, all time points, but only 1 observation per (i,j)-pair)
  a.ij <- matrix(a.i, nrow = design$n.clusters, ncol = design$total.time, byrow = FALSE)
  ## random slopes [MATRIX]
  ## (for all clusters, all time points, but only 1 observation per (i,j)-pair)
  rX.ij <- as.integer(X.ij>0) * matrix(r.i, nrow = design$n.clusters, ncol = design$total.time, byrow = FALSE)
  ## random time effects [MATRIX]
  g.ij <- matrix(g.i,nrow=design$n.clusters, ncol=design$total.time, byrow=FALSE)
  ##########
  #Combine fixed and random effects on link scale
  ##########
  ## link(mu.ijk) [MATRIX]
  ## (for all clusters, all time points, but only 1 observation per (i,j)-pair)
  link.mu.ij <- mu0 + beta.ij + thetaX.ij + a.ij + rX.ij + g.ij
  ## time [MATRIX]
  ## (for all clusters, all time points, but only 1 observation per (i,j)-pair)
  ## time.between provides the spacing between time points
  if (is.null(time.lab)) {
    time.lab <- 1:design$total.time
    time.ij <- matrix(rep(c(time.lab), design$n.clusters), design$n.clusters, design$total.time, byrow = TRUE)
  }
  else if (length(time.lab) == design$total.time) {
    time.ij <- matrix(rep(c(time.lab), design$n.clusters), design$n.clusters, design$total.time, byrow = TRUE)
  }
  else {
    stop("The length of the specified 'time.lab' vector is not equal to the total time points determined based on the SW design from specifying 'cluters', 'extra.time', and 'all.ctl.time0'. Ignore any results that follow and re-specify the vector 'time.lab' if you want to specify the label at each time point; specify 'NULL' to have the default labeling 1, 2, ... to the total number of time points.")
  }
  ## cluster [MATRIX]
  ## (for all clusters, all time points, but only 1 observation per (i,j)-pair)
  cluster.ij <- matrix(rep(c(1:design$n.clusters), each = design$total.time), design$n.clusters, design$total.time, byrow = TRUE)
  ##########
  ##Generate random data for each individual
  ##########
  #Note: some features were removed when updating to version 3, and the potential combinations of links and families were reduced.  The code is left here from the previous version, commented out, in case we want to add features back in or users want to see what old code was using.
  #
  ## creating response, tx, time, and cluster variables
  ## for specified observations 'n'
  ##   - scalar (same 'n' for each cluster AND each time point)
  ##   + vector (each element of vector is unique 'n[k]' for all time points within each cluster)
  ##   + matrix (each element of matrix is unique 'n[k,l]' for each cluster and each time points)
  #
  #####n a scalar
  if (is.vector(n) & length(n) == 1) {
    ## link(mu.ijk) [VECTOR]
    ## (for all clusters, all time points, AND all observations)
    tmplink.mu.ijk <- rep(as.vector(t(link.mu.ij)), each = n)
    if (zeta>0) {
      addto.mu.ijk <- NULL
      tmpIndivid <- NULL
      base <- 0
      for (i in 1:design$n.clusters) {
        addto.mu.ijk <- c(addto.mu.ijk,rep(rnorm(n,0,zeta),design$total.time))
        tmpIndivid <- c(tmpIndivid,rep(base+(1:n),design$total.time))
        base <- base+n
      }
      tmplink.mu.ijk <- tmplink.mu.ijk + addto.mu.ijk
    }
    miss = is.na(tmplink.mu.ijk)
    tmplink.mu.ijk = tmplink.mu.ijk[!miss]
    tmpTx <- rep(as.vector(t(X.ij)), each = n)[!miss]
    tmpTime <- rep(as.vector(t(time.ij)), each = n)[!miss]
    tmpTimeOnTx <- rep(as.vector(t(timeOnTx.ij)), each = n)[!miss]
    tmpCluster <- rep(as.vector(t(cluster.ij)), each = n)[!miss]
    if (zeta>0) tmpIndivid <- tmpIndivid[!miss]
    ## generate Y.ijk from specified distn [VECTOR]
    if (distn == "gaussian") {
      ## use inverse-link function [VECTOR]
      if (link == "identity") {
        mu.ijk <- tmplink.mu.ijk
      }
      else if (link == "logit") {
        stop("This combination of distribution and link is not supported.  See documentation.")
      }
      else if (link == "log") {
        stop("This combination of distribution and link is not supported.  Instead, use identity link and set log.gaussian=TRUE for lognormal data. See documentation.")
      }
      ## generating response variable (Y.ijk)
      Y.ijk <- rnorm(length(mu.ijk), mu.ijk, sd = sigma)
      if(log.gaussian == TRUE){
        Y.ijk <- exp(Y.ijk)
      }
    }
    else if (distn == "binomial") {
      ## use inverse-link function [VECTOR]
      if (link == "identity") {
        mu.ijk <- tmplink.mu.ijk
        mu.ijk.lessthan0 <- length(mu.ijk[mu.ijk < 0])
        mu.ijk.greaterthan1 <- length(mu.ijk[mu.ijk >1])
        mu.ijk[mu.ijk < 0] <- 0
        mu.ijk[mu.ijk > 1] <- 1
        if (mu.ijk.lessthan0 > 0 | mu.ijk.greaterthan1 >0)
          if (silent==FALSE) warning(paste("When using an identity link function with family=binomial, the probability of outcome for person k in cluster i at time j, denoted by mu.ijk in swSim(), cannot be less than 0 or greater than 1. Out of the ",
                        length(mu.ijk), " probabilities computed, ",
                        mu.ijk.lessthan0, " (", mu.ijk.lessthan0/length(mu.ijk) *
                          100, "%) were less than 0 [ACTION TAKEN BY swSim: forced these to be 0] and ",
                        mu.ijk.greaterthan1, " (", mu.ijk.greaterthan1/length(mu.ijk) *
                          100, "%) were greater than 1 [ACTION TAKEN BY swSim: forced these to be 1].",
                        sep = ""))
      }
      else if (link == "logit") {
        mu.ijk <- expit(tmplink.mu.ijk)
      }
      else if (link == "log") {
        mu.ijk <- exp(tmplink.mu.ijk)
        mu.ijk.greaterthan1 <- length(mu.ijk[mu.ijk >1])#Edit: fixed this from <0
        mu.ijk[mu.ijk > 1] <- 1
        if (mu.ijk.greaterthan1 > 1)
          if (silent==FALSE) warning(paste("When using a log link function with family=binomial, the mean of outcome for person k in cluster i at time j, denoted by mu.ijk in swSim(), cannot be greater than 1. Out of the ",
                        length(mu.ijk), " means computed, ", mu.ijk.greaterthan1,
                        " (", mu.ijk.greaterthan1/length(mu.ijk) *
                          100, "%) were greater than 1 [ACTION TAKEN BY swSim: forced these to be 1].",
                        sep = ""))
      }
      ## generating response variable (Y.ijk)
      Y.ijk <- rbinom(length(mu.ijk),1, mu.ijk)
      }
    else if (distn == "poisson") {
      ## use inverse-link function [VECTOR]
      if (link == "identity") {#Edit: removed restriction lambda<=1
        mu.ijk <- tmplink.mu.ijk
        mu.ijk.lessthan0 <- length(mu.ijk[mu.ijk < 0])
        mu.ijk[mu.ijk < 0] <- 0
        if (mu.ijk.lessthan0 > 0 )
          if (silent==FALSE) warning(paste("When using an identity link function with family=poisson, the expected count for person k in cluster i at time j, denoted by mu.ijk in swSim(), cannot be less than 0. Out of the ",
                        length(mu.ijk), " probabilities computed, ",
                        mu.ijk.lessthan0, " (", mu.ijk.lessthan0/length(mu.ijk) *
                          100, "%) were less than 0 [ACTION TAKEN BY swSim: forced these to be 0].",
                        sep = ""))
      }
      else if (link == "logit") {
        stop("This combination of distribution and link is not supported.  See documentation.")
      }
      else if (link == "log") {
        mu.ijk <- exp(tmplink.mu.ijk)
      }
      ## generating response variable (Y.ijk)
      Y.ijk <- rpois(length(mu.ijk), mu.ijk)
    }
  }
  #####n a vector or matrix
  else if ((is.vector(n) & length(n) > 1) | is.matrix(n)) {
    ## create nMat [MATRIX]
    if ((is.vector(n) & length(n) > 1)) {
      nMat <- matrix(rep(n, each = design$total.time), ncol = design$total.time, byrow = TRUE)
    }
    else if (is.matrix(n)) {
      nMat <- n
    }
    nMat[is.na(link.mu.ij)] <- 0
    ## convert nMat to nMat2Vec [VECTOR]
    nMat2Vec <- as.vector(t(nMat))
    tmplink.mu.ijk <- NULL
    tmpTx <- NULL
    tmpTime <- NULL
    tmpTimeOnTx <- NULL
    tmpCluster <- NULL
    addto.mu.ijk <- NULL
    tmpIndivid <- NULL
    clusterN=NULL
    base <- 0
    for (k in 1:length(nMat2Vec)) {
      ## link(mu.ijk) [VECTOR]## (for all clusters, all time points, AND all observations)
      if (zeta>0) {
        if (is.null(addto.mu.ijk) & nMat2Vec[k]>0)  {
# for closed cohort fix size of this cluster at first non-zero n for this cluster        
          clusterN = nMat2Vec[k]
          addto.mu.ijk = rnorm(clusterN,0,zeta)
        }
        if (nMat2Vec[k]>0){
          tmplink.mu.ijk <- c(tmplink.mu.ijk, rep(as.vector(t(link.mu.ij))[k],each = clusterN)+addto.mu.ijk)
          tmpIndivid <- c(tmpIndivid,base+(1:clusterN))
        }
        if (k%%design$total.time==0) {
# set up for new cluster at next k
          base = base+clusterN
          addto.mu.ijk = NULL
          clusterN = NULL
        }
      } else {
        tmplink.mu.ijk <- c(tmplink.mu.ijk, rep(as.vector(t(link.mu.ij))[k],each = nMat2Vec[k]))
      }
      tmpTx <- c(tmpTx, rep(as.vector(t(X.ij))[k], each = nMat2Vec[k]))
      tmpTime <- c(tmpTime, rep(as.vector(t(time.ij))[k],each = nMat2Vec[k]))
      tmpTimeOnTx <- c(tmpTimeOnTx, rep(as.vector(t(timeOnTx.ij))[k],each = nMat2Vec[k]))
      tmpCluster <- c(tmpCluster, rep(as.vector(t(cluster.ij))[k],each = nMat2Vec[k]))
    }
    ## generate Y.ijk from specified distn [VECTOR]
    if (distn == "gaussian") {
      if (link == "identity") {
        mu.ijk <- tmplink.mu.ijk
      }
      else if (link == "logit") {
        stop("This combination of distribution and link is not supported.  See documentation.")
      }
      else if (link == "log") {
        stop("This combination of distribution and link is not supported. Instead, use identity link and set log.gaussian=TRUE for lognormal data.  See documentation.")
      }
      Y.ijk <- rnorm(length(mu.ijk), mu.ijk, sd = sigma)
      if(log.gaussian == TRUE){
        Y.ijk <- exp(Y.ijk)
      }
    }
    else if (distn == "binomial") {
      if (link == "identity") {
        mu.ijk <- tmplink.mu.ijk
        mu.ijk.lessthan0 <- length(mu.ijk[mu.ijk < 0])
        mu.ijk.greaterthan1 <- length(mu.ijk[mu.ijk >1])
        mu.ijk[mu.ijk < 0] <- 0
        mu.ijk[mu.ijk > 1] <- 1
        if (mu.ijk.lessthan0 > 0 | mu.ijk.greaterthan1 >0)
          if (silent==FALSE) warning(paste("When using an identity link function with family=binomial, the probability of outcome for person k in cluster i at time j, denoted by mu.ijk in swSim(), cannot be less than 0 or greater than 1. Out of the ",
                        length(mu.ijk), " probabilities computed, ",
                        mu.ijk.lessthan0, " (", mu.ijk.lessthan0/length(mu.ijk) *
                          100, "%) were less than 0 [ACTION TAKEN BY swSim: forced these to be 0] and ",
                        mu.ijk.greaterthan1, " (", mu.ijk.greaterthan1/length(mu.ijk) *
                          100, "%) were greater than 1 [ACTION TAKEN BY swSim: forced these to be 1].",
                        sep = ""))
      }
      else if (link == "logit") {
        mu.ijk <- expit(tmplink.mu.ijk)
      }
      else if (link == "log") {
        mu.ijk <- exp(tmplink.mu.ijk)
        mu.ijk.greaterthan1 <- length(mu.ijk[mu.ijk >1])#Edit: fixed from < 0
        mu.ijk[mu.ijk > 1] <- 1
        if (mu.ijk.greaterthan1 > 1)
          if (silent==FALSE) warning(paste("When using a log link function with family=binomial, the mean of outcome for person k in cluster i at time j, denoted by mu.ijk in swSim(), cannot be greater than 1. Out of the ",
                        length(mu.ijk), " means computed, ", mu.ijk.greaterthan1,
                        " (", mu.ijk.greaterthan1/length(mu.ijk) *
                          100, "%) were greater than 1 [ACTION TAKEN BY swSim: forced these to be 1].",
                        sep = ""))
      }
      Y.ijk <- rbinom(length(mu.ijk), 1, mu.ijk)
    }
    else if (distn == "poisson") {
      if (link == "identity") {
        mu.ijk <- tmplink.mu.ijk
        mu.ijk.lessthan0 <- length(mu.ijk[mu.ijk < 0])
        mu.ijk[mu.ijk < 0] <- 0
        if (mu.ijk.lessthan0 > 0 )
          if (silent==FALSE) warning(paste("When using an identity link function with family=poisson, the expected count for person k in cluster i at time j, denoted by mu.ijk in swSim(), cannot be less than 0. Out of the ",
                        length(mu.ijk), " probabilities computed, ",
                        mu.ijk.lessthan0, " (", mu.ijk.lessthan0/length(mu.ijk) *
                          100, "%) were less than 0 [ACTION TAKEN BY swSim: forced these to be 0].",
                        sep = ""))
      }
      else if (link == "logit") {
        stop("This combination of distribution and link is not supported.  See documentation.")
      }
      else if (link == "log") {
        mu.ijk <- exp(tmplink.mu.ijk)
      }
      Y.ijk <- rpois(length(mu.ijk), mu.ijk)
    }
  }
  ##########
  ##Arrange data for export
  ##########
  ## response variable [VECTOR] (for each cluster, each time point, and each observation)
  response.var <- Y.ijk
  ## treatment variable [VECTOR] (for each cluster, each time point, and each observation)
  tx.var <- tmpTx
  ## time-on-treatment variable [VECTOR] (for each cluster, each time point, and each observation)
  timeOnTx.var <- tmpTimeOnTx
  ## time variable [VECTOR] (for each cluster, each time point, and each observation)
  time.var <- tmpTime
  ## cluster variable [VECTOR] (for each cluster, each time point, and each observation)
  cluster.var <- tmpCluster
  ## swData [DATAFRAME] (result of function)
  if (zeta > 0) {
  ## For closed cohort, individual id variable ]VECTOR]
    individ.var <- tmpIndivid
    swData <- data.frame(response.var, tx.var, timeOnTx.var, time.var, cluster.var, individ.var)
  }
  else {
    swData <- data.frame(response.var, tx.var, timeOnTx.var, time.var, cluster.var)
  }
  ## returning swData [DATAFRAME]
  swData
}


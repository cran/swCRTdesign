swSim <- function (design, family, log.gaussian=FALSE, n, mu0, mu1, time.effect, sigma, tau=NULL,
                         eta=NULL, rho=NULL, gamma=NULL, icc=NULL, cac=NULL, time.lab = NULL, retTimeOnTx = FALSE)
{
  #updated 8/5/2019 for v. 3.0
  #Changes made from the previous version which fixed errors that could have affected code are indicated by 'Edit:'

  #####
  #Warnings, unpack family, translate ICC/CAC
  #####
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

  #Basic input checks
  if(! all(n%%1 == 0)){
    stop("n (either scalar, vector, or matrix) must consist only of integers.")
  }

  #Checks to make sure people are using random effects OR ICC/CAC.
  param.icc.all <- !is.null(icc) & !is.null(cac)
  param.re.all <- !is.null(tau) & !is.null(eta) & !is.null(rho) & !is.null(gamma)

  param.icc.any <- !is.null(icc) | !is.null(cac)
  param.re.any <- !is.null(tau) | !is.null(eta) | !is.null(rho) | !is.null(gamma)

  if(param.icc.all == FALSE & param.re.all == FALSE){
    stop("Either enter values for both ICC and CAC, or all of tau, eta, gamma, and rho.  Note for users familiar with version 2.2.0: rho had a default value of 0, so you may need to add 'rho=0' to pre-existing code.")
  }
  if(param.re.any == TRUE & param.icc.any == TRUE){
    stop("The two parameterizations (random effects and ICC/CAC) are mututally exclusive.  Either enter values for both ICC and CAC, or all of tau, eta, gamma, and rho.")
  }

  #If using ICC and CAC, translate to random effects
  if (param.icc.all == TRUE){
    #Check range restrictions
    if(icc < 0 | icc > 1 | cac < 0 | cac > 1){
      stop("Both the ICC and CAC must be between 0 and 1.")
    }
    if(icc == 1){
      stop("There are multiple combinations of random effects that can make the ICC be 1; if you believe this is a realistic scenario, use the random effect parameterization.")
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
      #Note: for some links, it may be appropriate to try to extend the ICC/CAC definitions to non-Gaussian distributions.  Because of the links here, we're taking the option out for now.
      #mubar.temp <- (mu1+mu0)/2
      #sigmasq.temp <- mubar.temp*(1-mubar.temp)
      stop("For a Binomial outcome, please use the random effects parameterization instead of ICC/CAC.")
    }
    if (distn == 'poisson'){
      #mubar.temp <- (mu1+mu0)/2
      #sigmasq.temp <- mubar.temp
      stop("For a Poisson outcome, please use the random effects parameterization instead of ICC/CAC.")
    }
    if(cac == 1){
      gamma <- 0
      tau <- sqrt(sigmasq.temp*icc/(1-icc))
    }else{
      gamma <- sqrt(icc*sigmasq.temp*(1-cac)/(1-icc))
      tau <- sqrt(gamma^2*cac/(1-cac))
    }
  }

  #Basic restrictions on newly defined variance components
  if(rho < -1 | rho > 1){
    stop("Rho must be a numeral between -1 and 1.")
  }
  if(tau < 0 | eta < 0 | gamma < 0){
    stop("Tau, eta, and gamma must be non-negative numerals.")
  }
  if((tau == 0 | eta == 0) & rho != 0){
    stop("If either tau or eta is zero, rho must be zero, since fixed effects cannot be correlated.")
  }
  if(distn == 'gaussian' && sigma == 0 & gamma == 0){
    stop("For a non-deterministic Gaussian outcome, at least one of sigma and gamma needs to be non-zero.")
  }
  if (distn != "gaussian" & missing(sigma) == FALSE){
    warning("sigma is only used when family is gaussian")
  }
  if (distn != "gaussian" & log.gaussian == TRUE){
    stop("The log.gaussian=TRUE option must be paired with a Gaussian family.")
  }
  if (length(tau) > 1 | length(eta) > 1) {
    stop("Function cannot simulate stepped wedge design data for tau-vector or eta-vector; tau and eta must be scalars.")
  }
  if (any(rowSums(design$swDsn) == 0)) {
    warning("For the specified total number of clusters (I), total number of time periods (J), and number of cluster repetitions (I.rep), the specified stepped wedge design has at least one cluster which does not crossover from control(0) to treatment(1) arm.")
  }
  if (is.vector(n) & length(n) > 1 & length(n) != design$n.clusters){
    stop("The number of clusters in 'design' (design$n.clusters) and 'n' (length(n)) do not match.")
  }
  if (is.matrix(n) && ((nrow(n) != design$n.clusters)|(ncol(n) != design$total.time))){
    stop("The number of clusters and/or time steps in 'design' (design$n.clusters and design$total.time) and 'n' (number of rows and columns, respectively) do not match.")
  }
  if (length(n)>1){
    warning("When sample sizes are not uniform, power depends on order of clusters (see documentation).")
  }
  #####
  link <- family$link
  p0 <- mu0
  p1 <- mu1
  theta <- p1 - p0## treatment effect [SCALAR]
  CV.p0 <- tau/p0## standard deviation of random intercepts [SCALAR]
  CV.theta <- eta/abs(theta)## standard deviation of random slopes [SCALAR]
  X.ij <- design$swDsn## SW CRT design [MATRIX]
  if (length(time.effect) == 1) {
    time.effectVec <- rep(time.effect, design$total.time)
  }
  else if (length(time.effect) == design$total.time) {
    time.effectVec <- time.effect
  }
  else {
    warning("Invalid time.effects length. Either specify a scalar fixed time effect (i.e., the same fixed time effect at each time point), or specify a vector of fixed time effects for each time point. It is best to ignore any results which follow as a result of this warning message, and correctly assign value(s) to the 'time.effect' function argument of swSim().")
  }
  ## time effect [MATRIX]
  ## (for all clusters, all time points, but only 1 observation per (i,j)-pair)
  beta.ij <- matrix(rep(time.effectVec, design$n.clusters),design$n.clusters, design$total.time, byrow = TRUE)
  ## treatment effect [MATRIX]
  ## (for all clusters, all time points, but only 1 observation per (i,j)-pair)
  thetaX.ij <- X.ij * theta
  timeOnTx.ij <- swDsn(clusters = design$clusters, tx.effect = (1:design$total.time))$swDsn
  ## covariance matrix for random intercepts and random treatment effects [MATRIX] (can do gamma separately because assumed no correlation)
  sigMat <- matrix(c(tau^2, rho * tau * eta, rho * tau * eta, eta^2), 2, 2)
  ## set random seed; default is NULL, which does nothing.
  #Edit: removed this function from version 3.0
  #set.seed(seed)
  #####
  #Generate random effects at the coarsest level
  #####
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
  rX.ij <- X.ij * matrix(r.i, nrow = design$n.clusters, ncol = design$total.time, byrow = FALSE)
  ## random time effects [MATRIX]
  g.ij <- matrix(g.i,nrow=design$n.clusters, ncol=design$total.time, byrow=FALSE)
  #####
  #Combine fixed and random effects on link scale
  #####
  ## link(mu.ijk) [MATRIX]
  ## (for all clusters, all time points, but only 1 observation per (i,j)-pair)
  link.mu.ij <- mu0 + beta.ij + thetaX.ij + a.ij + rX.ij + g.ij
  ## expit [FUNCTION]
  expit <- function(x) {
    1/(1 + exp(-x))
  }
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
    warning("The length of the specified 'time.lab' vector is not equal to the total time points determined based on the SW design from specifying 'cluters', 'extra.time', and 'all.ctl.time0'. Ignore any results that follow and re-specify the vector 'time.lab' if you want to specify the label at each time point; specify 'NULL' to have the default labeling 1, 2, ... to the total number of time points.")
  }
  ## cluster [MATRIX]
  ## (for all clusters, all time points, but only 1 observation per (i,j)-pair)
  cluster.ij <- matrix(rep(c(1:design$n.clusters), each = design$total.time), design$n.clusters, design$total.time, byrow = TRUE)
  #####
  ##Generate random data for each individual
  #####
  #Note: some features were removed when updating to version 3, and the potential combinations of links and families were reduced.  The code is left here from the previous version, commented out, in case we want to add features back in or users want to see what old code was using.

  ## creating response, tx, time, and cluster variables
  ## for specified observations 'n'
  ##   - scalar (same 'n' for each cluster AND each time point)
  ##   + vector (each element of vector is unique 'n[k]' for all time points within each cluster)
  ##   + matrix (each element of matrix is unique 'n[k,l]' for each cluster and each time points)

  #####n a scalar
  if (is.vector(n) & length(n) == 1) {
    ## link(mu.ijk) [VECTOR]
    ## (for all clusters, all time points, AND all observations)
    tmplink.mu.ijk <- rep(as.vector(t(link.mu.ij)), each = n)
    tmpTx <- rep(as.vector(t(X.ij)), each = n)
    tmpTime <- rep(as.vector(t(time.ij)), each = n)
    tmpTimeOnTx <- rep(as.vector(t(timeOnTx.ij)), each = n)
    tmpCluster <- rep(as.vector(t(cluster.ij)), each = n)
    ## generate Y.ijk from specified distn [VECTOR]
    if (distn == "gaussian") {
      ## use inverse-link function [VECTOR]
      if (link == "identity") {
        mu.ijk <- tmplink.mu.ijk
      }
      else if (link == "logit") {
        #mu.ijk <- expit(tmplink.mu.ijk)
        #mu.ijk.lessthan0 <- length(mu.ijk[mu.ijk < 0])
        #mu.ijk.greaterthan1 <- length(mu.ijk[mu.ijk >
                                               #1])
        #mu.ijk[mu.ijk < 0] <- 0
        #mu.ijk[mu.ijk > 1] <- 1
        #if (mu.ijk.lessthan0 > 0 | mu.ijk.greaterthan1 >
            #0)
          #warning(paste("When using a logit link function with family=gaussian, the probability of outcome for person k in cluster i at time j, denoted by mu.ijk in swSim(), cannot be less than 0 or greater than 1. Out of the ",
                        #length(mu.ijk), " probabilities computed, ",
                        #mu.ijk.lessthan0, " (", mu.ijk.lessthan0/length(mu.ijk) *
                          #100, "%) were less than 0 [ACTION TAKEN BY swSim: forced these to be 0] and ",
                        #mu.ijk.greaterthan1, " (", mu.ijk.greaterthan1/length(mu.ijk) *
                          #100, "%) were greater than 1 [ACTION TAKEN BY swSim: forced these to be 1].",
                        #sep = ""))
        stop("As of version 3.0, this combination of distribution and link is no longer supported.  See documentation.")
      }
      else if (link == "log") {
        #mu.ijk <- exp(tmplink.mu.ijk)
        #mu.ijk.lessthan0 <- length(mu.ijk[mu.ijk < 0])
        #mu.ijk[mu.ijk < 0] <- 0
        #if (mu.ijk.lessthan0 > 0)
          #warning(paste("When using a log link function with family=gaussian, the mean of outcome for person k in cluster i at time j, denoted by mu.ijk in swSim(), cannot be less than 0. Out of the ",
                        #length(mu.ijk), " means computed, ", mu.ijk.lessthan0,
                        #" (", mu.ijk.lessthan0/length(mu.ijk) * 100,
                        #"%) were less than 0 [ACTION TAKEN BY swSim: forced these to be 0].",
                        #sep = ""))
        stop("As of version 3.0, this combination of distribution and link is no longer supported.  See documentation.")
      }
      ## generating response variable (Y.ijk)
      Y.ijk <- rnorm(n * design$n.clusters * design$total.time, mu.ijk, sd = sigma)
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
          warning(paste("When using an identity link function with family=binomial, the probability of outcome for person k in cluster i at time j, denoted by mu.ijk in swSim(), cannot be less than 0 or greater than 1. Out of the ",
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
          warning(paste("When using a log link function with family=binomial, the mean of outcome for person k in cluster i at time j, denoted by mu.ijk in swSim(), cannot be greater than 1. Out of the ",
                        length(mu.ijk), " means computed, ", mu.ijk.greaterthan1,
                        " (", mu.ijk.greaterthan1/length(mu.ijk) *
                          100, "%) were greater than 1 [ACTION TAKEN BY swSim: forced these to be 1].",
                        sep = ""))
      }
      ## generating response variable (Y.ijk)
      Y.ijk <- rbinom(n * design$n.clusters * design$total.time,1, mu.ijk)
    }
    else if (distn == "poisson") {
      ## use inverse-link function [VECTOR]
      if (link == "identity") {#Edit: removed restriction lambda<=1
        mu.ijk <- tmplink.mu.ijk
        mu.ijk.lessthan0 <- length(mu.ijk[mu.ijk < 0])
        #mu.ijk.greaterthan1 <- length(mu.ijk[mu.ijk >
                                               #1])
        mu.ijk[mu.ijk < 0] <- 0
        #mu.ijk[mu.ijk > 1] <- 1
        if (mu.ijk.lessthan0 > 0 )
          warning(paste("When using an identity link function with family=poisson, the expected count for person k in cluster i at time j, denoted by mu.ijk in swSim(), cannot be less than 0. Out of the ",
                        length(mu.ijk), " probabilities computed, ",
                        mu.ijk.lessthan0, " (", mu.ijk.lessthan0/length(mu.ijk) *
                          100, "%) were less than 0 [ACTION TAKEN BY swSim: forced these to be 0].",
                        sep = ""))
      }
      else if (link == "logit") {
        #mu.ijk <- expit(tmplink.mu.ijk)
        #mu.ijk.lessthan0 <- length(mu.ijk[mu.ijk < 0])
        #mu.ijk.greaterthan1 <- length(mu.ijk[mu.ijk >
                                               #1])
        #mu.ijk[mu.ijk < 0] <- 0
        #mu.ijk[mu.ijk > 1] <- 1
        #if (mu.ijk.lessthan0 > 0 | mu.ijk.greaterthan1 >
            #0)
          #warning(paste("When using a logit link function with family=poisson, the probability of outcome for person k in cluster i at time j, denoted by mu.ijk in swSim(), cannot be less than 0 or greater than 1. Out of the ",
                        #length(mu.ijk), " probabilities computed, ",
                        #mu.ijk.lessthan0, " (", mu.ijk.lessthan0/length(mu.ijk) *
                          #100, "%) were less than 0 [ACTION TAKEN BY swSim: forced these to be 0] and ",
                        #mu.ijk.greaterthan1, " (", mu.ijk.greaterthan1/length(mu.ijk) *
                          #100, "%) were greater than 1 [ACTION TAKEN BY swSim: forced these to be 1].",
                        #sep = ""))
        stop("As of version 3.0, this combination of distribution and link is no longer supported.  See documentation.")
      }
      else if (link == "log") {
        mu.ijk <- exp(tmplink.mu.ijk)
      }
      ## generating response variable (Y.ijk)
      Y.ijk <- rpois(n * design$n.clusters * design$total.time, mu.ijk)
    }
  }
  else if ((is.vector(n) & length(n) > 1) | is.matrix(n)) {#####n a vector or matrix
    ## create nMat [MATRIX]
    if ((is.vector(n) & length(n) > 1)) {
      nMat <- matrix(rep(n, each = design$total.time), ncol = design$total.time, byrow = TRUE)
    }
    else if (is.matrix(n)) {
      nMat <- n
    }
    ## convert nMat to nMat2Vec [VECTOR]
    nMat2Vec <- as.vector(t(nMat))
    tmplink.mu.ijk <- NULL
    tmpTx <- NULL
    tmpTime <- NULL
    tmpTimeOnTx <- NULL
    tmpCluster <- NULL
    for (k in 1:length(nMat2Vec)) {
      ## link(mu.ijk) [VECTOR]## (for all clusters, all time points, AND all observations)
      tmplink.mu.ijk <- c(tmplink.mu.ijk, rep(as.vector(t(link.mu.ij))[k],each = nMat2Vec[k]))
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
        #mu.ijk <- expit(tmplink.mu.ijk)
        #mu.ijk.lessthan0 <- length(mu.ijk[mu.ijk < 0])
        #mu.ijk.greaterthan1 <- length(mu.ijk[mu.ijk >
                                               #1])
        #mu.ijk[mu.ijk < 0] <- 0
        #mu.ijk[mu.ijk > 1] <- 1
        #if (mu.ijk.lessthan0 > 0 | mu.ijk.greaterthan1 >
            #0)
          #warning(paste("When using a logit link function with family=gaussian, the probability of outcome for person k in cluster i at time j, denoted by mu.ijk in swSim(), cannot be less than 0 or greater than 1. Out of the ",
                        #length(mu.ijk), " probabilities computed, ",
                        #mu.ijk.lessthan0, " (", mu.ijk.lessthan0/length(mu.ijk) *
                          #100, "%) were less than 0 [ACTION TAKEN BY swSim: forced these to be 0] and ",
                        #mu.ijk.greaterthan1, " (", mu.ijk.greaterthan1/length(mu.ijk) *
                          #100, "%) were greater than 1 [ACTION TAKEN BY swSim: forced these to be 1].",
                        #sep = ""))
        stop("As of version 3.0, this combination of distribution and link is no longer supported.  See documentation.")
      }
      else if (link == "log") {
        #mu.ijk <- exp(tmplink.mu.ijk)
        #mu.ijk.lessthan0 <- length(mu.ijk[mu.ijk < 0])
        #mu.ijk[mu.ijk < 0] <- 0
        #if (mu.ijk.lessthan0 > 0)
          #warning(paste("When using a log link function with family=gaussian, the mean of outcome for person k in cluster i at time j, denoted by mu.ijk in swSim(), cannot be less than 0. Out of the ",
                        #length(mu.ijk), " means computed, ", mu.ijk.lessthan0,
                        #" (", mu.ijk.lessthan0/length(mu.ijk) * 100,
                        #"%) were less than 0 [ACTION TAKEN BY swSim: forced these to be 0].",
                        #sep = ""))
        stop("As of version 3.0, this combination of distribution and link is no longer supported.  See documentation.")
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
          warning(paste("When using an identity link function with family=binomial, the probability of outcome for person k in cluster i at time j, denoted by mu.ijk in swSim(), cannot be less than 0 or greater than 1. Out of the ",
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
          warning(paste("When using a log link function with family=binomial, the mean of outcome for person k in cluster i at time j, denoted by mu.ijk in swSim(), cannot be greater than 1. Out of the ",
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
        #mu.ijk.greaterthan1 <- length(mu.ijk[mu.ijk >
                                               #1])
        mu.ijk[mu.ijk < 0] <- 0
        #mu.ijk[mu.ijk > 1] <- 1
        if (mu.ijk.lessthan0 > 0 )
          warning(paste("When using an identity link function with family=poisson, the expected count for person k in cluster i at time j, denoted by mu.ijk in swSim(), cannot be less than 0. Out of the ",
                        length(mu.ijk), " probabilities computed, ",
                        mu.ijk.lessthan0, " (", mu.ijk.lessthan0/length(mu.ijk) *
                          100, "%) were less than 0 [ACTION TAKEN BY swSim: forced these to be 0].",
                        sep = ""))
      }
      else if (link == "logit") {
        #mu.ijk <- expit(tmplink.mu.ijk)
        #mu.ijk.lessthan0 <- length(mu.ijk[mu.ijk < 0])
        #mu.ijk.greaterthan1 <- length(mu.ijk[mu.ijk >
                                               #1])
        #mu.ijk[mu.ijk < 0] <- 0
        #mu.ijk[mu.ijk > 1] <- 1
        #if (mu.ijk.lessthan0 > 0 | mu.ijk.greaterthan1 >
            #0)
          #warning(paste("When using a logit link function with family=poisson, the probability of outcome for person k in cluster i at time j, denoted by mu.ijk in swSim(), cannot be less than 0 or greater than 1. Out of the ",
                        #length(mu.ijk), " probabilities computed, ",
                        #mu.ijk.lessthan0, " (", mu.ijk.lessthan0/length(mu.ijk) *
                          #100, "%) were less than 0 [ACTION TAKEN BY swSim: forced these to be 0] and ",
                        #mu.ijk.greaterthan1, " (", mu.ijk.greaterthan1/length(mu.ijk) *
                          #100, "%) were greater than 1 [ACTION TAKEN BY swSim: forced these to be 1].",
                        #sep = ""))
        stop("As of version 3.0, this combination of distribution and link is no longer supported.  See documentation.")
      }
      else if (link == "log") {
        mu.ijk <- exp(tmplink.mu.ijk)
      }
      Y.ijk <- rpois(length(mu.ijk), mu.ijk)
    }
  }
  #####
  ##Arrange data for export
  #####
  ## response variable [VECTOR] (for each cluster, each time point, and each observation)
  response.var <- Y.ijk
  ## treatment variable [VECTOR] (for each cluster, each time point, and each observation)
  tx.var <- tmpTx
  ## time-on-treatment variable [VECTOR] (for each cluster, each time point, and each observation)
  tx.var <- ifelse(tmpTx > 0, 1, 0)
  timeOnTx.var <- tmpTimeOnTx
  ## time variable [VECTOR] (for each cluster, each time point, and each observation)
  time.var <- tmpTime
  ## cluster variable [VECTOR] (for each cluster, each time point, and each observation)
  cluster.var <- tmpCluster
  ## swData [DATAFRAME] (result of function)
  if (!retTimeOnTx) {
    swData <- data.frame(response.var, tx.var, time.var, cluster.var)
  }
  else if (retTimeOnTx) {
    swData <- data.frame(response.var, tx.var, timeOnTx.var, time.var, cluster.var)
  }
  ## returning swData [DATAFRAME]
  swData
}


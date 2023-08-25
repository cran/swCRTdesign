swPwr <- function (design, distn, n, mu0, mu1, H = NULL, sigma, tau=0, eta=0, rho=0, gamma=0, zeta=0, icc=0, cac=1, iac=0,
                         alpha = 0.05, retDATA = FALSE, silent = FALSE)
{
  #Version: 4/28/2023, v. 4.0, Jim Hughes
  # Note: zeta = 0 and iac = 0 means function defaults to a cross-sectional design; either zeta >0 or iac>0 means a closed cohort design
  ##### Helper functions #######
  replaceNAwithZero = function(x){
    ifelse(is.na(x),0,x)
  }
  ##########
  #Warnings
  ##########
  #Keep this warning around for a couple versions; if users didn't specify argument names, could change results with no indication that anything has happened.
  #if(silent == FALSE){
  #warning("The order of variance component arguments has changed for swPwr (in version 2.2.0, it was tau, eta, rho, sigma); please modify existing code if necessary. ")
  #}
  #Basic input checks
  if(alpha <= 0 | alpha >= 1){
    stop("Alpha must be strictly between 0 and 1.")
  }
  if(! all(n%%1 == 0)){
    stop("n (either scalar, vector, or matrix) should consist only of integers.")
  }
  if(!is.null(H)){
    if (design$nTxLev > 1){
    stop("Do not use H when you have specfied a design with multiple treatment levels. If you want an ETI model based estimator, specify a 0/1 intervention in swDsn or swDsn1 and use H to define the estimator.")
    }
    if (any(design$tx.effect.frac!=1)){
      stop("Do not specify fractional treatment effects in the design AND an ETI estimator")
    }
    maxET = max(as.vector(apply(replaceNAwithZero(design$swDsn),1,cumsum)))
    if (length(H) == 1) {
      H = rep(1/maxET,maxET)
      if (silent==FALSE) warning("ETI estimator with equal weighting for all exposure times assumed")
    }
    if (length(H) != maxET){
      stop("Length of H must be equal to number of exposure times. For this design, number of exposure times is ",maxET)
    }
	  if (sum(H) != 1) {
	    H = H/sum(H)
      if (silent==FALSE) warning("H has been renormalized to sum to 1.0")
	  }
  }
  if(design$nTxLev!=length(mu1)){
    stop("Length of mu1 must correspond to number of treatment levels specified in swDsn1")
  }
  if (distn == 'gaussian' & missing(sigma)) stop("If distribution is gaussian, a value for sigma must be entered")
  #Checks to make sure people are using random effects OR ICC/CAC/IAC:
  param.icc.any <- !missing(icc) | !missing(cac) | !missing(iac)
  param.re.any <- !missing(tau) | !missing(eta) | !missing(rho) | !missing(gamma) | !missing(zeta)
  if(param.re.any == TRUE & param.icc.any == TRUE){
    stop("The two parameterizations (random effects and ICC/CAC/IAC) are mutually exclusive.  Either enter values for ICC, CAC (and IAC, if cohort design), or tau, eta, gamma, rho (and zeta, if cohort design).")
  }
  ##########
  #If using ICC/CAC/IAC, translate to random effects
  ##########
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
      mubar.temp <- (mean(mu1)+mu0)/2
      sigmasq.temp <- mubar.temp*(1-mubar.temp)
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
  ##########
  #More warnings for variance components
  ##########
  #Basic restrictions on variance components
  if(rho < -1 | rho > 1){
    stop("rho must be a numeral between -1 and 1.")
  }
  if(tau < 0 | eta < 0 | gamma < 0 | zeta < 0){
    stop("tau, eta, gamma and zeta must be non-negative numerals.")
  }
  if((tau == 0 | eta == 0) & rho != 0){
    stop("If either tau or eta is zero, rho must be zero, since fixed effects cannot be correlated.")
  }
  if(distn == 'gaussian' && sigma == 0 & gamma == 0 & zeta == 0){
    stop("For a non-deterministic Gaussian outcome, at least one of sigma, zeta or gamma needs to be non-zero.")
  }
  if (length(tau) > 1 | length(eta) > 1 | length(zeta) > 1)
    stop("Function cannot compute stepped wedge design Power for tau-vector, zeta-vector or eta-vector; tau, zeta and eta must be scalars.")  ##########
  ##########
  obs.per.cluster.per.time <- n
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
  theta <- (mu1 - mu0)#treatment effect
  muBar <- (mu0 + mean(mu1))/2
  #The definition of sigSq is clear for a gaussian distribution; for a binomial distribution, we use a stand-in.
  if (distn == "gaussian"){
    sigSq <- sigma^2
  }else if (distn == "binomial") {
    if (!missing(sigma)& silent==FALSE) warning("sigma is not used when distn=binomial")
    sigSq <- muBar * (1 - muBar)
    sigma <- NA
    if ((tau^2 + eta^2 + gamma^2) > sigSq)
      stop("tau^2 + eta^2 + gamma^2 must be less than muBar*(1-muBar) when distn=binomial")
  }
#########
  I.rep <- design$clusters
  I <- design$n.clusters
  J <- design$total.time
  K <- design$nTxLev
  swDesign <- design$swDsn
  swDesignUnique <- design$swDsn.unique.clusters
  if (any(rowSums(swDesign,na.rm=TRUE) == 0)){#Note: this warning doesn't catch all cases, particularly for designs with transition periods
    warning("For the specified total number of clusters (I), total number of time periods (J), and number of cluster repetitions (I.rep), the specified stepped wedge design has at least one cluster which does not crossover from control to treatment arm.")
  }
#########
  ## Constructing the Treatment/Intervention Indicator Vector (X.ij) 
  if (is.null(H)) { 
    if (K==1) {
  # IT model with possible fractional treatment effects
      X.ij <- matrix(t(swDesign),I*J,1)
      colnames(X.ij) <- "treatment.1"      
    } else {
  # multi-level IT model (no fractional treatment effects)
      TxLev = design$TxLev
      X.ij <- matrix(0,I*J,K)
      colnames(X.ij) <- paste0("treatment.",TxLev)  
      for (k in 1:K){
        X.ij[,k] <- as.numeric(as.vector(t(swDesign)) == TxLev[k])
      }
    }
  } else {
  # ETI model (no fractional treatment effects)
  #assume no NA's in the middle of intervention periods (at end of intervention periods is okay)
	  tmp = as.vector(apply(replaceNAwithZero(swDesign),1,cumsum))
	  X.ij <- matrix(0,length(tmp),maxET)
	  for (j in 1:maxET) {X.ij[,j] = as.integer(tmp==j) }
  }  
  ## Constructing the Design Matrix (Xmat)
  beta.blk <- rbind(diag(1, J - 1, J - 1), 0)
  Xmat.blk <- matrix(rep(as.vector(t(cbind(1, beta.blk))),I), ncol = J, byrow = TRUE)
  Xmat <- cbind(Xmat.blk, X.ij)
  ##########
  ## Constructing the Covariance Matrix (Wmat); depends on configuration of n
  ##########
  #Make nMat, a matrix with a sample size for every cluster (rows) and time period (columns).
  if (length(n) == 1){
    Wmat.blk <- (tau^2+zeta^2/n) + diag(sigSq/n+gamma^2, J)#matrix V_i
    Wmat.partial <- kronecker(diag(1, I), Wmat.blk)#matrix V
    #above, I is the number of times that Wmat.blk gets repeated along the diagonal of a matrix filled with 0's otherwise
  }else{
    #Took this nMat creation from swSim: each row is for a specific cluster
    if ((is.vector(n) & length(n) > 1)) {
      #n is a vector with one entry per cluster
      if (length(n) != design$n.clusters){
        stop("The number of clusters in 'design' (design$n.clusters) and 'n' (length(n)) do not match.")
      }
      nMat <- matrix(rep(n, each = design$total.time),
                     ncol = design$total.time, byrow = TRUE)#Turns vector into COLUMNS of nMat
    }else if (is.matrix(n)) {
      #n is a matrix with one row per cluster and one column per time point
      if ((nrow(n) != design$n.clusters)|(ncol(n) != design$total.time)){
        stop("The number of clusters and/or time steps in 'design' (design$n.clusters and design$total.time) and 'n' (number of rows and columns, respectively) do not match.")
      }
      nMat <- n#
    }
    Wmat.partial <- matrix(0,nrow=I*J,ncol=I*J)#Wmat.blk doesn't get used anywhere else, so just setting up the Wmat.partial
    for(i in 1:I){#for each of the I total clusters...
#remember that sample size must be constant within a cluster if zeta > 0
      Wmat.partial.i <- (tau^2 + zeta^2/sqrt(t(nMat[i,,drop=F])%*%nMat[i,])) + diag(sigSq/nMat[i,]) + diag(gamma^2,length(nMat[i,]))#Wmat.partial for cluster i
      kronecker.diag.i <- rep(0,I)
      kronecker.diag.i[i] <- 1
      addition.i <- kronecker(diag(kronecker.diag.i), Wmat.partial.i)#Block diagonal matrix with the entries for block i
      #For a design with transition periods, need to clean out NaN's produced by multiplying 0*Inf (really want them to be 0)
      if(is.matrix(n)&any(n==0)){
        addition.i[is.nan(addition.i)] <- 0
      }
      Wmat.partial <- Wmat.partial + addition.i#For each cluster, we fill in another section of the block diagonal matrix
    }
  }
  #
  #Make a matrix with the random effects that vary by treatment (eta, rho) in the correct locations, and add it to the other random effects.
  Xij.Xil.ARRAY <- array(NA, c(J, J, length(I.rep)))
  for (i.indx in 1:length(I.rep)) {
    Xi.jl <- as.numeric(swDesignUnique[i.indx, ] > 0)
    Xij.Xil.ARRAY[, , i.indx] <- (Xi.jl %o% Xi.jl) * eta^2 + outer(Xi.jl, Xi.jl, "+") * rho * eta * tau
  }
  Xij.Xil.LIST_blk <- sapply(1:length(I.rep), function(x) NULL)
  for (i.indx in 1:length(I.rep)) {
    Xij.Xil.LIST_blk[[i.indx]] <- kronecker(diag(1, I.rep[i.indx]), Xij.Xil.ARRAY[, , i.indx])
  }
  W.eta <- blkDiag(Xij.Xil.LIST_blk)
  Wmat <- Wmat.partial + as.matrix(W.eta)#Covariance matrix
 
  #Wmat has entries of Inf in rows/columns that correspond to timexcluster periods with no observations.
  #Since these contain no information for the study, we can just remove them, and also adjust Xmat
  NAtimes <- is.na(as.vector(t(swDesign)))
  if((is.matrix(n)&any(n==0))|sum(NAtimes)>0){
    indices <- (rowSums(replaceNAwithZero(Wmat))==0) | NAtimes | (is.matrix(n) & as.vector(t(n)==0))
    Wmat <- Wmat[!indices,!indices]
    Xmat <- Xmat[!indices,]
  }
  ##########
  #Use design matrix and covariance matrix to calculate power
  ##########

  if (is.null(H)) { 
  # IT model
    np <- ncol(Xmat)
    var <- solve(t(Xmat) %*% solve(Wmat) %*% Xmat)[(np-K+1):np,(np-K+1):np,drop=F]
    var.theta.WLS <- diag(var)
	} else {
  # ETI model
    eti = (ncol(Xmat)-maxET+1):ncol(Xmat)
    var <- solve(t(Xmat) %*% solve(Wmat) %*% Xmat)[eti,eti]
    var.theta.WLS <- t(as.matrix(H))%*%var%*%as.matrix(H)
    }
  pwrWLS <- pnorm(abs(theta)/sqrt(var.theta.WLS) - qnorm(1 -alpha/2)) + pnorm(-abs(theta)/sqrt(var.theta.WLS) - qnorm(1 -alpha/2))
  rslt <- pwrWLS
  ##########
  ## Closed-form approach/solution
  ##########
  ##   eta==0 and gamma==0 and zeta==0 (i.e., *NO* random treatment or time or individual)
  # if (eta == 0 & gamma == 0 & zeta == 0 & length(n) == 1 & is.null(H)) {#we can only calculate closed form when n is an integer
    ## Closed-form formula
  #  X <- swDesign
  #  U <- sum(X)
  #  W <- sum(colSums(X)^2)
  #  V <- sum(rowSums(X)^2)
    ## 1. Obtain Variance Estimate of theta from Closed-form (var.theta.CLOSED)
    ## 2. Calculate Power for theta from Closed-form (pwrCLOSED)
    ## 3. Storing/Appending Resulting Power for theta from Closed-form (rslt)
  #  sigSq <- sigSq/n
  #  var.theta.CLOSED <- I * sigSq * (sigSq + J * tau^2)/((I *U - W) * sigSq + (U^2 + I * J * U - J * W - I * V) *tau^2)
  #  pwrCLOSED <- pnorm(abs(theta)/sqrt(var.theta.CLOSED) -qnorm(1 - alpha/2)) + pnorm(-abs(theta)/sqrt(var.theta.CLOSED) -qnorm(1 - alpha/2))
    ## From Excel Spreadsheet:
    ##
    ## varTheta1 <- (I*sigSq)*(sigSq+(J*tau^2))
    ## varTheta2 <- ((I*U)-W)*sigSq
    ## varTheta3 <- ((U^2)+(I*J*U)-(J*W)-(I*W))*tau^2
    ## varTheta <- varTheta1 / (varTheta2 + varTheta3)
    ##
    ## pwr1 <- sqrt(theta^2 / varTheta)
    ## pwr2 <- pwr1 - qnorm(1-alpha/2)
    ## pnorm( pwr2 )
 # }else {
 #   pwrCLOSED <- NA
 # }
  ##########
  ## Returning Resulting Power(s) for fixed theta of the specified SW design
  ##########
  if (retDATA)
    rslt <- list(design = design, n = n, mu0 = mu0, mu1 = mu1, H=H,
                 tau = tau, eta = eta, rho = rho, sigma = sigma, gamma=gamma, zeta=zeta, alpha = alpha,
                 Xmat = Xmat, Wmat = Wmat, var=var, var.theta.WLS = var.theta.WLS,
                 pwrWLS = pwrWLS)
  rslt
}

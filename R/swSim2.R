swSim2 <- function (design, family, log.gaussian=FALSE, n, mu0, mu1, time.effect, sigma, tau=0,
                         eta=0, rho=0, gamma=0, zeta=0, icc=0, cac=1, iac=0, ar=1, time.lab = NULL, 
                         silent = FALSE, nocheck=FALSE)
{
  #
  # Last update: 2/11/25, v. 4.1, Jim Hughes
  # In this file, 'random slopes' and 'random treatment effects' are used interchangeably
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
      stop("Length of mu1 must correspond to number of treatment levels specified in swDsn")
    }
    }
  if (distn == 'gaussian' & missing(sigma)) stop("If distribution is gaussian, a value for sigma must be entered")
  #Checks to make sure people are using random effects OR ICC/CAC.
    param.icc.any <- !missing(icc) | !missing(cac) | !missing(iac)
    param.re.any <- !missing(tau) | !missing(eta) | !missing(rho) | !missing(gamma) | !missing(zeta)
    if(param.re.any == TRUE & param.icc.any == TRUE){
      stop("The two parameterizations (random effects and ICC/CAC/IAC) are mutually exclusive.  Either enter values for ICC, CAC (and IAC, if cohort design), or tau, eta, gamma, rho (and zeta, if cohort design).")
    }
  if (ar!=1 & (zeta!=0 | iac!=0)) {
      stop("autoregressive structure (ar!=1) only allowed with cross-sectional sampling")
  }
  if (ar!=1 & (eta>0 | gamma>0)) {
      stop("random treatment effect and tandom time effect not allowed with autoregressive structure (ar!=1)")
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
  if(ar < -1 | ar > 1){
      stop("ar must be a numeral between -1 and 1.")
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
  if (length(tau) > 1 | length(eta) > 1 | length(zeta) > 1 | length(ar) > 1) {
    stop("Function cannot simulate stepped wedge design data for tau-vector, zeta-vector, eta-vector or ar-vector; tau, zeta, eta and ar must be scalars.")
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
  X.ij <- design$swDsn## SW CRT design [MATRIX]
  #Make matrix of time on treatment (Note: time on any treatment, regardless of fractional treatment effect)
  timeOnTx.ij.pre <- replaceNAwithZero(design$swDsn)
  timeOnTx.ij <- t(apply(timeOnTx.ij.pre>0, 1, cumsum))
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
  #####n a scalar
  if (is.vector(n) & length(n) == 1) {
    ## (for all clusters, all time points, AND all observations)
    if (zeta>0) {
      tmpIndivid <- NULL
      base <- 0
      for (i in 1:design$n.clusters) {
        tmpIndivid <- c(tmpIndivid,rep(base+(1:n),design$total.time))
        base <- base+n
      }
    }
    miss = is.na(as.vector(t(X.ij)))
    tmpTx <- rep(as.vector(t(X.ij)), each = n)[!miss]
    tmpTime <- rep(as.vector(t(time.ij)), each = n)[!miss]
    tmpTimeOnTx <- rep(as.vector(t(timeOnTx.ij)), each = n)[!miss]
    tmpCluster <- rep(as.vector(t(cluster.ij)), each = n)[!miss]
    if (zeta>0) tmpIndivid <- tmpIndivid[!miss]
  #####n a vector or matrix
  else if ((is.vector(n) & length(n) > 1) | is.matrix(n)) {
    ## create nMat [MATRIX]
    if ((is.vector(n) & length(n) > 1)) {
      nMat <- matrix(rep(n, each = design$total.time), ncol = design$total.time, byrow = TRUE)
    }
    else if (is.matrix(n)) {
      nMat <- n
    }
    nMat[is.na(X.ij)] <- 0
    ## convert nMat to nMat2Vec [VECTOR]
    nMat2Vec <- as.vector(t(nMat))
    tmpTx <- NULL
    tmpTime <- NULL
    tmpTimeOnTx <- NULL
    tmpCluster <- NULL
    tmpIndivid <- NULL
    clusterN=NULL
    base <- 0
    for (k in 1:length(nMat2Vec)) {
      ## link(mu.ijk) [VECTOR]## (for all clusters, all time points, AND all observations)
      if (zeta>0) {
        if (is.null(clusterN) & nMat2Vec[k]>0)  {
# for closed cohort fix size of this cluster at first non-zero n for this cluster        
          clusterN = nMat2Vec[k]
        }
        if (nMat2Vec[k]>0){
          tmpIndivid <- c(tmpIndivid,base+(1:clusterN))
        }
        if (k%%design$total.time==0) {
# set up for new cluster at next k
          base = base+clusterN
          clusterN = NULL
        }
      } 
      tmpTx <- c(tmpTx, rep(as.vector(t(X.ij))[k], each = nMat2Vec[k]))
      tmpTime <- c(tmpTime, rep(as.vector(t(time.ij))[k],each = nMat2Vec[k]))
      tmpTimeOnTx <- c(tmpTimeOnTx, rep(as.vector(t(timeOnTx.ij))[k],each = nMat2Vec[k]))
      tmpCluster <- c(tmpCluster, rep(as.vector(t(cluster.ij))[k],each = nMat2Vec[k]))
    }
  }
  ##########
  ##Arrange data for export
  ##########
  ## treatment variable [VECTOR] (for each cluster, each time point, and each observation)
  tx.var <- tmpTx
  ## time-on-treatment variable [VECTOR] (for each cluster, each time point, and each observation)
  timeOnTx.var <- tmpTimeOnTx
  ## time variable [VECTOR] (for each cluster, each time point, and each observation)
  time.var <- tmpTime
  ## cluster variable [VECTOR] (for each cluster, each time point, and each observation)
  cluster.var <- tmpCluster
   ## Make a dummy response variable for simulate_new (will be replaced with simulated response)
  response.var = runif(length(cluster.var),0,1) 
  ## swData [DATAFRAME] (result of function)
  if (zeta > 0) {
  ## For closed cohort, individual id variable ]VECTOR]
    individ.var <- tmpIndivid
  }
  ############
  # Set up parameters for simulate_new
  # Fixed effects
  if (length(time.effect) == 1) {
    time.effectVec <- rep(time.effect, design$total.time)
  }
  else if (length(time.effect) == design$total.time) {
    time.effectVec <- time.effect
  } else {
    stop("Invalid time.effects length. Either specify a scalar fixed time effect (i.e., the same fixed time effect at each time point), or specify a vector of fixed time effects for each time point. It is best to ignore any results which follow as a result of this warning message, and correctly assign value(s) to the 'time.effect' function argument of swSim().")
  }   
  theta <- mu1 - mu0## treatment effect  
  fixedeffects = c(mu0+time.effectVec[1],time.effectVec[-1],theta)
  ###########
  # Set up dataframe and model
  ftimeontx.var=as.factor(timeOnTx.var)
  ftime.var = as.factor(time.var)
  fcluster.var = as.factor(cluster.var)
  ftx.var = as.factor(tx.var)}
  itx.var = as.integer(tx.var>0)
  if (zeta>0) {
      findivid.var = as.factor(individ.var)  
      junk = data.frame(response.var,ftime.var,ftx.var,itx.var,ftimeontx.var,fcluster.var,findivid.var)
    } else {
      junk = data.frame(response.var,ftime.var,ftx.var,itx.var,ftimeontx.var,fcluster.var)  
    }
  itot = ifelse(length(theta)==1,0,1)
  if (ar==1){
    # No autoregressive
    model = as.formula(paste0(" ~ "," ftime.var",
                            ifelse(itot==0,"+ ftx.var","+ ftimeontx.var"),
                            ifelse(tau>0&(eta==0|(eta>0&rho==0))," + (1|fcluster.var)",""),
                            ifelse(gamma>0," + (1|fcluster.var:ftime.var)",""),
                            ifelse(eta>0&rho!=0," + (itx.var | fcluster.var)",""),
                            ifelse(eta>0&rho==0," + (0 + itx.var | fcluster.var)",""),
                            ifelse(zeta>0,"+ (1|findivid.var)",""))) 
    varcomps = NULL
    if (tau>0&(eta==0|(eta>0&rho==0))) varcomps=c(varcomps,log(tau))
    if (gamma>0&ar==1) varcomps=c(varcomps,log(gamma))
    if (eta>0&rho!=0) varcomps=c(varcomps,log(tau),log(eta),rho/sqrt(1-rho^2))
    if (eta>0&rho==0) varcomps=c(varcomps,log(eta))
    if (zeta>0) varcomps=c(varcomps,log(zeta))
  } else {
    # Autoregressive
    model = as.formula(paste0(" ~ "," ftime.var",
                              ifelse(itot==0,"+ ftx.var","+ ftimeontx.var"),
                              ifelse(tau>0," + ar1(0 + ftime.var|fcluster.var)","")))
#                              ,ifelse(gamma>0," + (1|fcluster.var:ftime.var)",""))) # this only contributes to diagonal of covariance matrix so no need to apply ar1 to this term
    # Random effects
    varcomps = NULL
    if (tau>0) varcomps=c(varcomps,log(tau))
    varcomps = c(varcomps,ar/sqrt(1-ar^2))
#    if (gamma>0) varcomps=c(varcomps,log(gamma))
  }
  if (distn=="gaussian") {
    params = list(beta=fixedeffects,betadisp=log(sigma),theta=varcomps)
    } else {
    params = list(beta=fixedeffects,theta=varcomps)
    }
  ##########
  sim_obj = simulate_new(model,nsim=1,seed=NULL,family=family,newdata=junk,newparams=params)
  junk$response.var = sim_obj[[1]]
   ## returning swData [DATAFRAME]
  ## swData [DATAFRAME] (result of function)
  if (zeta > 0) {
    ## For closed cohort, individual id variable ]VECTOR]
    swData <- data.frame(response.var=junk$response.var, tx.var=junk$itx.var, timeOnTx.var=junk$ftimeontx.var, time.var=junk$ftime.var, cluster.var=junk$fcluster.var, individ.var=junk$findivid.var)
  }
  else {
    swData <- data.frame(response.var=junk$response.var, tx.var=junk$itx.var, timeOnTx.var=junk$ftimeontx.var, time.var=junk$ftime.var, cluster.var=junk$fcluster.var)
  }
  swData
}
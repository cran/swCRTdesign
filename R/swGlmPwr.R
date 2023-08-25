swGlmPwr <- function(design,distn,n,fixed.intercept,fixed.treatment.effect,fixed.time.effect,H = NULL,tau=0,eta=0,rho=0,gamma=0,zeta=0,
                     alpha=0.05,retDATA=FALSE,silent=FALSE)
{
  #Version: 4/28/2023, v. 4.0, Jim Hughes
  # Note: zeta = 0 means function defaults to a cross-sectional design; either zeta >0 or iac>0 means a closed cohort design
### Help functions:
  logit <- function(x){log(x/(1 - x))}
  expit <- function(x){exp(x)/(1 + exp(x))}
  replaceNAwithZero = function(x){
    ifelse(is.na(x),0,x)
  }
### Checks and warnings  
  if (alpha <= 0 | alpha >= 1) {
    stop("Alpha must be strictly between 0 and 1.")
  }
  if (!all(n%%1 == 0)) {
    if (silent==FALSE) warning("n (either scalar, vector, or matrix) should consist only of integers.")
  }
  if(is.null(H)) {
    if(design$nTxLev!=length(fixed.treatment.effect)){
      stop("Length of fixed.treatment.effect must correspond to number of treatment levels specified in swDsn/swDsn1")
    }
  } else {
    if (design$nTxLev > 1){
      stop("Do not use H when you have specfied a design with multiple treatment levels. If you want an ETI model based estimator, specify a 0/1 intervention in swDsn and use H to define the estimator.")
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
    if (length(fixed.treatment.effect)!=length(H)) {
      stop("ETI model - length of fixed.treatment.effect must be same as length of H")
    }
    if (sum(H) != 1) {
      H = H/sum(H)
      if (silent==FALSE) warning("H has been renormalized to sum to 1.0")
    }
  }
  if (rho < -1 | rho > 1) {
    stop("The correlation between the random cluster effect and the random treatment effect must be a numeral between -1 and 1.")
  }
  
  if (tau < 0 | eta < 0 | gamma < 0 | zeta < 0) {
    stop("Standard deviations of the random effects must be non-negative numerals.")
  }
  
  if(distn =="binomial")
  { 
    phi <- 1
    a <- 1
    gprime <- function(x){1/(x*(1-x))}
    v <- function(x){x*(1-x)}
    h <- expit
    g <- logit
    print("Outcome type: binomial")
  } 
  if(distn == "poisson")
  {
    phi <- 1
    #take 1 for now
    a <- 1 
    gprime <- function(x){1/x}
    v <- function(x){x}
    h <- function(x){exp(x)}
    g <- function(x){log(x)}
    print("Outcome type: poisson")
  } 
  if(distn != "binomial" && distn != "poisson")
  {
    stop("Valid outcome distributions are 'binomial' and 'poisson'")
  }
  #check inputs of fixed and random effect components  
  if(is.null(fixed.intercept) | is.null(fixed.time.effect) | is.null(fixed.treatment.effect))
  {
    stop("Parameters for fixed effects (intercept, treatment effect, and time effect) must be specified.")
  }
  if(is.null(tau) | is.null(gamma) | is.null(eta) | is.null(zeta))
  {
    if (silent==FALSE) warning("Standard deviations of all random effects are needed. If a random effect does not exist, its standard deviation should be set to 0. The default standard deviations are set to 0. ")
    if(is.null(tau)) {tau <- 0}
    if(is.null(gamma)) {gamma <- 0}
    if(is.null(eta)) {eta <- 0}
    if(is.null(zeta)) {zeta <- 0}
  }
  if (tau==0 & (gamma>0 | eta>0 | zeta>0)){
    stop("If gamma, eta or zeta greater than 0 then tau must be greater than 0")
  }
  if (is.null(rho)){
    rho <- 0 
    if (eta*gamma != 0){
      if (silent==FALSE) warning("The correlation between the random cluster effect and the random treatment effect (rho) is set to 0.")
    }
  }
  if (eta*gamma == 0 & (!is.null(rho)) & (rho!=0)){
    if (silent==FALSE) warning("The correlation between the random cluster effect and the random treatment effect (rho) is not used.")
  }

### Key functions 
  varfun <- function(w,clustersize,Ntrt,treatment.temp,Ntime,time.temp.mtrx,time.temp.mtrx.z,tau,eta,gamma,rho,zeta)
  {
    # Calculations for all clusters in a given sequence
    # Build (most of) D
    dimD <- 1*as.integer(tau>0)+Ntime*as.integer(gamma>0)+Ntrt*as.integer(eta>0)+1*as.integer(zeta>0)
    if (dimD>0){
     D <- matrix(0,dimD,dimD)
     ind=0
     if (tau>0){
       D[1,1] = tau*tau
       ind=1
     }
     if (gamma>0){
       D[(ind+1):(ind+Ntime),(ind+1):(ind+Ntime)] = diag(gamma*gamma,Ntime)
       ind=ind+Ntime
     }
     if (eta>0){
       D[(ind+1):(ind+Ntrt), (ind+1):(ind+Ntrt)] = diag(eta*eta,Ntrt)
       D[1,(ind+1):(ind+Ntrt)] = D[(ind+1):(ind+Ntrt),1] = tau*eta*rho
       ind=ind+Ntrt
     }
     if (zeta>0){
       D[(ind+1),(ind+1)] = zeta*zeta
     }
     D.inv = solve(D)
    }
    # 

    info.seq <- 0
    # loop over all clusters in this sequence
    for(i in 1:nrow(clustersize))
    {
      xTx <- xTz <- zTz <- 0
      for(j in 1:Ntime)
      {
      if (clustersize[i,j]>0){
# this removes all NA periods also since they were given size = 0
        kx.j <-  as.matrix(c(time.temp.mtrx[j,],treatment.temp[j,]))
        xTx = xTx + clustersize[i,j]*(kx.j%*%t(kx.j))/w[j]

        if (dimD>0){
            # Build other matricies 
              temp.j <- 1   # assume tau>0 if dimD >0
              if (gamma>0) temp.j <- c(temp.j,time.temp.mtrx.z[j,])
              if (eta>0) temp.j <- c(temp.j,as.integer(treatment.temp[j,]>0))
              if (zeta>0) temp.j <- c(temp.j,1/sqrt(clustersize[i,j])) 
              kz.j <-  as.matrix(temp.j)
              xTz <- xTz + clustersize[i,j]*(kx.j %*% t(kz.j))/w[j]
              zTz <- zTz + clustersize[i,j]*(kz.j %*% t(kz.j))/w[j]
              }
        }
       }

       if (dimD==0){
          increase <- xTx
          } else {       
          increase <- (xTx - xTz%*%solve(D.inv + zTz)%*%t(xTz))
          } 
      
       info.seq <- info.seq + increase
    }
    info.seq
  }
  
  power <- function(var.alt,var.null,betap,alpha)
  {
    p1 <- pnorm((betap-qnorm(1-alpha/2)*sqrt(var.null))/sqrt(var.alt))
    p2 <- pnorm((betap+qnorm(1-alpha/2)*sqrt(var.null))/sqrt(var.alt))
    return(p1 + 1 - p2)
  }

### Setup
  #design related
  Nseq <- design$n.waves
  cldistninseq <- design$clusters
  design.mtrx <- design$swDsn.unique.clusters
  Ntime <- design$total.time
  time.mt <- matrix(rep(1:Ntime,each=Nseq),ncol=Ntime)
  Ntrt = length(fixed.treatment.effect)
  TxLev = design$TxLev

  #check if n is correctly specified  
  if (length(n) == 1) 
  {
    size.matrix <- matrix(n,nrow = sum(cldistninseq),ncol = Ntime)*(!is.na(design$swDsn))
  }
  else 
  {
    if ((is.vector(n) & length(n) > 1)) 
    {
      if (length(n) != design$n.clusters) 
      {
        stop("The number of clusters in 'design' (design$n.clusters) and 'n' (length(n)) do not match.")
      }
      size.matrix <- matrix(rep(n,Ntime),ncol = Ntime)*(!is.na(design$swDsn))
    }
    else if (is.matrix(n)) 
    {
      if (zeta > 0){
        for (i in 1:nrow(n)){
          if (var(n[i,n[i,]>0])>0) stop("For closed cohort design (zeta > 0) sample size must be constant across time for each cluster (exception: sample size can be 0 to denote time periods which will not be included in the analysis)")
        }
      }     
      if ((nrow(n) != design$n.clusters) | (ncol(n) != design$total.time)) 
      {
        stop("The number of clusters and/or time steps in 'design' (design$n.clusters and design$total.time) and 'n' (number of rows and columns, respectively) do not match.")
      }
      size.matrix <- n*(!is.na(design$swDsn))
    }
  }
  
  ##fixed.effects
  #a vector of fixed effect under the alternative. c(fixed.intercept, fixed.time effect, fixed.treatment effect)
  if(length(fixed.time.effect)==1)
  {
    fixed.time.effect <- rep(fixed.time.effect,Ntime - 1)
  }
  else
  {
    if((is.vector(fixed.time.effect)) & length(fixed.time.effect) > 1)
    {
      if(length(fixed.time.effect) != Ntime - 1)
      {
        stop(paste0("Time effect must be specified for time 2 through time", Ntime, ", as it is coded as dummy variables with time 1 being the reference group."))
      }
    }
  }
  fixed.effects <- c(fixed.intercept,fixed.time.effect,fixed.treatment.effect)
  beta0.vector <- beta1.vector <- fixed.effects
  beta0.vector[(Ntime+1):(Ntime+Ntrt)] <- 0  #Null

  ### Calculate power by looping through the sequences  
  #approximation
  Info1 <- Info0 <- 0
  for(seq in 1:Nseq)
  {
    if(seq ==1)
    {
      cluster.index <- 1: cldistninseq[seq]
    }
    if(seq > 1)
    {
      cluster.index <- (sum(cldistninseq[1:(seq-1)])+1):sum(cldistninseq[1:(seq)])
    }
    
    clustersize <- matrix(size.matrix[cluster.index,],ncol=Ntime)
    if (is.null(H)){
      if (Ntrt==1) {
  # IT model with possible fractional treatment effects
        treatment.temp <- as.matrix(design.mtrx[seq,])
        colnames(treatment.temp) <- "treatment.1"
      } else {
  # multi-level IT model (no fractional treatment effects)
        treatment.temp <- matrix(0,length(design.mtrx[seq,]),Ntrt)
        colnames(treatment.temp) <- paste0("treatment.",TxLev)  
        for (k in 1:Ntrt){
          treatment.temp[,k] <- as.numeric(design.mtrx[seq,] == TxLev[k])
        }
      }
    } else {
      # ETI model (no fractional treatment effects)
      # assume no NA's in the middle of intervention periods (at end of intervention periods is okay)
      tmp = cumsum(replaceNAwithZero(design.mtrx[seq,]))
      treatment.temp <- matrix(0,length(tmp),Ntrt)
      for (k in 1:Ntrt){ treatment.temp[,k] = as.integer(tmp==k) }
      colnames(treatment.temp) <- paste0("treatment.",1:Ntrt)
    }
# note that treatment.temp.z = treatment.temp, so use treatment.temp for both    
    time.temp <- time.mt[seq,]
    time.temp.mtrx <- model.matrix(~factor(time.temp))
    time.temp.mtrx.z <- diag(1,nrow = Ntime,ncol = Ntime)
    
    X.seq.i  <- cbind(time.temp.mtrx, treatment.temp)
    mu.null.i <- h(X.seq.i  %*% beta0.vector)
    mu.alt.i <- h(X.seq.i  %*% beta1.vector)
    
    w.null.i <- (phi*a*v(mu.null.i)*(gprime(mu.null.i))^2)
    w.alt.i <- (phi*a*v(mu.alt.i)*(gprime(mu.alt.i))^2)
   
    Info0 <- Info0 +  varfun(w.null.i,clustersize,Ntrt,treatment.temp,Ntime,time.temp.mtrx,time.temp.mtrx.z,tau,eta,gamma,rho,zeta)
    Info1 <- Info1 +  varfun(w.alt.i,clustersize,Ntrt,treatment.temp,Ntime,time.temp.mtrx,time.temp.mtrx.z,tau,eta,gamma,rho,zeta)
   }
  if (is.null(H)){
    #IT model
    beta_x1 = fixed.treatment.effect
    var.theta.alt <- diag(solve(Info1)[(Ntime+1):(Ntime+Ntrt),(Ntime+1):(Ntime+Ntrt),drop=FALSE])
    var.theta.null <- diag(solve(Info0)[(Ntime+1):(Ntime+Ntrt),(Ntime+1):(Ntime+Ntrt),drop=FALSE])
  } else {
    #ETI model
    beta_x1 = sum(H*fixed.treatment.effect)
    H = as.matrix(H,Ntrt,1)
    var.theta.alt <- t(H)%*%solve(Info1)[(Ntime+1):(Ntime+Ntrt),(Ntime+1):(Ntime+Ntrt)]%*%H
    var.theta.null <- t(H)%*%solve(Info0)[(Ntime+1):(Ntime+Ntrt),(Ntime+1):(Ntime+Ntrt)]%*%H
  }

  pwrGLM <- power(var.theta.alt,var.theta.null,beta_x1,alpha)
  
  
  if(retDATA==TRUE)
  {
    res.list <- list(design, distn,n,fixed.intercept,fixed.treatment.effect,fixed.time.effect,H,tau,eta,rho,gamma,alpha,
                     var.theta.null,var.theta.alt,pwrGLM)
    names(res.list) <- c("design","distn","n","fixed.intercept","fixed.treatment.effect","fixed.time.effect","H",
                         "tau","eta","rho","gamma","alpha",
                         "var.theta.null","var.theta.alt","pwrGLM")
    return(res.list)
  }
  
  if(retDATA==FALSE)
  {
    return(pwrGLM)
  }
  
}
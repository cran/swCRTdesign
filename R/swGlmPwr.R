swGlmPwr <- function(design,distn,n,fixed.intercept,fixed.treatment.effect,fixed.time.effect,tau=0,eta=0,rho=0,gamma=0,alpha=0.05,retDATA=FALSE)
{
  #help functions:
  logit <- function(x){log(x/(1 - x))}
  expit <- function(x){exp(x)/(1 + exp(x))}
  if (alpha <= 0 | alpha >= 1) {
    stop("Alpha must be strictly between 0 and 1.")
  }
  if (!all(n%%1 == 0)) {
    warning("n (either scalar, vector, or matrix) should consist only of integers.")
  }
  if (rho < -1 | rho > 1) {
    stop("The correlation between the random cluster effect and the random treatment effect must be a numeral between -1 and 1.")
  }
  if (tau < 0 | eta < 0 | eta < 0) {
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
    stop("Parameters for fixed effects (intercept, treatment effect, and time effect) 
           must be specified.")
  }
  if(is.null(tau) | is.null(gamma) | is.null(eta))
  {
    warning("Standard deviations of all random effects are needed. If a random effect does not exist, its standard deviation should be set to 0. The default standard deviations are set to 0. ")
    if(is.null(tau)) {tau <- 0}
    if(is.null(gamma)) {gamma <- 0}
    if(is.null(eta)) {eta <- 0}
  }
  if(is.null(rho))
  {
    rho <- 0 
    if(eta*gamma != 0)
    {
      warning("The correlation between the random cluster effect and the random treatment effect (rho) is set to 0.")
    }
  }
  if(eta*gamma == 0 & (!is.null(rho)) & (rho!=0))
  {
    warning("The correlation between the random cluster effect and the random treatment effect (rho) is not used.")
  }
  
  power <- function(var.alt,var.null,betap,alpha)
  {
    p1 <- pnorm((betap-qnorm(1-alpha/2)*sqrt(var.null))/sqrt(var.alt))
    p2 <- pnorm((betap+qnorm(1-alpha/2)*sqrt(var.null))/sqrt(var.alt))
    return(p1 + 1 - p2)
  }
  
  #design related
  Nseq <- design$n.waves
  cldistninseq <- design$clusters
  design.mtrx <- design$swDsn.unique.clusters
  Ntime <- design$total.time
  time.mt <- matrix(rep(1:Ntime,each=Nseq),ncol=Ntime)
  
  #check if n is correctly specified  
  if (length(n) == 1) 
  {
    size.matrix <- matrix(n,nrow = sum(cldistninseq),ncol = Ntime)
  }
  else 
  {
    if ((is.vector(n) & length(n) > 1)) 
    {
      if (length(n) != design$n.clusters) 
      {
        stop("The number of clusters in 'design' (design$n.clusters) and 'n' (length(n)) do not match.")
      }
      size.matrix <- matrix(rep(n,Ntime),ncol = Ntime)
    }
    else if (is.matrix(n)) 
    {
      if ((nrow(n) != design$n.clusters) | (ncol(n) != design$total.time)) 
      {
        stop("The number of clusters and/or time steps in 'design' (design$n.clusters and design$total.time) and 'n' (number of rows and columns, respectively) do not match.")
      }
      size.matrix <- n
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
  fixed.effects <- c(fixed.intercept, fixed.time.effect,fixed.treatment.effect)
  beta0.vector <- beta1.vector <- fixed.effects
  beta_x1 <- beta1.vector[length(beta0.vector)]
  beta0.vector[length(beta0.vector)] <- 0
  
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
    treatment.temp <- as.matrix(design.mtrx[seq,])
    colnames(treatment.temp) <- "treatment.i"
    time.temp <- time.mt[seq,]
    time.temp.mtrx <- model.matrix(~factor(time.temp))
    time.temp.mtrx.z <- diag(1,nrow = Ntime,ncol = Ntime)
    
    X.seq.i  <- cbind(time.temp.mtrx, treatment.temp)
    mu.null.i <- h(X.seq.i  %*% beta0.vector)
    mu.alt.i <- h(X.seq.i  %*% beta1.vector)
    
    w.null.i <- (phi*a*v(mu.null.i)*(gprime(mu.null.i))^2)
    w.alt.i <- (phi*a*v(mu.alt.i)*(gprime(mu.alt.i))^2)
    
    varfun <- function(type,clustersize,treatment.temp,Ntime,time.temp.mtrx,tau,eta,gamma,rho)
    {
      if(type == "null")
      {
        w <- w.null.i 
      }
      if(type == "alt")
      {
        w <- w.alt.i
      }
      info.seq <- 0
      for(i in 1:nrow(clustersize))
      {
        tr <- treatment.temp
        LxTLx <- LxTtr <- trTtr <-  0
        for(j in 1:Ntime)
        {
          kx.j <-  as.matrix(time.temp.mtrx[j,])
          LxTLx <-  LxTLx + clustersize[i,j]*(kx.j%*%t(kx.j))/w[j]
          LxTtr <- LxTtr + clustersize[i,j]*(kx.j*tr[j])/w[j]
          trTtr <- trTtr + clustersize[i,j]*(tr[j])^2/w[j]
        }
        
        xTx <- rbind(cbind(LxTLx,LxTtr),cbind(t(LxTtr),trTtr))
        
        if(tau==0)
        {
          stop("The cluster random effect cannot be 0.")
        }
        #cluster + treatment 
        if(eta!=0 & gamma==0)
        {
          dimD <- 2
          D.inv <- matrix(0,dimD,dimD)
          D.inv <- (tau*tau*(1 - rho^2))^(-1)*matrix(c(eta^2,-tau*eta*rho,-tau*eta*rho,tau^2),2,2)
          
          LxTLz <- LzTtr <- LzTLz <- 0
          for(j in 1:Ntime)
          {
            kx.j <-  as.matrix(time.temp.mtrx[j,])
            kz.j <-  1
            LxTLz <-  LxTLz + clustersize[i,j]*(kx.j)/w[j]
            LzTtr <- LzTtr + clustersize[i,j]*(tr[j])/w[j]
            LzTLz <- LzTLz + clustersize[i,j]/w[j]
          }
          
          zTz <- rbind(cbind(LzTLz,LzTtr),cbind(t(LzTtr),trTtr))
          xTz <- rbind(cbind(LxTLz,LxTtr),cbind(t(LzTtr),trTtr))
          increase <- (xTx - xTz%*%solve(D.inv + zTz)%*%t(xTz))
        }
        #cluster + time + treatment
        if(eta!=0 & gamma!=0)
        {
          dimD <- 2+Ntime
          D.inv <- matrix(0,dimD,dimD)
          D.inv[1:2,1:2] <- (tau*tau*(1 - rho^2))^(-1)*matrix(c(eta^2,-tau*eta*rho,-tau*eta*rho,tau^2),2,2)
          Dtime.inv <- diag(rep(gamma^(-2),Ntime))
          D.inv[3:dimD,3:dimD] <- Dtime.inv
          
          LxTLz <- LzTtr <- LzTLz <- 0
          for(j in 1:Ntime)
          {
            kx.j <-  as.matrix(time.temp.mtrx[j,])
            kz.j <-  as.matrix(c(1,time.temp.mtrx.z[j,]))
            LxTLz <-  LxTLz + clustersize[i,j]*(kx.j %*% t(kz.j))/w[j]
            LzTtr <- LzTtr + clustersize[i,j]*(kz.j*tr[j])/w[j]
            LzTLz <- LzTLz + clustersize[i,j]*(kz.j %*% t(kz.j))/w[j]
          }
          
          zTz <- rbind(cbind(LzTLz,LzTtr),cbind(t(LzTtr),trTtr))
          xTz <- rbind(cbind(LxTLz,LxTtr),cbind(t(LzTtr),trTtr))
          increase <- (xTx - xTz%*%solve(D.inv + zTz)%*%t(xTz))
        }
        #cluster + time
        if(eta==0 & gamma!=0)
        {
          dimD <- 1 + Ntime
          D.inv <- matrix(0,dimD,dimD)
          D.inv[1,1] <- tau^(-2)
          Dtime.inv <- diag(rep(gamma^(-2),Ntime))
          D.inv[2:dimD,2:dimD] <- Dtime.inv
          
          zTz <- LxTLz <- LzTtr <- 0
          for(j in 1:Ntime)
          {
            kx.j <-  as.matrix(time.temp.mtrx[j,])
            kz.j <-  as.matrix(c(1,time.temp.mtrx.z[j,]))
            zTz <-  zTz + clustersize[i,j]*(kz.j %*% t(kz.j))/w[j]
            LxTLz <-  LxTLz + clustersize[i,j]*(kx.j %*% t(kz.j))/w[j]
            LzTtr <-  LzTtr + clustersize[i,j]*(kz.j * tr[j])/w[j]
          }
          xTz <- rbind(LxTLz,t(LzTtr))
          increase <- (xTx - xTz%*%solve(D.inv + zTz)%*%t(xTz))
        }
        #only random cluster effect
        if(eta==0 & gamma==0)
        {
          D.inv <- tau^(-2)
          
          zTz <- LxTLz <-LzTtr <-  0
          for(j in 1:Ntime)
          {
            kx.j <-  as.matrix(time.temp.mtrx[j,])
            kz.j <-  1
            zTz <-  zTz + clustersize[i,j]/w[j]
            LxTLz <-  LxTLz + clustersize[i,j]*(kx.j)/w[j]
            LzTtr <-  LzTtr + clustersize[i,j]*(tr[j])/w[j]
          }
          xTz <- rbind(LxTLz,t(LzTtr))
          increase <- (xTx - xTz%*%t(xTz)*(D.inv + zTz)^(-1))
        }
        info.seq <- info.seq + increase
      }
      info.seq
    }
    
    Info0 <- Info0 +   varfun("null",clustersize,treatment.temp,Ntime,time.temp.mtrx,tau,eta,gamma,rho)
    Info1 <- Info1 +   varfun("alt",clustersize,treatment.temp,Ntime,time.temp.mtrx,tau,eta,gamma,rho)
  }
  var.theta.alt <- solve(Info1)[Ntime+1,Ntime+1]
  var.theta.null <- solve(Info0)[Ntime+1,Ntime+1]
  
  pwrGLM <- power(var.theta.alt,var.theta.null,beta_x1,alpha)
  
  
  
  if(retDATA==TRUE)
  {
    res.list <- list(design, distn,n,fixed.intercept,fixed.treatment.effect,fixed.time.effect,tau,eta,rho,gamma,alpha,
                     var.theta.null,var.theta.alt,pwrGLM)
    names(res.list) <- c("design","distn","n","fixed.intercept","fixed.treatment.effect","fixed.time.effect",
                         "tau","eta","rho","gamma","alpha",
                         "var.theta.null","var.theta.alt","pwrGLM")
    return(res.list)
  }
  
  if(retDATA==FALSE)
  {
    return(pwrGLM)
  }
  
}
swPwr <-
function( design, distn, n, mu0, mu1, tau, eta, rho=0, sigma, alpha=0.05, retDATA=FALSE ) {
##
## (v2.0 2015.04.14)
##
## design = design object as from swDsn (must include at least the following components: 
##          clusters, n.clusters, total.time, swDsn, swDsn.unique.clusters)
## distn = distribution ("gaussian", "binomial" allowed)
##     n = number of observations per cluster per time period
##    mu0 = mean in control group
##    mu1 = mean in treatment group##   tau = standard deviation of random intercept
##   eta = standard deviation of random treatment effect
##   rho = correlation between tau and eta
## sigma = residual stadnard deviation (not needed if distn is binomial or possson)
## alpha = statistical significance level
##
    ####
    obs.per.cluster.per.time <- n
    ####
    theta <- (mu1 - mu0)
    muBar <- (mu0 + mu1) / 2
#**#
if(distn=="gaussian") sigSq <- sigma^2
else if(distn=="binomial") {
    if( !missing(sigma) ) warning("sigma is not used when distn=binomial")
    sigSq <- muBar * (1 - muBar); sigma <- NA
    if((tau^2+eta^2)>sigSq) stop("tau^2 + eta^2 must be less than muBar*(1-muBar) when distn=binomial")
#}else if(distn=="poisson") {
#    sigSq <- muBar; sigma <- NA
#    if((tau^2+eta^2)>sigSq) stop("tau^2 + eta^2 must be less than muBar when distn=binomial") 
}
#**#
	####
	## Warning message
	####
	if( length(tau) > 1  |  length(eta) > 1 ) warning("Function cannot compute stepped wedge design Power for tau-vector or eta-vector; tau and eta must be scalars. Ignore any results that follow.")
	####
##
I.rep <- design$clusters
I <- design$n.clusters
J <- design$total.time  
swDesign <- design$swDsn
swDesignUnique <- design$swDsn.unique.clusters
##
	####
	## Warning message:
	####
	if( any(rowSums(swDesign) == 0) ) warning("For the specified total number of clusters (I), total number of time periods (J), and number of cluster repetitions (I.rep), the specified stepped wedge design has at least one cluster which does not crossover from control(0) to treatment(1) arm.")
	####
	## Constructing the Treatment/Intervention Indicator Vector (X.ij)
	####
	X.ij <- as.vector( t(swDesign) )
	####
	## Constructing the Design Matrix (Xmat)
	####
	beta.blk <- rbind(diag(1,J-1,J-1), 0)
	Xmat.blk <- matrix( rep( as.vector(t(cbind(1,beta.blk))), I ), ncol=J, byrow=TRUE )	
	Xmat <- cbind(Xmat.blk, X.ij)
	####
	## Constructing the Covariance Matrix (Wmat)
	####
  Wmat.blk <- tau^2 + diag(sigSq/n, J)
	Wmat.partial <- kronecker( diag(1,I), Wmat.blk )
    #
    Xij.Xil.ARRAY <- array(NA,c(J,J,length(I.rep))) 
    for(i.indx in 1:length(I.rep)) {
				Xi.jl <- swDesignUnique[i.indx,]
		    Xij.Xil.ARRAY[,,i.indx] <- (Xi.jl %o% Xi.jl)*eta^2 + outer(Xi.jl, Xi.jl, "+")*rho*eta*tau 
    }
    Xij.Xil.LIST_blk <- sapply(1:length(I.rep), function(x) NULL) 
    for( i.indx in 1:length(I.rep) ) {
	    Xij.Xil.LIST_blk[[i.indx]] <- kronecker( diag(1,I.rep[i.indx]), Xij.Xil.ARRAY[,,i.indx] )
    }
    W.eta <- blkDiag(Xij.Xil.LIST_blk)
    Wmat <- Wmat.partial + as.matrix(W.eta)    
	var.theta.WLS <- solve( t(Xmat)%*%solve(Wmat)%*%Xmat )["X.ij", "X.ij"]
	pwrWLS <- pnorm( abs(theta)/sqrt( var.theta.WLS ) - qnorm(1-alpha/2) ) + pnorm( -abs(theta)/sqrt( var.theta.WLS ) - qnorm(1-alpha/2) )
	rslt <- pwrWLS
	####
	## Closed-form approach/solution
	##   eta==0 (i.e., *NO* random slopes)
	#### 
	if(eta==0){
		####
		## Closed-form formula
		####
		X <- swDesign ###round( upper.tri(matrix(1,I,J)))
		U <- sum(X)
		W <- sum( colSums(X)^2 )
		V <- sum( rowSums(X)^2 )
		####
		## 1. Obtain Variance Estimate of theta from Closed-form (var.theta.CLOSED)
		## 2. Calculate Power for theta from Closed-form (pwrCLOSED)
		## 3. Storing/Appending Resulting Power for theta from Closed-form (rslt)
		####
    sigSq <- sigSq / n
		var.theta.CLOSED <- I*sigSq*(sigSq + J*tau^2) / ( (I*U - W)*sigSq + (U^2 + I*J*U - J*W - I*V)*tau^2 )
		pwrCLOSED <- pnorm( abs(theta)/sqrt(var.theta.CLOSED) - qnorm(1-alpha/2) ) + pnorm( -abs(theta)/sqrt(var.theta.CLOSED) - qnorm(1-alpha/2) )
		####
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
		########
	} else {
    pwrCLOSED <- NA
	}
	####
	## Returning Resulting Power(s) for fixed theta of the specified SW design
	####
	if(retDATA) rslt <- list(design=design, n=n, mu0=mu0, mu1=mu1, tau=tau, 
                           eta=eta, rho=rho, sigma=sigma, alpha=alpha, 
                           Xmat=Xmat, Wmat=Wmat, 
                           var.theta.WLS=var.theta.WLS, pwrWLS=pwrWLS, 
                           pwrCLOSED=pwrCLOSED )
	rslt
}

swSimPwr <- function (design, family, n, mu0, mu1, time.effect=0, H = NULL, sigma, tau=0,
                      eta=0, rho=0, gamma=0, zeta=0, icc=0, cac=1, iac=0, alpha = 0.05, nsim = 500, retDATA = FALSE, silent = FALSE, counter=TRUE)
{
#
#Last update: 5/22/23, v. 4.0, Jim Hughes
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
  ## taken from beginning of glmer() in lme4 package
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
  #Basic input checks (that are not part of swSim already)
  if(alpha <= 0 | alpha >= 1){
    stop("Alpha must be strictly between 0 and 1.")
  }
  if(!is.null(H)){
    if (design$nTxLev > 1){
      stop("Do not use H when you have specfied a design with multiple treatment levels. If you want an ETI model based estimator, specify a 0/1 intervention in swDsn and use H to define the estimator.")
    }
    if (sum(H) != 1) {
      H = H/sum(H)
      if (silent==FALSE) warning("H has been renormalized to sum to 1.0")
    }
    Hext = matrix(c(0,H,rep(0,design$total.time-1)),ncol=1)
  }
  #Basic input checks that are normally part of swSim, but put them here and 
  #set nocheck=TRUE in swSim call so they are only run once
  log.gaussian=FALSE
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
  frac = (any(design$tx.effect.frac < 1))
  if ((!is.null(H)|design$nTxLev>1)&frac) stop("Do not specify time-varying treatment effect or multiple treatment effects with fractional treatment effects")
  #############
  model = as.formula(paste0("response.var ~ ",
                            ifelse(is.null(H),"ftx.var","ftimeontx.var"),
                            " + ftime.var",
                            ifelse(tau>0&(eta==0|(eta>0&rho==0))," + (1|fcluster.var)",""),
                            ifelse(gamma>0," + (1|clustime)",""),
                            ifelse(eta>0&rho!=0," + (itx.var | fcluster.var)",""),
                            ifelse(eta>0&rho==0," + (0 + itx.var | fcluster.var)",""),
                            ifelse(zeta>0,"+ (1|findivid.var)",""))) 
  #############
  power=NULL
  pval=est=se=df=matrix(NA,nsim,design$nTxLev)
  nvcomp = (tau>0)+(eta>0)+(rho!=0)+(gamma>0)+(zeta>0)+as.numeric(distn == 'gaussian')
  varcomp = matrix(NA,nsim,nvcomp)
  code=rep(0,nsim)
  msg=rep("",nsim)
  for (i in 1:nsim){
    if (counter) cat(i)
    data = swSim(design=design,family=family,n=n,mu0=mu0,mu1=mu1,
                 time.effect=time.effect, sigma=sigma, tau=tau,
                 eta=eta, rho=rho, gamma=gamma, zeta=zeta, 
                 time.lab = NULL, silent=silent,nocheck=TRUE)
    #
    data$ftime.var = as.factor(data$time.var)
    data$fcluster.var = as.factor(data$cluster.var)
    data$clustime = interaction(data$fcluster,data$ftime)
    if (zeta>0) data$findivid.var = as.factor(data$individ.var)
    if (frac) {data$ftx.var = data$tx.var} else {data$ftx.var = as.factor(data$tx.var)}
    data$itx.var = as.integer(data$tx.var>0)
    data$ftimeontx.var = as.factor(data$timeOnTx.var)
    #
    if (distn == 'gaussian'){
      rslt = suppressMessages(lmerTest::lmer(model,data=data))
      if (length(rslt@optinfo$conv$lme4$messages)>0)  msg[i] = paste0("[1] ",rslt@optinfo$conv$lme4$messages[1])
      if (length(rslt@optinfo$conv$lme4$messages)>1) {
        for (k in 2:length(rslt@optinfo$conv$lme4$messages)){
          msg[i] = paste0(msg[i]," [",k,"] ",rslt@optinfo$conv$lme4$messages[k])
        }}
      if (is.null(H)){
        for (j in 1:design$nTxLev){
          if (frac) {varname="ftx.var"} else {varname=paste0("ftx.var",j)}
          est[i,j] = summary(rslt)$coef[varname,"Estimate"]
          se[i,j] = summary(rslt)$coef[varname,"Std. Error"]
          df[i,j] = summary(rslt)$coef[varname,"df"]
          pval[i,j] = summary(rslt)$coef[varname,"Pr(>|t|)"]
        }
      } else {
        contrast = lmerTest::contest(rslt,as.vector(Hext))
        est[i,1] = sum(lme4::fixef(rslt)*Hext)
        se[i,1] = as.numeric(sqrt(t(Hext)%*%vcov(rslt)%*%Hext))
        df[i,1] = as.numeric(contrast["DenDF"])
        pval[i,1] = as.numeric(contrast["Pr(>F)"])
      }
    } else {
    if (distn =="binomial" | distn == "poisson") {
       rslt = suppressMessages(lme4::glmer(model,data=data,family=family))
       if (length(rslt@optinfo$conv$lme4$messages)>0)  msg[i] = paste0("[1] ",rslt@optinfo$conv$lme4$messages[1])
       if (length(rslt@optinfo$conv$lme4$messages)>1) {
         for (k in 2:length(rslt@optinfo$conv$lme4$messages)){
           msg[i] = paste0(msg[i]," [",k,"] ",rslt@optinfo$conv$lme4$messages[k])
         }}
       if (is.null(H)){
         for (j in 1:design$nTxLev){
           if (frac) {varname="ftx.var"} else {varname=paste0("ftx.var",j)}
           est[i,j] = summary(rslt)$coef[varname,"Estimate"]
           se[i,j] = summary(rslt)$coef[varname,"Std. Error"]
           pval[i,j] = summary(rslt)$coef[varname,"Pr(>|z|)"]
         }
       } else {
         est[i,1] = sum(lme4::fixef(rslt)*Hext)
         se[i,1] = as.numeric(sqrt(t(Hext)%*%vcov(rslt)%*%Hext))
         z = est[i,1]/se[i,1]
         pval[i,1] = 2*(1-pnorm(z))
       }
    } else {
    stop("Distribution ",distn," not recognized")
    }}
    varcomp[i,] = as.data.frame(lme4::VarCorr(rslt))$sdcor
    if (!is.null(rslt@optinfo$conv$lme4$code)) code[i]=rslt@optinfo$conv$lme4$code
  }
  colnames(varcomp) = apply(as.data.frame(lme4::VarCorr(rslt))[,1:3],1,paste,collapse="-")
  colnames(varcomp) = gsub("-NA","",colnames(varcomp))
  colnames(varcomp) = sub("Residual","sigma",colnames(varcomp),fixed=TRUE)
  colnames(varcomp) = sub("fcluster.var-(Intercept)-itx.var","rho",colnames(varcomp),fixed=TRUE)
  colnames(varcomp) = sub("clustime-(Intercept)","gamma",colnames(varcomp),fixed=TRUE)
  colnames(varcomp) = sub("clustime-(Intercept)","gamma",colnames(varcomp),fixed=TRUE)
  colnames(varcomp) = sub("fcluster.var-itx.var","eta",colnames(varcomp),fixed=TRUE)
  colnames(varcomp) = sub("fcluster.var.1-itx.var","eta",colnames(varcomp),fixed=TRUE)
  colnames(varcomp) = sub("fcluster.var-(Intercept)","tau",colnames(varcomp),fixed=TRUE)
  colnames(varcomp) = sub("fcluster.var.1-(Intercept)","tau",colnames(varcomp),fixed=TRUE)
  for (j in 1:design$nTxLev){
    power[j] = sum(pval[,j]<=alpha)/nsim
  }
  swDesign=design$swDsn
  TxLev = design$TxLev
  names(power) <- paste0("treatment.",TxLev) 
  cat("\n",sum(code!=0)," of ",nsim," simulations failed to converge (see code)\n")
  cat(sum(msg!="")," of ",nsim," simulations generated warning messages (see msg)\n")
  rslt = power
  if (retDATA) rslt = list(power=power,nsim=nsim,est=est,se=se,sdcor=varcomp,df=df,pvalue=pval,code=code,msg=msg)
  rslt
}
  
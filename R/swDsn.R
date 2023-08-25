swDsn <- function (clusters, tx.effect.frac = 1, extra.ctl.time = 0, extra.trt.time = 0, all.ctl.time0 = TRUE, swBlk, extra.time=NULL)
{
  #Last update: 6/15/2023, v. 4.0, Jim Hughes
  ##########
  #Warnings
  ##########
  #
  if (missing(swBlk)){
  # Original form with extra.ctl.time , extra.trt.time , all.ctl.time0 arguments
  if (!missing(extra.time)) {
    extra.trt.time=extra.time
    warning("extra.time interpreted as extra.trt.time")
  }
  if(extra.ctl.time < 0 | extra.ctl.time%%1 != 0){
    stop("extra.ctl.time must be a non-negative integer.")
  }
  if(extra.trt.time < 0 | extra.trt.time%%1 != 0){
    stop("extra.trt.time must be a non-negative integer.")
  }
  if(! all(clusters%%1 == 0)){
    stop("'clusters' must be a vector of integers.")
  }
  if(! all(tx.effect.frac >= 0 & tx.effect.frac <= 1)){
      #Removed condition below because I eliminated within-package use of tx.effect.frac.  If it were added back in, would recommend adding an argument for within-package use.
      #if(length(tx.effect.frac) != (length(clusters)+extra.time+as.numeric(all.ctl.time0)) || !all(tx.effect.frac == (1:(length(clusters)+extra.time+as.numeric(all.ctl.time0))))){#The extra condition was because swSim uses swDsn's fractional treatment effect to build a matrix of time on treatment.
      warning("Typically, the fractional treatment effect should be between 0 and 1.  Note that the 'fractional treatment effect' is different from the 'treatment effect'; see documentation.")
  }#}  
   ##########
  #Make design matrix
  ##########
  #Note: 'wave' and 'sequence' have the same meaning; code tends to use wave, while documentation tends to use sequence.
  clusters.per.wave <- clusters[clusters != 0]
  total.clusters <- sum(clusters.per.wave)
  waves <- length(clusters.per.wave)
  lclusters <- length(clusters)
  total.time <- (lclusters + 1 + extra.ctl.time + extra.trt.time)#The extra 1 is for all.ctl.time0; will be adjusted later if set to FALSE
  nTxLev = 1
  TxLev = as.vector(1)
  #swBlk <- round(upper.tri(matrix(1, length(clusters), total.time)))#Preliminary design matrix with one row per sequence
  swBlk <- cbind(matrix(0,lclusters,extra.ctl.time),round(upper.tri(matrix(1, lclusters, lclusters))),matrix(1,lclusters,extra.trt.time+1))
  swPreDesign.MAT <- cbind(clusters, swBlk)
  #Expand swBlk to represent individual clusters instead of sequences
  swPreDesign.LIST <- apply(swPreDesign.MAT, 1, function(z) {
    rep(z[-1], z[1])
  })
  #Design matrix with one row per cluster, assuming all clusters start on control
  swDsn <- matrix(unlist(swPreDesign.LIST), ncol = total.time,
                  byrow = TRUE)
  ##########
  #Adjust design matrix for non-standard designs
  #Future work: some of this could be easily incorporated into the section above instead of changing things afterwards
  ##########
  #For designs with no all-control time, remove that baseline time from matrices and total.time.
  if (!all.ctl.time0) {
    swBlk <- swBlk[, -1]
    swDsn <- swDsn[, -1]
    #Following was removed because extra time doesn't affect waves, and the clusters.per.wave change was incorrect; errors may have affected versions <= 2.1.
    #clusters.per.wave <- clusters.per.wave[-c(waves)]
    #total.clusters <- sum(clusters.per.wave)
    #waves <- length(clusters.per.wave)
    total.time <- dim(swDsn)[2]
  }
  #For designs with a treatment effect, replace the design matrices with ones that have the fractional treatment effect instead of just 0s and 1s
  if (!is.null(tx.effect.frac)) {
    swTxEffectPreDesign.LIST <- apply(swDsn, 1, function(z) {
      Nctl <- sum(z == 0)
      rowWithTxEffect <- c(rep(0, Nctl), tx.effect.frac, rep(1,
                                                             length(z)))[1:length(z)]
      rowWithTxEffect
    })
    swTxEffectDsn <- matrix(unlist(swTxEffectPreDesign.LIST),
                            ncol = total.time, byrow = TRUE)
    swDsn <- swTxEffectDsn
    swBlk <- unique(swDsn)
  }
  } else {
  # New form with swBlk argument
    if (any(tx.effect.frac!=1)) stop("tx.effect.frac not available with swBlk argument")
    if (!missing(extra.ctl.time) | !missing(extra.trt.time) | !missing(all.ctl.time0) | !missing(extra.time)) stop("Use swBlk argument OR extra.trt.time, extra.ctl.time, all.ctl.time0 arguments, not both")
    ##########
    #Warnings
    ##########
    if(! all(clusters%%1 == 0) | any(clusters<=0)){
      stop("'clusters' must be a vector of non-zero positive integers.")
    }
    if (!is.numeric(swBlk) | !is.matrix(swBlk) | ! all(swBlk%%1 == 0 | is.na(swBlk))){
      stop("swBlk must be a matrix of integers (NA's allowed)")
    }
    if (length(clusters) != dim(swBlk)[1]) {
      stop("length of 'clusters' must be equal to number of rows of 'swBlk'")
    }
    # drop any column of swBlk that is all NA
    swBlk = swBlk[,!(apply(is.na(swBlk),2,sum)==dim(swBlk)[1])]
    ##########
    #Make design matrix
    ##########
    #Note: 'wave' and 'sequence' have the same meaning; code tends to use wave, while documentation tends to use sequence.
    clusters.per.wave <- clusters
    total.clusters <- sum(clusters.per.wave)
    waves <- dim(swBlk)[1]
    total.time <- dim(swBlk)[2]
    nTxLev = length(unique(swBlk[!is.na(swBlk)])) - 1
    TxLev = unique(swBlk[!is.na(swBlk)&swBlk!=0])
    swPreDesign.MAT <- cbind(clusters, swBlk)
    #Expand swBlk to represent individual clusters instead of sequences
    swPreDesign.LIST <- apply(swPreDesign.MAT, 1, function(z) {
      rep(z[-1], z[1])
    })
    #Design matrix with one row per cluster
    swDsn <- matrix(unlist(swPreDesign.LIST), ncol = total.time,
                    byrow = TRUE)
    }
  ##########
  #Return results
  ##########
  list(swDsn = swDsn, swDsn.unique.clusters = swBlk, n.waves = waves,
       clusters = clusters.per.wave, n.clusters = total.clusters,
       tx.effect.frac = tx.effect.frac, total.time = total.time, nTxLev = nTxLev, TxLev = TxLev)
}

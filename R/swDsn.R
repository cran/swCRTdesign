swDsn <-
function (clusters, tx.effect = 1, extra.time = 0, all.ctl.time0 = TRUE) 
{
## (v.2014.08.06)
	##
    clusters.per.wave <- clusters
    total.clusters <- sum(clusters.per.wave)
    waves <- length(clusters.per.wave)
    ##
    total.time <- (waves + 1 + extra.time)
    ##
    swBlk <- round(upper.tri(matrix(1, waves, total.time)))
    swPreDesign.MAT <- cbind(clusters.per.wave, swBlk)
    swPreDesign.LIST <- apply(swPreDesign.MAT, 1, function(z) {
        rep(z[-1], z[1])
    })
    swDsn <- matrix(unlist(swPreDesign.LIST), ncol = total.time, 
        byrow = TRUE)
    ##
    if (!all.ctl.time0) {
        swBlk <- swBlk[, -1]
        swDsn <- swDsn[, -1]
        clusters.per.wave <- clusters.per.wave[-c(waves)]
        total.clusters <- sum(clusters.per.wave)
        waves <- length(clusters.per.wave)
        total.time <- dim(swDsn)[2]
    }
    ##
    if( !is.null(tx.effect) ) {
    	swTxEffectPreDesign.LIST <- apply( swDsn, 1, function(z){ 
    		Nctl <- sum(z==0)
    		rowWithTxEffect <- c( rep(0,Nctl), tx.effect, 
    						rep(1,length(z)) )[1:length(z)]
    		rowWithTxEffect
    	 })
    	swTxEffectDsn <- matrix(unlist(swTxEffectPreDesign.LIST), ncol = total.time, 
        byrow = TRUE)
        swDsn <- swTxEffectDsn
        ## fixing unique swDsn
        swBlk <- unique(swDsn)
    }
    ##
    list(swDsn = swDsn, swDsn.unique.clusters = swBlk, n.waves = waves, 
        clusters = clusters.per.wave, n.clusters = total.clusters, 
        tx.effect = tx.effect, total.time = total.time, extra.time = extra.time)
}

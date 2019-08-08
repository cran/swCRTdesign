swPlot <- function(response.var, tx.var, time.var, cluster.var, data, choose.mfrow=NULL, by.wave=TRUE, combined.plot=TRUE, choose.xlab="Time", choose.main=NULL, choose.pch=NULL, choose.cex=1, choose.tx.col=NULL, choose.ncol=2, choose.tx.pos="topright", choose.legend.pos="right")
{

## (v.2014.04.29); updated 8/5/2019 for v. 3.0
swSummaryFcnCall <- swSummary(response.var=substitute(response.var), tx.var=substitute(tx.var), time.var=substitute(time.var), cluster.var=substitute(cluster.var), data=data, type="mean", digits=16, fcn.Call=TRUE)
swDsnTmp <- swSummaryFcnCall$swDsn
##
if( is.null(choose.tx.col) ) {
	## This line of could should be the same as the following line; however, I created 'swDsnTmp' for most of the remaining code ----> 	uniqueTx <- sort( unique(as.vector(swSummaryFcnCall$swDsn)), decreasing=FALSE )
	##
	uniqueTx <- sort(unique( as.vector(swDsnTmp) ))
	##
	colorTx0 <- "blue"
	colorTx1 <- "red"
	colorTx0.1 <- c("green3", "orange", "brown")
	colorTx1.plus <- c("purple", "lavender")
	##
	Tx0 <- which(uniqueTx > -1e-8 & uniqueTx < 1e-8)
	Tx1 <- which(uniqueTx > 1 - 1e-8 & uniqueTx < 1 + 1e-8)
	Tx0.1 <- which(uniqueTx > 0 & uniqueTx < 1)
	Tx1.plus <- which(uniqueTx > 1)
	##
	uniqueTxColors <- rep(NA, length(uniqueTx))
	uniqueTxColors[Tx0] <- colorTx0
	uniqueTxColors[Tx1] <- colorTx1
	uniqueTxColors[Tx0.1] <- colorTx0.1[1:length(Tx0.1)]
	uniqueTxColors[Tx1.plus] <- colorTx1.plus[1:length(Tx1.plus)]
}
###*********************************************
###
### (by.wave == TRUE)
###
###*********************************************
if( by.wave==TRUE ) {
	swMeanResponse.Wave <- swSummaryFcnCall$response.wave
	swMeanResponse.Cluster <- swSummaryFcnCall$response.cluster
	sw.xMin.Wave <- min( swSummaryFcnCall$time.at.each.wave )
	sw.xMax.Wave <- max( swSummaryFcnCall$time.at.each.wave )
	sw.yMin.Wave <- min( swMeanResponse.Wave )
	sw.yMax.Wave <- max( swMeanResponse.Wave )
#*# want to figure out how to specify xlim and ylim, maybe ylim=NULL as default
	###*********************************************
	###
	### (by.wave == TRUE & combined.plot == TRUE)
	###
	###*********************************************
	if( combined.plot==TRUE ) {
		par( mfrow=c(1,1) )
#*#choose.main <- s...............
		plot(0,0, type="n", xlim=c(sw.xMin.Wave, sw.xMax.Wave), ylim=c(sw.yMin.Wave, sw.yMax.Wave), xaxt="n", xlab=choose.xlab, ylab="Mean Response", main=choose.main)
		axis(1, at=c(swSummaryFcnCall$time.at.each.wave), labels=c(swSummaryFcnCall$time.at.each.wave))
		for( i in 1:dim(swMeanResponse.Wave)[1] ){
			lines( swSummaryFcnCall$time.at.each.wave, swMeanResponse.Wave[i,], lty=1, col="grey" )
		}
		if( is.null(choose.pch) ) choose.pch <- c( 0:(swSummaryFcnCall$n.waves - 1) )
#*# need warning message if choose.pch-vector is NOT of length equal to total number of waves, then need to fix issue.....
###**************************************************************************************************************************
###**************************************************************************************************************************
### ADDITION (wave): Allows different Tx types.......(e.g., 0, .8, .9, 1, 2)
###
		for( j in 1:length(uniqueTx) ){
			swDsnTmp.j <- swDsnTmp
			swDsnTmp.j[ (swDsnTmp.j < uniqueTx[j] - 1e-8) | (swDsnTmp > uniqueTx[j] + 1e-8) ] <- NA
			swDsnTmp.j[ !is.na(swDsnTmp.j) ] <- 1
			unique(swDsnTmp)
			swPlotTx.j.Wave <- unique(cbind(swDsnTmp.j, swDsnTmp))[,1:dim(swDsnTmp.j)[2]] * swMeanResponse.Wave
			if( sum(swDsnTmp.j, na.rm=TRUE)==0 ) swPlotTx.j.Wave <- (1 - unique(cbind(swDsnTmp.j, swDsnTmp))[,1:dim(swDsnTmp.j)[2]]) * swMeanResponse.Wave
			for( i in 1:dim( unique(swDsnTmp) )[1] ){
				points( swSummaryFcnCall$time.at.each.wave, swPlotTx.j.Wave[i,], pch=choose.pch[i], cex=choose.cex, col=uniqueTxColors[j] )
			}
		}
###**************************************************************************************************************************
###**************************************************************************************************************************
##
		legend( choose.tx.pos, c("Tx:", uniqueTx), text.col=c("black", uniqueTxColors), ncol=(1 + length(uniqueTx)), x.intersp=0.01, xjust=1, y.intersp=0.5, cex=0.9)
		##
		if(!is.null(choose.legend.pos) && choose.legend.pos == "mouseclick"){
		  cat("Click on plot to add a legend\n")
		  #utils::flush.console()#If there are difficulties with printing the message at the appropriate time, this may be turned on.
		  legend( locator(1), legend=1:swSummaryFcnCall$n.waves, col="black", pch=choose.pch, ncol=choose.ncol, title="Wave" )
		}
		else if(!is.null(choose.legend.pos)) legend( choose.legend.pos, legend=1:swSummaryFcnCall$n.waves, col="black", pch=choose.pch, ncol=choose.ncol, title="Wave" )
	}
	###*********************************************
	###
	### (by.wave == TRUE & combined.plot == FALSE)
	###
	###*********************************************
	else if( combined.plot==FALSE ){
		##
		if( is.null(choose.mfrow) ) {
			if( swSummaryFcnCall$n.waves==1 ){ choose.mfrow <- c(1,1) }
			else if( swSummaryFcnCall$n.waves==2 ){ choose.mfrow <- c(1,2) }
			else if( swSummaryFcnCall$n.waves==3 | swSummaryFcnCall$n.waves==4 ){ choose.mfrow <- c(2,2) }
			else if( swSummaryFcnCall$n.waves==5 | swSummaryFcnCall$n.waves==6 ){ choose.mfrow <- c(2,3) }
			else if( swSummaryFcnCall$n.waves==7 | swSummaryFcnCall$n.waves==8 ){ choose.mfrow <- c(2,4) }
			else if( swSummaryFcnCall$n.waves==9 ){ choose.mfrow <- c(3,3) }
			else if( swSummaryFcnCall$n.waves>=10 | swSummaryFcnCall$n.waves<=12 ){ choose.mfrow <- c(3,4) }
		}
		par( mfrow=choose.mfrow )
		##
		if( is.null(choose.pch) ) choose.pch <- c( 0:(swSummaryFcnCall$n.waves - 1) )
#*# need warning message if choose.pch-vector is NOT of length equal to total number of waves, then need to fix issue.....
		##
		for( plot.indx in 1:swSummaryFcnCall$n.waves ){
			plot(0,0, type="n", xlim=c(sw.xMin.Wave, sw.xMax.Wave), ylim=c(sw.yMin.Wave, sw.yMax.Wave), xaxt="n", xlab=choose.xlab, ylab="Mean Response", main=paste("Wave ", plot.indx, sep=""))
			axis(1, at=c(swSummaryFcnCall$time.at.each.wave), labels=c(swSummaryFcnCall$time.at.each.wave))
			##
			lines( swSummaryFcnCall$time.at.each.wave, swMeanResponse.Wave[plot.indx,], lty=1, col="grey" )
			##
			###**************************************************************************************************************************
			###**************************************************************************************************************************
			### ADDITION (wave): Allows different Tx types.......(e.g., 0, .8, .9, 1, 2)
			###
			for( j in 1:length(uniqueTx) ){
				swDsnTmp.j <- swDsnTmp
				swDsnTmp.j[ (swDsnTmp.j < uniqueTx[j] - 1e-8) | (swDsnTmp > uniqueTx[j] + 1e-8) ] <- NA
				swDsnTmp.j[ !is.na(swDsnTmp.j) ] <- 1
				unique(swDsnTmp)
				swPlotTx.j.Wave <- unique(cbind(swDsnTmp.j, swDsnTmp))[,1:dim(swDsnTmp.j)[2]] * swMeanResponse.Wave
				points( swSummaryFcnCall$time.at.each.wave, swPlotTx.j.Wave[plot.indx,], pch=choose.pch[plot.indx], cex=choose.cex, col=uniqueTxColors[j] )
			}
			###**************************************************************************************************************************
			###**************************************************************************************************************************
			##
			select.uniqueTx <- sort(unique( unique(swDsnTmp)[plot.indx,] ))
			select.uniqueTxColors <- uniqueTxColors[1:length(select.uniqueTx)]
			legend( choose.tx.pos, c("Tx:", select.uniqueTx), text.col=c("black", select.uniqueTxColors), ncol=(1 + length(select.uniqueTx)), x.intersp=0.01, xjust=1, y.intersp=0.5, cex=0.9)
		}
	}
}
###*********************************************
###
### (by.wave == FALSE)
###
###*********************************************
else if( by.wave == FALSE ){
	swMeanResponse.Wave <- swSummaryFcnCall$response.wave
	swMeanResponse.Cluster <- swSummaryFcnCall$response.cluster
	sw.xMin.Cluster <- min( swSummaryFcnCall$time.at.each.wave ) - 0.5
	sw.xMax.Cluster <- max( swSummaryFcnCall$time.at.each.wave ) + 0.5
	sw.yMin.Cluster <- min( swMeanResponse.Cluster )
	sw.yMax.Cluster <- max( swMeanResponse.Cluster )
#*# want to figure out how to specify xlim and ylim, maybe ylim=NULL as default
	###*********************************************
	###
	### (by.wave == FALSE & combined.plot == TRUE)
	###
	###*********************************************
	if( combined.plot==TRUE ) {
		par( mfrow=c(1,1) )
		plot(0,0, type="n", xlim=c(sw.xMin.Cluster, sw.xMax.Cluster), ylim=c(sw.yMin.Cluster, sw.yMax.Cluster), xaxt="n", xlab=choose.xlab, ylab="Mean Response", main=choose.main)
		axis(1, at=c(swSummaryFcnCall$time.at.each.wave), labels=c(swSummaryFcnCall$time.at.each.wave))
		##
		for( i in 1:dim(swMeanResponse.Cluster)[1] ){
			lines( swSummaryFcnCall$time.at.each.wave, swMeanResponse.Cluster[i,], lty=1, col="grey" )
		}
		##
#*# need to modify warning message(s) for choose.col????
		##
###**************************************************************************************************************************
###**************************************************************************************************************************
### ADDITION (cluster): Allows different Tx types.......(e.g., 0, .8, .9, 1, 2)
###
		for( j in 1:length(uniqueTx) ){
			swDsnTmp.j <- swDsnTmp
			swDsnTmp.j[ (swDsnTmp.j < uniqueTx[j] - 1e-8) | (swDsnTmp > uniqueTx[j] + 1e-8) ] <- NA
			swDsnTmp.j[ !is.na(swDsnTmp.j) ] <- 1
			swPlotTx.j.Cluster <- cbind(swDsnTmp.j, swDsnTmp)[,1:dim(swDsnTmp.j)[2]] * swMeanResponse.Cluster
			for( i in 1:dim( swDsnTmp )[1] ){
				text( swSummaryFcnCall$time.at.each.wave, swPlotTx.j.Cluster[i,], labels=as.character(i), cex=choose.cex, col=uniqueTxColors[j] )
			}
		}
###**************************************************************************************************************************
###**************************************************************************************************************************
		##
		custom.legend.vals <- cbind( cumsum(swSummaryFcnCall$clusters) - swSummaryFcnCall$clusters + 1, cumsum(swSummaryFcnCall$clusters), 1:swSummaryFcnCall$n.waves )
		custom.legend <- apply( custom.legend.vals, 1, function(z){ paste(z[1], ":", z[2], " (", z[3], ")", sep="") })
		##
		legend( choose.tx.pos, c("Tx:", uniqueTx), text.col=c("black", uniqueTxColors), ncol=(1 + length(uniqueTx)), x.intersp=0.01, xjust=1, y.intersp=0.5, cex=0.9)
		##
		if(!is.null(choose.legend.pos) && choose.legend.pos == "mouseclick"){
		  cat("Click on plot to add a legend\n")
		  #utils::flush.console()
		  legend( locator(1), legend=custom.legend, col="black", ncol=choose.ncol, title="Cluster (Wave)" )
		}
		else if(!is.null(choose.legend.pos)) legend( choose.legend.pos, legend=custom.legend, col="black", ncol=choose.ncol, title="Cluster (Wave)" )
	}
	###*********************************************
	###
	### (by.wave == FALSE & combined.plot == FALSE)
	###
	###*********************************************
	else if( combined.plot==FALSE ){
		##
		if( is.null(choose.mfrow) ) {
			if( swSummaryFcnCall$n.waves==1 ){ choose.mfrow <- c(1,1) }
			else if( swSummaryFcnCall$n.waves==2 ){ choose.mfrow <- c(1,2) }
			else if( swSummaryFcnCall$n.waves==3 | swSummaryFcnCall$n.waves==4 ){ choose.mfrow <- c(2,2) }
			else if( swSummaryFcnCall$n.waves==5 | swSummaryFcnCall$n.waves==6 ){ choose.mfrow <- c(2,3) }
			else if( swSummaryFcnCall$n.waves==7 | swSummaryFcnCall$n.waves==8 ){ choose.mfrow <- c(2,4) }
			else if( swSummaryFcnCall$n.waves==9 ){ choose.mfrow <- c(3,3) }
			else if( swSummaryFcnCall$n.waves>=10 | swSummaryFcnCall$n.waves<=12 ){ choose.mfrow <- c(3,4) }
		}
		par( mfrow=choose.mfrow )
		##
#*# need to modify warning message(s) for choose.col????
		##
		custom.main.vals <- cbind( cumsum(swSummaryFcnCall$clusters) - swSummaryFcnCall$clusters + 1, cumsum(swSummaryFcnCall$clusters), 1:swSummaryFcnCall$n.waves )
		##
		for( plot.indx in 1:swSummaryFcnCall$n.waves ){
			##
			if( custom.main.vals[plot.indx, 2] - custom.main.vals[plot.indx, 1] < 2 ) {
				choose.main <- paste( "Wave ", custom.main.vals[plot.indx, 3], ", Cluster ", custom.main.vals[plot.indx, 1], "-", custom.main.vals[plot.indx, 2], sep="" )
			}else { choose.main <- paste( "Wave ", custom.main.vals[plot.indx, 3], ", Clusters ", custom.main.vals[plot.indx, 1], "-", custom.main.vals[plot.indx, 2], sep="" ) }
			##
			plot(0,0, type="n", xlim=c(sw.xMin.Cluster, sw.xMax.Cluster), ylim=c(sw.yMin.Cluster - (sw.yMin.Cluster/10), sw.yMax.Cluster + (sw.yMin.Cluster/10)), xaxt="n", xlab=choose.xlab, ylab="Mean Response", main=choose.main)
			axis(1, at=c(swSummaryFcnCall$time.at.each.wave), labels=c(swSummaryFcnCall$time.at.each.wave))
			##
			for( i in custom.main.vals[plot.indx, 1]:custom.main.vals[plot.indx, 2] ) {
				lines( swSummaryFcnCall$time.at.each.wave, swMeanResponse.Cluster[i,], lty=1, col="grey" )
			}
			##
			###**************************************************************************************************************************
			###**************************************************************************************************************************
			### ADDITION (cluster): Allows different Tx types.......(e.g., 0, .8, .9, 1, 2)
			###
			for( j in 1:length(uniqueTx) ){
				swDsnTmp.j <- swDsnTmp
				swDsnTmp.j[ (swDsnTmp.j < uniqueTx[j] - 1e-8) | (swDsnTmp > uniqueTx[j] + 1e-8) ] <- NA
				swDsnTmp.j[ !is.na(swDsnTmp.j) ] <- 1
				swPlotTx.j.Cluster <- cbind(swDsnTmp.j, swDsnTmp)[,1:dim(swDsnTmp.j)[2]] * swMeanResponse.Cluster
				for( i in custom.main.vals[plot.indx, 1]:custom.main.vals[plot.indx, 2] ) {
					text( swSummaryFcnCall$time.at.each.wave, swPlotTx.j.Cluster[i,], labels=as.character(i), cex=choose.cex, col=uniqueTxColors[j] )
				}
				## should not need to return: swPlotTx.j.Cluster
			}
			###**************************************************************************************************************************
			###**************************************************************************************************************************
			##
			select.uniqueTx <- sort(unique( unique(swDsnTmp)[plot.indx,] ))
			select.uniqueTxColors <- uniqueTxColors[1:length(select.uniqueTx)]
			legend( choose.tx.pos, c("Tx:", select.uniqueTx), text.col=c("black", select.uniqueTxColors), ncol=(1 + length(select.uniqueTx)), x.intersp=0.01, xjust=1, y.intersp=0.5, cex=0.9)
		}
	}
}
## Function ending with next brace.
}

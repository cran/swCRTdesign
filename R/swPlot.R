swPlot <- function(response.var, tx.var, time.var, cluster.var, data, 
                   choose.mfrow=NULL, by.wave=TRUE, combined.plot=TRUE, choose.xlab="Time", 
                   choose.main=NULL, choose.pch=NULL, choose.cex=1, choose.tx.col=NULL, 
                   choose.tx.lty = c(2,1), choose.ncol=2, xlim=NULL, ylim=NULL,choose.tx.pos="topright", choose.legend.pos="right")
{
#Last update: Jim Hughes 5/1/2023, v. 4.0
##########
#Use swSummary to process data
##########
#Future extension: some warning messages making sure user inputs are consistent with expectations
swSummaryFcnCall <- swSummary(response.var=substitute(response.var), tx.var=substitute(tx.var), time.var=substitute(time.var), cluster.var=substitute(cluster.var), data=data, type="mean", digits=16, fcn.Call=TRUE)
swDsnTmp <- swSummaryFcnCall$swDsn
##
##########
#Set up colors for plots
##########
#For traditional data (one Tx level) only the colorTx0 and colorTx1 will be used.
	uniqueTx <- sort(unique( as.vector(swDsnTmp) ))
if( is.null(choose.tx.col) ) {
	uniqueTxColors <- c("blue","red","green3", "orange", "brown","purple", "lavender")
} else { 
  if (length(choose.tx.col)!=length(uniqueTx)) stop("Length of colors vector shorter than number of intervention states; colors will be recycled")
  uniqueTxColors <- choose.tx.col
}
###*********************************************
###
### (by.wave == TRUE)
###
###*********************************************
if( by.wave==TRUE ) {
	swMeanResponse.Wave <- swSummaryFcnCall$response.wave
	swMeanResponse.Cluster <- swSummaryFcnCall$response.cluster
	if (is.null(xlim)){
	  sw.xMin.Wave <- min( swSummaryFcnCall$time.at.each.wave )
	  sw.xMax.Wave <- max( swSummaryFcnCall$time.at.each.wave )
	} else {
	  if (length(xlim)!=2) stop("Length of xlim must be 2 if provided")
	  sw.xMin.Wave <- xlim[1]
	  sw.xMax.Wave <- xlim[2]
	}
	if (is.null(ylim)){
	  sw.yMin.Wave <- min( swMeanResponse.Wave, na.rm=TRUE )
  	sw.yMax.Wave <- max( swMeanResponse.Wave, na.rm=TRUE )
	} else {
	  if (length(ylim)!=2) stop("Length of ylim must be 2 if provided")
	  sw.yMin.Wave <- ylim[1]
	  sw.yMax.Wave <- ylim[2]
	}
	###*********************************************
	###
	### (by.wave == TRUE & combined.plot == TRUE)
	###
	###*********************************************
	if( combined.plot==TRUE ) {
		par( mfrow=c(1,1) )
	  #Set up empty plot
		plot(0,0, type="n", xlim=c(sw.xMin.Wave, sw.xMax.Wave), ylim=c(sw.yMin.Wave, sw.yMax.Wave), xaxt="n", xlab=choose.xlab, ylab="Mean Response", main=choose.main)
		axis(1, at=c(swSummaryFcnCall$time.at.each.wave), labels=c(swSummaryFcnCall$time.at.each.wave))
		#If only 2 unique treatments, plot control lines as dash (default lty=2) and treatment lines as solid (default lty=1). First plot all lines as dashed (lty=2) [here]. Then later draw solid lines for treatment times (lty=1) [after].
		#If more than 2 unique treatments, plot all lines as solid (default lty=1)
		if( length(uniqueTx) == 2 ){
		  for (i in 1:dim(swMeanResponse.Wave)[1]) {
		    lines(swSummaryFcnCall$time.at.each.wave, swMeanResponse.Wave[i,], lty = choose.tx.lty[1], col = "grey")
		  }
		}else if( length(uniqueTx) > 2 ){
		  for (i in 1:dim(swMeanResponse.Wave)[1]) {
		    lines(swSummaryFcnCall$time.at.each.wave, swMeanResponse.Wave[i,], lty = 1, col = "grey")
		  }
		}
		if( is.null(choose.pch) ) choose.pch <- c( 0:(swSummaryFcnCall$n.waves - 1) )
		#Potential future addition: warning message if choose.pch-vector is NOT of length equal to total number of waves, then need to fix issue
		#For each of the unique treatment conditions (usually just control and treatment), plot points for each cluster, colored by treatment and with shape determined by cluster
		for( j in 1:length(uniqueTx) ){
			swDsnTmp.j <- swDsnTmp
			swDsnTmp.j[ (swDsnTmp.j < uniqueTx[j] - 1e-8) | (swDsnTmp > uniqueTx[j] + 1e-8) ] <- NA
			swDsnTmp.j[ !is.na(swDsnTmp.j) ] <- 1
			unique(swDsnTmp)
			swPlotTx.j.Wave <- unique(cbind(swDsnTmp.j, swDsnTmp))[,1:dim(swDsnTmp.j)[2]] * swMeanResponse.Wave
			if( sum(swDsnTmp.j, na.rm=TRUE)==0 ) swPlotTx.j.Wave <- (1 - unique(cbind(swDsnTmp.j, swDsnTmp))[,1:dim(swDsnTmp.j)[2]]) * swMeanResponse.Wave
			for( i in 1:dim( unique(swDsnTmp) )[1] ){
				points( swSummaryFcnCall$time.at.each.wave, swPlotTx.j.Wave[i,], pch=choose.pch[i], cex=choose.cex, col=uniqueTxColors[j] )
			  #If only 2 unique treatments (j=2), plot control lines as dash (lty=2) and treatment lines as solid (lty=1). First plot all lines as dashed (lty=2) [already done above]. Then later draw solid lines for treatment times (lty=1) [here].
			  #Draw solid white line first, to cover lines drawn before.  (Not necessary if trt is lty 1)
			  if( j==2 ){
			    lines(swSummaryFcnCall$time.at.each.wave, swPlotTx.j.Wave[i, ], lty = 1, col = "white")
			    lines(swSummaryFcnCall$time.at.each.wave, swPlotTx.j.Wave[i, ], lty = choose.tx.lty[2], col = "grey")
			  }
			  }
		}
		#Make legend
		legend( choose.tx.pos, c("Tx:", uniqueTx), text.col=c("black", uniqueTxColors), ncol=(1 + length(uniqueTx)), x.intersp=0.01, xjust=1, y.intersp=0.5, cex=0.9)
		#Mouseclick option: used to be default, but that was causing issues, so changed it
		if(!is.null(choose.legend.pos) && choose.legend.pos == "mouseclick"){
		  cat("Click on plot to add a legend\n")#Print a message to the user
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
		#Pick plot grid layout (Future extension: make this a little more automated so it can accept any number of waves as input)
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
		if( is.null(choose.pch) ) choose.pch <- c( 0:(swSummaryFcnCall$n.waves - 1) )
    #Future extension: need warning message if choose.pch-vector is NOT of length equal to total number of waves, then need to fix issue
		for( plot.indx in 1:swSummaryFcnCall$n.waves ){
		  #Set up empty plot
		  plot(0,0, type="n", xlim=c(sw.xMin.Wave, sw.xMax.Wave), ylim=c(sw.yMin.Wave, sw.yMax.Wave), xaxt="n", xlab=choose.xlab, ylab="Mean Response", main=paste("Wave ", plot.indx, sep=""))
			axis(1, at=c(swSummaryFcnCall$time.at.each.wave), labels=c(swSummaryFcnCall$time.at.each.wave))
			#If only 2 unique treatments, plot control lines as dash (default lty=2) and treatment lines as solid (default lty=1). First plot all lines as dashed (lty=2) [here]. Then later draw solid lines for treatment times (lty=1) [after].
			#If more than 2 unique treatments, plot all lines as solid (default lty=1)
			#Don't need for i here because each plot only has one wave.
			if( length(uniqueTx) == 2 ){
			  lines(swSummaryFcnCall$time.at.each.wave, swMeanResponse.Wave[plot.indx, ], lty = choose.tx.lty[1], col = "grey")
			}else if( length(uniqueTx) > 2 ){
			  lines(swSummaryFcnCall$time.at.each.wave, swMeanResponse.Wave[plot.indx, ], lty = 1, col = "grey")
			}
			#For each of the unique treatment conditions (usually just control and treatment), plot points for each cluster, colored by treatment and with shape determined by cluster
			for( j in 1:length(uniqueTx) ){
				swDsnTmp.j <- swDsnTmp
				swDsnTmp.j[ (swDsnTmp.j < uniqueTx[j] - 1e-8) | (swDsnTmp > uniqueTx[j] + 1e-8) ] <- NA
				swDsnTmp.j[ !is.na(swDsnTmp.j) ] <- 1
				unique(swDsnTmp)
				swPlotTx.j.Wave <- unique(cbind(swDsnTmp.j, swDsnTmp))[,1:dim(swDsnTmp.j)[2]] * swMeanResponse.Wave
				points( swSummaryFcnCall$time.at.each.wave, swPlotTx.j.Wave[plot.indx,], pch=choose.pch[plot.indx], cex=choose.cex, col=uniqueTxColors[j] )
				#If only 2 unique treatments (j=2), plot control lines as dash (lty=2) and treatment lines as solid (lty=1). First plot all lines as dashed (lty=2) [already done above]. Then later draw solid lines for treatment times (lty=1) [here].
				#Don't need for i here because each plot only has one wave.
				if( j==2 ){
				  lines(swSummaryFcnCall$time.at.each.wave, swPlotTx.j.Wave[plot.indx, ], lty = 1, col = "white")
				  lines(swSummaryFcnCall$time.at.each.wave, swPlotTx.j.Wave[plot.indx, ], lty = choose.tx.lty[2], col = "grey")
				}
				}
			#Make legend
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
	if (is.null(xlim)){
	  sw.xMin.Cluster <- min( swSummaryFcnCall$time.at.each.wave ) - 0.5
	  sw.xMax.Cluster <- max( swSummaryFcnCall$time.at.each.wave ) + 0.5
	} else {
	  sw.xMin.Cluster <- xlim[1]
	  sw.xMax.Cluster <- xlim[2]
	}
	if (is.null(ylim)){
	  sw.yMin.Cluster <- min( swMeanResponse.Cluster, na.rm=TRUE )
	  sw.yMax.Cluster <- max( swMeanResponse.Cluster, na.rm=TRUE )
	} else {
	  sw.yMin.Cluster <- ylim[1]
	  sw.yMax.Cluster <- ylim[2]
	}
	###*********************************************
	###
	### (by.wave == FALSE & combined.plot == TRUE)
	###
	###*********************************************
	if( combined.plot==TRUE ) {
		par( mfrow=c(1,1) )
	  #Set up empty plot
		plot(0,0, type="n", xlim=c(sw.xMin.Cluster, sw.xMax.Cluster), ylim=c(sw.yMin.Cluster, sw.yMax.Cluster), xaxt="n", xlab=choose.xlab, ylab="Mean Response", main=choose.main)
		axis(1, at=c(swSummaryFcnCall$time.at.each.wave), labels=c(swSummaryFcnCall$time.at.each.wave))
		#If only 2 unique treatments, plot control lines as dash (lty=2) and treatment lines as solid (lty=1). First plot all lines as dashed (lty=2) [here]. Then later draw solid lines for treatment times (lty=1) [after].
		#If more than 2 unique treatments, plot all lines as solid (lty=1)
		if( length(uniqueTx) == 2 ){
		  for (i in 1:dim(swMeanResponse.Cluster)[1]) {
		    lines(swSummaryFcnCall$time.at.each.wave, swMeanResponse.Cluster[i, ], lty = choose.tx.lty[1], col = "grey")
		  }
		}else if( length(uniqueTx) > 2 ){
		  for (i in 1:dim(swMeanResponse.Wave)[1]) {
		    lines(swSummaryFcnCall$time.at.each.wave, swMeanResponse.Cluster[i, ], lty = 1, col = "grey")
		  }
		}
		#For each of the unique treatment conditions (usually just control and treatment), plot points for each cluster, colored by treatment and with shape determined by cluster
		for( j in 1:length(uniqueTx) ){
			swDsnTmp.j <- swDsnTmp
			swDsnTmp.j[ (swDsnTmp.j < uniqueTx[j] - 1e-8) | (swDsnTmp > uniqueTx[j] + 1e-8) ] <- NA
			swDsnTmp.j[ !is.na(swDsnTmp.j) ] <- 1
			swPlotTx.j.Cluster <- cbind(swDsnTmp.j, swDsnTmp)[,1:dim(swDsnTmp.j)[2]] * swMeanResponse.Cluster
			for( i in 1:dim( swDsnTmp )[1] ){
				text( swSummaryFcnCall$time.at.each.wave, swPlotTx.j.Cluster[i,], labels=as.character(i), cex=choose.cex, col=uniqueTxColors[j] )
			  #If only 2 unique treatments (j=2), plot control lines as dash (lty=2) and treatment lines as solid (lty=1). First plot all lines as dashed (lty=2) [already done above]. Then later draw solid lines for treatment times (lty=1) [here].
			  if( j==2 ){
			    lines(swSummaryFcnCall$time.at.each.wave, swPlotTx.j.Cluster[i, ], lty = 1, col = "white")
			    lines(swSummaryFcnCall$time.at.each.wave, swPlotTx.j.Cluster[i, ], lty = choose.tx.lty[2], col = "grey")
			  }
			  }
		}
		#Make legend
		custom.legend.vals <- cbind( cumsum(swSummaryFcnCall$clusters) - swSummaryFcnCall$clusters + 1, cumsum(swSummaryFcnCall$clusters), 1:swSummaryFcnCall$n.waves )
		custom.legend <- apply( custom.legend.vals, 1, function(z){ paste(z[1], ":", z[2], " (", z[3], ")", sep="") })
		legend( choose.tx.pos, c("Tx:", uniqueTx), text.col=c("black", uniqueTxColors), ncol=(1 + length(uniqueTx)), x.intersp=0.01, xjust=1, y.intersp=0.5, cex=0.9)
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
		#Future extension: extend this better to larger numbers of waves
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
		custom.main.vals <- cbind( cumsum(swSummaryFcnCall$clusters) - swSummaryFcnCall$clusters + 1, cumsum(swSummaryFcnCall$clusters), 1:swSummaryFcnCall$n.waves )
		for( plot.indx in 1:swSummaryFcnCall$n.waves ){
			if( custom.main.vals[plot.indx, 2] - custom.main.vals[plot.indx, 1] < 2 ) {
				choose.main <- paste( "Wave ", custom.main.vals[plot.indx, 3], ", Cluster ", custom.main.vals[plot.indx, 1], "-", custom.main.vals[plot.indx, 2], sep="" )
			}else { choose.main <- paste( "Wave ", custom.main.vals[plot.indx, 3], ", Clusters ", custom.main.vals[plot.indx, 1], "-", custom.main.vals[plot.indx, 2], sep="" ) }
		  #Set up empty plot
			plot(0,0, type="n", xlim=c(sw.xMin.Cluster, sw.xMax.Cluster), ylim=c(sw.yMin.Cluster - (sw.yMin.Cluster/10), sw.yMax.Cluster + (sw.yMin.Cluster/10)), xaxt="n", xlab=choose.xlab, ylab="Mean Response", main=choose.main)
			axis(1, at=c(swSummaryFcnCall$time.at.each.wave), labels=c(swSummaryFcnCall$time.at.each.wave))
			#If only 2 unique treatments, plot control lines as dash (lty=2) and treatment lines as solid (lty=1). First plot all lines as dashed (lty=2) [here]. Then later draw solid lines for treatment times (lty=1) [after].
			#If more than 2 unique treatments, plot all lines as solid (lty=1)
			if( length(uniqueTx) == 2 ){
			  for (i in custom.main.vals[plot.indx, 1]:custom.main.vals[plot.indx, 2]) {
			    lines(swSummaryFcnCall$time.at.each.wave, swMeanResponse.Cluster[i, ], lty = choose.tx.lty[1], col = "grey")
			  }
			}else if( length(uniqueTx) > 2 ){
			  for (i in custom.main.vals[plot.indx, 1]:custom.main.vals[plot.indx, 2]) {
			    lines(swSummaryFcnCall$time.at.each.wave, swMeanResponse.Cluster[i, ], lty = 1, col = "grey")
			  }
			}
			#For each of the unique treatment conditions (usually just control and treatment), plot points for each cluster, colored by treatment and with shape determined by cluster
			for( j in 1:length(uniqueTx) ){
				swDsnTmp.j <- swDsnTmp
				swDsnTmp.j[ (swDsnTmp.j < uniqueTx[j] - 1e-8) | (swDsnTmp > uniqueTx[j] + 1e-8) ] <- NA
				swDsnTmp.j[ !is.na(swDsnTmp.j) ] <- 1
				swPlotTx.j.Cluster <- cbind(swDsnTmp.j, swDsnTmp)[,1:dim(swDsnTmp.j)[2]] * swMeanResponse.Cluster
				for( i in custom.main.vals[plot.indx, 1]:custom.main.vals[plot.indx, 2] ) {
					text( swSummaryFcnCall$time.at.each.wave, swPlotTx.j.Cluster[i,], labels=as.character(i), cex=choose.cex, col=uniqueTxColors[j] )
				  #If only 2 unique treatments (j=2), plot control lines as dash (lty=2) and treatment lines as solid (lty=1). First plot all lines as dashed (lty=2) [already done above]. Then later draw solid lines for treatment times (lty=1) [here].
				  if( j==2 ){
				    lines(swSummaryFcnCall$time.at.each.wave, swPlotTx.j.Cluster[i, ], lty = 1, col = "white")
				    lines(swSummaryFcnCall$time.at.each.wave, swPlotTx.j.Cluster[i, ], lty = choose.tx.lty[2], col = "grey")
				  }
				  }
				## should not need to return: swPlotTx.j.Cluster
			}
			#Make legend
			select.uniqueTx <- sort(unique( unique(swDsnTmp)[plot.indx,] ))
			select.uniqueTxColors <- uniqueTxColors[1:length(select.uniqueTx)]
			legend( choose.tx.pos, c("Tx:", select.uniqueTx), text.col=c("black", select.uniqueTxColors), ncol=(1 + length(select.uniqueTx)), x.intersp=0.01, xjust=1, y.intersp=0.5, cex=0.9)
		}
	}
}
## Function ending with next brace.
}

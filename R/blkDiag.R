blkDiag <-
function(z) {
##
## (v.2014.03.04)
##
    ##
    twoMatListBlkDiagFCN <- function(zTwoMatList) {
        i <- 1
        zList <- zTwoMatList
        zTop <- cbind(  zList[[i]], matrix( 0, dim(zList[[i]])[1], dim(zList[[i+1]])[1] )  )
        zBot <- cbind(  matrix( 0, dim(zList[[i+1]])[1], dim(zList[[i]])[1] ), zList[[i+1]]  )
        zTwoBlkDiag <- rbind( zTop, zBot )
        zTwoBlkDiag
    }
    ##
	if( is.array(z) ) {
    	zArray <- z
		for( i in 1:(length(zArray[1,1,])-1) ) {
        	if(i==1) {
    	        zMat1 <- zArray[,,1]
	            zMat2 <- zArray[,,2]
          	}
            else {
            	zMat1 <- zArrayBlkDiag
                zMat2 <- zArray[,,i+1]
			}
            zArrayBlkDiag <- twoMatListBlkDiagFCN( zTwoMatList=list(zMat1, zMat2) )
		}
        zBlkDiag <- zArrayBlkDiag
	}
    else if( is.list(z) ) {
		zList <- z
        for( i in 1:(length(zList)-1) ) {
        	if(i==1) {
            	zMat1 <- zList[[1]]
                zMat2 <- zList[[2]]
           	}
            else {
            	zMat1 <- zListBlkDiag
                zMat2 <- zList[[i+1]]
            }
            zListBlkDiag <- twoMatListBlkDiagFCN( zTwoMatList=list(zMat1, zMat2) )
		}
        zBlkDiag <- zListBlkDiag
	}
    zBlkDiag
}

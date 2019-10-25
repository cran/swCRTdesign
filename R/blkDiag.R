blkDiag <-
function(z) {
##
## (v.2014.03.04)
##
##########
#twoMatListBlkDiagFCN function: take two blocks and return a block diagonal matrix
##########
    twoMatListBlkDiagFCN <- function(zTwoMatList) {
        i <- 1
        zList <- zTwoMatList
        #Take the first block, and add a block of 0's on its right-hand side
        zTop <- cbind(  zList[[i]], matrix( 0, dim(zList[[i]])[1], dim(zList[[i+1]])[1] )  )
        #Take the second block, and add a block of 0's on its left-hand side
        zBot <- cbind(  matrix( 0, dim(zList[[i+1]])[1], dim(zList[[i]])[1] ), zList[[i+1]]  )
        zTwoBlkDiag <- rbind( zTop, zBot )
        zTwoBlkDiag#Block diagonal matrix
    }
##########
#Blocks are an array:
##########
	if( is.array(z) ) {
    	zArray <- z
    #Taking blocks one at a time, assemble into a block diagonal matrix
		for( i in 1:(length(zArray[1,1,])-1) ) {
        	if(i==1) {
    	        zMat1 <- zArray[,,1]#First block
	            zMat2 <- zArray[,,2]#Second block
          	}
            else {
            	zMat1 <- zArrayBlkDiag
                zMat2 <- zArray[,,i+1]#Next block to be added
			}
            zArrayBlkDiag <- twoMatListBlkDiagFCN( zTwoMatList=list(zMat1, zMat2) )#Adds the next block to the matrix
		}
        zBlkDiag <- zArrayBlkDiag
	}
##########
#Blocks are a list:
##########
    else if( is.list(z) ) {
		zList <- z
		    #Taking blocks one at a time, assemble into a block diagonal matrix
        for( i in 1:(length(zList)-1) ) {
        	if(i==1) {
            	zMat1 <- zList[[1]]#First block
                zMat2 <- zList[[2]]#Second block
           	}
            else {
            	zMat1 <- zListBlkDiag
                zMat2 <- zList[[i+1]]#Next block to be added
            }
            zListBlkDiag <- twoMatListBlkDiagFCN( zTwoMatList=list(zMat1, zMat2) )#Adds the next block to the matrix
		}
        zBlkDiag <- zListBlkDiag
	}
##########
#Return block diagonal matrix:
##########
  zBlkDiag
}

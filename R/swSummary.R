swSummary <-
function (response.var, tx.var, time.var, cluster.var, data=NULL,
    type = "mean", digits = 16, fcn.Call = FALSE)
{
  #Last update: unknown; added comments 10/23/2019, v. 3.1
  ##########
  #Unpack data
  ##########
    if (!is.null(data) & !fcn.Call) {
        Response <- eval(parse(text = paste("data$", substitute(response.var),
            sep = "")))
        Tx <- eval(parse(text = paste("data$", substitute(tx.var),
            sep = "")))
        Time <- eval(parse(text = paste("data$", substitute(time.var),
            sep = "")))
        Cluster <- eval(parse(text = paste("data$", substitute(cluster.var),
            sep = "")))
    }
    else if (!is.null(data) & fcn.Call) {
        Response <- eval(parse(text = paste("data$", response.var,
            sep = "")))
        Tx <- eval(parse(text = paste("data$", tx.var, sep = "")))
        Time <- eval(parse(text = paste("data$", time.var, sep = "")))
        Cluster <- eval(parse(text = paste("data$", cluster.var,
            sep = "")))
    }
    else if (is.null(data)) {
        Response <- response.var
        Tx <- tx.var
        Time <- time.var
        Cluster <- cluster.var
    }
    df <- data.frame(cbind(Response, Tx, Cluster, Time))
    names(df) <- c("response", "tx", "cluster", "time")
    ##########
    ##########
    tmpCluster <- sort(unique(df$cluster))#List of unique cluster IDs
    tmpTime <- sort(unique(df$time))#List of unique time IDs
    tmpWave <- NULL#Initialize list of wave IDs (for each cluster)
    for (clusterID in tmpCluster) {
        tmpPreWave <- NULL
        for (timeID in tmpTime) {
            tmpPreWave <- c(tmpPreWave, 1 - mean(df[df$cluster ==
                clusterID & df$time == timeID, "tx"]))
        }
        tmpWave <- c(tmpWave, sum(tmpPreWave))
    }
tmpWaveMat <- cbind(1:length(unique(tmpWave)), as.vector(table(tmpWave)) )
tmpWave <- as.vector( unlist( apply(tmpWaveMat, 1, function(z) rep(z[1], each=z[2])) ) )#For each cluster, list of wave IDs
    tmpClusterWave <- cbind(tmpCluster, tmpWave)
    df$wave <- df$cluster#Initialize wave column, and fill in below with a wave ID for each row
    for (clusterID in tmpCluster) {
        df[df$cluster == clusterID & df$wave == clusterID, ]$wave <- tmpClusterWave[clusterID,
            "tmpWave"]
    }
    ##########
    #Collect info for each wave
    ##########
    tmpTx.Wave.mean <- NULL
    tmpResponse.Wave.mean <- NULL
    tmpResponse.Wave.sum <- NULL
    tmpResponse.Wave.n <- NULL
    for (waveID in unique(tmpWave)) {
        tmpTx.Wave.Time.mean <- NULL
        tmpResponse.Wave.Time.mean <- NULL
        tmpResponse.Wave.Time.sum <- NULL
        tmpResponse.Wave.Time.n <- NULL
        for (timeID in tmpTime) {
            tmpTx.Wave.Time.Vec <- df[df$wave == waveID & df$time ==
                timeID, "tx"]#Treatment indicators for people from this wave at this time
            tmpTx.Wave.Time.mean <- c(tmpTx.Wave.Time.mean, mean(tmpTx.Wave.Time.Vec,
                na.rm = TRUE))#Append treatment status at this wave and time to tmpTx.Wave.Time.mean
            tmpResponse.Wave.Time.Vec <- df[df$wave == waveID &
                df$time == timeID, "response"]#Responses for people from this wave at this time
            tmpResponse.Wave.Time.mean <- c(tmpResponse.Wave.Time.mean,
                as.numeric(sprintf(fmt = paste("%8.", digits,
                  "f", sep = ""), mean(tmpResponse.Wave.Time.Vec,
                  na.rm = TRUE))))#Append average response at this wave and time
            tmpResponse.Wave.Time.sum <- c(tmpResponse.Wave.Time.sum,
                sum(tmpResponse.Wave.Time.Vec, na.rm = TRUE))#Append sum of responses at this wave and time; not sure why this doesn't also have the rounding option (future extension?)
            tmpResponse.Wave.Time.n <- c(tmpResponse.Wave.Time.n,
                sum(!is.na(tmpResponse.Wave.Time.Vec)))#Append number of recorded responses at this wave and time
        }
        tmpTx.Wave.mean <- rbind(tmpTx.Wave.mean, tmpTx.Wave.Time.mean)
        tmpResponse.Wave.mean <- rbind(tmpResponse.Wave.mean,
            tmpResponse.Wave.Time.mean)
        tmpResponse.Wave.sum <- rbind(tmpResponse.Wave.sum, tmpResponse.Wave.Time.sum)
        tmpResponse.Wave.n <- rbind(tmpResponse.Wave.n, tmpResponse.Wave.Time.n)
    }
    #Add labels for waves and times:
    rownames(tmpTx.Wave.mean) <- paste("(wave ", unique(tmpWave),
        ")", sep = "")
    colnames(tmpTx.Wave.mean) <- paste("(time ", tmpTime, ")",
        sep = "")
    rownames(tmpResponse.Wave.mean) <- paste("(wave ", unique(tmpWave),
        ")", sep = "")
    colnames(tmpResponse.Wave.mean) <- paste("(time ", tmpTime,
        ")", sep = "")
    rownames(tmpResponse.Wave.sum) <- paste("(wave ", unique(tmpWave),
        ")", sep = "")
    colnames(tmpResponse.Wave.sum) <- paste("(time ", tmpTime,
        ")", sep = "")
    rownames(tmpResponse.Wave.n) <- paste("(wave ", unique(tmpWave),
        ")", sep = "")
    colnames(tmpResponse.Wave.n) <- paste("(time ", tmpTime,
        ")", sep = "")
    ##########
    #Collect info for each cluster
    ##########
    tmpTx.Cluster.mean <- NULL
    tmpResponse.Cluster.mean <- NULL
    tmpResponse.Cluster.sum <- NULL
    tmpResponse.Cluster.n <- NULL
    for (clusterID in tmpCluster) {
        tmpTx.Cluster.Time.mean <- NULL
        tmpResponse.Cluster.Time.mean <- NULL
        tmpResponse.Cluster.Time.sum <- NULL
        tmpResponse.Cluster.Time.n <- NULL
        for (timeID in tmpTime) {
            tmpTx.Cluster.Time.Vec <- df[df$cluster == clusterID &
                df$time == timeID, "tx"]#Treatment indicators for people from this cluster at this time
            tmpTx.Cluster.Time.mean <- c(tmpTx.Cluster.Time.mean,
                mean(tmpTx.Cluster.Time.Vec, na.rm = TRUE))#Append treatment status at this cluster and time
            tmpResponse.Cluster.Time.Vec <- df[df$cluster ==
                clusterID & df$time == timeID, "response"]#Responses for people from this cluster at this time
            tmpResponse.Cluster.Time.mean <- c(tmpResponse.Cluster.Time.mean,
                as.numeric(sprintf(fmt = paste("%8.", digits,
                  "f", sep = ""), mean(tmpResponse.Cluster.Time.Vec,
                  na.rm = TRUE))))#Append average response at this cluster and time
            tmpResponse.Cluster.Time.sum <- c(tmpResponse.Cluster.Time.sum,
                as.numeric(sprintf(fmt = paste("%8.", digits,
                  "f", sep = ""), sum(tmpResponse.Cluster.Time.Vec,
                  na.rm = TRUE))))#Append sum of responses at this cluster and time
            tmpResponse.Cluster.Time.n <- c(tmpResponse.Cluster.Time.n,
                sum(!is.na(tmpResponse.Cluster.Time.Vec)))#Append number of recorded responses at this cluster and time
        }
        tmpTx.Cluster.mean <- rbind(tmpTx.Cluster.mean, tmpTx.Cluster.Time.mean)
        tmpResponse.Cluster.mean <- rbind(tmpResponse.Cluster.mean,
            tmpResponse.Cluster.Time.mean)
        tmpResponse.Cluster.sum <- rbind(tmpResponse.Cluster.sum,
            tmpResponse.Cluster.Time.sum)
        tmpResponse.Cluster.n <- rbind(tmpResponse.Cluster.n,
            tmpResponse.Cluster.Time.n)
    }
    #Add labels for clusters and times:
    rownames(tmpTx.Cluster.mean) <- paste("(cluster ", tmpCluster,
        ")", sep = "")
    colnames(tmpTx.Cluster.mean) <- paste("(time ", tmpTime,
        ")", sep = "")
    rownames(tmpResponse.Cluster.mean) <- paste("(cluster ",
        tmpCluster, ")", sep = "")
    colnames(tmpResponse.Cluster.mean) <- paste("(time ", tmpTime,
        ")", sep = "")
    rownames(tmpResponse.Cluster.sum) <- paste("(cluster ", tmpCluster,
        ")", sep = "")
    colnames(tmpResponse.Cluster.sum) <- paste("(time ", tmpTime,
        ")", sep = "")
    rownames(tmpResponse.Cluster.n) <- paste("(cluster ", tmpCluster,
        ")", sep = "")
    colnames(tmpResponse.Cluster.n) <- paste("(time ", tmpTime,
        ")", sep = "")
    ##########
    #Return correct type of result
    ##########
    if (type == "mean") {
        rslt <- list(swDsn = tmpTx.Cluster.mean, swDsn.unique.clusters = tmpTx.Wave.mean,
            n.waves = length(unique(tmpWave)), clusters = as.vector(table(tmpWave)),
            n.clusters = length(unique(tmpCluster)), time.at.each.wave = unique(tmpTime),
            response.cluster = tmpResponse.Cluster.mean, response.wave = tmpResponse.Wave.mean)
    }
    else if (type == "sum") {
        rslt <- list(swDsn = tmpTx.Cluster.mean, swDsn.unique.clusters = tmpTx.Wave.mean,
            n.waves = length(unique(tmpWave)), clusters = as.vector(table(tmpWave)),
            n.clusters = length(unique(tmpCluster)), time.at.each.wave = unique(tmpTime),
            response.cluster = tmpResponse.Cluster.sum, response.wave = tmpResponse.Wave.sum)
    }
    else if (type == "n") {
        rslt <- list(type = type, swDsn = tmpTx.Cluster.mean,
            swDsn.unique.clusters = tmpTx.Wave.mean, n.waves = length(unique(tmpWave)),
            clusters = as.vector(table(tmpWave)), n.clusters = length(unique(tmpCluster)),
            time.at.each.wave = unique(tmpTime), response.cluster = tmpResponse.Cluster.n,
            response.wave = tmpResponse.Wave.n)
    }
    rslt
}

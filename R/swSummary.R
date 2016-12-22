swSummary <-
function (response.var, tx.var, time.var, cluster.var, data, 
    type = "mean", digits = 16, fcn.Call = FALSE) 
{
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
    tmpCluster <- sort(unique(df$cluster))
    tmpTime <- sort(unique(df$time))
    tmpWave <- NULL
    for (clusterID in tmpCluster) {
        tmpPreWave <- NULL
        for (timeID in tmpTime) {
            tmpPreWave <- c(tmpPreWave, 1 - mean(df[df$cluster == 
                clusterID & df$time == timeID, "tx"]))
        }
        tmpWave <- c(tmpWave, sum(tmpPreWave))
    }
tmpWaveMat <- cbind(1:length(unique(tmpWave)), as.vector(table(tmpWave)) )
tmpWave <- unlist( apply(tmpWaveMat, 1, function(z) rep(z[1], each=z[2])) )
    tmpClusterWave <- cbind(tmpCluster, tmpWave)
    df$wave <- df$cluster
    for (clusterID in tmpCluster) {
        df[df$cluster == clusterID & df$wave == clusterID, ]$wave <- tmpClusterWave[clusterID, 
            "tmpWave"]
    }
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
                timeID, "tx"]
            tmpTx.Wave.Time.mean <- c(tmpTx.Wave.Time.mean, mean(tmpTx.Wave.Time.Vec, 
                na.rm = TRUE))
            tmpResponse.Wave.Time.Vec <- df[df$wave == waveID & 
                df$time == timeID, "response"]
            tmpResponse.Wave.Time.mean <- c(tmpResponse.Wave.Time.mean, 
                as.numeric(sprintf(fmt = paste("%8.", digits, 
                  "f", sep = ""), mean(tmpResponse.Wave.Time.Vec, 
                  na.rm = TRUE))))
            tmpResponse.Wave.Time.sum <- c(tmpResponse.Wave.Time.sum, 
                sum(tmpResponse.Wave.Time.Vec, na.rm = TRUE))
            tmpResponse.Wave.Time.n <- c(tmpResponse.Wave.Time.n, 
                sum(!is.na(tmpResponse.Wave.Time.Vec)))
        }
        tmpTx.Wave.mean <- rbind(tmpTx.Wave.mean, tmpTx.Wave.Time.mean)
        tmpResponse.Wave.mean <- rbind(tmpResponse.Wave.mean, 
            tmpResponse.Wave.Time.mean)
        tmpResponse.Wave.sum <- rbind(tmpResponse.Wave.sum, tmpResponse.Wave.Time.sum)
        tmpResponse.Wave.n <- rbind(tmpResponse.Wave.n, tmpResponse.Wave.Time.n)
    }
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
                df$time == timeID, "tx"]
            tmpTx.Cluster.Time.mean <- c(tmpTx.Cluster.Time.mean, 
                mean(tmpTx.Cluster.Time.Vec, na.rm = TRUE))
            tmpResponse.Cluster.Time.Vec <- df[df$cluster == 
                clusterID & df$time == timeID, "response"]
            tmpResponse.Cluster.Time.mean <- c(tmpResponse.Cluster.Time.mean, 
                as.numeric(sprintf(fmt = paste("%8.", digits, 
                  "f", sep = ""), mean(tmpResponse.Cluster.Time.Vec, 
                  na.rm = TRUE))))
            tmpResponse.Cluster.Time.sum <- c(tmpResponse.Cluster.Time.sum, 
                as.numeric(sprintf(fmt = paste("%8.", digits, 
                  "f", sep = ""), sum(tmpResponse.Cluster.Time.Vec, 
                  na.rm = TRUE))))
            tmpResponse.Cluster.Time.n <- c(tmpResponse.Cluster.Time.n, 
                sum(!is.na(tmpResponse.Cluster.Time.Vec)))
        }
        tmpTx.Cluster.mean <- rbind(tmpTx.Cluster.mean, tmpTx.Cluster.Time.mean)
        tmpResponse.Cluster.mean <- rbind(tmpResponse.Cluster.mean, 
            tmpResponse.Cluster.Time.mean)
        tmpResponse.Cluster.sum <- rbind(tmpResponse.Cluster.sum, 
            tmpResponse.Cluster.Time.sum)
        tmpResponse.Cluster.n <- rbind(tmpResponse.Cluster.n, 
            tmpResponse.Cluster.Time.n)
    }
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

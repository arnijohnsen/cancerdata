# Load data files (if data files have already been loaded they are not
#  loaded again)
cat("Loading data files\n")
if(!exists("linkedProbesGenes")){
  load("../Rdata/BRCA/info/linkedProbesGenes.Rdata")
}
if(!exists("normalMethyl")){
  load("../Rdata/BRCA/data/BRCA-NMP.Rdata")
}
if(!exists("cancerMethyl")){
  load("../Rdata/BRCA/data/BRCA-CMP.Rdata")
}
if(!exists("normalRnaseq")){
  load("../Rdata/BRCA/data/BRCA-NEA.Rdata")
}
if(!exists("cancerRnaseq")){
  load("../Rdata/BRCA/data/BRCA-CEA.Rdata")
}

# Select samples which are in both expression and methylation data sets
normalSamples <- intersect(rownames(normalMethyl), rownames(normalRnaseq))
cancerSamples <- intersect(rownames(cancerMethyl), rownames(cancerRnaseq))

# Ask user for gene or probe names to plot
search <- readline("Enter gene or probe name: ")
search <- paste("^", search, "$", sep="")
# Create index list with all matching data
n <- union(grep(search, linkedProbesGenes$probes), grep(search, linkedProbesGenes$genes))
eps <- 1e-10
if (length(n) == 0){
  cat("No hits\n")
} else {
  cat(length(n), "hits\n")
  for (i in 1:length(n)){
    m <- n[i];
    # Select probe and sample
    probe <- as.character(linkedProbesGenes$probes[m])
    gene  <- as.character(linkedProbesGenes$genes[m])
    
    # Split data in normal/cancer x(methylation)/y(expression)
    cat("n =", m, "\n")
    xNdata <- normalMethyl[normalSamples, probe] + eps
    yNdata <- normalRnaseq[normalSamples, gene] + eps
    xCdata <- cancerMethyl[cancerSamples, probe] + eps
    yCdata <- cancerRnaseq[cancerSamples, gene] + eps

    # Get plot limits for regular plot
    ymax <- max(c(yNdata, yCdata))

    # Perform correlation tests
    test <- cor.test(c(xNdata,xCdata), c(yNdata, yCdata))
    logtest <- cor.test(log(c(xNdata,xCdata)), log(c(yNdata, yCdata)))
    par(mfrow=c(1,2))
    # Plot regular data
    plot(xCdata, yCdata, col="red", pch=20, xlim=c(0,1), ylim=c(0,ymax),
         xlab="Promoter methylation", ylab="RNA expression", 
	 sub=paste("r=", format(test$estimate, digits=3), 
	       ", R^2=", format(test$estimate^2, digits=3), 
	         ", p=", format(test$p.value, digits=3), sep=""), 
	 main=paste(probe, gene, sep="-"))
    points(xNdata, yNdata, col="blue", pch=20)

    # Get plot limits for log plot
    xlogmin <- min(log(c(xNdata, xCdata)), na.rm=T)
    xlogmax <- max(log(c(xNdata, xCdata)), na.rm=T)
    ylogmin <- min(log(c(yNdata, yCdata)), na.rm=T)
    ylogmax <- max(log(c(yNdata, yCdata)), na.rm=T)
    
    # Plot log transformed data
    plot(log(xCdata), log(yCdata), col="red", pch=20, 
         xlim=c(xlogmin,xlogmax), ylim=c(ylogmin,ylogmax),
         xlab="log Promoter methylation", ylab="log RNA expression", 
	 sub=paste("r=", format(logtest$estimate, digits=3), 
	       ", R^2=", format(logtest$estimate^2, digits=3), 
	         ", p=", format(logtest$p.value, digits=3), sep=""), 
	 main=paste(probe, gene, sep="-"))
    points(log(xNdata), log(yNdata), col="blue", pch=20)

    # Ask user for input to plot next plot
    cat ("Press [enter] for next plot")
    line <- readline()
  }
}

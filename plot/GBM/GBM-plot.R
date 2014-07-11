# Load data files (if data files have already been loaded they are not
#  loaded again)
cat("Loading data files\n")
if(!exists("linkedProbesGenes")){
  load("../Rdata/GBM/info/linkedProbesGenes.Rdata")
}
if(!exists("cancerMethyl")){
  load("../Rdata/GBM/data/GBM-CMP.Rdata")
}
if(!exists("cancerRnaseq")){
  load("../Rdata/GBM/data/GBM-CEA.Rdata")
}

# Select samples which are in both expression and methylation data sets
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
    cat("n =", m, "\n")
    
    # Split data in x(methylation)/y(expression)
    xCdata <- cancerMethyl[cancerSamples, probe] + eps
    yCdata <- cancerRnaseq[cancerSamples, gene] + eps

    # Get plot limits for regular plot
    ymax <- max(yCdata)

    # Perform correlation tests
    test <- cor.test(xCdata, yCdata)
    logtest <- cor.test(log(xCdata), log(yCdata))
    par(mfrow=c(1,2))
    # Plot regular data
    plot(xCdata, yCdata, col="red", pch=20, xlim=c(0,1), ylim=c(0,ymax),
         xlab="Promoter methylation", ylab="RNA expression", 
	 sub=paste("r=", format(test$estimate, digits=3), 
	       ", R^2=", format(test$estimate^2, digits=3), 
	         ", p=", format(test$p.value, digits=3), sep=""), 
	 main=paste(probe, gene, sep="-"))

    # Get plot limits for log plot
    xlogmin <- min(log(xCdata), na.rm=T)
    xlogmax <- max(log(xCdata), na.rm=T)
    ylogmin <- min(log(yCdata), na.rm=T)
    ylogmax <- max(log(yCdata), na.rm=T)
    
    # Plot log transformed data
    plot(log(xCdata), log(yCdata), col="red", pch=20, 
         xlim=c(xlogmin,xlogmax), ylim=c(ylogmin,ylogmax),
         xlab="log Promoter methylation", ylab="log RNA expression", 
	 sub=paste("r=", format(logtest$estimate, digits=3), 
	       ", R^2=", format(logtest$estimate^2, digits=3), 
	         ", p=", format(logtest$p.value, digits=3), sep=""), 
	 main=paste(probe, gene, sep="-"))

    # Ask user for input to plot next plot
    cat ("Press [enter] for next plot")
    line <- readline()
  }
}

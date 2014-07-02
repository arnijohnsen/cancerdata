# Set wdir
# setwd("~/cancerdata/")

# Load data files
cat("Loading data files\n")
if(!exists("linkedProbesGenes")){
  load("../Rdata/linkedProbesGenes.Rdata")
}
if(!exists("normalMethyl")){
  load("../Rdata/normalMethylPromoters.Rdata")
}
if(!exists("cancerMethyl")){
  load("../Rdata/cancerMethylPromoters.Rdata")
}
if(!exists("normalRnaseq")){
  load("../Rdata/normalRnaseqAllGenes.Rdata")
}
if(!exists("cancerRnaseq")){
  load("../Rdata/cancerRnaseqAllGenes.Rdata")
}

normalSamples <- intersect(rownames(normalMethyl), rownames(normalRnaseq))
cancerSamples <- intersect(rownames(cancerMethyl), rownames(cancerRnaseq))

search <- readline("Enter gene or probe name: ")
search <- paste("^", search, "$", sep="")
n <- union(grep(search, linkedProbesGenes$probes), grep(search, linkedProbesGenes$genes))

if (length(n) == 0){
  cat("No hits\n")
} else {
  cat(length(n), "hits\n")
  for (i in 1:length(n)){
    m <- n[i];
    probe <- as.character(linkedProbesGenes$probes[m])
    gene  <- as.character(linkedProbesGenes$genes[m])
    cat("n =", m, "\n")
    xNdata <- normalMethyl[normalSamples, probe]
    yNdata <- normalRnaseq[normalSamples, gene]
    xCdata <- cancerMethyl[cancerSamples, probe]
    yCdata <- cancerRnaseq[cancerSamples, gene]
    #ymax <- max(quantile(yNdata, probs=prob, na.rm=TRUE), quantile(yCdata, probs=prob, na.rm=TRUE))
    ymax <- max(c(yNdata, yCdata))
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
    cat ("Press [enter] for next plot")
    line <- readline()
  }
}

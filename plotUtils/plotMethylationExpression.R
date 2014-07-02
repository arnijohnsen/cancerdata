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

prob <- 0.99
search <- readline("Enter gene or probe name: ")

n <- union(grep(search, linkedProbesGenes$probes), pmatch(search, linkedProbesGenes$genes))

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
    plot(xCdata, yCdata, col="red", pch=20, xlim=c(0,1), ylim=c(0,ymax),
         xlab="Promoter methylation", ylab="RNA expression", main=paste(probe, gene, sep="-"))
    points(xNdata, yNdata, col="blue", pch=20)
    cat ("Press [enter] for next plot")
    line <- readline()
  }
}

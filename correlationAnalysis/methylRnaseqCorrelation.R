# Load data files
cat("Loading data files\n")
if(!exists("linkedProbesGenes")){
  load("../Rdata/BRCA/linkedProbesGenes.Rdata")
}
if(!exists("normalMethyl")){
  load("../Rdata/BRCA/normalMethylPromoters.Rdata")
}
if(!exists("cancerMethyl")){
  load("../Rdata/BRCA/cancerMethylPromoters.Rdata")
}
if(!exists("normalRnaseq")){
  load("../Rdata/BRCA/normalRnaseqAllGenes.Rdata")
}
if(!exists("cancerRnaseq")){
  load("../Rdata/BRCA/cancerRnaseqAllGenes.Rdata")
}

cat("Resizing data frames\n")
if(!exists("allMethyl") || !exists("allRnaseq")){
normalSamples <- intersect(rownames(normalMethyl), rownames(normalRnaseq))
cancerSamples <- intersect(rownames(cancerMethyl), rownames(cancerRnaseq))

allMethyl <- rbind(normalMethyl[normalSamples,], 
                   cancerMethyl[cancerSamples,])
allRnaseq <- rbind(normalRnaseq[normalSamples,],
                   cancerRnaseq[cancerSamples,])
}

# Make smaller for testing
n <- dim(allMethyl)[2]
eps <- 1e-10
result <- data.frame(r2 = numeric(0), 
                     intEst  = numeric(0), sloEst  = numeric(0), 
		     intStd  = numeric(0), sloStd  = numeric(0), 
		     intTval = numeric(0), sloTval = numeric(0), 
		     intPval = numeric(0), sloPval = numeric(0), 
		     can99q  = numeric(0), can95q  = numeric(0), 
		     canMean = numeric(0), canSd   = numeric(0), 
		     norMean = numeric(0), norSd   = numeric(0), 
		     norExpr = numeric(0), canExpr = numeric(0), 
		     nonMethQuant = numeric(0))

pb <- txtProgressBar(min=1, max=n, style=3)
cat("Running loop for lm\n")
for(i in 1:n){
  setTxtProgressBar(pb, i)
  # lm model
  probe <- as.character(linkedProbesGenes$probes[i])
  gene  <- as.character(linkedProbesGenes$genes[i])
  x <- allMethyl[,probe]
  y <- allRnaseq[,gene] + eps
  fit <- lm(log(y) ~ log(x))
  sum <- summary(fit)
  # 99th quantile in cancer, mean/sd in cancer, mean/sd in normal
  xc <- x[73:599]
  yc <- y[73:599]
  result[i,] <- c(sum$r.squared, 
                  sum$coefficients, 
		  quantile(x[73:559], 0.99, na.rm=T), 
		  quantile(x[73:559], 0.95, na.rm=T), 
		  mean(x[73:559], na.rm=T), 
		  sd(x[73:599], na.rm=T), 
		  mean(x[1:72], na.rm=T), 
		  sd(x[1:72], na.rm=T), 
		  mean(y[1:72], na.rm=T),
		  mean(y[73:599], na.rm=T), 
		  quantile(yc[xc<0.2], 0.10, na.rm=T))

}
cat("\n")
rownames(result) <- paste(linkedProbesGenes$probes[1:n], linkedProbesGenes$genes[1:n], sep="-")
result$sloPadj <- p.adjust(result$sloPval, method="BH")
result$intPadj <- p.adjust(result$intPval, method="BH")
result$r <- sqrt(result$r2)*sign(result$sloEst)


# Filter data
save(result, "../Rdata/BRCA/lmMethylRnaseq.Rdata")
quit(save="no")

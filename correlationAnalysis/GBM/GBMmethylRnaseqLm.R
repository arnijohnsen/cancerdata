# Load data files
library(doParallel)
registerDoParallel(2)
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

cat("Resizing data frames\n")
if(!exists("allMethyl") || !exists("allRnaseq")){
cancerSamples <- intersect(rownames(cancerMethyl), rownames(cancerRnaseq))

allMethyl <- cancerMethyl[cancerSamples,]
allRnaseq <- cancerRnaseq[cancerSamples,]
}
rm(cancerMethyl, cancerRnaseq)
gc()

# Make smaller for testing
n <- dim(allMethyl)[2]
#n <- 256
eps <- 1e-10
pb <- txtProgressBar(min=1, max=n, style=3)
prt <- proc.time()
cat("Running loop for lm\n")

result <- data.frame(r2 = numeric(0), 
                     sloEst  = numeric(0), 
		     sloPval = numeric(0), 
		     can99q  = numeric(0),
		     nonMethQuant = numeric(0))
for(i in 1:n){
  setTxtProgressBar(pb, i)
  # lm model
  probe <- as.character(linkedProbesGenes$probes[i])
  gene  <- as.character(linkedProbesGenes$genes[i])
  x <- allMethyl[,probe]
  y <- allRnaseq[,gene] + eps
  fit <- lm(log(y) ~ log(x), weights=x^3)
  sum <- summary(fit)
  # 99th quantile in cancer, mean/sd in cancer, mean/sd in normal
  result[i,] <- c(sum$r.squared,                      # R^2
                  sum$coefficients[c(2,8)],           # Slope and P-value
		  quantile(x, 0.99, na.rm=T),        # 99th methyl quant. of cancer
		  quantile(y[x<0.2], 0.10, na.rm=T))# 10th expr. quant of non-methylated cancer

}
print(proc.time() - prt)
cat("\n")
colnames(result) <- c("r2", "sloEst", "sloPval", "can99q", "nonMethQuant")
rownames(result) <- paste(linkedProbesGenes$probes[1:n], linkedProbesGenes$genes[1:n], sep="-")
result$sloPadj <- p.adjust(result$sloPval, method="BH")
result$r <- sqrt(result$r2)*sign(result$sloEst)

# Filter data
save(result, file="../Rdata/GBM/info/lmMethylRnaseq.Rdata")
quit(save="no")

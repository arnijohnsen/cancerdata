# Load data files
cat("Loading data files\n")
load("../Rdata/BRCA/info/linkedProbesGenes.Rdata")
load("../Rdata/BRCA/data/BRCA-NMP.Rdata")
load("../Rdata/BRCA/data/BRCA-CMP.Rdata")
load("../Rdata/BRCA/data/BRCA-NEA.Rdata")
load("../Rdata/BRCA/data/BRCA-CEA.Rdata")

# Resize data frames
cat("Resizing data frames\n")
# Use only samples iwhich have both methylation and expression data
normalSamples <- intersect(rownames(normalMethyl), rownames(normalRnaseq))
cancerSamples <- intersect(rownames(cancerMethyl), rownames(cancerRnaseq))

# Combine cancer and normal data in one frame, using selected samples
allMethyl <- rbind(normalMethyl[normalSamples,], cancerMethyl[cancerSamples,])
allRnaseq <- rbind(normalRnaseq[normalSamples,], cancerRnaseq[cancerSamples,])

# Remove other not-used data.frames
rm(normalMethyl, cancerMethyl, normalRnaseq, cancerRnaseq)
gc()

n <- dim(allMethyl)[2]
eps <- 1e-10
pb <- txtProgressBar(min=1, max=n, style=3)
prt <- proc.time()

cat("Running loop for lm\n")
# Create data frame with result
result <- data.frame(r2 = numeric(0), 
                     sloEst  = numeric(0), 
		     sloPval = numeric(0), 
		     can99q  = numeric(0),
		     nonMethQuant = numeric(0))
# Run loop, in each iteration a log~log linear model is computed for
# a pair of gene-probe
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
  xc <- x[73:599]
  yc <- y[73:599]
  result[i,] <- c(sum$r.squared,                      # R^2
                  sum$coefficients[c(2,8)],           # Slope and P-value
		  quantile(xc, 0.99, na.rm=T),        # 99th methyl quant. of cancer
		  quantile(yc[xc<0.2], 0.10, na.rm=T))# 10th expr. quant of non-methylated cancer

}
print(proc.time() - prt)
cat("\n")

# Add row names (ex. cg12345678-ABCD1)
rownames(result) <- paste(linkedProbesGenes$probes[1:n], linkedProbesGenes$genes[1:n], sep="-")
# Create column with adjusted p-values and r
result$sloPadj <- p.adjust(result$sloPval, method="BH")
result$r <- sqrt(result$r2)*sign(result$sloEst)

# Save result and exit
save(result, file="../Rdata/BRCA/calc/BRCA-lmMethylRnaseq.Rdata")
quit(save="no")

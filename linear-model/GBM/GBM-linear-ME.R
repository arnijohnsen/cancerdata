# Load data files
cat("Loading data files\n")
load("../Rdata/GBM/info/GBM-linked-probes-genes.Rdata")
load("../Rdata/GBM/data/GBM-CMP.Rdata")
load("../Rdata/GBM/data/GBM-CEA.Rdata")

#size <- dim(GBM.linked.probes.genes)[1]
# Resize data frames
cat("Resizing data frames ..\n")
# Use only samples which have both methylation and expression data
cancer.samples <- intersect(rownames(GBM.CMP), rownames(GBM.CEA))

chunk.size <- 1000
total.size <- dim(GBM.linked.probes.genes)[1]
chunks <- ceiling(total.size/chunk.size)
res.all <- c()
pb <- txtProgressBar(min=1, max=chunks, style=3)
for(i in 1:chunks){
  setTxtProgressBar(pb, i)
  index <- (1+chunk.size*(i-1)):(min(chunk.size*i, total.size))
  # Combine cancer data in one frame, using selected samples
  GBM.AMP <- GBM.CMP[cancer.samples, as.character(GBM.linked.probes.genes$probes)[index]]
  GBM.AEA <- GBM.CEA[cancer.samples, as.character(GBM.linked.probes.genes$genes)[index]]
  colnames(GBM.AMP) <- paste(GBM.linked.probes.genes$probe[index],
                              GBM.linked.probes.genes$genes[index], sep=".")
  colnames(GBM.AEA) <- paste(GBM.linked.probes.genes$probe[index],
                              GBM.linked.probes.genes$genes[index], sep=".")
  GBM.all.data <- rbind(GBM.AMP, GBM.AEA)

  n <- length(cancer.samples)
  res <- apply(GBM.all.data, 2, function(x) {
    tmp <- lm(x[(n+1):(2*n)] ~ x[1:n], weights=x[1:n]^3)
    sum <- summary(tmp)
    c(sum$r.squared, sum$coefficients[c(2,8)])
  })
  res.all <- cbind(res.all, res)
}

cat("Creating data frame..\n")
rownames(res.all) <- c("r.squared", "slope", "p.raw")
GBM.linear.ME <- as.data.frame(t(res.all))
cat("Adjusting p-values and computing r\n")
GBM.linear.ME$p.adj <- p.adjust(GBM.linear.ME$p.raw, method="BH")
GBM.linear.ME$r <- sqrt(GBM.linear.ME$r.squared)*sign(GBM.linear.ME$slope)

save(GBM.linear.ME, file="../Rdata/GBM/calc/GBM-linear-ME.Rdata")
quit(save="no")

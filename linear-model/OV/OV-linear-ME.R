# Load data files
cat("Loading data files\n")
load("../Rdata/COAD/info/COAD-linked-probes-genes.Rdata")
load("../Rdata/COAD/data/COAD-CMP.Rdata")
load("../Rdata/COAD/data/COAD-CEA.Rdata")

#size <- dim(COAD.linked.probes.genes)[1]
# Resize data frames
cat("Resizing data frames ..\n")
# Use only samples iwhich have both methylation and expression data
cancer.samples <- intersect(rownames(COAD.CMP), rownames(COAD.CEA))

chunk.size <- 1000
total.size <- dim(COAD.linked.probes.genes)[1]
chunks <- ceiling(total.size/chunk.size)
res.all <- c()
pb <- txtProgressBar(min=1, max=chunks, style=3)
for(i in 1:chunks){
  setTxtProgressBar(pb, i)
  index <- (1+chunk.size*(i-1)):(min(chunk.size*i, total.size))
  # Combine cancer and normal data in one frame, using selected samples
  COAD.AMP <- COAD.CMP[cancer.samples, as.character(COAD.linked.probes.genes$probes)[index]]
  COAD.AEA <- COAD.CEA[cancer.samples, as.character(COAD.linked.probes.genes$genes)[index]]
  colnames(COAD.AMP) <- paste(COAD.linked.probes.genes$probe[index], 
                              COAD.linked.probes.genes$genes[index], sep=".")
  colnames(COAD.AEA) <- paste(COAD.linked.probes.genes$probe[index], 
                              COAD.linked.probes.genes$genes[index], sep=".")
  COAD.all.data <- rbind(COAD.AMP, COAD.AEA)

  n <- length(cancer.samples)
  res <- apply(COAD.all.data, 2, function(x) { 
    a <- x[1:n]
    b <- x[(n+1):(2*n)] + 1e-10
    tmp <- lm(log(b) ~ log(a), weights=a^3)
    sum <- summary(tmp)
    c(sum$r.squared, sum$coefficients[c(2,8)])
  })
  res.all <- cbind(res.all, res)
}

cat("Creating data frame..\n")
rownames(res.all) <- c("r.squared", "slope", "p.raw")
COAD.linear.ME <- as.data.frame(t(res.all))
cat("Adjusting p-values and computing r\n")
COAD.linear.ME$p.adj <- p.adjust(COAD.linear.ME$p.raw, method="BH")
COAD.linear.ME$r <- sqrt(COAD.linear.ME$r.squared)*sign(COAD.linear.ME$slope)

save(COAD.linear.ME, file="../Rdata/COAD/calc/COAD-linear-ME.Rdata")
quit(save="no")

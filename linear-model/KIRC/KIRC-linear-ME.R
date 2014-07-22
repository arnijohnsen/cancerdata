# Load data files
cat("Loading data files\n")
load("../Rdata/KIRC/info/KIRC-linked-probes-genes.Rdata")
load("../Rdata/KIRC/data/KIRC-CMP.Rdata")
load("../Rdata/KIRC/data/KIRC-CEA.Rdata")

#size <- dim(KIRC.linked.probes.genes)[1]
# Resize data frames
cat("Resizing data frames ..\n")
# Use only samples which have both methylation and expression data
cancer.samples <- intersect(rownames(KIRC.CMP), rownames(KIRC.CEA))

chunk.size <- 1000
total.size <- dim(KIRC.linked.probes.genes)[1]
chunks <- ceiling(total.size/chunk.size)
res.all <- c()
pb <- txtProgressBar(min=1, max=chunks, style=3)
for(i in 1:chunks){
  setTxtProgressBar(pb, i)
  index <- (1+chunk.size*(i-1)):(min(chunk.size*i, total.size))
  # Combine cancer data in one frame, using selected samples
  KIRC.AMP <- KIRC.CMP[cancer.samples, as.character(KIRC.linked.probes.genes$probes)[index]]
  KIRC.AEA <- KIRC.CEA[cancer.samples, as.character(KIRC.linked.probes.genes$genes)[index]]
  colnames(KIRC.AMP) <- paste(KIRC.linked.probes.genes$probe[index],
                              KIRC.linked.probes.genes$genes[index], sep=".")
  colnames(KIRC.AEA) <- paste(KIRC.linked.probes.genes$probe[index],
                              KIRC.linked.probes.genes$genes[index], sep=".")
  KIRC.all.data <- rbind(KIRC.AMP, KIRC.AEA)

  n <- length(cancer.samples)
  res <- apply(KIRC.all.data, 2, function(x) {
    tmp <- lm(x[(n+1):(2*n)] ~ x[1:n], weights=x[1:n]^3)
    sum <- summary(tmp)
    c(sum$r.squared, sum$coefficients[c(2,8)])
  })
  res.all <- cbind(res.all, res)
}

cat("Creating data frame..\n")
rownames(res.all) <- c("r.squared", "slope", "p.raw")
KIRC.linear.ME <- as.data.frame(t(res.all))
cat("Adjusting p-values and computing r\n")
KIRC.linear.ME$p.adj <- p.adjust(KIRC.linear.ME$p.raw, method="BH")
KIRC.linear.ME$r <- sqrt(KIRC.linear.ME$r.squared)*sign(KIRC.linear.ME$slope)

save(KIRC.linear.ME, file="../Rdata/KIRC/calc/KIRC-linear-ME.Rdata")
quit(save="no")

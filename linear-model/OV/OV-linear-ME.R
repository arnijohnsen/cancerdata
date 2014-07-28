# Load data files
cat("Loading data files\n")
load("../Rdata/OV/info/OV-linked-probes-genes.Rdata")
load("../Rdata/OV/data/OV-CMP.Rdata")
load("../Rdata/OV/data/OV-CEA.Rdata")

#size <- dim(OV.linked.probes.genes)[1]
# Resize data frames
cat("Resizing data frames ..\n")
# Use only samples which have both methylation and expression data
cancer.samples <- intersect(rownames(OV.CMP), rownames(OV.CEA))

chunk.size <- 1000
total.size <- dim(OV.linked.probes.genes)[1]
chunks <- ceiling(total.size/chunk.size)
res.all <- c()
pb <- txtProgressBar(min=1, max=chunks, style=3)
for(i in 1:chunks){
  setTxtProgressBar(pb, i)
  index <- (1+chunk.size*(i-1)):(min(chunk.size*i, total.size))
  # Combine cancer data in one frame, using selected samples
  OV.AMP <- OV.CMP[cancer.samples, as.character(OV.linked.probes.genes$probes)[index]]
  OV.AEA <- OV.CEA[cancer.samples, as.character(OV.linked.probes.genes$genes)[index]]
  colnames(OV.AMP) <- paste(OV.linked.probes.genes$probe[index],
                              OV.linked.probes.genes$genes[index], sep=".")
  colnames(OV.AEA) <- paste(OV.linked.probes.genes$probe[index],
                              OV.linked.probes.genes$genes[index], sep=".")
  OV.all.data <- rbind(OV.AMP, OV.AEA)

  n <- length(cancer.samples)
  res <- apply(OV.all.data, 2, function(x) {
    tmp <- lm(log(x[(n+1):(2*n)]) ~ log(x[1:n]), weights=x[1:n]^3)
    sum <- summary(tmp)
    c(sum$r.squared, sum$coefficients[c(2,8)])
  })
  res.all <- cbind(res.all, res)
}

cat("Creating data frame..\n")
rownames(res.all) <- c("r.squared", "slope", "p.raw")
OV.linear.ME <- as.data.frame(t(res.all))
cat("Adjusting p-values and computing r\n")
OV.linear.ME$p.adj <- p.adjust(OV.linear.ME$p.raw, method="BH")
OV.linear.ME$r <- sqrt(OV.linear.ME$r.squared)*sign(OV.linear.ME$slope)

save(OV.linear.ME, file="../Rdata/OV/calc/OV-linear-ME.Rdata")
quit(save="no")

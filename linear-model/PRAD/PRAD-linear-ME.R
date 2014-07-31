# Load data files
cat("Loading data files\n")
load("../Rdata/PRAD/info/PRAD-linked-probes-genes.Rdata")
load("../Rdata/PRAD/data/PRAD-NMP.Rdata")
load("../Rdata/PRAD/data/PRAD-CMP.Rdata")
load("../Rdata/PRAD/data/PRAD-NEA.Rdata")
load("../Rdata/PRAD/data/PRAD-CEA.Rdata")

#size <- dim(PRAD.linked.probes.genes)[1]
# Resize data frames
cat("Resizing data frames ..\n")
# Use only samples iwhich have both methylation and expression data
normal.samples <- intersect(rownames(PRAD.NMP), rownames(PRAD.NEA))
cancer.samples <- intersect(rownames(PRAD.CMP), rownames(PRAD.CEA))

chunk.size <- 1000
total.size <- dim(PRAD.linked.probes.genes)[1]
chunks <- ceiling(total.size/chunk.size)
res.all <- c()
pb <- txtProgressBar(min=1, max=chunks, style=3)
for(i in 1:chunks){
  setTxtProgressBar(pb, i)
  index <- (1+chunk.size*(i-1)):(min(chunk.size*i, total.size))
  # Combine cancer and normal data in one frame, using selected samples
  PRAD.AMP <- rbind(PRAD.NMP[normal.samples, as.character(PRAD.linked.probes.genes$probes)[index]], 
                    PRAD.CMP[cancer.samples, as.character(PRAD.linked.probes.genes$probes)[index]])
  PRAD.AEA <- rbind(PRAD.NEA[normal.samples, as.character(PRAD.linked.probes.genes$genes)[index]], 
                    PRAD.CEA[cancer.samples, as.character(PRAD.linked.probes.genes$genes)[index]]) 
  colnames(PRAD.AMP) <- paste(PRAD.linked.probes.genes$probe[index], 
                              PRAD.linked.probes.genes$genes[index], sep=".")
  colnames(PRAD.AEA) <- paste(PRAD.linked.probes.genes$probe[index], 
                              PRAD.linked.probes.genes$genes[index], sep=".")
  PRAD.all.data <- rbind(PRAD.AMP, PRAD.AEA)

  n <- length(normal.samples) + length(cancer.samples)
  res <- apply(PRAD.all.data, 2, function(x) { 
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
PRAD.linear.ME <- as.data.frame(t(res.all))
cat("Adjusting p-values and computing r\n")
PRAD.linear.ME$p.adj <- p.adjust(PRAD.linear.ME$p.raw, method="BH")
PRAD.linear.ME$r <- sqrt(PRAD.linear.ME$r.squared)*sign(PRAD.linear.ME$slope)

save(PRAD.linear.ME, file="../Rdata/PRAD/calc/PRAD-linear-ME.Rdata")
quit(save="no")

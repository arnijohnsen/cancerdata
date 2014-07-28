# Load data files
cat("Loading data files\n")
load("../Rdata/BRCA/info/BRCA-linked-probes-genes.Rdata")
load("../Rdata/BRCA/data/BRCA-NMP.Rdata")
load("../Rdata/BRCA/data/BRCA-CMP.Rdata")
load("../Rdata/BRCA/data/BRCA-NEA.Rdata")
load("../Rdata/BRCA/data/BRCA-CEA.Rdata")

#size <- dim(BRCA.linked.probes.genes)[1]
# Resize data frames
cat("Resizing data frames ..\n")
# Use only samples iwhich have both methylation and expression data
normal.samples <- intersect(rownames(BRCA.NMP), rownames(BRCA.NEA))
cancer.samples <- intersect(rownames(BRCA.CMP), rownames(BRCA.CEA))

chunk.size <- 1000
total.size <- dim(BRCA.linked.probes.genes)[1]
chunks <- ceiling(total.size/chunk.size)
res.all <- c()
pb <- txtProgressBar(min=1, max=chunks, style=3)
for(i in 1:chunks){
  setTxtProgressBar(pb, i)
  index <- (1+chunk.size*(i-1)):(min(chunk.size*i, total.size))
  # Combine cancer and normal data in one frame, using selected samples
  BRCA.AMP <- rbind(BRCA.NMP[normal.samples, as.character(BRCA.linked.probes.genes$probes)[index]], 
                    BRCA.CMP[cancer.samples, as.character(BRCA.linked.probes.genes$probes)[index]])
  BRCA.AEA <- rbind(BRCA.NEA[normal.samples, as.character(BRCA.linked.probes.genes$genes)[index]], 
                    BRCA.CEA[cancer.samples, as.character(BRCA.linked.probes.genes$genes)[index]]) 
  colnames(BRCA.AMP) <- paste(BRCA.linked.probes.genes$probe[index], 
                              BRCA.linked.probes.genes$genes[index], sep=".")
  colnames(BRCA.AEA) <- paste(BRCA.linked.probes.genes$probe[index], 
                              BRCA.linked.probes.genes$genes[index], sep=".")
  BRCA.all.data <- rbind(BRCA.AMP, BRCA.AEA)

  n <- length(normal.samples) + length(cancer.samples)
  res <- apply(BRCA.all.data, 2, function(x) { 
    tmp <- lm(log(x[(n+1):(2*n)]) ~ log(x[1:n]), weights=x[1:n]^3)
    sum <- summary(tmp)
    c(sum$r.squared, sum$coefficients[c(2,8)])
  })
  res.all <- cbind(res.all, res)
}

cat("Creating data frame..\n")
rownames(res.all) <- c("r.squared", "slope", "p.raw")
BRCA.linear.ME <- as.data.frame(t(res.all))
cat("Adjusting p-values and computing r\n")
BRCA.linear.ME$p.adj <- p.adjust(BRCA.linear.ME$p.raw, method="BH")
BRCA.linear.ME$r <- sqrt(BRCA.linear.ME$r.squared)*sign(BRCA.linear.ME$slope)

save(BRCA.linear.ME, file="../Rdata/BRCA/calc/BRCA-linear-ME.Rdata")
quit(save="no")

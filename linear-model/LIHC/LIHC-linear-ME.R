# Load data files
cat("Loading data files\n")
load("../Rdata/LIHC/info/LIHC-linked-probes-genes.Rdata")
load("../Rdata/LIHC/data/LIHC-NMP.Rdata")
load("../Rdata/LIHC/data/LIHC-CMP.Rdata")
load("../Rdata/LIHC/data/LIHC-NEA.Rdata")
load("../Rdata/LIHC/data/LIHC-CEA.Rdata")

# Fix one troublemaker probe
LIHC.linked.probes.genes <- LIHC.linked.probes.genes[-53182,]

#size <- dim(LIHC.linked.probes.genes)[1]
# Resize data frames
cat("Resizing data frames ..\n")
# Use only samples iwhich have both methylation and expression data
normal.samples <- intersect(rownames(LIHC.NMP), rownames(LIHC.NEA))
cancer.samples <- intersect(rownames(LIHC.CMP), rownames(LIHC.CEA))

chunk.size <- 1000
total.size <- dim(LIHC.linked.probes.genes)[1]
chunks <- ceiling(total.size/chunk.size)
res.all <- c()
pb <- txtProgressBar(min=1, max=chunks, style=3)
for(i in 1:chunks){
  setTxtProgressBar(pb, i)
  index <- (1+chunk.size*(i-1)):(min(chunk.size*i, total.size))
  # Combine cancer and normal data in one frame, using selected samples
  LIHC.AMP <- rbind(LIHC.NMP[normal.samples, as.character(LIHC.linked.probes.genes$probes)[index]], 
                    LIHC.CMP[cancer.samples, as.character(LIHC.linked.probes.genes$probes)[index]])
  LIHC.AEA <- rbind(LIHC.NEA[normal.samples, as.character(LIHC.linked.probes.genes$genes)[index]], 
                    LIHC.CEA[cancer.samples, as.character(LIHC.linked.probes.genes$genes)[index]]) 
  colnames(LIHC.AMP) <- paste(LIHC.linked.probes.genes$probe[index], 
                              LIHC.linked.probes.genes$genes[index], sep=".")
  colnames(LIHC.AEA) <- paste(LIHC.linked.probes.genes$probe[index], 
                              LIHC.linked.probes.genes$genes[index], sep=".")
  LIHC.all.data <- rbind(LIHC.AMP, LIHC.AEA)

  n <- length(normal.samples) + length(cancer.samples)
  res <- apply(LIHC.all.data, 2, function(x) { 
    tmp <- lm(log(x[(n+1):(2*n)]) ~ log(x[1:n]), weights=x[1:n]^3)
    sum <- summary(tmp)
    c(sum$r.squared, sum$coefficients[c(2,8)])
  })
  res.all <- cbind(res.all, res)
}

cat("Creating data frame..\n")
rownames(res.all) <- c("r.squared", "slope", "p.raw")
LIHC.linear.ME <- as.data.frame(t(res.all))
cat("Adjusting p-values and computing r\n")
LIHC.linear.ME$p.adj <- p.adjust(LIHC.linear.ME$p.raw, method="BH")
LIHC.linear.ME$r <- sqrt(LIHC.linear.ME$r.squared)*sign(LIHC.linear.ME$slope)

save(LIHC.linear.ME, file="../Rdata/LIHC/calc/LIHC-linear-ME.Rdata")
quit(save="no")

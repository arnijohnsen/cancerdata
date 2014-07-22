# Load data files
cat("Loading data files\n")
load("../Rdata/LUSC/info/LUSC-linked-probes-genes.Rdata")
load("../Rdata/LUSC/data/LUSC-NMP.Rdata")
load("../Rdata/LUSC/data/LUSC-CMP.Rdata")
load("../Rdata/LUSC/data/LUSC-NEA.Rdata")
load("../Rdata/LUSC/data/LUSC-CEA.Rdata")

#size <- dim(LUSC.linked.probes.genes)[1]
# Resize data frames
cat("Resizing data frames ..\n")
# Use only samples iwhich have both methylation and expression data
normal.samples <- intersect(rownames(LUSC.NMP), rownames(LUSC.NEA))
cancer.samples <- intersect(rownames(LUSC.CMP), rownames(LUSC.CEA))

chunk.size <- 1000
total.size <- dim(LUSC.linked.probes.genes)[1]
chunks <- ceiling(total.size/chunk.size)
res.all <- c()
pb <- txtProgressBar(min=1, max=chunks, style=3)
for(i in 1:chunks){
  setTxtProgressBar(pb, i)
  index <- (1+chunk.size*(i-1)):(min(chunk.size*i, total.size))
  # Combine cancer and normal data in one frame, using selected samples
  LUSC.AMP <- rbind(LUSC.NMP[normal.samples, as.character(LUSC.linked.probes.genes$probes)[index]], 
                    LUSC.CMP[cancer.samples, as.character(LUSC.linked.probes.genes$probes)[index]])
  LUSC.AEA <- rbind(LUSC.NEA[normal.samples, as.character(LUSC.linked.probes.genes$genes)[index]], 
                    LUSC.CEA[cancer.samples, as.character(LUSC.linked.probes.genes$genes)[index]]) 
  colnames(LUSC.AMP) <- paste(LUSC.linked.probes.genes$probe[index], 
                              LUSC.linked.probes.genes$genes[index], sep=".")
  colnames(LUSC.AEA) <- paste(LUSC.linked.probes.genes$probe[index], 
                              LUSC.linked.probes.genes$genes[index], sep=".")
  LUSC.all.data <- rbind(LUSC.AMP, LUSC.AEA)
  
  n <- length(normal.samples) + length(cancer.samples)
  res <- apply(LUSC.all.data, 2, function(x) { 
    tmp <- lm(x[(n+1):(2*n)] ~ x[1:n], weights=x[1:n]^3)
    sum <- summary(tmp)
    c(sum$r.squared, sum$coefficients[c(2,8)])
  })
  res.all <- cbind(res.all, res)
}

cat("Creating data frame..\n")
rownames(res.all) <- c("r.squared", "slope", "p.raw")
LUSC.linear.ME <- as.data.frame(t(res.all))
cat("Adjusting p-values and computing r\n")
LUSC.linear.ME$p.adj <- p.adjust(LUSC.linear.ME$p.raw, method="BH")
LUSC.linear.ME$r <- sqrt(LUSC.linear.ME$r.squared)*sign(LUSC.linear.ME$slope)

save(LUSC.linear.ME, file="../Rdata/LUSC/calc/LUSC-linear-ME.Rdata")
quit(save="no")

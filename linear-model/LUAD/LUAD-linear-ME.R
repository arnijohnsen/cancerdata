# Load data files
cat("Loading data files\n")
load("../Rdata/LUAD/info/LUAD-linked-probes-genes.Rdata")
load("../Rdata/LUAD/data/LUAD-NMP.Rdata")
load("../Rdata/LUAD/data/LUAD-CMP.Rdata")
load("../Rdata/LUAD/data/LUAD-NEA.Rdata")
load("../Rdata/LUAD/data/LUAD-CEA.Rdata")

#size <- dim(LUAD.linked.probes.genes)[1]
# Resize data frames
cat("Resizing data frames ..\n")
# Use only samples iwhich have both methylation and expression data
normal.samples <- intersect(rownames(LUAD.NMP), rownames(LUAD.NEA))
cancer.samples <- intersect(rownames(LUAD.CMP), rownames(LUAD.CEA))

chunk.size <- 1000
total.size <- dim(LUAD.linked.probes.genes)[1]
chunks <- ceiling(total.size/chunk.size)
res.all <- c()
pb <- txtProgressBar(min=1, max=chunks, style=3)
for(i in 1:chunks){
  setTxtProgressBar(pb, i)
  index <- (1+chunk.size*(i-1)):(min(chunk.size*i, total.size))
  # Combine cancer and normal data in one frame, using selected samples
  LUAD.AMP <- rbind(LUAD.NMP[normal.samples, as.character(LUAD.linked.probes.genes$probes)[index]], 
                    LUAD.CMP[cancer.samples, as.character(LUAD.linked.probes.genes$probes)[index]])
  LUAD.AEA <- rbind(LUAD.NEA[normal.samples, as.character(LUAD.linked.probes.genes$genes)[index]], 
                    LUAD.CEA[cancer.samples, as.character(LUAD.linked.probes.genes$genes)[index]]) 
  colnames(LUAD.AMP) <- paste(LUAD.linked.probes.genes$probe[index], 
                              LUAD.linked.probes.genes$genes[index], sep=".")
  colnames(LUAD.AEA) <- paste(LUAD.linked.probes.genes$probe[index], 
                              LUAD.linked.probes.genes$genes[index], sep=".")
  LUAD.all.data <- rbind(LUAD.AMP, LUAD.AEA)
  
  n <- length(normal.samples) + length(cancer.samples)
  res <- apply(LUAD.all.data, 2, function(x) { 
    tmp <- lm(x[(n+1):(2*n)] ~ x[1:n], weights=x[1:n]^3)
    sum <- summary(tmp)
    c(sum$r.squared, sum$coefficients[c(2,8)])
  })
  res.all <- cbind(res.all, res)
}

cat("Creating data frame..\n")
rownames(res.all) <- c("r.squared", "slope", "p.raw")
LUAD.linear.ME <- as.data.frame(t(res.all))
cat("Adjusting p-values and computing r\n")
LUAD.linear.ME$p.adj <- p.adjust(LUAD.linear.ME$p.raw, method="BH")
LUAD.linear.ME$r <- sqrt(LUAD.linear.ME$r.squared)*sign(LUAD.linear.ME$slope)

save(LUAD.linear.ME, file="../Rdata/LUAD/calc/LUAD-linear-ME.Rdata")
quit(save="no")

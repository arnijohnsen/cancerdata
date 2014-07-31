# Load data files
cat("Loading data files\n")
load("../Rdata/READ/info/READ-linked-probes-genes.Rdata")
load("../Rdata/READ/data/READ-NMP.Rdata")
load("../Rdata/READ/data/READ-CMP.Rdata")
load("../Rdata/READ/data/READ-NEA.Rdata")
load("../Rdata/READ/data/READ-CEA.Rdata")

#size <- dim(READ.linked.probes.genes)[1]
# Resize data frames
cat("Resizing data frames ..\n")
# Use only samples iwhich have both methylation and expression data
normal.samples <- intersect(rownames(READ.NMP), rownames(READ.NEA))
cancer.samples <- intersect(rownames(READ.CMP), rownames(READ.CEA))

chunk.size <- 1000
total.size <- dim(READ.linked.probes.genes)[1]
chunks <- ceiling(total.size/chunk.size)
res.all <- c()
pb <- txtProgressBar(min=1, max=chunks, style=3)
for(i in 1:chunks){
  setTxtProgressBar(pb, i)
  index <- (1+chunk.size*(i-1)):(min(chunk.size*i, total.size))
  # Combine cancer and normal data in one frame, using selected samples
  READ.AMP <- rbind(READ.NMP[normal.samples, as.character(READ.linked.probes.genes$probes)[index]], 
                    READ.CMP[cancer.samples, as.character(READ.linked.probes.genes$probes)[index]])
  READ.AEA <- rbind(READ.NEA[normal.samples, as.character(READ.linked.probes.genes$genes)[index]], 
                    READ.CEA[cancer.samples, as.character(READ.linked.probes.genes$genes)[index]]) 
  colnames(READ.AMP) <- paste(READ.linked.probes.genes$probe[index], 
                              READ.linked.probes.genes$genes[index], sep=".")
  colnames(READ.AEA) <- paste(READ.linked.probes.genes$probe[index], 
                              READ.linked.probes.genes$genes[index], sep=".")
  READ.all.data <- rbind(READ.AMP, READ.AEA)

  n <- length(normal.samples) + length(cancer.samples)
  res <- apply(READ.all.data, 2, function(x) { 
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
READ.linear.ME <- as.data.frame(t(res.all))
cat("Adjusting p-values and computing r\n")
READ.linear.ME$p.adj <- p.adjust(READ.linear.ME$p.raw, method="BH")
READ.linear.ME$r <- sqrt(READ.linear.ME$r.squared)*sign(READ.linear.ME$slope)

save(READ.linear.ME, file="../Rdata/READ/calc/READ-linear-ME.Rdata")
quit(save="no")

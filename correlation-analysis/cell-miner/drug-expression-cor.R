library(psych)

cat("Reading drug data..\n")
load("../Rdata/cell-miner/data/cell-miner-drug.Rdata")
cat("Reading expression data..\n")
load("../Rdata/cell-miner/data/cell-miner-expression.Rdata")

cat("Using only common samples..\n")
common.samples <- intersect(rownames(cell.miner.drug),
                            rownames(cell.miner.expression))

cell.miner.drug       <- cell.miner.drug[common.samples,]
cell.miner.expression <- cell.miner.expression[common.samples,]

options("warn"=-1)

cat("Staring computation..\n")
max.drug   <- 21101
chunk.size <- 100
total.size <- dim(cell.miner.expression)[2]
chunks     <- ceiling(total.size/chunk.size)

cell.miner.cor <- data.frame(row=numeric(0), col=numeric(0), r=numeric(0), n=numeric(0), p=numeric(0))
pb <- txtProgressBar(min=1, max=chunks, style=3)
for(i in 1:chunks){
  setTxtProgressBar(pb,i)

  index   <- (1+chunk.size*(i-1)):(min(chunk.size*i, total.size))
  r       <- cor(cell.miner.drug[,1:max.drug], cell.miner.expression[,index], use="p")
  n       <- count.pairwise(cell.miner.drug[,1:max.drug], cell.miner.expression[,index])
  t       <- r*sqrt(n-2)/sqrt(1-r^2)
  p.lower <- pt(t, n-2, lower.tail=T)
  p.upper <- pt(t, n-2, lower.tail=F)
  p.min   <- pmin(p.lower, p.upper)
  p.adj   <- matrix(p.adjust(p.min, method="BH", n=21101*54100), nrow=dim(p.min)[1], ncol=dim(p.min)[2])
  idx     <- (p.adj < 0.01)&!is.na(p.adj)

  tmp <- data.frame(which(idx, arr.ind=T), r=r[idx], n=n[idx], p.adj=p.adj[idx])
  tmp$col <- tmp$col + (i-1)*chunk.size
  rownames(tmp) <- NULL
  cell.miner.cor <- rbind(cell.miner.cor, tmp)
}
cat("\n")

save(cell.miner.cor, file="../Rdata/cell-miner/calc/cell-miner-cor.Rdata")
quit(save="no")

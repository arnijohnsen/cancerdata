# Load files
load("../Rdata/GBM/calc/GBM-linear-ME.Rdata")
load("../Rdata/GBM/calc/GBM-nm-expr.Rdata")
load("../Rdata/GBM/calc/GBM-quantiles.Rdata")

GBM.statistics <- cbind(GBM.linear.ME, GBM.nm.expr, GBM.quantiles)
colnames(GBM.statistics)[6] <- "nm.expr"
save(GBM.statistics, file="../Rdata/GBM/calc/GBM-statistics.Rdata")
quit(save="no")
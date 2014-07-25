# Load files
load("../Rdata/PRAD/calc/PRAD-linear-ME.Rdata")
load("../Rdata/PRAD/calc/PRAD-nm-expr.Rdata")
load("../Rdata/PRAD/calc/PRAD-quantiles.Rdata")

PRAD.statistics <- cbind(PRAD.linear.ME, PRAD.nm.expr, PRAD.quantiles)
colnames(PRAD.statistics)[6] <- "nm.expr"
save(PRAD.statistics, file="../Rdata/PRAD/calc/PRAD-statistics.Rdata")
quit(save="no")
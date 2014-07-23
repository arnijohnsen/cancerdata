# Load files
load("../Rdata/LUSC/calc/LUSC-linear-ME.Rdata")
load("../Rdata/LUSC/calc/LUSC-nm-expr.Rdata")
load("../Rdata/LUSC/calc/LUSC-quantiles.Rdata")

LUSC.statistics <- cbind(LUSC.linear.ME, LUSC.nm.expr, LUSC.quantiles)
colnames(LUSC.statistics)[6] <- "nm.expr"
save(LUSC.statistics, file="../Rdata/LUSC/calc/LUSC-statistics.Rdata")
quit(save="no")

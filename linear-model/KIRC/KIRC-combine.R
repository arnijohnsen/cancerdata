# Load files
load("../Rdata/KIRC/calc/KIRC-linear-ME.Rdata")
load("../Rdata/KIRC/calc/KIRC-nm-expr.Rdata")
load("../Rdata/KIRC/calc/KIRC-quantiles.Rdata")

KIRC.statistics <- cbind(KIRC.linear.ME, KIRC.nm.expr, KIRC.quantiles)
colnames(KIRC.statistics)[6] <- "nm.expr"
save(KIRC.statistics, file="../Rdata/KIRC/calc/KIRC-statistics.Rdata")
quit(save="no")

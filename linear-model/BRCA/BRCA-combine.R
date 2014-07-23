# Load files
load("../Rdata/BRCA/calc/BRCA-linear-ME.Rdata")
load("../Rdata/BRCA/calc/BRCA-nm-expr.Rdata")
load("../Rdata/BRCA/calc/BRCA-quantiles.Rdata")

BRCA.statistics <- cbind(BRCA.linear.ME, BRCA.nm.expr, BRCA.quantiles)
colnames(BRCA.statistics)[6] <- "nm.expr"
save(BRCA.statistics, file="../Rdata/BRCA/calc/BRCA-statistics.Rdata")
quit(save="no")

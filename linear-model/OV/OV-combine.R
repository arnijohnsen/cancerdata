# Load files
load("../Rdata/OV/calc/OV-linear-ME.Rdata")
load("../Rdata/OV/calc/OV-nm-expr.Rdata")
load("../Rdata/OV/calc/OV-quantiles.Rdata")
load("../Rdata/OV/calc/OV-diff-fold.Rdata")

OV.statistics <- cbind(OV.linear.ME, OV.nm.expr, OV.quantiles, OV.diff.fold)
colnames(OV.statistics)[6] <- "nm.expr"
save(OV.statistics, file="../Rdata/OV/calc/OV-statistics.Rdata")
quit(save="no")

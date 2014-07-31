# Load files
load("../Rdata/LUAD/calc/LUAD-linear-ME.Rdata")
load("../Rdata/LUAD/calc/LUAD-nm-expr.Rdata")
load("../Rdata/LUAD/calc/LUAD-quantiles.Rdata")
load("../Rdata/LUAD/calc/LUAD-diff-fold.Rdata")

LUAD.statistics <- cbind(LUAD.linear.ME, LUAD.nm.expr, LUAD.quantiles, LUAD.diff.fold)
colnames(LUAD.statistics)[6] <- "nm.expr"
save(LUAD.statistics, file="../Rdata/LUAD/calc/LUAD-statistics.Rdata")
quit(save="no")
